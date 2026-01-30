import polars as pl
import os
import pysam

## Local dependencies
from common import snv_filter, ct_filter, natural_sort_variants, read_variants

def get_neighborhood_sequence(struct_series: pl.Series, fasta_path: str) -> pl.Series:
    """
    Extracts neighborhood sequences for variants using a single fetch per variant.

    This function is designed to be used within a Polars `map_batches` context.
    It constructs a sequence consisting of:
    [20bp Left Flank] + [ALT Allele] + [20bp Right Flank].

    Note:
        This function instantiates a `pysam.FastaFile` object locally within the 
        function scope. This ensures thread safety when running in parallel 
        (e.g., via Polars threading), avoiding race conditions associated with 
        sharing a single pysam pointer.

    Args:
        struct_series (pl.Series): A Polars Series of Structs containing the fields:
            - "Chr" (str): Chromosome name.
            - "Pos" (int): 1-based position.
            - "Ref" (str): Reference allele.
            - "Alt" (str): Alternative allele.
        fasta_path (str): Absolute path to the reference genome FASTA file.

    Returns:
        pl.Series: A Series of strings containing the constructed neighborhood sequences.
                   Returns an empty string or raises an error if fetching fails.

    Raises:
        RuntimeError: If Pysam fails to fetch the sequence for a specific coordinate.
    """
    # Deconstruct columns to Python lists for speed
    chroms = struct_series.struct.field("Chr").to_list()
    poss = struct_series.struct.field("Pos").to_list()
    refs = struct_series.struct.field("Ref").to_list()
    alts = struct_series.struct.field("Alt").to_list()

    results = []

    # Open FASTA handle locally (thread-safe)
    with pysam.FastaFile(fasta_path) as genome:
        
        for chrom, pos, ref, alt in zip(chroms, poss, refs, alts):
            
            ref_len = len(ref)
            
            # --- Coordinate Calculation ---
            # We want to fetch: [20bp Left] + [Ref Sequence] + [20bp Right]
            # Pysam is 0-based, half-open [start, end)
            
            # Start: (Pos - 1) moves to 0-based index of the variant start. 
            # Subtract 20 to get the start of the left flank.
            # End: (Pos - 1) + ref_len is the end of the Ref. 
            # Add 20 to get the end of the right flank.
            fetch_start = (pos - 1 - 20)
            fetch_end = (pos - 1) + ref_len + 20

            # --- Single Fetch ---
            try:
                # Fetches "LeftFlank + Ref + RightFlank" in one contiguous block
                wide_context = genome.fetch(chrom, fetch_start, fetch_end).upper()
            except (ValueError, IndexError):
                raise RuntimeError(f"Pysam failed to fetch {chrom}:{fetch_start}-{fetch_end}")

            # --- String Assembly ---
            # 1. Take the first 20 chars (Left Flank)
            # 2. Add the ALT allele
            # 3. Skip the REF part (start at 20 + len(ref)) and take the rest (Right Flank)
            full_seq = wide_context[:20] + alt + wide_context[20 + ref_len:]
            results.append(full_seq.upper())

    return pl.Series(results)


def get_read_length(bam_path: str, max_reads: int = 100000) -> int:
    """
    Calculates the mode read length from a BAM file.

    Iterates through the first `max_reads` mapped reads in the BAM file
    to determine the most common query length.

    Args:
        bam_path (str): Path to the BAM file.
        max_reads (int, optional): The maximum number of reads to inspect. 
            Defaults to 100,000.

    Returns:
        int: The mode (most frequent) read length found in the sample.

    Raises:
        ValueError: If no mapped reads are found within the first `max_reads`.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        lengths = []
        
        for i, read in enumerate(bam):
            if read.is_unmapped:
                continue
            if i >= max_reads:
                break
            lengths.append(read.query_length)
            
    if not lengths:
        raise ValueError("Read Length could not be retrieved from BAM")
        
    return pl.Series(lengths).mode()[0]


def annotate_simple_repeats(mut_info: pl.DataFrame, simple_repeats_df: pl.DataFrame) -> pl.DataFrame:
    """
    Annotates variants indicating if they fall within a simple repeat region.

    Performs a backward search using `join_asof` to find the nearest preceding 
    repeat region start, then checks if the variant position falls within that 
    region's end coordinate.

    Args:
        mut_info (pl.DataFrame): DataFrame containing variants. Must have "Chr" and "Pos" columns.
        simple_repeats_df (pl.DataFrame): DataFrame containing repeat regions. 
            Must have "chrom", "start", and "end" columns.

    Returns:
        pl.DataFrame: The input `mut_info` DataFrame with an additional column:
            - "SimpleRepeat_TRF": "Y" if inside a repeat, "N" otherwise.
    """
    
    # Sort both DataFrames
    # join_asof requires the 'on' column (Pos/start) to be sorted within the 'by' groups (Chr).
    # Sorting by ["Chr", "Pos"] achieves this and keeps the data human-readable.
    variants_sorted = natural_sort_variants(mut_info, chr_col="Chr", pos_col="Pos")
    repeats_sorted = simple_repeats_df.sort(["chrom", "start"])

    # Perform the ASOF Join
    # We find the repeat segment that started most recently before (or exactly at) the variant position.
    joined = variants_sorted.join_asof(
        repeats_sorted,
        left_on="Pos",
        right_on="start",
        by_left="Chr",
        by_right="chrom",
        strategy="backward" # Look for the closest start position <= variant position
    )

    # Check overlap
    # We found the closest start. Now, check if the variant actually inside that segment,
    # i.e., Variant Pos <= Repeat End?
    final_df = joined.with_columns(
        pl.when(pl.col("Pos") <= pl.col("end"))
        .then(pl.lit("Y"))
        .otherwise(pl.lit("N"))
        .alias("SimpleRepeat_TRF")
    ).drop(["start", "end"]) # Clean up BED columns

    return final_df


def make_mut_info(
    vcf_path: str, 
    ref_fasta_path: str, 
    sample_name: str, 
    simple_repeats_bed: pl.DataFrame | None = None, 
    ct_only: bool = False, 
    snv_only: bool = False
) -> pl.DataFrame:
    """
    Generates a mutation info DataFrame for MicroSEC input from a VCF file.

    This function orchestrates the following steps:
    1. Reads variants from the VCF.
    2. Filters based on SNV status or C>T changes (optional).
    3. Filters variants close to chromosome boundaries.
    4. Fetches neighborhood sequences (Left+Alt+Right).
    5. Classifies mutation types (SNV, Ins, Del).
    6. Annotates simple repeats (optional).

    Args:
        vcf_path (str): Path to the input VCF file.
        ref_fasta_path (str): Path to the reference genome FASTA.
        sample_name (str): Name of the sample (added as a column).
        simple_repeats_bed (pl.DataFrame | None, optional): DataFrame of simple repeats 
            for annotation. Defaults to None.
        ct_only (bool, optional): If True, filters for C>T (or G>A) mutations only. 
            Defaults to False.
        snv_only (bool, optional): If True, filters for Single Nucleotide Variants only. 
            Defaults to False.

    Returns:
        pl.DataFrame: A DataFrame containing the following columns:
            - Sample, Mut_type, Chr, Pos, Ref, Alt, SimpleRepeat_TRF, Neighborhood_sequence
    """
    
    mut_info = read_variants(vcf_path)

    if snv_only:
        mut_info = snv_filter(mut_info)
    if ct_only:
        mut_info = ct_filter(mut_info)

    mut_info = mut_info.rename({
        "chrom": "Chr",
        "pos": "Pos",
        "ref": "Ref",
        "alt": "Alt"
    })
 
    mut_info = mut_info.with_columns(pl.lit(sample_name).alias("Sample"))

    # Filter by chromosome boundaries
    # We open the genome here just to get lengths, which is fast and safe (read-only metadata)
    with pysam.FastaFile(ref_fasta_path) as genome:
        chrom_lengths = pl.DataFrame({
            "Chr": genome.references,
            "chr_len_ref": genome.lengths
        })

    mut_info = (
        mut_info.join(chrom_lengths, on="Chr", how="inner")
        .filter(
            (pl.col("Pos") > 200) & 
            (pl.col("Pos") < (pl.col("chr_len_ref") - 200))
        )
        .drop("chr_len_ref")
    )

    # Optimized Batch Processing
    # pass the PATH to the fasta, not the object.
    mut_info = mut_info.with_columns(
        pl.struct(["Chr", "Pos", "Ref", "Alt"])
        .map_batches(
            lambda batch: get_neighborhood_sequence(batch, ref_fasta_path), 
            return_dtype=pl.Utf8
        )
        .alias("Neighborhood_sequence")
    )
    
    # Mutation Type Logic
    snv_mask = mut_info["Ref"].str.len_chars() == mut_info["Alt"].str.len_chars()
    deletion_mask = (mut_info["Ref"].str.len_chars() > 1) & (mut_info["Alt"].str.len_chars() == 1)
    insertion_mask = (mut_info["Ref"].str.len_chars() == 1) & (mut_info["Alt"].str.len_chars() > 1)
    
    mut_info = (
        mut_info.with_columns([
            pl.when(snv_mask)
            .then(mut_info["Ref"].str.len_chars().cast(pl.Utf8) + "-snv")
            .when(deletion_mask)
            .then((mut_info["Ref"].str.len_chars() - 1).cast(pl.Utf8) + "-del")
            .when(insertion_mask)
            .then((mut_info["Alt"].str.len_chars() - 1).cast(pl.Utf8) + "-ins")
            .otherwise(pl.lit("-"))
            .alias("Mut_type")
        ])
    )

    if simple_repeats_bed is not None:
        mut_info = annotate_simple_repeats(mut_info, simple_repeats_bed)
    else:
        mut_info = mut_info.with_columns(pl.lit("-").alias("SimpleRepeat_TRF"))

    mut_info = mut_info.select([
        "Sample",
        "Mut_type",
        "Chr",
        "Pos",
        "Ref",
        "Alt",
        "SimpleRepeat_TRF",
        "Neighborhood_sequence"
    ])

    return mut_info