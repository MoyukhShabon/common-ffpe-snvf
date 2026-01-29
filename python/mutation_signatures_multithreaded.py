"""mutation_signatures.py

This module provides functions for processing variant call files (VCF/SNV),
extracting mutation contexts, and generating mutation signature matrices (e.g., 96-channel).

It utilizes Polars for high-performance data manipulation and Pysam for 
efficient FASTA querying. Thread-local storage is used to handle Pysam's 
non-thread-safe C-backend during parallel Polars execution.

Dependencies:
    - polars
    - pysam
    - numpy
    - seaborn
    - matplotlib

Author: Moyukh Shabon Khan
"""

import math
import threading
from typing import Any, Optional, Union, Dict

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import pysam
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.ticker import MaxNLocator

# local dependency
from common import snv_filter

# Thread-local storage container
# This ensures that every thread created by Polars gets its own
# independent Pysam file handle, preventing race conditions.
_thread_local_data = threading.local()


def extract_trinucleotide_context(row: Dict[str, Any], fasta_path: str) -> str:
    """
    Extracts the 5' and 3' flanking bases for a given variant position.

    This function is designed to be used within a Polars `map_elements` context.
    It manages a thread-local Pysam FastaFile instance to ensure thread safety.

    Args:
        row (Dict[str, Any]): A dictionary containing 'chrom', 'pos', and 'ref'.
        fasta_path (str): Path to the reference genome FASTA file.

    Returns:
        str: The trinucleotide context (e.g., "ACG"). Returns "N{ref}N" on failure.
    """
    # Check if this specific thread already has the genome open
    if not hasattr(_thread_local_data, "genome"):
        _thread_local_data.genome = pysam.FastaFile(fasta_path)

    genome = _thread_local_data.genome

    chrom = row["chrom"]
    pos = row["pos"]
    ref = row["ref"]

    try:
        # Pysam coordinates are 0-based.
        # If pos is 1-based from VCF, we need to be careful.
        # Assuming 'pos' passed here is 1-based standard VCF:
        # Left flank (pos-2) to (pos-1) gets the base before.
        
        chrom_length = genome.get_reference_length(chrom)

        if pos <= 1:
            left_flank = "N"
        else:
            # fetch(region, start, end) -> [start, end)
            left_flank = genome.fetch(chrom, pos - 2, pos - 1)

        if pos >= chrom_length:
            right_flank = "N"
        else:
            right_flank = genome.fetch(chrom, pos, pos + 1)

    except (ValueError, KeyError):
        return f"N{ref}N"

    return (left_flank + ref + right_flank).upper()


def canonicalize_mutation(mutation: str) -> str:
    """
    Maps a mutation to its pyrimidine-centric canonical form (e.g., G>T becomes C>A).

    Args:
        mutation (str): The raw mutation string (e.g., "G>T").

    Returns:
        str: The canonical mutation string.
    """
    complement_to_canonical = {
        "G>T": "C>A",
        "G>C": "C>G",
        "G>A": "C>T",
        "A>T": "T>A",
        "A>G": "T>C",
        "A>C": "T>G"
    }
    return complement_to_canonical.get(mutation, mutation)


def reverse_complement_seq(seq: str) -> str:
    """
    Generates the reverse complement of a DNA sequence.

    Args:
        seq (str): Input DNA sequence.

    Returns:
        str: Reverse complemented sequence.
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement.get(base, base) for base in reversed(seq))


def canonicalize_context(context: str) -> str:
    """
    Converts a trinucleotide context to the canonical pyrimidine reference.

    If the central base is a purine (A or G), the entire trinucleotide
    is reverse complemented.

    Args:
        context (str): The trinucleotide context (e.g., "ATG").

    Returns:
        str: The canonical context.

    Raises:
        TypeError: If input is not a string.
        ValueError: If length is less than 3 or even.
    """
    if not isinstance(context, str):
        raise TypeError("Input context must be a string")
    if len(context) < 3:
        raise ValueError("Input context string length must be at least 3")
    if len(context) % 2 == 0:
        raise ValueError("Input context string length must be odd to have a middle base")

    mid_idx = len(context) // 2
    ref_base = context[mid_idx]

    # If the reference base is Purine (G or A), we must Reverse Complement the WHOLE context
    # to match the Pyrimidine (C or T) mutation conversion.
    if ref_base in ["G", "A"]:
        return reverse_complement_seq(context)

    return context


def annotate_variant_context(
    snv: pl.DataFrame, 
    ref_fasta_path: str, 
    sample_name: Optional[str] = None, 
    c96: bool = True
) -> pl.DataFrame:
    """
    Annotates a SNV DataFrame with mutation types and trinucleotide contexts.

    Args:
        snv (pl.DataFrame): Input DataFrame containing variants.
        ref_fasta_path (str): Path to the reference genome FASTA.
        sample_name (Optional[str]): If provided, adds a 'sample_name' column.
        c96 (bool): If True, converts mutations/contexts to 96-channel canonical format.

    Returns:
        pl.DataFrame: The annotated DataFrame.
    """
    # Apply external filter
    snv = snv.pipe(snv_filter)

    if sample_name:
        snv = snv.with_columns(pl.lit(sample_name).alias("sample_name"))

    # 1. Obtain mutation type i.e. SBS change e.g C>T, G>A
    snv = snv.with_columns([
        pl.concat_str([pl.col("ref"), pl.lit(">"), pl.col("alt")]).alias("mutation_type"),
    ])

    # 2. Obtain mutation trinucleotide context via Thread-Safe Map
    snv = snv.with_columns([
        pl.struct(["chrom", "pos", "ref"])
        .map_elements(
            lambda row: extract_trinucleotide_context(row, ref_fasta_path), 
            return_dtype=pl.String
        )
        .alias("trinucleotide"),
    ])

    if c96:
        snv = (
            snv.with_columns([
                pl.col("mutation_type")
                .map_elements(canonicalize_mutation, return_dtype=pl.String)
                .alias("mutation_type_collapsed"),
                
                pl.col("trinucleotide")
                .map_elements(canonicalize_context, return_dtype=pl.String)
                .alias("trinuclotide_collapsed")
            ])
            .drop(["mutation_type", "trinucleotide"])
            .rename({
                "mutation_type_collapsed": "mutation_type",
                "trinuclotide_collapsed": "trinucleotide"
            })
        )

    return snv


def build_mutation_count_matrix(snv_96c: pl.DataFrame) -> pl.DataFrame:
    """
    Aggregates annotated variants into a 96-channel count matrix.

    Args:
        snv_96c (pl.DataFrame): DataFrame output from `annotate_variant_context`.

    Returns:
        pl.DataFrame: A matrix where rows are (MutationType, Trinucleotide) and 
                      columns are counts (or sample names if multiple samples exist).
    """
    mutation_types = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    bases = ["A", "C", "G", "T"]

    # Generate all 96 combinations
    channels = []
    for mut in mutation_types:
        ref_base = mut[0]
        for left in bases:
            for right in bases:
                trinuc = left + ref_base + right
                channels.append((mut, trinuc))

    channels_df = pl.DataFrame(channels, schema=["mutation_type", "trinucleotide"], orient="row")

    if "sample_name" in snv_96c.columns:
        mut_matrix_96c = (
            snv_96c
            .group_by(["mutation_type", "trinucleotide", "sample_name"])
            .len()
            .pivot(
                values="len",
                index=["mutation_type", "trinucleotide"],
                on="sample_name"
            )
            .fill_null(0)
        )
        # Sort sample columns alphabetically
        index_cols = ["mutation_type", "trinucleotide"]
        sample_cols = sorted([col for col in mut_matrix_96c.columns if col not in index_cols])
        mut_matrix_96c = mut_matrix_96c.select(index_cols + sample_cols).sort(index_cols)

    else:
        mut_matrix_96c = (
            snv_96c
            .group_by(["mutation_type", "trinucleotide"])
            .len()
            .rename({"len": "count"})
            .fill_null(0)
        )

    # Left join with the full 96-channel template to ensure no missing channels
    mut_matrix_96c = (
        channels_df
        .join(mut_matrix_96c, on=["mutation_type", "trinucleotide"], how="left")
        .fill_null(0)
        .sort(["mutation_type", "trinucleotide"])
    )

    return mut_matrix_96c


def plot_sbs96_signature(
    sig: np.ndarray, 
    label: str = "", 
    name: str = "", 
    file: Optional[str] = None, 
    norm: bool = False,
    width: int = 10, 
    height: int = 2, 
    bar_width: float = 1, 
    xticks_label: bool = False, 
    grid: float = 0.2, 
    s: int = 10, 
    ylim: Optional[float] = None, 
    ax: Optional[Axes] = None
) -> None:
    """
    Plots a standard 96-channel Single Base Substitution (SBS) signature.

    Args:
        sig (np.ndarray): Array of 96 counts or frequencies.
        label (str): Label for the Y-axis side strip.
        name (str): Title/Annotation inside the plot.
        file (Optional[str]): Path to save the file. If None, shows plot.
        norm (bool): If True, normalizes the signature to sum to 1.
        width (int): Figure width.
        height (int): Figure height.
        bar_width (float): Width of bars.
        xticks_label (bool): Whether to show x-axis labels (trinucleotides).
        grid (float): Grid line width.
        s (int): Base font size scale.
        ylim (Optional[float]): Manual Y-axis limit.
        ax (Optional[Axes]): Matplotlib Axes object. If None, creates new figure.
    """
    channel = 96
    col_set = ['deepskyblue', 'black', 'red', 'lightgrey', 'yellowgreen', 'pink']
    col_list = []
    for i in range(len(col_set)):
        col_list += [col_set[i]] * 16

    # Determine if standalone or subplot
    if ax is None:
        fig = plt.figure(figsize=(width, height))
        ax = plt.gca()
        is_standalone = True
    else:
        is_standalone = False

    # Set theme
    sns.set_theme(style="whitegrid", color_codes=True,
                  rc={"grid.linewidth": grid, 'grid.color': '.7', 'ytick.major.size': 2,
                      'axes.edgecolor': '.3', 'axes.linewidth': 1.35})

    channel6 = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    channel96 = [
        # 1. C>A
        'ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'CCC', 'CCG', 'CCT',
        'GCA', 'GCC', 'GCG', 'GCT', 'TCA', 'TCC', 'TCG', 'TCT',
        # 2. C>G
        'ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'CCC', 'CCG', 'CCT',
        'GCA', 'GCC', 'GCG', 'GCT', 'TCA', 'TCC', 'TCG', 'TCT',
        # 3. C>T
        'ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'CCC', 'CCG', 'CCT',
        'GCA', 'GCC', 'GCG', 'GCT', 'TCA', 'TCC', 'TCG', 'TCT',
        # 4. T>A
        'ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'CTC', 'CTG', 'CTT',
        'GTA', 'GTC', 'GTG', 'GTT', 'TTA', 'TTC', 'TTG', 'TTT',
        # 5. T>C
        'ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'CTC', 'CTG', 'CTT',
        'GTA', 'GTC', 'GTG', 'GTT', 'TTA', 'TTC', 'TTG', 'TTT',
        # 6. T>G
        'ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'CTC', 'CTG', 'CTT',
        'GTA', 'GTC', 'GTG', 'GTT', 'TTA', 'TTC', 'TTG', 'TTT'
    ]

    # Plot Normalized Version
    if norm:
        normed_sig = sig / np.sum(sig)
        ax.bar(range(channel), normed_sig, width=bar_width, color=col_list)
        ax.set_xticks(range(channel))
        if xticks_label:
            ax.set_xticklabels(channel96, rotation=90, ha="center", va="center", size=7)
        else:
            ax.set_xticklabels([])

        ax.set_ylim(0, np.max(normed_sig) * 1.15)
        ax.annotate(name, (90 - len(name), np.max(sig) * 0.95), size=s)
        ax.set_ylabel("Frequency")

    # Plot Count Version
    else:
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.bar(range(channel), sig, width=bar_width, color=col_list)

        if xticks_label:
            ax.set_xticks(range(channel))
            ax.set_xticklabels(channel96, rotation=90, ha="center", va="center", size=s * 1.15)
        else:
            ax.set_xticks([])
            ax.set_xticklabels([])

        if ylim is not None:
            ax.set_ylim(0, math.ceil(ylim * 1.15))
        else:
            ax.set_ylim(0, math.ceil(np.max(sig) * 1.15))

        if np.round(np.sum(sig)) != 1:
            ax.annotate(
                f"Total Count : {np.sum(sig):,}\nC>T Count : {np.sum(sig[32:48]):,}",
                xy=(0.02, 0.95),
                xycoords='axes fraction',
                size=s * 1.5,
                verticalalignment='top'
            )
        ax.set_ylabel("Number of\nSBSs", size=s * 2)
        ax.annotate(name, (90 - len(name), np.max(sig) * 0.95), size=s * 1.15)

    ax.tick_params(axis='y', labelsize=s * 1.15)

    # Plot Top Colored Strips
    text_col = ["w", "w", "w", "black", "black", "black"]
    for i in range(6):
        left, width_rect = 0 + 1 / 6 * i + 0.001, 1 / 6 - 0.002
        bottom, height_rect = 1.003, 0.15
        right = left + width_rect
        top = bottom + height_rect

        p = plt.Rectangle((left, bottom), width_rect, height_rect, fill=True, color=col_set[i])
        p.set_transform(ax.transAxes)
        p.set_clip_on(False)
        ax.add_patch(p)
        ax.text(0.5 * (left + right), 0.5 * (bottom + top), channel6[i],
                color=text_col[i], weight='bold', size=s * 1.8,
                horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes)

    # Plot Side Label
    if label != "":
        left, width_rect = 1.003, 0.05
        bottom, height_rect = 0, 1
        right = left + width_rect
        top = bottom + height_rect

        p = plt.Rectangle((left, bottom), width_rect, height_rect, fill=True, color="silver", alpha=0.3)
        p.set_transform(ax.transAxes)
        p.set_clip_on(False)
        ax.add_patch(p)
        ax.text(0.505 * (left + right), 0.5 * (bottom + top), label, color="black", size=s * 2,
                horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes, rotation=90)

    ax.margins(x=0.002, y=0.002)

    if is_standalone:
        plt.tight_layout()
        if file:
            plt.savefig(file, bbox_inches="tight", dpi=300)
            plt.close()
        else:
            plt.show()
            plt.close()

