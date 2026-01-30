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

def extract_trinucleotide_context(struct_series: pl.Series, fasta_path: str) -> pl.Series:
	"""
	Extracts the 5' and 3' flanking bases for a given variant position.

	Optimized batch processing for trinucleotide context extraction.


	Args:
		row (Dict[str, Any]): A dictionary containing 'chrom', 'pos', and 'ref'.
		fasta_path (str): Path to the reference genome FASTA file.

	Returns:
		str: The trinucleotide context (e.g., "ACG"). Returns "N{ref}N" on failure.
	
	"""
	#  Extract columns to Python lists
	chroms = struct_series.struct.field("chrom").to_list()
	poss = struct_series.struct.field("pos").to_list()
	refs = struct_series.struct.field("ref").to_list()
	
	results = []
	
	# Open FASTA handle locally. 
	# Since Polars runs map_batches in parallel threads, this creates 
	# a unique file handle per thread/batch. No race conditions.
	with pysam.FastaFile(fasta_path) as genome:
		
		
		for chrom, pos, ref in zip(chroms, poss, refs):
			try:
				# Validation: Check chromosome existence/bounds
				chrom_len = genome.get_reference_length(chrom)
				
				# 0-based logic for pysam (VCF is 1-based)
				# Left flank: VCF pos - 2
				# Right flank: VCF pos
				if pos <= 1 or pos >= chrom_len:
					results.append(f"N{ref}N")
					continue

				# Fetch 3 bases at once: (pos-2) to (pos+1)
				# region is [start, end)
				seq = genome.fetch(chrom, pos - 2, pos + 1)
				
				# Sanity check: ensure the middle base matches ref
				if seq[1].upper() != ref.upper():
					raise ValueError("Middle Base of Trinucleotide does not match Reference Base")

				results.append(seq.upper())

			except (ValueError, KeyError, IndexError):
				results.append(f"N{ref}N")

	return pl.Series(results)


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

	# Obtain mutation type i.e. SBS change e.g C>T, G>A
	snv = snv.with_columns([
		pl.concat_str([pl.col("ref"), pl.lit(">"), pl.col("alt")]).alias("mutation_type"),
	])

	# Obtain mutation trinucleotide context via batch processing, 
	# Passing the struct of columns ["chrom", "pos", "ref"] to map_batches.
	snv = snv.with_columns([
		pl.struct(["chrom", "pos", "ref"])
		.map_batches(
			lambda batch: extract_trinucleotide_context(batch, ref_fasta_path), 
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


def build_mutation_count_matrix(snv_96c: pl.DataFrame, sample_name_col: str = "sample_name") -> pl.DataFrame:
	"""
	Aggregates annotated variants into a 96-channel count matrix.

	Args:
		snv_96c (pl.DataFrame): DataFrame output from `annotate_variant_context`.
		sample_name_col (str): Name if column with sample_name annotation, if exists.

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

	if sample_name_col in snv_96c.columns:
		mut_matrix_96c = (
			snv_96c
			.group_by(["mutation_type", "trinucleotide", sample_name_col])
			.len()
			.pivot(
				values="len",
				index=["mutation_type", "trinucleotide"],
				on=sample_name_col
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
	base_fontsize: int = 10,
	ylim: Optional[float] = None, 
	ax: Optional[Axes] = None
) -> None:
	"""
	Plots a standard 96-channel Single Base Substitution (SBS) signature.
	"""
	
	# Font Configuration ---
	fonts = {
		"xtick": base_fontsize * 0.7,      # Trinucleotide labels (very small)
		"ytick": base_fontsize * 1.2,      # Y-axis numbers
		"title": base_fontsize * 1.2,      # Inner plot annotation
		"stats": base_fontsize * 1.5,      # "Total Count" text
		"strip": base_fontsize * 1.8,      # Top colored strip labels (C>A, etc)
		"axis_label": base_fontsize * 2.0, # Y-axis label ("Frequency")
		"side_label": base_fontsize * 2.2  # side-label ("Frequency")
	}

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

	# Set theme and specific rcParams for this plot
	sns.set_theme(style="whitegrid", color_codes=True,
				  rc={"grid.linewidth": grid, 'grid.color': '.7', 'ytick.major.size': 2,
					  'axes.edgecolor': '.3', 'axes.linewidth': 1.35})

	channel6 = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
	bases = "ACGT"
	muts = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
	channel96 = []
	for m in muts:
		ref = m[0]
		for b1 in bases:
			for b2 in bases:
				channel96.append(f"{b1}{ref}{b2}")

	# --- Plotting Logic ---

	# Plot Normalized Version
	if norm:
		normed_sig = sig / np.sum(sig)
		ax.bar(range(channel), normed_sig, width=bar_width, color=col_list)
		ax.set_xticks(range(channel))
		
		if xticks_label:
			ax.set_xticklabels(channel96, rotation=90, ha="center", va="center", 
							   fontsize=fonts["xtick"])
		else:
			ax.set_xticklabels([])

		ax.set_ylim(0, np.max(normed_sig) * 1.15)
		ax.annotate(name, (90 - len(name), np.max(sig) * 0.95), fontsize=fonts["title"])
		ax.set_ylabel("Frequency", fontsize=fonts["axis_label"])

	# Plot Count Version
	else:
		ax.yaxis.set_major_locator(MaxNLocator(integer=True))
		ax.bar(range(channel), sig, width=bar_width, color=col_list)

		if xticks_label:
			ax.set_xticks(range(channel))
			ax.set_xticklabels(channel96, rotation=90, ha="center", va="center", 
							   fontsize=fonts["xtick"])
		else:
			ax.set_xticks([])
			ax.set_xticklabels([])

		if ylim is not None:
			ax.set_ylim(0, math.ceil(ylim * 1.15))
		else:
			ax.set_ylim(0, math.ceil(np.max(sig) * 1.15))

		if np.round(np.sum(sig)) != 1:
			stats_text = f"Total Count : {np.sum(sig):,}\nC>T Count : {np.sum(sig[32:48]):,}"
			ax.annotate(
				stats_text,
				xy=(0.02, 0.95),
				xycoords='axes fraction',
				fontsize=fonts["stats"],
				verticalalignment='top'
			)
		
		ax.set_ylabel("Number of\nSBSs", fontsize=fonts["axis_label"])
		ax.annotate(name, (90 - len(name), np.max(sig) * 0.95), fontsize=fonts["title"])

	ax.tick_params(axis='y', labelsize=fonts["ytick"])

	# --- Annotations (Top Strips) ---
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
		
		# Using fonts["strip"] here
		ax.text(0.5 * (left + right), 0.5 * (bottom + top), channel6[i],
				color=text_col[i], weight='bold', fontsize=fonts["strip"],
				horizontalalignment='center', verticalalignment='center',
				transform=ax.transAxes)

	# --- Side Label ---
	if label != "":
		left, width_rect = 1.003, 0.05
		bottom, height_rect = 0, 1
		right = left + width_rect
		top = bottom + height_rect

		p = plt.Rectangle((left, bottom), width_rect, height_rect, fill=True, color="silver", alpha=0.3)
		p.set_transform(ax.transAxes)
		p.set_clip_on(False)
		ax.add_patch(p)
		
		# Using fonts["axis_label"] here
		ax.text(0.505 * (left + right), 0.5 * (bottom + top), label, color="black", 
				fontsize=fonts["axis_label"],
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
			
