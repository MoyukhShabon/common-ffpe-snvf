"""mutation-signatures.py

This module provides functions for processing variant call files (VCF/SNV),
extracting mutation contexts, and generating mutation signature matrices (e.g., 96-channel).

Functions:
	- import_formatted_vcf: Import and standardize VCF/SNV files.
	- get_context: Extract trinucleotide context for a variant.
	- collapse_mut_96c: Collapse mutation types to canonical 96-channel format.
	- collapse_context_96c: Collapse trinucleotide contexts.
	- variant_mut_profile: Generate mutation profiles with context.
	- create_96c_mutation_matrix: Create a 96-channel mutation count matrix.

Dependencies:
	- polars
	- pysam
	- numpy
	- seaborn
	- matplotlib

Author: Moyukh Shabon Khan
"""

import polars as pl
import pysam
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from common import snv_filter

def import_formatted_vcf(vcf_path, sample_name, snv_only = True, normal_chr=True) -> pl.DataFrame:
	
	"""
	Import a VCF or SNV file and return a Polars DataFrame with standardized columns.

	Parameters
	----------
	vcf_path : str
		Path to the VCF or SNV file.
	snv_only : bool, optional
		If True, filter to include only single nucleotide variants (default: True).

	Returns
	-------
	pl.DataFrame
		DataFrame with columns: chrom, pos, ref, alt, sample_id.
		Multi-allelic variants are exploded into separate rows.
	"""
	
	allowed_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
	
	# print(sample_name)
	
	if vcf_path.endswith(".vcf") or vcf_path.endswith(".vcf.gz"):
		vcf = (
			pl.read_csv(vcf_path, separator="\t", comment_prefix="##", columns=["#CHROM", "POS", "REF", "ALT"])
			# .filter(pl.col("#CHROM").is_in(allowed_chroms))
			.rename({
						"#CHROM" : "chrom",
						"POS" : "pos", 
						"REF": "ref", 
						"ALT" : "alt"
			})
		)
	elif (vcf_path.endswith(".snv") | vcf_path.endswith(".tsv")):
		vcf = pl.read_csv(vcf_path, separator="\t", infer_schema_length=10000)
		vcf = vcf.rename({col: col.lower() for col in vcf.columns})
		
	
	vcf = (
		vcf.with_columns([
			pl.lit(sample_name).alias("sample_id"),
			pl.col("ref").str.to_uppercase(),
			pl.col("alt").str.to_uppercase()
		])
		# .filter(pl.col("chrom").is_in(allowed_chroms))
		.with_columns(
			pl.col("alt").str.split(",")
		)
		.explode("alt")
	)
	
	if normal_chr:
		vcf = vcf.filter(pl.col("chrom").is_in(allowed_chroms))
 
	if snv_only:
		vcf = vcf.filter((pl.col("alt").str.len_chars() == 1) & (pl.col("ref").str.len_chars() == 1))
	
	return vcf



def get_context(row, genome):
	chrom = row["chrom"]
	pos = row["pos"]
	ref = row["ref"]

	# Get chromosome length
	chrom_length = genome.get_reference_length(chrom)

	# Handle 1-based to 0-based conversion and edge cases
	# If left context is out of bounds, use 'N'
	if pos <= 1:
		left_flank = "N"
	else:
		left_flank = genome.fetch(chrom, pos - 2, pos - 1)

	# If right context is out of bounds, use 'N'
	if pos >= chrom_length:
		right_flank = "N"
	else:
		right_flank = genome.fetch(chrom, pos, pos + 1)

	return (left_flank + ref + right_flank).upper()


def collapse_mut_96c(x):
	complement_to_canonical = {
		"G>T": "C>A",
		"G>C": "C>G",
		"G>A": "C>T",
		"A>T": "T>A",
		"A>G": "T>C",
		"A>C": "T>G"
	}
	return complement_to_canonical.get(x, x)
	
def collapse_context_96c(x):
	if not isinstance(x, str):
		raise TypeError("Input context must be a string")
	if len(x) < 3:
		raise ValueError("Input context string length must be at least 3")
	if len(x) % 2 == 0:
		raise ValueError("Input context string length must be odd to have a middle base")
	mid = len(x) // 2
	base_map = {"G": "C", "A": "T"}
	new_mid = base_map.get(x[mid], x[mid])
	return x[:mid] + new_mid + x[mid+1:]


def variants_mut_profile(snv, genome, sample_name = None, c96=True):

	snv = snv.pipe(snv_filter)

	if sample_name:
		snv = snv.with_columns(pl.lit(sample_name).alias("sample_name"))

	## obtain mutation type
	snv = snv.with_columns([
		pl.concat_str([pl.col("ref"), pl.lit(">"), pl.col("alt")]).alias("mutation_type"),
	])
	
	## Obtain mutation trinucleotide context
	snv = snv.with_columns([
		pl.struct(["chrom", "pos", "ref"]).map_elements(lambda row: get_context(row, genome), return_dtype=pl.String).alias("trinucleotide"),
	])
	
	if c96:
		snv = (
			snv.with_columns([
				pl.col("mutation_type").map_elements(collapse_mut_96c, return_dtype=pl.String).alias("mutation_type_collapsed"),
				pl.col("trinucleotide").map_elements(collapse_context_96c, return_dtype=pl.String).alias("trinuclotide_collapsed")
			])
			.drop(["mutation_type", "trinucleotide"])
			.rename({
				"mutation_type_collapsed": "mutation_type",
				"trinuclotide_collapsed": "trinucleotide"
			})
		)

	return snv


def create_96c_mutation_counts(snv_96c):
	mutation_types = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
	bases = ["A", "C", "G", "T"]

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
		# Sort sample columns
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
	
	mut_matrix_96c = channels_df.join(mut_matrix_96c, on=["mutation_type", "trinucleotide"], how="left").fill_null(0).sort(["mutation_type", "trinucleotide"])
	
	return mut_matrix_96c


	
def SBS96_plot(sig, label="", name="", file=None, norm=False,
			   width=10, height=2, bar_width=1, 
			   xticks_label=False, grid=0.2, s=10, ylim=None, ax=None):
	"""
	Modified to accept an 'ax' argument for subplotting.
	"""
	
	channel = 96
	col_set = ['deepskyblue','black','red','lightgrey','yellowgreen','pink']
	col_list = []
	for i in range(len(col_set)):
		col_list += [col_set[i]] * 16
	
	# LOGIC CHANGE: If no ax is provided, create a new figure (standalone mode)
	# If ax IS provided, use it (subplot mode)
	if ax is None:
		fig = plt.figure(figsize=(width, height))
		ax = plt.gca()
		is_standalone = True
	else:
		is_standalone = False

	# Set theme (Note: this sets global theme, might affect other plots if not careful)
	sns.set_theme(style="whitegrid", color_codes=True, 
				  rc={"grid.linewidth": grid, 'grid.color': '.7', 'ytick.major.size': 2,
					  'axes.edgecolor': '.3', 'axes.linewidth': 1.35,})
	

	channel6 = ['C>A','C>G','C>T','T>A','T>C','T>G']
	channel96 = ['ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'CCC', 'CCG', 'CCT', 'GCA',
			   'GCC', 'GCG', 'GCT', 'TCA', 'TCC', 'TCG', 'TCT', 'ACA', 'ACC',
			   'ACG', 'ACT', 'CCA', 'CCC', 'CCG', 'CCT', 'GCA', 'GCC', 'GCG',
			   'GCT', 'TCA', 'TCC', 'TCG', 'TCT', 'ACA', 'ACC', 'ACG', 'ACT',
			   'CCA', 'CCC', 'CCG', 'CCT', 'GCA', 'GCC', 'GCG', 'GCT', 'TCA',
			   'TCC', 'TCG', 'TCT', 'ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'CTC',
			   'CTG', 'CTT', 'GTA', 'GTC', 'GTG', 'GTT', 'TTA', 'TTC', 'TTG',
			   'TTT', 'ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'CTC', 'CTG', 'CTT',
			   'GTA', 'GTC', 'GTG', 'GTT', 'TTA', 'TTC', 'TTG', 'TTT']
	
	## plot the normalized version:
	if norm:
		normed_sig = sig / np.sum(sig)
		ax.bar(range(channel), normed_sig, width=bar_width, color=col_list)
		# Handle xticks
		ax.set_xticks(range(channel))
		if xticks_label:
			ax.set_xticklabels(channel96, rotation=90, ha="center", va="center", size=7)
		else:
			ax.set_xticklabels([])
			
		ax.set_ylim(0, np.max(normed_sig) * 1.15)
		ax.annotate(name, (90 - len(name), np.max(sig) * 0.95), size=s)
		ax.set_ylabel("Frequency")

	## plot the original version:
	else:
		ax.yaxis.set_major_locator(MaxNLocator(integer=True))
		ax.bar(range(channel), sig, width=bar_width, color=col_list)
		
		# Handle xticks
		if xticks_label:
			ax.set_xticks(range(channel))
			ax.set_xticklabels(channel96, rotation=90, ha="center", va="center", size=7)
		else:
			ax.set_xticks([])
			ax.set_xticklabels([])

		if ylim is not None:
			ax.set_ylim(0, int(ylim * 1.15))
		else:
			ax.set_ylim(0, int(np.max(sig) * 1.15))
		
		if np.round(np.sum(sig)) != 1:
			ax.annotate(
				f"Total Count : {np.sum(sig):,}\nC>T Count : {np.sum(sig[32:48]):,}", 
				xy=(0.02, 0.95), 
				xycoords='axes fraction', 
				size=s*1.5,
				verticalalignment='top'
			)
		ax.set_ylabel("Number of\nSBSs", size=s*2)
		ax.annotate(name, (90 - len(name), np.max(sig) * 0.95), size=s)
	
	ax.tick_params(axis='y', labelsize=s)
	
	## plot the bar annotation (Top colored strips):
	text_col = ["w","w","w","black","black","black"]
	for i in range(6):       
		left, width = 0 + 1/6 * i + 0.001, 1/6 - 0.002
		bottom, height = 1.003, 0.15
		right = left + width
		top = bottom + height
		
		# Use ax.add_patch instead of plt.gca().add_patch
		p = plt.Rectangle((left, bottom), width, height, fill=True, color=col_set[i])
		p.set_transform(ax.transAxes)
		p.set_clip_on(False)
		ax.add_patch(p)
		ax.text(0.5 * (left + right), 0.5 * (bottom + top), channel6[i], 
				color=text_col[i], weight='bold', size=s*1.8,
				horizontalalignment='center', verticalalignment='center', 
				transform=ax.transAxes)
	
	## plot the name annotation (Side label)
	if label != "":
		left, width = 1.003, 0.05
		bottom, height = 0, 1
		right = left + width
		top = bottom + height
		
		p = plt.Rectangle((left, bottom), width, height, fill=True, color="silver", alpha=0.3)
		p.set_transform(ax.transAxes)
		p.set_clip_on(False)
		ax.add_patch(p)
		ax.text(0.505 * (left + right), 0.5 * (bottom + top), label, color="black", size=s*2,
				horizontalalignment='center', verticalalignment='center',
				transform=ax.transAxes, rotation=90)

	ax.margins(x=0.002, y=0.002)
	
	# Only save/show if we are NOT in subplot mode
	if is_standalone:
		plt.tight_layout()
		if file:
			plt.savefig(file, bbox_inches="tight", dpi=300)
			plt.close()
		else:
			plt.show()
			plt.close()
	
## A polars mask for filtering C:G>T:A mutations. Invertable using ~
## Usage df.filter(ct_mask) or df.filter(~ct_mask). Where df: pl.DataFrame
ct_mask = (((pl.col("ref") == "C") & (pl.col("alt") == "T")) | ((pl.col("ref") == "G") & (pl.col("alt") == "A")))

