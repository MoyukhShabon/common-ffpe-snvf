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
import math
from matplotlib.ticker import MaxNLocator
from common import snv_filter
import threading

def get_context(row, genome, lock):
    chrom = row["chrom"]
    pos = row["pos"]
    ref = row["ref"]

    # Fetching sequence is NOT thread-safe in pysam
    with lock:
        try:
            chrom_length = genome.get_reference_length(chrom)
            
            # Handle 1-based to 0-based conversion and edge cases
            if pos <= 1:
                left_flank = "N"
            else:
                # fetch is 0-based, half-open [start, end)
                left_flank = genome.fetch(chrom, pos - 2, pos - 1)

            if pos >= chrom_length:
                right_flank = "N"
            else:
                right_flank = genome.fetch(chrom, pos, pos + 1)
        except ValueError:
            # Handle cases where chrom might not exist in fasta
            return "N" + ref + "N"

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

## Bug Identified: The whole sequence should be complemented.
# def collapse_context_96c(x):
# 	if not isinstance(x, str):
# 		raise TypeError("Input context must be a string")
# 	if len(x) < 3:
# 		raise ValueError("Input context string length must be at least 3")
# 	if len(x) % 2 == 0:
# 		raise ValueError("Input context string length must be odd to have a middle base")
# 	mid = len(x) // 2
# 	base_map = {"G": "C", "A": "T"}
# 	new_mid = base_map.get(x[mid], x[mid])
# 	return x[:mid] + new_mid + x[mid+1:]

def reverse_complement_seq(seq):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement.get(base, base) for base in reversed(seq))

def collapse_context_96c(x):
    if not isinstance(x, str):
        raise TypeError("Input context must be a string")
    if len(x) < 3:
        raise ValueError("Input context string length must be at least 3")
    if len(x) % 2 == 0:
        raise ValueError("Input context string length must be odd to have a middle base")
    
    mid_idx = len(x) // 2
    ref_base = x[mid_idx]
    
    # If the reference base is Purine (G or A), we must Reverse Complement the WHOLE context
    # to match the Pyrimidine (C or T) mutation conversion.
    if ref_base in ["G", "A"]:
        return reverse_complement_seq(x)
    
    return x


def variants_mut_profile(snv, genome, sample_name = None, c96=True):

	genome_lock = threading.Lock()

	snv = snv.pipe(snv_filter)

	if sample_name:
		snv = snv.with_columns(pl.lit(sample_name).alias("sample_name"))

	## obtain mutation type i.e. SBS change e.g C>T, G>A rtc
	snv = snv.with_columns([
		pl.concat_str([pl.col("ref"), pl.lit(">"), pl.col("alt")]).alias("mutation_type"),
	])
	
	## Obtain mutation trinucleotide context
	snv = snv.with_columns([
		pl.struct(["chrom", "pos", "ref"]).map_elements(lambda row: get_context(row, genome, genome_lock), return_dtype=pl.String).alias("trinucleotide"),
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
	channel96 =     channel96 = [
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
			ax.set_xticklabels(channel96, rotation=90, ha="center", va="center", size=s*1.15)
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
				size=s*1.5,
				verticalalignment='top'
			)
		ax.set_ylabel("Number of\nSBSs", size=s*2)
		ax.annotate(name, (90 - len(name), np.max(sig) * 0.95), size=s*1.15)
	
	ax.tick_params(axis='y', labelsize=s*1.15)
	
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

