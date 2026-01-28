import polars as pl
import numpy as np
from scipy.stats import false_discovery_control
from sklearn.metrics import confusion_matrix


def read_variants(path:str, columns: list = ["#CHROM", "POS", "REF", "ALT"]) -> pl.DataFrame:
	"""
	Reads Variants from a VCF file into a Polars DataFrame.
	By default, only reads essential columns.
	Splits multi-allelic variants into separate rows.
	Can be extended to read additional columns as needed or
	read custom variant files with similar structure.
	
	:param path: Path to VCF file
	:type path: str
	:param columns: List of column names to read from the VCF file. Defaults to ["#CHROM", "POS", "REF", "ALT"].
	:type columns: list
	:return: DataFrame containing the variants
	:rtype: DataFrame
	"""
	variants = (
		pl.read_csv(path, separator="\t", comment_prefix="##", infer_schema_length=1000, columns=columns)
		.rename(lambda x: x.lstrip("#").lower())
		.with_columns(pl.col("alt").str.split(","))
		.explode("alt")
	)
	return variants

def snv_filter(variants: pl.DataFrame) -> pl.DataFrame:
	return variants.filter(
		(pl.col("ref").str.len_chars() == 1) &
		(pl.col("alt").str.len_chars() == 1)
	)

def ct_filter(variants: pl.DataFrame, invert=False) -> pl.DataFrame:
	mask = (
		((pl.col("ref").str.to_uppercase() == "C") & (pl.col("alt").str.to_uppercase() == "T")) |
		((pl.col("ref").str.to_uppercase() == "G") & (pl.col("alt").str.to_uppercase() == "A"))
	)
	if invert:
		return variants.filter(~mask)
	else:
		return variants.filter(mask)


def natural_sort_variants(df: pl.DataFrame, chr_col: str = "chrom", pos_col: str = "pos") -> pl.DataFrame:
	"""
	Docstring for natural_sort_variants
	
	:param df: Polars DataFrame containing variants
	:type df: pl.DataFrame
	:param chr_col: Name of the chromosome column
	:type chr_col: str
	:param pos_col: Name of the position column
	:type pos_col: str
	:return: DataFrame containing the variants sorted naturally by chromosome and position
	:rtype: DataFrame
	"""
	return (
		df.with_columns(
			# Create a temporary column 'chr_rank' for sorting
			pl.col(chr_col)
			.str.replace("chr", "") # Remove 'chr' prefix
			.str.replace("X", "23") # Handle Sex chromosomes
			.str.replace("Y", "24")
			.str.replace("M", "25") # Handle Mitochondria if present
			.str.replace("MT", "26")
			.cast(pl.Int32, strict=False) # Convert to Integer (strict=False turns unknown contigs to null)
			.fill_null(999) # Put weird contigs at the end
			.alias("chr_rank")
		)
		.sort(["chr_rank", pos_col]) # Sort
		.drop("chr_rank")
	)


def adaptive_fdr_cut(df: pl.DataFrame, fp_cut: float, score_col: str = "q") -> pl.DataFrame:
	"""
	Takes a DataFrame containing a 'q' column, calculates the adaptive cutoff metric,
	and adds a boolean 'pred' column using Polars' built-in rank method. Applies a 
	cutoff where the expected number of false positives in the selected set is less 
	than a specific number defined by fp_cut
	
	:param df: DataFrame containing model (MOBSNVF) scores
	:type df: pl.DataFrame
	:param fp_cut: False positive cutoff threshold
	:type fp_cut: float
	:param score_col: Name of the score column
	:type score_col: str
	:return: DataFrame with a boolean 'pred' column indicating predictions
	:rtype: DataFrame
	"""
	return df.with_columns(
		((pl.col(score_col).rank(method="ordinal") * pl.col(score_col)) < fp_cut).alias("pred")
	)


def fdr_cut_pred(df: pl.DataFrame, score_col: str, fp_cut: float = 0.5) -> pl.DataFrame:
	"""
	Labels model predictions based on adaptive FDR cut using SciPy for BH correction.
	
	:param df: Polars DataFrame containing model (MOBSNVF) scores
	:type df: pl.DataFrame
	:param score_col: Name of the score column
	:type score_col: str
	:param fp_cut: False positive cutoff threshold
	:type fp_cut: float
	:return: DataFrame with a boolean 'pred' column indicating predictions
	:rtype: DataFrame
	"""
	
	# Split Complete Cases (C>T) and Nulls
	df_ct = df.filter(pl.col(score_col).is_not_null())
	df_nct = df.filter(pl.col(score_col).is_null())
	
	# Handle Non-C>T mutations (Null scores)
	# Non-C>T mutations are assumed to be real mutations in this context, hence we set pred = True
	df_nct = df_nct.with_columns([
		pl.lit(None, dtype=pl.Float64).alias(score_col),
		pl.lit(None, dtype=pl.Float64).alias("q"),
		pl.lit(True).alias("pred")
	])
		
	# Handle C>T mutations (Scores exist)
	if df_ct.height > 0:

		# Define machine epsilon for float type (smallest positive float with which 1.0 + eps != 1.0)
		machine_eps = np.finfo(float).eps
		
		df_ct = (
			df_ct
			## Substitute zeros with machine epsilon 
			.with_columns(
				pl.when(pl.col(score_col) == 0)
				.then(machine_eps)
				.otherwise(pl.col(score_col))
				.alias(score_col)
			)
			## Calculate q-values using BH correction from SciPy
			.with_columns(
				pl.col(score_col).map_batches(
					lambda x: false_discovery_control(x, method='bh'), 
					return_dtype=pl.Float64
				).alias("q")
			)
			## Apply adaptive FDR cut
			.pipe(lambda df: adaptive_fdr_cut(df, fp_cut))		
		)
	

	# 4. Combine and sort
	final_df = pl.concat([df_ct, df_nct], how="vertical").pipe(natural_sort_variants)
	
	return final_df


def get_eval_metrics(truth: pl.Series, pred: pl.Series) -> dict:
	"""
	Calculates evaluation metrics based on truth and predicted labels.
	
	:param truth: Series containing true labels
	:type truth: pl.Series
	:param pred: Series containing predicted labels
	:type pred: pl.Series
	:return: Dictionary containing evaluation metrics
	:rtype: dict
	"""

	cm = confusion_matrix(truth, pred)
	TN, FP, FN, TP = cm.ravel()

	eval_metrics = {
		"precision": TP / (TP + FP),
		"recall": TP / (TP + FN),
		"specificity": TN / (TN + FP),
		"TP": TP,
		"TN": TN,
		"FP": FP,
		"FN": FN,
		"confusion_matrix": cm
	}

	return eval_metrics
