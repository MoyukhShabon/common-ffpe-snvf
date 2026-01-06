import polars as pl
import numpy as np
from scipy.stats import false_discovery_control
from sklearn.metrics import confusion_matrix


def read_variants(path:str, columns: list = ["#CHROM", "POS", "REF", "ALT", "FILTER"]) -> pl.DataFrame:
	variants = (
		pl.read_csv(path, separator="\t", comment_prefix="##", infer_schema_length=1000, columns=columns)
		.rename(lambda x: x.lstrip("#").lower())
		.with_columns(pl.col("alt").str.split(","))
		.explode("alt")
	)
	return variants


def natural_sort_variants(df: pl.DataFrame, chr_col: str = "chrom", pos_col: str = "pos") -> pl.DataFrame:
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


def adaptive_fdr_cut(df: pl.DataFrame, fp_cut: float) -> pl.DataFrame:
	"""
	Takes a DataFrame containing a 'q' column, calculates the adaptive cutoff metric,
	and adds a boolean 'pred' column using Polars' built-in rank method. Applies a 
	cutoff where the expected number of false positives in the selected set is less 
	than a specific number defined by fp_cut
	"""
	return df.with_columns(
		((pl.col("q").rank(method="ordinal") * pl.col("q")) < fp_cut).alias("pred")
	)



def fdr_cut_pred(df: pl.DataFrame, score_col: str, fp_cut: float = 0.5) -> pl.DataFrame:
	"""
	Labels model predictions based on adaptive FDR cut using SciPy for BH correction.
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
