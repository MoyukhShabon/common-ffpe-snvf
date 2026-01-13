import polars as pl
import numpy as np
from scipy.stats import false_discovery_control 
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay 
from tqdm import tqdm
import glob
import os


# %% [markdown]
# ## Helper Functions


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
		.sort(["chr_rank", pos_col]) # Sort by Rank, then Position
		.drop("chr_rank") # Remove the temp column
	)


def predict(df, score_col, threshold):
	return df.with_columns(
		## Lower scores from MOBSNVF represents real muatations
		(pl.col(score_col) < threshold).alias("pred")
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


def adaptive_fdr_pred(df: pl.DataFrame, score_col: str, cutoff: float = 0.5, ct_only: bool = True) -> pl.DataFrame:
	"""
	Applies FDR Benjamini-Hochberg correction on the specified score column and
	generates predictions from the the q-values using an adaptive false positive cutoff.
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
		

	if df_ct.height == 0:
		return df_ct if ct_only else df_nct

	# Handle C>T mutations (Scores exist)
	## Define machine epsilon for float type (smallest positive float with which 1.0 + eps != 1.0)
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
		## -x-Apply adaptive FDR cut based on -x-
		## Predict based on q value cutoff specified
		.pipe(lambda df: predict(df, "q", cutoff))
	)
	
	if ct_only:
		return df_ct

	# Combine and sort
	final_df = pl.concat([df_ct, df_nct], how="vertical").pipe(natural_sort_variants)
	
	return final_df



def get_eval_metrics(truth: pl.Series, pred: pl.Series) -> dict:

	cm = confusion_matrix(truth, pred)
	TN, FP, FN, TP = cm.ravel()

	eval_metrics = {
		"precision": TP / (TP + FP),
		"recall": TP / (TP + FN),
		"specificity": TN / (TN + FP),
		"1-PPV": (1-(TP / (TP + FP))),
		"TP": TP,
		"TN": TN,
		"FP": FP,
		"FN": FN,
		# "confusion_matrix": cm
	}

	return eval_metrics

def post_process(df):
	return (
		df
		.with_columns(
			pl.col("cutoff").log().alias("log(q-cut)"),
			pl.col("1-PPV").log().alias("log(1-PPV)"),
		)
		.filter(pl.col("log(1-PPV)").is_finite())
	)


score_label_path = glob.glob("score_labels/*.tsv")

thresholds = np.geomspace(1e-16, 1.0, num=100)

all_metrics = []
for path in score_label_path:
	
	dataset = "_".join(os.path.basename(path).split("_")[:2])
	print(f"Processing dataset: {dataset}")
	
	df = pl.read_csv(path, separator="\t", infer_schema_length=1000)
	df = df.with_columns(dataset = pl.lit(dataset))

	samples = df["sample_name"].unique()
	
	dataset_metrics = []
	for sample in (samples):
		s_df = df.filter(pl.col("sample_name") == sample)
		if (s_df.filter(pl.col("truth")).height == 0) | (s_df.filter(~pl.col("truth")).height == 0):
			continue

		for thresh in thresholds:
			s_pred = adaptive_fdr_pred(s_df, score_col="FOBP", cutoff=thresh)
			eval_metrics = get_eval_metrics(s_pred["truth"], s_pred["pred"])
			eval_metrics = {"dataset": dataset, "sample_name": sample} | {"cutoff": thresh} | eval_metrics

			dataset_metrics.append(eval_metrics)

	pl.DataFrame(dataset_metrics).pipe(post_process).write_csv(f"{dataset}-cutoff_metrics.tsv", separator="\t")
	all_metrics.extend(dataset_metrics)

pl.DataFrame(all_metrics).pipe(post_process).write_csv("combined.tsv", separator = "\t")

