#!/usr/bin/env Rscript

# This script contains functions for evaluation
## This is a work in progress. The functions still needs to be improved for generalizability across multiple projects and datasets

library(io)
library(precrec)

#### A list of Chromosomes 1-22, X and Y
standard_chromosomes <- paste0("chr", c(1:22, "X", "Y"))

## Read a VCF File
## @param path string | VCF path
## @param columns vector | (optional) column names to keep in lowercase
## @return  dataframe | containing vcf columns
read_vcf <- function(path, columns = NULL) {
	all_lines <- readLines(path)
	filtered_lines <- grep("^##", all_lines, value = TRUE, invert = TRUE)

	vcf <- read.delim(
		text = filtered_lines,
		sep = "\t",
		check.names = FALSE
	)

	names(vcf)[1] <- sub("^#", "", names(vcf)[1])
	names(vcf) <- tolower(names(vcf))

	if (!is.null(columns)) {
		vcf <- vcf[, columns]
	}

	vcf
}

# Filter C>T variants
# @params d	data.frame of variants
# @params ref_col	string containing Reference allele column name
# @params alt_col	string containing Alternate allele column name
# @returns dataframe with only C:G>T:A mutations
ct_filter <- function(d, ref_col = "ref", alt_col = "alt") {
	ct_mask = ((d[["ref"]] == "C") & (d[["alt"]] == "T")) | ((d[["ref"]] == "G") & (d[["alt"]] == "A"))
	d[ct_mask,]
}


# Create variant ID based on chrom, pos, ref, and alt
# @param d  data.frame of variants
# @return data.frame of variant with ID annotation
add_id <- function(d) {
	d$snv <- with(d, paste(chrom, pos, ref, alt, sep="_"))
	d
}


# Annotate called variant with whether each variant is in the ground truth
# variant call set
# @param  d  data.frame of called variants to be annotated
# @param  truth  data.frame of ground truth variants
# @returns data.frame with annotated ground truth label
annotate_truth <- function(d, truth) {
	d$truth <- d$snv %in% truth$snv;
	d
}

#### Create a set of SNVs across multiple samples
# @params paths list of strings with vcf paths
# @return data.frame of variant union set
snv_union <- function(paths) {

	combined <- data.frame()

	for (path in paths){
		# Read the VCF data
		sample <- read_vcf(path, columns = c("chrom", "pos", "ref", "alt"))
		# Select only the required columns using base R
		combined <- rbind(combined, sample)
	}

	# only keep one instance of each variant
	union <- unique(combined)

	return(union)
}


#### Function to classify variants based on a set of Q values and a False Positive cutoff threshold
## @params q vector | vector of FDR adjusted p-values i.e. q-values
## @returns a boolean vector
adaptive_fdr_cut <- function(q, fp.cut) {
	n <- length(q);
	# adaptive threshold chosen to expect less than one false positive
	idx <- order(q);
	under <- which((1:n) * q[idx] < fp.cut);
	if (length(under) > 0) {
		top <- under[length(under)];
		pred <- logical(n);
		##### select the top significant results
		pred[idx[1:top]] <- TRUE;
		pred
	} else {
		##### none passes
		rep(FALSE, n)
	}
}


#### Function to label the models prediction based on the adaptive_fdr_cut function
# @params df (data.frame) | with model's scores
# @params score (string) | column name of score column
# @params fp-cut (numeric) | the false positive cut off value. Default = 0.5
# @return data.frame of model's score with annotated prediction based on fp-cut
fdr_cut_pred <- function(df, score, fp.cut=0.5) {
	
	##### Split C>T and non C>T mutations into two dataframes using base R
	df.ct <- df[complete.cases(df[[score]]), ]
	df.nct <- df[!complete.cases(df[[score]]), ]

	##### For non C>T, mutate columns using direct assignment
	df.nct$score <- NA
	df.nct$q <- NA
	df.nct$pred <- TRUE

	##### Predict artifacts based on q-values and adaptive_fdr_cut function
	# The original dataframe is modified by direct assignment
	df.ct$score <- ifelse(df.ct[[score]] == 0, .Machine$double.eps, df.ct[[score]])
	df.ct$q <- p.adjust(df.ct$score, "BH")
	## TRUE being real muatation and FALSE being artifacts
	df.ct$pred <- ifelse(adaptive_fdr_cut(df.ct$q, fp.cut), TRUE, FALSE)

	##### Combine the C>T and non C>T dataframes back into one single dataframe
	df <- rbind(df.ct, df.nct)
	# Sort
	df <- df[order(df$chrom, df$pos), ]
	df
}


#### Calculate True/False Positives/Negatives from prediction and ground truth
# @param data.frame with models prediction and ground truth columns
# @return list with numeric FP, TP, FN, and TP
calc.confusion.matrix <- function(df, pred_col = "pred", truth_col = "truth") {
	TP <- nrow(df[df[[pred_col]] & df[[truth_col]], ])
	FP <- nrow(df[df[[pred_col]] & !df[[truth_col]], ])
	FN <- nrow(df[!df[[pred_col]] & df[[truth_col]], ])
	TN <- nrow(df[!df[[pred_col]] & !df[[truth_col]], ])
	
	list(
		TP = TP,
		FP = FP,
		FN = FN,
		TN = TN
	)
}


#### Calculate evaluation metrics based on confusion matrix
# @param data.frame with models prediction and ground truth columns
# @return list with sensitivity, specificity, precision, and recall
calc_eval_metrics <- function(df, pred_col = "pred", truth_col = "truth") {

	confusion_matrix <- calc.confusion.matrix(df, pred_col, truth_col)

	precision <- confusion_matrix$TP / (confusion_matrix$TP + confusion_matrix$FP)
	recall <- confusion_matrix$TP / (confusion_matrix$TP + confusion_matrix$FN)
	sensitivity <- recall  # sensitivity is same as recall
	specificity <- confusion_matrix$TN / (confusion_matrix$TN + confusion_matrix$FP)
	
	list(
		precision = precision,
		recall = recall,
		sensitivity = sensitivity,
		specificity = specificity
	)
}

## FIX ME remove after refactoring is complete
## Function to create multimodel AUROC and AUPRC table based on a dataframe containing ROC and PRC
## Column names for ROC df must be fpr and tpr
## Column names for ROC df must be precision and recall
get_multimodel_auroc_auprc <- function(multimodel_roc_df,  multimodel_prc_df, sample_name, model_names = NULL, model_col = "model") {

	if(is.null(model_names)) {
		model_names <- unique(multimodel_roc_df$model)
	}

	auc_list <- list()

	for (i in seq_along(model_names)) {
		model_name <- model_names[i]

		# Filter using base R subsetting instead of dplyr::filter
		model_roc <- multimodel_roc_df[multimodel_roc_df[[model_col]] == model_name, ]
		model_prc <- multimodel_prc_df[multimodel_prc_df[[model_col]] == model_name, ]

		model_auroc <- calculate.auc(model_roc, x_col = "fpr", y_col = "tpr")
		model_auprc <- calculate.auc(model_prc, x_col = "recall", y_col = "precision")

		auc_list[[i]] <- data.frame(
			sample_id = sample_name,
			model = model_name,
			auroc = model_auroc,
			auprc = model_auprc
		)

	}

	# Combine list of data frames using base R
	aucs <- do.call(rbind, auc_list)

	return(aucs)
}

############################

# @params sample_name	unique identifier for the ffpe sample
# @params filter_name	name of the ffpe filter used in the path structure
# @params ffpe_snvf.dir	root directory containing the ffpe filter outputs
read_snv <- function(sample_name, filter_name, ffpe_snvf.dir) {
	path <- file.path(ffpe_snvf.dir, filter_name, sample_name, sprintf("%s.%s.snv", sample_name, filter_name))
	read.delim(path)
}


# @params sample_name	unique identifier for the ffpe sample
# @params vcf.dir	root directory for VCFs
# @params filter_name	Unused. Present for compatibility with process_sample()
# @params read_vcf_f	function to read VCFs. By default set to "read_vcf" which is defined in eval.R
read_gatk_snv <- function(sample_name, filter_name, vcf.dir, read_vcf_f = read_vcf){
	d <- read_vcf(file.path(vcf.dir, sample_name, sprintf("%s.vcf.gz", sample_name)), columns = c("chrom", "pos", "ref", "alt", "filter"))
	d
}

# Construct the Ground Truth SNVs
## To do this we make a union set for all the FF samples matched to this FFPE sample
## This dataset has matched FFPE and FF from the same patient.
## All samples within the same tissue type are matched
## Hence, we select the frozen variants from the same tissue type as the FFPE sample
## @params annot_d		data.frame containing sample annotation for frozen samples
## @params tissue_type		string describing the type of tissue
## @return data.frame of ground truth variants
construct_ground_truth <- function(annot_d, tissue, vcf.dir){

	annot_d <- annot_d[annot_d$tissue_type == tissue, ]

	## Obtain the sample names for the frozen samples within the same tissue
	samples <- annot_d$sample_name
	paths <- file.path(vcf.dir, samples, sprintf("%s.vcf.gz", samples))
	truths <- snv_union(paths)
	truths <- add_id(truths)
	truths
}

# @param d  data.frame of variant annotation by mobsnvf
# @param truths  data.frame of ground-truth variants
# @return data.frame of variants with id and ground truth annotation
preprocess_mobsnvf <- function(d, truths) {
	# mobsnvf sets FOBP to NA for variants that are not C>T
	d <- d[!is.na(d$FOBP), ];
	# lower score signifies real mutation:  
	# hence, we flip scores to make higher score signify a real mutation
	d$score <- -d$FOBP;
	d <- add_id(d);
	d <- annotate_truth(d, truths)
	d
}

# @param d  data.frame of variant annotation by vafsnvf
# @param truths  data.frame of ground-truth variants
# @return data.frame of variants with id and ground truth annotation
preprocess_vafsnvf <- function(d, truths) {
	d <- d[!is.na(d$VAFF), ]
	# vafsnvf sets VAFF to NA for variants that are not C>T
	d$score <- d$VAFF;
	d <- add_id(d);
	d <- annotate_truth(d, truths)
	d
}

# @param d  data.frame of variant annotation by sobdetector
# @param truths  data.frame of ground-truth variants
# @return data.frame of variants with id and ground truth annotation
preprocess_sobdetector <- function(d, truths) {
	# SOBDetector output column explanations:
	# 		artiStatus: Binary classification made by SOBDetector. Values are "snv" or "artifact"
	# 		SOB: This is the strand oreintation bias score column which ranges from 0 and 1. Exception values: "." or NaN. 
	# variants ignored by SOBdetector have score of "."
	d <- d[!(d$SOB == "."), ]
	# now, it's safe to convert to numeric
	d$SOB <- as.numeric(d$SOB)
	# SOBdetector score = 0 indicates that it is not artifact
	# variants classified by SOBdetect as "snv" have score of NaN
	d$SOB <- ifelse(is.nan(d$SOB), 0, d$SOB)
	### Precaution to remove the any NA Scores if present
	d <- d[!is.na(d$SOB), ]
	# in this set the TRUE label indicates real mutation and FALSE indicates artifacts
	# Hence, scores needs to be adjusted so that higher score represents a real mutation.
	d$score <- -d$SOB;
	# Keep only C>T variants
	d <- ct_filter(d)
	d <- add_id(d);
	d <- annotate_truth(d, truths)
	d
}

# @param d  data.frame of variant annotation by GATK Orientation Bias Mixture Model
# @param truths  data.frame of ground-truth variants
# @return data.frame of variants with id and ground truth annotation
preprocess_gatk_obmm <- function(d, truths){
	### GATK Orientation Bias mixture model makes binary classification. This is casted into scores 0 and 1
	d$score <- ifelse(grepl("orientation", d$filter), 0, 1)
	d <- d[!is.na(d$score), ]
	# Keep only C>T variants
	d <- ct_filter(d)
	d <- add_id(d)
	d <- annotate_truth(d, truths)
	d
}


## Function to evaluate an FFPE filter's performance
# @param d       data.frame of variant annotation
# @param model_names	string containing name to model to be evaluated
# @param sample_name 	string containing sample name being evaluated
# @return list with components:
# 	- eval	precrec eval object
#	- auc	data.frame containing AUROC and AUPRC of each model
#	- roc	data.frame containing ROC curve coordinated
#	- prc	data.frame containing PRC curve coordinated
evaluate_filter <- function(d, model_name, sample_name) {;
	eval_obj <- with(d, evalmod(scores = score, labels = truth, modnames = model_name));

	# Format AUC table
	eval_auc <- auc(eval_obj)
	eval_auc$dsids <- NULL
	eval_auc <- reshape(eval_auc, idvar = "modnames", timevar = "curvetypes", direction = "wide")
	colnames(eval_auc) <- c("model", "auroc", "auprc")
	eval_auc$sample_name <- sample_name
	eval_auc <- eval_auc[, c("sample_name", "model", "auroc", "auprc")]

	# Collect ROC PRC coordinates
	## Columns -> (x, y, modname, type):
	## x and y are coordinates for making ROC and PRC plots,
	## plot type values are ROC or PRC
	roc_prc <- as.data.frame(eval_obj)
	roc_prc$sample_name <- sample_name
	colnames(roc_prc)[colnames(roc_prc) == "modname"] <- "model"
	roc_prc <- roc_prc[, c("sample_name", "model", "type", "x", "y")]
	# Stratify ROC PRC
	roc <- roc_prc[roc_prc$type == "ROC", ]
	prc <- roc_prc[roc_prc$type == "PRC", ]

	list(
		eval = eval_obj,
		auc = eval_auc,
		roc = roc,
		prc = prc
	)
}


# Remark: evalmod will be replaced with another function in the future!

# Function to obtain necessary metadata
# @params annot_d		data.frame containing sample annotation
# @params index		the index of the sample to process
# @return list with "sample_name" and "tissue" type to be evaluated as string
set_up <- function(annot_d, index) {
	list(
		sample_name = annot_d[index, "sample_name"],
		tissue = annot_d[index, "tissue_type"]
	)
}


# Higher-order wrapper function to process a single FFPE sample with a given filter
# All functions are passed in as arguments to make this a reusable higher-order function.
#
# @param read_f        function(sample_name, filter_name, snvf_dir) -> data.frame
#                      Reads the filter-specific variant output for a sample
#
# @param preprocess_f  function(d, truths) -> data.frame
#                      constructs the ground truth variants for the FFPE sample
#
# @param truth_f       function(gt_annot_d, tissue, gt_vcf_dir) -> data.frame
#                      Constructs the ground-truth variant set for the FFPE sample.
#
# @param evaluate_f    function(d, filter_name) -> list
#                      evaluates the preprocessed variants (e.g. computes ROC/PRC/AUC)
#
# @param sample_name   character. Unique identifier for the FFPE sample (used by read_f and for naming)
# @param tissue        character. Tissue type string used by truth_f to select matching frozen samples
# @param filter_name   character. Name of the filter (used to build file paths or as model name in evaluate_f)
# @param gt_annot_d    data.frame. Annotation table used by truth_f to select matched samples (e.g. frozen_tumoral)
# @param snvf_dir      character. Root directory containing filter outputs (passed to read_f)
# @param gt_vcf_dir    character. Root directory containing ground-truth VCFs (passed to truth_f)
#
# @return list with components:
#  - d   : data.frame. The preprocessed variant table annotated with "score", "id", and "truth" columns.
#  - res : list. Evaluation result object returned by evaluate_f (e.g. contains eval, auc, roc, prc).
#  - gt  : data.frame. Ground-truth variants constructed by truth_f for this sample.
process_sample <- function(read_f, truth_f, preprocess_f, evaluate_f, sample_name, tissue, filter_name, gt_annot_d, snvf_dir, gt_vcf_dir) {
	d <- read_f(sample_name, filter_name, snvf_dir);
	gt <- truth_f(gt_annot_d, tissue, gt_vcf_dir)
	d <- preprocess_f(d, gt);
	res <- evaluate_f(d, filter_name, sample_name);
	
	list(
		d = d,
		res = res
	)
}

# Function to write results for a processed sample
# @params processed_obj 	R-object, processed object returned by function: process_sample()
# @params outdir_root		string, root output directory for evaluation
# @params sample_name		string, sample name for the processed sample
# @params model_name		name of ffpe snv filter being evaluated
write_sample_eval <- function(processed_obj, outdir_root, sample_name, model_name){
	
	# Saving the variant set with scores and ground truth labels for each sample
	out_dir <- file.path(score_truth_outdir, sample_name)
	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
	qwrite(processed_obj$d, file.path(out_dir, sprintf("%s_%s-scores_truths.tsv", sample_name, model_name)))
	
	# Create output directory for roc prc and auc evaluation
	out_dir <- file.path(eval_outdir, sample_name)
	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

	qwrite(processed_obj$res$eval, file.path(out_dir, sprintf("%s_%s_precrec_eval.rds", sample_name, model_name)))
	qwrite(processed_obj$res$auc, file.path(out_dir, sprintf("%s_%s_auc_table.tsv", sample_name, model_name)))
	qwrite(processed_obj$res$roc, file.path(out_dir, sprintf("%s_%s_roc_coordinates.tsv", sample_name, model_name)))
	qwrite(processed_obj$res$prc, file.path(out_dir, sprintf("%s_%s_prc_coordinates.tsv", sample_name, model_name)))
}

# Function to write results for overall evaluation
# @params score_truth_d		data.frame, model scores with ground truth annotation
# @params processed_obj 	R-object, processed object returned by function: process_sample()
# @params outdir_root		string, root output directory for evaluation
# @params sample_name		string, sample name for the processed sample
# @params model_name		name of ffpe snv filter being evaluated
write_overall_eval <- function(score_truth_d, result_obj, score_truth_outdir, eval_outdir, sample_name, model_name){
	
	dir.create(score_truth_outdir, recursive = TRUE, showWarnings = FALSE)
	qwrite(score_truth_d, file.path(score_truth_outdir, sprintf("%s_%s-scores_truths.tsv", sample_name, model_name)))
	
	dir.create(eval_outdir, recursive = TRUE, showWarnings = FALSE)
	qwrite(result_obj$eval, file.path(eval_outdir, sprintf("%s_%s_precrec_eval.rds", sample_name, model_name)))
	qwrite(result_obj$auc, file.path(eval_outdir, sprintf("%s_%s_auc_table.tsv", sample_name, model_name)))
	qwrite(result_obj$roc, file.path(eval_outdir, sprintf("%s_%s_roc_coordinates.tsv", sample_name, model_name)))
	qwrite(result_obj$prc, file.path(eval_outdir, sprintf("%s_%s_prc_coordinates.tsv", sample_name, model_name)))
}
