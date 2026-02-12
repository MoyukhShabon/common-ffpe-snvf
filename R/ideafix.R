#!/usr/bin/env Rscript
library(ideafix)
library(io)
library(argparser)
library(tidyr)
library(dplyr)
library(tictoc)

#' Read a VCF File
#'
#' @param path [string] The path to the VCF file.
#' @param columns [vector] of column names to keep in lowercase (optional).
#' @param split_multiallelic [logical] Splits multiallelic sties into biallelic (optional).
#'
#' @return [data.frame] Table containing VCF columns
read_vcf <- function(path, columns = NULL, split_multiallelic = TRUE) {
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

	if (split_multiallelic) {
		vcf <- separate_rows(vcf, alt, sep = ",")
	}

	return(vcf)
}

#' Filter C>T variants
#' @param d	[data.frame] of variants
#' @param ref_col	[string] containing Reference allele column name
#' @param alt_col	[string] containing Alternate allele column name
#' @returns [data.frame] with only C:G>T:A mutations
ct_filter <- function(d, ref_col = "ref", alt_col = "alt", invert = FALSE) {
	ct_mask = ((d[["ref"]] == "C") & (d[["alt"]] == "T")) | ((d[["ref"]] == "G") & (d[["alt"]] == "A"))
	if (invert){
		d[!ct_mask,]
	} else {
		d[ct_mask,]
	}
}

#' Ideafix only runs on variants with VAF < 0.3 and Filter == PASS. 
#' Essentially, calling them VAF>0.3 non-demaination (real mutations) and Filter != PASS as artifacts.
#' Thus, the filtered variants are restored with a score of 0 i.e highest confidence non-deaminations.
#' @param result [data.frame] This is the result from ideafix's classify_variants function.
#' @param vcf [string] path to the vcf file for the sample being analyzed
#' @param apply_ct_filter (logical) if TRUE (default) retains only C:G>T:A variants
#' @return [data.frame] Table with the VAF filtered variants restored
postprocess <- function(result, vcf, apply_ct_filter=TRUE) {
	tic("Postprocessing")
	# changes ideafix result col names to lowercase (chrom	pos	ref	alt	deam_score	deamination)
	names(result) <- tolower(names(result))

	# Reads all variants from vcf and applies a C>T filter if specified
	all_vars <- read_vcf(vcf, c("chrom", "pos", "ref", "alt", "filter"))
	if (apply_ct_filter) {
    	all_vars <- ct_filter(all_vars)
  	}

	# Merges the result table with all variants in the to restore variants filtered out by the VAF<0.30 filter
	if (nrow(result) != 0){
		result_all_var <- merge(x = all_vars, y = result, by = c("chrom", "pos", "ref", "alt"), all.x = TRUE)
	} else {
		all_vars$deam_score <- NA
		all_vars$deamination <- NA
		result_all_var <- all_vars
	}

	result_all_var$deamination <- as.character(result_all_var$deamination)

	# Score is deamination probability. Therefore hard filtered variants are set to 0 i.e non-deamination
	result_all_var$deam_score <- ifelse(
		!is.na(result_all_var$deam_score),
		result_all_var$deam_score,
		ifelse(result_all_var$filter == "PASS", 0, 1)
	)

	result_all_var$deamination <- ifelse(
		!is.na(result_all_var$deamination),
		result_all_var$deamination,
		ifelse(result_all_var$filter == "PASS", "vaf>0.3", "non-PASS")
	)
	toc(log=TRUE)
	return(result_all_var)
}

# Parse Command-Line Arguments
p <- arg_parser("Perform FFPE SNVF with IdeaFix")
p <- add_argument(p, "--vcf", help = "path to the VCF file")
p <- add_argument(p, "--ref", help = "path to the reference genome")
p <- add_argument(p, "--outdir", help = "Output directory for the results", default = ".")
p <- add_argument(p, "--sample_name", help = "Name of Sample being filtered", default = NA)
p <- add_argument(p, "--run-rf", help = "Name of Sample being filtered", default = FALSE)
argv <- parse_args(p)

tic("Total Execution")
tic("Preamble")
# Assign Variables
vcf_path <- argv$vcf
ref_genome <- argv$ref
sample_name <- ifelse(is.na(argv$sample_name), unlist(strsplit(basename(vcf_path), "\\."))[1], argv$sample_name)

outdir <- file.path(argv$outdir, "ideafix", sample_name)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

toc(log=TRUE)
tic("Feature Extraction")
message(cat(sprintf("\nProcessing %s \n", sample_name)))
message(cat("\tGetting descriptors from VCF \n"))
descriptors <- get_descriptors(vcf_filename = vcf_path, fasta_filename = ref_genome)
## An error is thrown if inf values exist in the descriptors They are removed
descriptors <- filter(descriptors, if_all(where(is.numeric), ~ !is.infinite(.)))
toc(log=TRUE)

qwrite(descriptors, file.path(outdir, sprintf("%s.ideafix.descriptors.tsv", sample_name)))

## Ideafix imposes an AF <= 0.3 and Filter == PASS. If no variants pass these then an empty descriptor data.frame is returned.
if(nrow(descriptors) == 0){
	message("No variant descriptors pass the filters")
	results <- postprocess(data.frame(), vcf_path)
		qwrite(results, file.path(outdir, sprintf("%s.ideafix.tsv", sample_name)))
	if (argv$run_rf){
		qwrite(results, file.path(outdir, sprintf("%s.ideafix-rf.tsv", sample_name)))
	}

	quit(save = "no")
}

tic("Inference - XGBoost")
message(cat("\n\tRunning IdeaFix XGBoost Model \n"))
predictions_xgboost <- classify_variants(variant_descriptors = descriptors, algorithm = "XGBoost")
toc(log=TRUE)
predictions_xgboost <- postprocess(predictions_xgboost, vcf_path)
message(cat(sprintf("\n\tWriting outputs to: %s/ \n", outdir)))
qwrite(predictions_xgboost, file.path(outdir, sprintf("%s.ideafix.tsv", sample_name)))


if (argv$run_rf){
	tic("Inference - XGBoost")
	message(cat("\n\tRunning IdeaFix RF Model \n"))
	predictions_rf <- classify_variants(variant_descriptors = descriptors, algorithm = "RF")
	predictions_rf <- postprocess(predictions_rf, vcf_path)
	qwrite(predictions_rf, file.path(outdir, sprintf("%s.ideafix-rf.tsv", sample_name)))
	toc(log=TRUE)
}

message(cat("\n\tComplete"))

toc(log=TRUE)

# Process Timing Logs ---
raw_log <- tic.log(format = FALSE)
timings_df <- do.call(rbind, lapply(raw_log, function(x) {
	data.frame(
		section = x$msg,
		seconds = x$toc - x$tic
	)
}))
timings_df$minutes <- round(timings_df$seconds / 60, 4)
timings_df$seconds <- round(timings_df$seconds, 4)

qwrite(timings_df, file.path(outdir, sprintf("%s.runtime.tsv", sample_name)))