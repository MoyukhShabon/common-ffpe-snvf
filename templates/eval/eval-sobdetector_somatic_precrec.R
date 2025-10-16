library(io)
library(precrec)

source("../../common-ffpe-snvf/R/eval.R")

#####################################################################################

# Function for performing per sample evaluation
# @params: ffpe_tumor_annot (data.frame) | annotation table describing FFPE tumor samples
# @params: ff_tumor_annot (data.frame) |  annotation table describing Fresh Frozen tumor samples
# @params: case_id_col (string) | column name with unique identifier for each case/patient. Used for finding match FF for a FFPE sample
# @params: sample_name_col (string) | column with the sample name used in the FFPE
# @params: model_name (string) | name of model being evaluated and used in path
# @params: ffpe_snvf_dir (string) | directory with the results of ffpe snvf filters
# @params: outdir_root (string) | output directory for evaluation results
# @params: vcf_ext (string) | vcf extension (.vcf or .vcf.gz)

process_samples <- function(ffpe_tumor_annot, ff_tumor_annot, case_id_col, variant_caller_col, sample_name_col, model_name, ffpe_snvf_dir, outdir_root, vcf_ext = ".vcf.gz", ct_only=FALSE) {
	# Evaluate Filter per sample
	## Each sample is evaluated first due to the necessity of independently annotating the scores with ground truth

	message("Processing samples:")
	for (index in seq_len(nrow(ffpe_tumor_annot))){

		## Get FFPE sample metadata
		sample_name <- sprintf("%s",ffpe_tumor_annot[index, sample_name_col])
		case_id <- ffpe_tumor_annot[index, case_id_col]
		wtype <- ffpe_tumor_annot[index, variant_caller_col]

		message(sprintf("	%s", sample_name))

		## Getting matched FF metadata by matching patient ID
		matched_ff_metadata <- frozen_tumoral[(frozen_tumoral[[case_id_col]] == case_id) & (frozen_tumoral[[variant_caller_col]] == wtype), ]
		matched_ff_sample_name <- matched_ff_metadata[[sample_name_col]]
		matched_ff_paths <- file.path(vcf_dir, matched_ff_sample_name, sprintf("%s%s", matched_ff_sample_name, vcf_ext))

		sobdetector_path <- file.path(ffpe_snvf_dir, tolower(wtype), model_name, sample_name, sprintf("%s.%s.snv", sample_name, model_name))
		
		if(!file.exists(sobdetector_path)){
			message(sprintf("		Result table for %s does not exist at: %s", model_name, sobdetector_path))
			message(sprintf("		Skippping %s ...", sample_name))
			next
		}

		d <- read.delim(sobdetector_path)
		truth <- snv_union(matched_ff_paths)
		d <- preprocess_sobdetector(d, truth, ct_only)

		## Check if truth labels are not exclusively TRUE or FALSE in d (variant_score_truth table
		## Cases like these are skipped as evaluation is not supported by precrec
		if(nrow(d[d$truth, ]) == 0 | nrow(d[!d$truth, ]) == 0){
			message(sprintf("		truth labels consists of only one class for %s", sample_name))
			next
		}

		# Evaluate the filter's performance
		sobdetector_res <- evaluate_filter(d, model_name, sample_name)

		# write results
		write_sample_eval(d, sobdetector_res, outdir_root, sample_name, model_name)
	}
}


# Function for combining the ground truth annotated SNV score table
# @params: ffpe_tumoral_annot (data.frame) | annotation table describing FFPE tumor samples
# @params: sample_name_col (string) | column with the sample name used in the FFPE
# @params: model_name (string) | name of model being evaluated. The name used in path
# @params: score_truth_outdir (srring) | directory where the ground truth annotated scores and truth were saved
combine_snv_score_truth <- function(ffpe_tumoral_annot, sample_name_col, variant_caller_col, model_name, score_truth_outdir) {
	message("Combining all the per sample ground truth SNV score tables into one")
	do.call(
		rbind,
		lapply(seq_len(nrow(ffpe_tumoral_annot)), function(i) {
			sample_name <- ffpe_tumoral_annot[i, sample_name_col]
			wtype <- ffpe_tumoral_annot[i, variant_caller_col]
			message(sprintf("	%s", sample_name))
			path <- file.path(score_truth_outdir, sample_name, sprintf("%s_%s-scores_truths.tsv", sample_name, model_name))
			if (!file.exists(path)){
				message(sprintf("		Warning: %s was not was not found. SKIPPING", path))
			} else {
				d <- read.delim(path)
				d$sample_name <- sample_name
				d$workflow_type <- wtype
				d
			}
		})
	)
}

#################################################################################################

model_name <- "sobdetector"
message(sprintf("Evaluating %s: ", model_name))

#################################  TCGA  ########################################

# Setup Directories and lookup table for TCGA

dataset_id <- "TCGA"
message(sprintf("Dataset: %s", dataset_id))

# Directory for the Somatic VCFs
vcf_dir <- "../vcf-pass-filter"
ffpe_snvf_dir <- "../vcf-pass-filter-ffpe-snvf"
# Output directory
outdir_root <- "."
score_truth_outdir <- "model-scores_truths" 
eval_outdir <- "roc-prc-auc/precrec"


# Read Annotation Table
lookup_table <- read.delim(sprintf("../annot2/vcf_annotation.tsv"))
lookup_table$is_ffpe <- as.logical(lookup_table$is_ffpe)

# Stratify annotation table based on FFPE and FF Somatic Variants
ffpe_tumoral <- lookup_table[(lookup_table$is_ffpe), ]
frozen_tumoral <- lookup_table[!(lookup_table$is_ffpe), ]


# ------------------------------------------------------

# Perform per sample evaluation
process_samples(ffpe_tumor_annot = ffpe_tumoral, ff_tumor_annot = frozen_tumoral, case_id_col = "case_id", variant_caller_col = "workflow_type", sample_name_col = "vcf_fid", model_name, ffpe_snvf_dir, outdir_root)

# Combine ground truth annotated score tables into one
sobdetector_all_score_truth <- combine_snv_score_truth(ffpe_tumoral_annot = ffpe_tumoral, sample_name_col = "vcf_fid", variant_caller_col = "workflow_type", model_name, score_truth_outdir)


# Stratify each variant caller
# mutect2
message("Performing Evaluation across mutect2 samples")
sobdetector_mutect2_score_truth <- sobdetector_all_score_truth[grepl("mutect2", tolower(sobdetector_all_score_truth$workflow_type)), ]
sobdetector_mutect2_overall_res <- evaluate_filter(sobdetector_mutect2_score_truth, model_name, "all-mutect2-samples")
write_overall_eval(sobdetector_mutect2_score_truth, sobdetector_mutect2_overall_res, score_truth_outdir, eval_outdir, "all-mutect2-samples", model_name)

# varscan2
message("Performing Evaluation across varscan2 samples")
sobdetector_varscan2_score_truth <- sobdetector_all_score_truth[grepl("varscan2", tolower(sobdetector_all_score_truth$workflow_type)), ]
sobdetector_varscan2_overall_res <- evaluate_filter(sobdetector_varscan2_score_truth, model_name, "all-varscan2-samples")
write_overall_eval(sobdetector_varscan2_score_truth, sobdetector_varscan2_overall_res, score_truth_outdir, eval_outdir, "all-varscan2-samples", model_name)

# muse
message("Performing Evaluation across muse samples")
sobdetector_muse_score_truth <- sobdetector_all_score_truth[grepl("muse", tolower(sobdetector_all_score_truth$workflow_type)), ]
sobdetector_muse_overall_res <- evaluate_filter(sobdetector_muse_score_truth, model_name, "all-muse-samples")
write_overall_eval(sobdetector_muse_score_truth, sobdetector_muse_overall_res, score_truth_outdir, eval_outdir, "all-muse-samples", model_name)

# somaticsniper 
message("Performing Evaluation across somaticsniper samples")
sobdetector_somaticsniper_score_truth <- sobdetector_all_score_truth[grepl("somaticsniper", tolower(sobdetector_all_score_truth$workflow_type)), ]
sobdetector_somaticsniper_overall_res <- evaluate_filter(sobdetector_somaticsniper_score_truth, model_name, "all-somaticsniper-samples")
write_overall_eval(sobdetector_somaticsniper_score_truth, sobdetector_somaticsniper_overall_res, score_truth_outdir, eval_outdir, "all-somaticsniper-samples", model_name)


##############################################################################################

