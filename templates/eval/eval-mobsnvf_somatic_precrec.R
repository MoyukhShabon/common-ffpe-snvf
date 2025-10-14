library(io)
library(precrec)

source("../common-ffpe-snvf/R/eval.R")

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

process_samples <- function(ffpe_tumor_annot, ff_tumor_annot, case_id_col, sample_name_col, model_name, ffpe_snvf_dir, outdir_root, vcf_ext = ".vcf") {
	# Evaluate Filter per sample
	## Each sample is evaluated first due to the necessity of independently annotating the scores with ground truth

	message("Processing samples:")
	for (index in seq_len(nrow(ffpe_tumor_annot))){

		## Get FFPE sample metadata
		sample_name <- sprintf("%s",ffpe_tumor_annot[index, sample_name_col])
		case_id <- ffpe_tumor_annot[index, case_id_col]

		message(sprintf("	%s", sample_name))

		## Getting matched FF metadata by matching patient ID
		matched_ff_metadata <- frozen_tumoral[(frozen_tumoral[[case_id_col]] == case_id), ]
		matched_ff_sample_name <- matched_ff_metadata[[sample_name_col]]
		matched_ff_paths <- file.path(vcf_dir, matched_ff_sample_name, sprintf("%s%s", matched_ff_sample_name, vcf_ext))

		mobsnvf_path <- file.path(ffpe_snvf_dir, model_name, sample_name, sprintf("%s.%s.snv", sample_name, model_name))
		
		if(!file.exists(mobsnvf_path)){
			message(sprintf("		Result table for %s does not exist at: %s", model_name, mobsnvf_path))
			message(sprintf("		Skippping %s ...", sample_name))
			next
		}

		d <- read.delim(mobsnvf_path)
		truth <- snv_union(matched_ff_paths)
		d <- preprocess_mobsnvf(d, truth)

		## Check if truth labels are not exclusively TRUE or FALSE in d (variant_score_truth table
		## Cases like these are skipped as evaluation is not supported by precrec
		if(nrow(d[d$truth, ]) == 0 | nrow(d[!d$truth, ]) == 0){
			message(sprintf("		no true labels exist for %s", sample_name))
			next
		}

		# Evaluate the filter's performance
		mobsnvf_res <- evaluate_filter(d, model_name, sample_name)

		# write results
		write_sample_eval(d, mobsnvf_res, outdir_root, sample_name, model_name)
	}
}


# Function for combining the ground truth annotated SNV score table
# @params: ffpe_tumoral_annot (data.frame) | annotation table describing FFPE tumor samples
# @params: sample_name_col (string) | column with the sample name used in the FFPE
# @params: model_name (string) | name of model being evaluated. The name used in path
# @params: score_truth_outdir (srring) | directory where the ground truth annotated scores and truth were saved
combine_snv_score_truth <- function(ffpe_tumoral_annot, sample_name_col, model_name, score_truth_outdir) {
	message("Combining all the per sample ground truth SNV score tables into one")
	do.call(
		rbind,
		lapply(seq_len(nrow(ffpe_tumoral_annot)), function(i) {
			sample_name <- sprintf("%s",ffpe_tumoral_annot[i, sample_name_col])
			message(sprintf("	%s", sample_name))
			path <- file.path(score_truth_outdir, sample_name, sprintf("%s_%s-scores_truths.tsv", sample_name, model_name))
			if (!file.exists(path)){
				message(sprintf("		Warning: %s was not was not found. SKIPPING", path))
			} else {
				d <- read.delim(path)
				d$sample_name <- sample_name
				d
			}
		})
	)
}

#################################################################################################

model_name <- "mobsnvf"
message(sprintf("Evaluating %s: ", model_name))

#################################  ENA PRJEB44073  ##############################################

# Setup Directories and lookup table for ENA PRJEB44073

dataset_id <- "PRJEB44073"
message(sprintf("Dataset: %s", dataset_id))
# Directory for FFPE SNVF inputs
ffpe_snvf_dir <- sprintf("../ffpe-snvf/%s/vcf_filtered_pass_orientation", dataset_id)
# Directory for the Somatic VCFs
vcf_dir <- sprintf("../vcf/%s/vcf_filtered_pass_orientation", dataset_id)
# Output directory
outdir_root <- sprintf("%s/vcf_filtered_pass_orientation", dataset_id)
score_truth_outdir <- sprintf("%s/vcf_filtered_pass_orientation/model-scores_truths", dataset_id)
eval_outdir <- sprintf("%s/vcf_filtered_pass_orientation/roc-prc-auc/precrec", dataset_id)


# Read Annotation Table
lookup_table <- read.delim(sprintf("../annot/%s/sample-info_stage2.tsv", dataset_id))

# Stratify annotation table based on FFPE and FF Somatic Variants
ffpe_tumoral <- lookup_table[(lookup_table$preservation == "FFPE"), ]
frozen_tumoral <- lookup_table[(lookup_table$preservation == "Frozen"), ]


# ------------------------------------------------------

# Perform per sample evaluation
process_samples(ffpe_tumor_annot = ffpe_tumoral, ff_tumor_annot = frozen_tumoral, case_id_col = "case_id", sample_name_col = "sample_alias", model_name, ffpe_snvf_dir, outdir_root)

# Combine ground truth annotated score tables into one
mobsnvf_all_score_truth <- combine_snv_score_truth(ffpe_tumoral_annot = ffpe_tumoral, sample_name_col = "sample_alias", model_name, score_truth_outdir)

# Evaluate across all samples
message("Performing Evaluation across all samples")
mobsnvf_overall_res <- evaluate_filter(mobsnvf_all_score_truth, model_name, "all-samples")
write_overall_eval(mobsnvf_all_score_truth, mobsnvf_overall_res, score_truth_outdir, eval_outdir, "all-samples", model_name)

##############################################################################################

