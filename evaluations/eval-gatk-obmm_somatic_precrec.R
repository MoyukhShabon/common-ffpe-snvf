library(io)
library(precrec)
library(argparser)

# TODO move common functions into module
# TODO refactoring in progress

source("../R/eval.R")

# Setup Directories

## specific for each dataset

dataset_id <- "EGAD00001004066"

# Directory for FFPE SNVF inputs
ffpe_snvf.dir <- sprintf("../ffpe-snvf/%s/somatic_vcf", dataset_id)
# Directory for the Somatic VCFs
vcf.dir <- sprintf("../vcf/%s/somatic_vcf", dataset_id)
# Output directory
main.outdir <- sprintf("%s/somatic_vcf", dataset_id)
score_truth_outdir <- sprintf("%s/somatic_vcf/model-scores_truths/", dataset_id)
eval_outdir <- sprintf("%s/somatic_vcf/roc-prc-auc/precrec", dataset_id)


# Read Annotation Table
lookup_table <- read.delim(sprintf("../annot/%s/sample_annotations.tsv", dataset_id))

# Stratify annotation table based on FFPE and FF Somatic Variants
ffpe_tumoral <- lookup_table[(lookup_table$preservation == "FFPE" & lookup_table$sample_type == "Tumoral"), ]
frozen_tumoral <- lookup_table[(lookup_table$preservation == "Frozen" & lookup_table$sample_type == "Tumoral"), ]


#######################################

# Evaluate gatk_obmm
message("Evaluating gatk-obmm:")
model_name <- "gatk-obmm"
for (index in seq_len(nrow(ffpe_tumoral))){
	
	metadata <- set_up(ffpe_tumoral, index)
	message(sprintf("	%s", metadata$sample_name))
	
	gatk_obmm_processed <- process_sample(
		read_gatk_snv,
		construct_ground_truth,
		preprocess_gatk_obmm,
		evaluate_filter,
		sample_name = metadata$sample_name,
		tissue = metadata$tissue,
		filter_name = model_name,
		snvf_dir = vcf.dir,
		gt_vcf_dir = vcf.dir,
		gt_annot_d = frozen_tumoral
	)

	write_sample_eval(gatk_obmm_processed, main.outdir, metadata$sample_name, model_name)
}

# Overall Evaluation
## The scores annotated with ground truth is combined into a single dataframe
message("	performing Evaluation across all samples")
gatk_obmm_all_score_truth <- do.call(
	rbind,
	lapply(seq_len(nrow(ffpe_tumoral)), function(i) {
		meta <- set_up(ffpe_tumoral, i)
		path <- file.path(score_truth_outdir, meta$sample_name, sprintf("%s_%s-scores_truths.tsv", meta$sample_name, model_name))
		d <- read.delim(path)
		d$sample_name <- meta$sample_name
		d
	})
)


# Evaluate across all samples
gatk_obmm_overall_res <- evaluate_filter(gatk_obmm_all_score_truth, model_name, "all-samples")
write_overall_eval(gatk_obmm_all_score_truth, gatk_obmm_overall_res, score_truth_outdir, eval_outdir, "all-samples", model_name)


# Evaluate across colon samples
message("	performing Evaluation across all colon samples")
gatk_obmm_colon_score_truth <- gatk_obmm_all_score_truth[grepl("Colon", gatk_obmm_all_score_truth$sample_name), ]
gatk_obmm_colon_res <- evaluate_filter(gatk_obmm_colon_score_truth, model_name, "all-colon-samples")
write_overall_eval(gatk_obmm_colon_score_truth, gatk_obmm_colon_res, score_truth_outdir, eval_outdir, "all-colon-samples", model_name)


## Evaluate across liver samples
message("	performing Evaluation across all liver samples")
gatk_obmm_liver_score_truth <- gatk_obmm_all_score_truth[grepl("Liver", gatk_obmm_all_score_truth$sample_name), ]
gatk_obmm_liver_res <- evaluate_filter(gatk_obmm_liver_score_truth, model_name, "all-liver-samples")
write_overall_eval(gatk_obmm_liver_score_truth, gatk_obmm_liver_res, score_truth_outdir, eval_outdir, "all-liver-samples", model_name)

message("Done.")
