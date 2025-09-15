## This script is used to make ROC and PRC plots from the evaluation results aggregated using combine_results.R
source("../R/plot.R")

## Set output directory
outdir_root = "EGAD00001004066/somatic_vcf/plots"

## Make a vector with paths to all the precrec eval objects
roc_coord_paths <- sort(list.files("EGAD00001004066/somatic_vcf/roc-prc-auc/precrec", pattern = "all-models_roc_coordinates.tsv", recursive = TRUE, full.names = TRUE))
prc_coord_paths <- sort(list.files("EGAD00001004066/somatic_vcf/roc-prc-auc/precrec", pattern = "all-models_prc_coordinates.tsv", recursive = TRUE, full.names = TRUE))

## Plot display size parameters for debugging
# options(repr.plot.width = 6, repr.plot.height = 5)

## Create plots for each of the evaluated samples
message("Creating ROC PRC plot for:")

outdir = glue("{outdir_root}/roc_prc_plots")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(length(roc_coord_paths))) {

	sample_name <- match_return_sample_name(roc_coord_paths[i], prc_coord_paths[i])
	message(glue("\t {i}. {sample_name}"))

	# ## Read in the variant set to count number and annotate the number of SNVs in the plot
	# eval_snv_set <- qread(glue("EGAD00001004066/somatic_vcf/model-scores_truths/{sample_name}/{sample_name}_model-scores_truths.tsv"))
	# snv_count <- nrow(eval_snv_set)

	roc_coord <- qread(roc_coord_paths[i])
	prc_coord <- qread(prc_coord_paths[i])

	roc_plot <- ggplot(roc_coord, aes(x = x, y = y, color = model))  +
		geom_abline(linetype = "dashed", color = "lightgrey") +
		geom_line() +
		labs(
			title = "ROC",
			x = "False Positive Rate (1 - Specificity)",
			y = "True Positive Rate (Sensitivity)",
			color = "Models"
		) +
		coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
		theme_minimal() +
		theme(
			panel.grid.major = element_blank(), # Remove major grid lines
			panel.grid.minor = element_blank(), # Remove minor grid lines
			panel.background = element_blank(), # Optional: Remove panel background
			axis.line = element_line(color = "darkgrey"), # Optional: Add axis lines
			legend.position = "bottom",
			plot.title = element_text(hjust = 0.5) # Center the plot title
		)

	prc_plot <- ggplot(prc_coord, aes(x = x, y = y, color = model))  +
		geom_line() +
		labs(
			title = "PRC",
			x = "Recall",
			y = "Precision",
			color = "Models"
		) +
		coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
		theme_minimal() +
		theme(
			panel.grid.major = element_blank(), # Remove major grid lines
			panel.grid.minor = element_blank(), # Remove minor grid lines
			panel.background = element_blank(), # Optional: Remove panel background
			axis.line = element_line(color = "darkgrey"), # Optional: Add axis lines
			legend.position = "bottom",
			plot.title = element_text(hjust = 0.5) # Center the plot title
		)

	roc_prc_plot <- (roc_plot | prc_plot) +
		plot_annotation(
			title = "EGAD00001004066 Dataset",
			subtitle = sample_name,
			caption = "Only C>T SNVs Evaluated"
		) +
		plot_layout(widths = c(1, 1), guides = "collect") &
		theme(
			plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(t = 10, b = 10)),
			plot.subtitle = element_text(hjust = 0.5),
			plot.caption = element_text(hjust = 0),
			legend.position = "bottom",
			axis.title.x = element_text(size = 10),
			axis.title.y = element_text(size = 10)
		)

	qdraw(roc_prc_plot, glue("{outdir}/{sample_name}_roc_prc_plot.pdf"), width = 6, height = 5)

}
