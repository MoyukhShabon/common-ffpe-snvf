#!/usr/bin/env Rscript

# This script contains functions for making ROC and PRC plots
## This is a work in progress. The functions still needs to be improved for generalizability across multiple projects and datasets


### Load necessary libraries

library(io)
library(ggplot2)
library(dplyr)
library(glue)
library(patchwork)
library(grid)
library(hrbrthemes)
library(viridis)

### Plotting

# Function to check if the roc and prc coordinate table is from the same sample and returns sample name if they match
# @param roc_path (string) | path to roc coordinates table
# @param prc_path (string) | path to prc coordinates table
# @param roc_suffix (string) | default: "_all-models_roc_coordinates.tsv" | roc coordinates table suffix after sample name
# @param prc_suffix (string) | default: "_all-models_roc_coordinates.tsv" | prc coordinates table suffix after sample name
# @return boolean
match_return_sample_name <- function(roc_path, prc_path, roc_suffix = "_all-models_roc_coordinates.tsv", prc_suffix = "_all-models_roc_coordinates.tsv"){

	sample_name_roc <- gsub(roc_suffix, "", basename(roc_path))
	sample_name_prc <- gsub(prc_suffix, "", basename(roc_path))
	
	if (sample_name_roc != sample_name_prc) {
		stop(paste("Mismatch found at index", i, ": ROC sample '", 
			sample_name_roc, "' does not match PRC sample '", 
			sample_name_prc, "'. Halting execution.", sep=""))
	}

	sample_name_roc
}

## Function to make roc prc plot based on roc and prc coordinate
# @param roc_coord (data.frame) | dataframe of roc coordinates
# @param prc_coord (data.frame) | dataframe of prc coordinates
# @param title (string) | Plot title
# @param subtitle (string) | Plot subtitle
# @param x_col (name) | The name of the column to be used for the x-axis. Default: x
# @param y_col (name) | The name of the column to be used for the y-axis. Default: y
# @param model_col (name) | The name of the column for color grouping. Default: model
# @param caption (string) | Plot caption | default: "Only C>T SNVs Evaluated"
# @param legend_rows (integer) | The number of rows to wrap the legend into. Default: NULL (single row).
# @param individual_plots (boolean) | creates separate ROC and PRC plots alongside combined ROC PRC plot. Returns list of plot objects (roc_prc, roc, prc)
# @return ggplot2 object | containing both roc and prc plot
make_roc_prc_plot <- function(
	roc_coord,
	prc_coord,
	title,
	subtitle,
	caption = "Only C>T SNVs Evaluated",
	x_col = x,
	y_col = y,
	model_col = model,
	text_scale = 1,
	line_width = 0.5,
	legend_scale = 1,
	legend_rows = NULL,
	individual_plots = FALSE
	) {

	# ROC Plot
	roc_plot <- ggplot(roc_coord, aes(x = {{ x_col }}, y = {{ y_col }}, color = {{ model_col }})) +
		geom_abline(linetype = "dashed", color = "lightgrey") +
		geom_line(linewidth = line_width, alpha = 0.8) +
		coord_fixed() +
		labs(
			title = "ROC",
			x = "1 - Specificity",
			y = "Sensitivity",
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
			legend.title = element_blank(),
			legend.key.width = unit(line_width*legend_scale, "cm"),
			legend.text = element_text(size = 10*text_scale*legend_scale),
			axis.title.x = element_text(size = 10*text_scale),
			axis.title.y = element_text(size = 10*text_scale),
			axis.text = element_text(size = 8*text_scale),
			plot.title = element_text(size = 12*text_scale, face = "plain", hjust = 0.5) # Resize, Center the plot title
		)
		#  +
		# scale_color_hue(h = c(60, 420))

	# PRC Plot
	prc_plot <- ggplot(prc_coord, aes(x = {{ x_col }}, y = {{ y_col }}, color = {{ model_col }})) +
		geom_line(linewidth = line_width, alpha = 0.8) +
		coord_fixed() +
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
			legend.title = element_blank(),
			legend.key.width = unit(line_width*legend_scale, "cm"),
			legend.text = element_text(size = 10*text_scale*legend_scale),
			axis.title.x = element_text(size = 10*text_scale),
			axis.title.y = element_text(size = 10*text_scale),
			axis.text = element_text(size = 8*text_scale),
			plot.title = element_text(size = 12*text_scale, face = "plain", hjust = 0.5) # Resize, Center the plot title
		)
		#  +
		# scale_color_hue(h = c(60, 420))

	if (!is.null(legend_rows)) {
		legend_guide <- guides(color = guide_legend(nrow = legend_rows))
		roc_plot <- roc_plot + legend_guide
		prc_plot <- prc_plot + legend_guide
	}

	# Combining the plots using the 'patchwork' package
	roc_prc_plot <- (roc_plot + prc_plot) +
		plot_annotation(
			title = title,
			subtitle = subtitle,
			caption = caption
		) +
		plot_layout(widths = c(1, 1), guides = "collect") &
		theme(
			plot.title = element_text(hjust = 0.5, face = "bold", size = 12*text_scale, margin = margin(t = 10, b = 10)),
			plot.subtitle = element_text(hjust = 0.5),
			plot.caption = element_text(hjust = 0),
			legend.position = "bottom",
			legend.title = element_blank(),
			legend.key.width = unit(line_width*legend_scale, "cm"),
			legend.text = element_text(size = 10*text_scale*legend_scale),
			axis.title.x = element_text(size = 10*text_scale),
			axis.title.y = element_text(size = 10*text_scale),
			axis.text = element_text(size = 8*text_scale)
		)


	if(!individual_plots){
		return(roc_prc_plot)
	} else {
		return(
			list(
				roc_prc = roc_prc_plot,
				roc = roc_plot,
				prc = prc_plot
			)
		)
	}
}

#### Creates a text panel containing all the AUC metrics for each model
make_plot_auc_text <- function(multi.model.eval.object, model.names = c("mobsnvf", "vafsnvf", "sobdetector", "microsec")) {
	
	##### Get AUCs to include in plot
	all.model.aucs <- auc(multi.model.eval.object)
	
	##### Extract AUROC and AUPRC for each model in model.names
	auroc <- sapply(model.names, function(m) {
		all.model.aucs |> filter(modnames == m & curvetypes == "ROC") |> pull(aucs)
	})
	auprc <- sapply(model.names, function(m) {
		all.model.aucs |> filter(modnames == m & curvetypes == "PRC") |> pull(aucs)
	})

	##### Dynamically build AUROC and AUPRC text lines for each model
	auroc_lines <- paste0(model.names, "=", round(auroc[model.names], 3))
	auprc_lines <- paste0(model.names, "=", round(auprc[model.names], 3))

	##### Make AUCROC and AUPRC texts to include in the plots
	auc_text <- glue(
		"\nAUROC: \n{paste(auroc_lines, collapse = '\n')} \n\nAUPRC: \n{paste(auprc_lines, collapse = '\n')}"
	)

	##### Make plot text
	auc_grob <- textGrob(
		auc_text,
		x = 0, y = 1, just = c("left", "top"),
		gp = grid::gpar(fontsize = 8, fontfamily = "mono")
	)

	list(
		text = auc_text,
		text.plot.object = auc_grob
	)
}

#### Function to mark the cutoff point on the ROC and PRC plots based on different values and save them into separate plots
mark_mobsnvf_cutoff <- function(scores_truth_df, out_dir, cut_thresholds = NULL, title = NULL){

	if (is.null(cut_thresholds)){
		n <- c(1, 2, 3, 4, 5, 6, 7, 8, 10, 12)
		cut_thresholds <- 5 * 10^(-n) |> sort(decreasing = TRUE)
	}

	print(glue("Determining mobsnvf cutoff point across different values of fp-cut:"))

	for (threshold_value in cut_thresholds) {

		print(glue("\t{threshold_value}"))
		mobsnvf.cut.metrics <- scores_truth_df |> 
			select(chrom, pos, ref, alt, snv, truth, FOBP) |>
			# Restore the original uninverted score
			mutate(FOBP = -FOBP) |>
			fdr_cut_pred("FOBP", fp.cut = threshold_value) |>
			calc.eval.metrics()
		
		mobsnvf.cut.roc.plot <- all.model.roc.plot +
			geom_point(aes(x=mobsnvf.cut.metrics$fp_rate , y=mobsnvf.cut.metrics$tp_rate), color="red", shape=4)

		mobsnvf.cut.prc.plot <- all.model.prc.plot + 
			geom_point(aes(x=mobsnvf.cut.metrics$recall , y=mobsnvf.cut.metrics$precision), color="red", shape =4)


		mobsnvf.cut.roc.prc.plot <- make.roc.prc.plot(
			mobsnvf.cut.roc.plot, 
			mobsnvf.cut.prc.plot, 
			plot.auc.text,
			title = 'title',
			subtitle = subtitle,
			caption = glue("fp-cut: {threshold_value} \n{caption}")
		)
		
		dir.create(glue("{out_dir}/mobsnvf_cutoffs"), recursive = TRUE, showWarnings = FALSE)
		qdraw(mobsnvf.cut.roc.prc.plot, 
			glue("{out_dir}/mobsnvf_cutoffs/{sample_name}_cutoff-{threshold_value}_roc_prc.pdf"),
			width = 7, height = 5
		)
	}
}


#### Box plotting function
make_auc_boxplot <- function(
    df,
    auc_type_col,
    model_col = "model",
    scale = 1,
    text_scale = 1,
    width = 6,
    height = 6,
	title = NULL,
	subtitle = NULL,
    grids = FALSE
){
    # options(repr.plot.width = width * scale, repr.plot.height = height * scale)

	if(is.null(title)){
		title <- c(
			"auroc" = "Area under Receiver Operating Characteristic Curve",
			"auprc" = "Area under Precision-Recall Curve"
		)
	}

    auc_model_boxplot <- ggplot(df, aes(x = .data[[model_col]], y = .data[[auc_type_col]], fill = .data[[model_col]])) +
        geom_boxplot(width = 0.5, alpha = 0.9) +
        geom_jitter(color = "black", size = 0.5*scale, alpha = 0.7) +
        scale_fill_viridis(discrete = TRUE) +
        theme_ipsum(base_family = "sans") +
        theme(
            legend.position = "none",
            plot.title = element_text(size = 12 * text_scale, margin = margin(b = 10), hjust = 0.5), # Control title size, padding and position
            plot.subtitle = element_text(size = 10 * text_scale, margin = margin(b = 15), hjust = 0.5), # subtitle size, padding 
            axis.title.x = element_text(size = 10 * text_scale, margin = margin(t = 10), hjust = 0.5), # x label, padding 
            axis.title.y = element_text(size = 10 * text_scale, margin = margin(r = 10), hjust = 0.5), # y label, padding 
            axis.text.x = element_text(size = 7 * text_scale, color = "grey30"), # x tick labels
            axis.text.y = element_text(size = 7 * text_scale, color = "grey40"), # y tick labels
            # legend.text = element_text(size = 16 * text_scale), # legend text size
            # legend.title = element_text(size = 18 * text_scale, face = "bold"), # legend title size and properties
            # legend.key.size = unit(1.5, "cm"), # legend keys sizes
            # legend.key.width = unit(1.5, "cm"), # legend keys sizes
            strip.text = element_text(size = 9 * text_scale, hjust = 0.5),
            panel.grid.minor = element_line(color = "grey95", linewidth = 0.5), # minor grids properties
            panel.grid.major = element_line(color = "grey95", linewidth = 0.5), # major grids properties
        ) +
        labs(
            x = "Model", 
            y = toupper(auc_type_col),
            title = title,
            subtitle= paste0(if (is.null(subtitle)) "" else subtitle)
        )

    if (grids) {
        auc_model_boxplot <- auc_model_boxplot +
            scale_y_continuous(
                breaks = seq(0, 1, by = 0.1), # major y-tick granularity
                minor_breaks = seq(0, 1, by = 0.02) # minor y-tick granularity
            ) + 
            theme(
                panel.grid.minor = element_line(color = "grey90", linewidth = 0.5), # minor grids properties
                panel.grid.major = element_line(color = "grey85", linewidth = 0.5), # major grids properties
            )
    }

    return(auc_model_boxplot)
}


#### Function to make AUPRC and AUROC plots from multiple variant callers using the make_auc_boxplot function. 
make_all_auc_boxplots <- function(results, scale = 0.4, text_scale = 0.5, w = 14*scale, h = 12*scale, 
	variant_callers = c("all", "muse", "mutect2", "somaticsniper", "varscan2"), subtitle = NULL,
	tissue_types = c("all", "Liver", "Colon")
	) {
	auc_types <- c("auroc", "auprc")
	
	plots <- list()

	for (tissue in tissue_types){
		
		plots[[tissue]] <- list()

		if (tissue != "all"){
			tissue_df <- results |> filter(str_detect(sample_id, tissue))
		} 
		else {
			tissue_df <- results
		}

		for (auc_type in auc_types) {

			plots[[tissue]][[auc_type]] <- list()
			for (caller in variant_callers) {
				if (caller == "all") {
					df <- tissue_df
				} else {
					df <- tissue_df |> filter(variant_caller == caller)
				}
				if (is.null(subtitle)) {
					subtitle <- caller
				} 
				plots[[tissue]][[auc_type]][[caller]] <- make_auc_boxplot(
					df, auc_type, scale = scale, text_scale = text_scale, subtitle = subtitle
				)
			}
		}
	}
	return(plots)
}
