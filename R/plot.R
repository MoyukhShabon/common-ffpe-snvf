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
	individual_plots = TRUE,
	model_colors = c(
		"MOBSNVF"     = "#E41A1C", # Crimson Red      
		"VAFSNVF"     = "#FF7F00", # Standard Orange      
		"GATK-OBMM"   = "#FFC020", # Golden Yellow    
		"SOBDetector" = "#377EB8", # Solid Blue      
		"Ideafix"     = "#4DAF4A", # Solid Green      
		"FFPolish"    = "#984EA3"  # Solid Purple    
	),
	# High opacity for MOBSNVF and reference VAFSNVF, lower for others
	model_alphas = c(
		"MOBSNVF" = 0.9,
		"VAFSNVF" = 0.9,
		"GATK-OBMM" = 0.6,
		"SOBDetector" = 0.6,
		"MicroSEC" = 0.6,
		"Ideafix" = 0.6,
		"FFPolish" = 0.6
	)
) {

	# ROC Plot
	# Added alpha to the aesthetic mapping
	roc_plot <- ggplot(roc_coord, aes(x = {{ x_col }}, y = {{ y_col }}, color = {{ model_col }}, alpha = {{ model_col }})) +
		geom_abline(linetype = "dashed", color = "lightgrey") +
		geom_line(linewidth = line_width) +
		coord_fixed() +
		labs(
			title = "ROC",
			x = "1 - Specificity",
			y = "Sensitivity",
			color = "Models"
		) +
		coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
		scale_color_manual(values = model_colors) +
		scale_alpha_manual(values = model_alphas, guide = "none") + # guide="none" keeps legend colors 100% opaque
		theme_minimal() +
		theme(
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.line = element_line(color = "darkgrey"),
			axis.ticks = element_line(color = "darkgrey"),
			legend.position = "bottom",
			legend.title = element_blank(),
			legend.key.width = unit(line_width*legend_scale, "cm"),
			legend.text = element_text(size = 10*text_scale*legend_scale),
			axis.title.x = element_text(size = 10*text_scale),
			axis.title.y = element_text(size = 10*text_scale),
			axis.text = element_text(size = 8*text_scale),
			plot.title = element_text(size = 12*text_scale, face = "plain", hjust = 0.5)
		)

	# PRC Plot
	# Added alpha to the aesthetic mapping
	prc_plot <- ggplot(prc_coord, aes(x = {{ x_col }}, y = {{ y_col }}, color = {{ model_col }}, alpha = {{ model_col }})) +
		geom_line(linewidth = line_width) +
		coord_fixed() +
		labs(
			title = "PRC",
			x = "Recall",
			y = "Precision",
			color = "Models"
		) +
		coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
		scale_color_manual(values = model_colors) +
		scale_alpha_manual(values = model_alphas, guide = "none") + # guide="none" keeps legend colors 100% opaque
		theme_minimal() +
		theme(
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.line = element_line(color = "darkgrey"),
			axis.ticks = element_line(color = "darkgrey"),
			legend.position = "bottom",
			legend.title = element_blank(),
			legend.key.width = unit(line_width*legend_scale, "cm"),
			legend.text = element_text(size = 10*text_scale*legend_scale),
			axis.title.x = element_text(size = 10*text_scale),
			axis.title.y = element_text(size = 10*text_scale),
			axis.text = element_text(size = 8*text_scale),
			plot.title = element_text(size = 12*text_scale, face = "plain", hjust = 0.5)
		)

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
