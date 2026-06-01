#!/usr/bin/env Rscript

# This script contains functions for making ROC and PRC plots

### Load necessary libraries


library(ggplot2)
library(ggpubr)
library(arrow)
library(patchwork)
library(LTTB)
library(io)
# library(dplyr)
# library(glue)
# library(grid)
# library(hrbrthemes)
# library(viridis)

### Plotting

MODEL_COLORS <- c(
	"MOBSNVF"     = "#E41A1C", # Crimson Red      
	"VAFSNVF"     = "#FF7F00", # Standard Orange      
	"GATK-OBMM"   = "#FFC020", # Golden Yellow    
	"SOBDetector" = "#377EB8", # Solid Blue      
	"Ideafix"     = "#4DAF4A", # Solid Green      
	"FFPolish"    = "#984EA3",  # Solid Purple
	"FFPErase"    = "#a04300" # Solid Brown
)
# High opacity for MOBSNVF and reference VAFSNVF, lower for others
MODEL_ALPHAS <- c(
	"MOBSNVF" = 0.9,
	"VAFSNVF" = 0.9,
	"GATK-OBMM" = 0.6,
	"SOBDetector" = 0.6,
	"MicroSEC" = 0.6,
	"Ideafix" = 0.6,
	"FFPolish" = 0.6,
	"FFPErase" = 0.6
)

#' Function to check if the roc and prc coordinate table is from the same sample and returns sample name if they match
#' @param roc_path (string) | path to roc coordinates table
#' @param prc_path (string) | path to prc coordinates table
#' @param roc_suffix (string) | default: "_all-models_roc_coordinates.tsv" | roc coordinates table suffix after sample name
#' @param prc_suffix (string) | default: "_all-models_prc_coordinates.tsv" | prc coordinates table suffix after sample name
#' @return (string) name of the sample
#' To be used by **/eval{...}/plot.R
match_return_sample_name <- function(roc_path, prc_path, roc_suffix = "_all-models_roc_coordinates.tsv", prc_suffix = "_all-models_prc_coordinates.tsv"){

	sample_name_roc <- gsub(roc_suffix, "", basename(roc_path), fixed = TRUE)
	sample_name_prc <- gsub(prc_suffix, "", basename(prc_path), fixed = TRUE)
	
	if (sample_name_roc != sample_name_prc) {
		stop(paste("Mismatch found... '", 
			sample_name_roc, "' does not match PRC sample '", 
			sample_name_prc, "'. Halting execution.", sep=""))
	}

	sample_name_roc
}

#' Compute baseline precision (positive-class prevalence) from a sample's ground truth TSV.
#' @param dir `character(1)`. The directory in which model-scores_truth exisit.
#' @param sample_name `character(1)`. Sample identifier used to locate the `model-scores_truths/.../<sample_name>_mobsnvf-scores_truths.tsv` file.
#'   If it contains `"all-samples"`, the file is read from the top-level `model-scores_truths/` directory instead of a per-sample subfolder.
#' @return `numeric(1)` in `[0, 1]`: mean of the logical `truth` column.
#' To be used by **/eval{...}/plot.R
get_baseline_precision <- function(dir, sample_name) {
	.read_file <- function(path) {
		if (!file.exists(path)) {
			stop(sprintf(
				"Table for sample '%s' with ground truth not found at: %s",
				sample_name, path
			))
		}
		df <- read_delim_arrow(path, delim = "\t")
		if (!is.logical(df$truth)) {
			stop(sprintf(
				"Ground truth column is not logical (got %s) in file: %s",
				typeof(df$truth), path
			))
		}
		return(df)
	}

	path <- if (grepl("all-", sample_name) & grepl("-samples", sample_name)) {
		file.path(dir, "model-scores_truths", paste0(sample_name, "_mobsnvf-scores_truths.tsv"))
	} else {
		file.path(dir, "model-scores_truths", sample_name, paste0(sample_name, "_mobsnvf-scores_truths.tsv"))
	}

	df <- .read_file(path)
	baseline_precision <- mean(df$truth)
	return(baseline_precision)
}

#' Function to make roc prc plot based on roc and prc coordinate
#' @param roc_coord (data.frame) | dataframe of roc coordinates
#' @param prc_coord (data.frame) | dataframe of prc coordinates
#' @param baseline_precision (numeric) | If not NULL, draws a dotted horizontal line at this y-value on the PRC plot (typically the positive class prevalence). Default: NULL
#' @param title (string) | Plot title. Default: NULL
#' @param subtitle (string) | Plot subtitle. Default: NULL
#' @param caption (string) | Plot caption. Default: "Only C>T SNVs Evaluated"
#' @param x_col (name) | The name of the column to be used for the x-axis. Default: x
#' @param y_col (name) | The name of the column to be used for the y-axis. Default: y
#' @param model_col (name) | The name of the column for color grouping. Default: model
#' @param base_size (numeric) | Base font size in points; theme_pubr derives all other text sizes from this. Default: 20
#' @param (base_size * title_scale) (numeric) | Font size of the overall (patchwork) plot title. Default: 30
#' @param line_width (numeric) | Width of the ROC/PRC curves. Default: 2
#' @param subplot_legend_ncol (integer) | Number of legend columns for the standalone ROC / PRC subplots returned in the list. Default: 3
#' @param combined_legend_ncol (integer) | Number of legend columns for the combined (patchwork) ROC + PRC plot. Applied right before stitching. Default: 4
#' @param individual_plots (boolean) | creates separate ROC and PRC plots alongside combined ROC PRC plot. Returns list of plot objects (roc_prc, roc, prc)
#' @param downsample_threshold (integer) | Per-model row count above which LTTB downsampling is applied. Default: 5000. Set to Inf to disable.
#' @param downsample_n (integer) | Target number of points per model after LTTB downsampling. Default: 2000. Actual output will be at most n_bins + 2 where n_bins = downsample_n - 2.
#' @param model_colors (named character) | Named vector mapping model names to hex colors. Default: module-level `model_colors`.
#' @param model_alphas (named numeric) | Named vector mapping model names to alpha (opacity) values in [0, 1]. Default: module-level `model_alphas`.
#' @return ggplot2 object (or named list if individual_plots = TRUE) containing the ROC and PRC plots
make_roc_prc_plot <- function(
	roc_coord,
	prc_coord,
	x_col                = x,
	y_col                = y,
	baseline_precision   = NULL,
	title                = NULL,
	subtitle             = NULL,
	caption              = NULL,
	model_col            = model,
	base_size            = 20,
	title_scale          = 1.5,
	subtitle_scale       = 1.2,
	line_width           = 2,
	subplot_legend_ncol  = 3,
	combined_legend_ncol = 6,
	individual_plots     = TRUE,
	downsample_threshold = 5000,
	downsample_n         = 2000,
	model_colors         = MODEL_COLORS,
	model_alphas         = MODEL_ALPHAS
) {

	# ── LTTB Downsampling (per model) ─────────────────────────────────
	# Resolve tidy-eval column names to strings for programmatic use
	x_nm     <- rlang::as_name(rlang::enquo(x_col))
	y_nm     <- rlang::as_name(rlang::enquo(y_col))
	model_nm <- rlang::as_name(rlang::enquo(model_col))

	.apply_lttb_per_model <- function(df, coord_label) {
		models <- unique(df[[model_nm]])
		result <- lapply(models, function(m) {
			sub <- df[df[[model_nm]] == m, , drop = FALSE]
			if (nrow(sub) <= downsample_threshold) return(sub)

			message(
				sprintf(
					"[%s] Downsampling '%s' from %d to ~%d points via LTTB",
					coord_label, m, nrow(sub), downsample_n
				)
			)

			# LTTB::LTTB expects a 2-column numeric matrix sorted by column 1
			# and returns a matrix with at most n_bins + 2 rows
			mat <- data.frame(sub[[x_nm]], sub[[y_nm]])
			downsampled <- LTTB(mat, n_bins = downsample_n - 2)

			# Match downsampled points back to original rows to preserve all columns
			# LTTB selects existing points (no interpolation), so exact matching works
			key_orig <- paste(sub[[x_nm]], sub[[y_nm]], sep = "\t")
			key_ds   <- paste(downsampled[, 1], downsampled[, 2], sep = "\t")
			idx      <- match(key_ds, key_orig)

			sub[idx, , drop = FALSE]
		})
		do.call(rbind, result)
	}

	roc_coord <- .apply_lttb_per_model(roc_coord, "ROC")
	prc_coord <- .apply_lttb_per_model(prc_coord, "PRC")
	# ───────────────────────────────────────────────────────────────────

	# ── Shared style layers for both ROC and PRC plots ────────────────
	# All text sizes are derived from theme_pubr(base_size); no manual
	# axis/legend size overrides so the proportions stay as verified.
	# `legend_ncol` is parameterized so we can build subplots with one
	# layout and re-emit the combined patchwork plot with a different one.
	build_style <- function(legend_ncol) {
		list(
			geom_line(linewidth = line_width),
			coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1), expand = FALSE),
			scale_color_manual(values = model_colors),
			scale_alpha_manual(values = model_alphas),
			theme_pubr(base_size = base_size),
			guides(
				color = guide_legend(ncol = legend_ncol),
				alpha = guide_legend(ncol = legend_ncol)
			),
			labs(color = "", alpha = ""),
			theme(
				legend.position = "bottom",
				plot.title      = element_text(hjust = 0.5, face = "bold", margin = margin(t = 10, b = 10)),
				plot.subtitle   = element_text(hjust = 0.5, face = "bold")
			)
		)
	}

	# ROC Plot (uses subplot legend layout)
	roc_plot <- ggplot(
		roc_coord,
		aes(x = {{ x_col }}, y = {{ y_col }}, color = {{ model_col }}, alpha = {{ model_col }})
	) +
		geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
		build_style(subplot_legend_ncol) +
		labs(
			title = title, 
			subtitle = "ROC", 
			x = "False positive rate", 
			y = "True positive rate"
		)

	# PRC Plot (uses subplot legend layout)
	prc_plot <- ggplot(
		prc_coord,
		aes(x = {{ x_col }}, y = {{ y_col }}, color = {{ model_col }}, alpha = {{ model_col }})
	)

	# Optional baseline (positive-class prevalence) on the PRC plot
	if (!is.null(baseline_precision)) {
		if (!is.numeric(baseline_precision) || length(baseline_precision) != 1) {
			stop("`baseline_precision` must be a single numeric value.")
		}
		prc_plot <- prc_plot +
			geom_hline(
				yintercept = baseline_precision,
				linetype = "dashed",
				color = "grey"
			)
	}

	prc_plot <- prc_plot +
		build_style(subplot_legend_ncol) +
		labs(
			title = title, 
			subtitle = "PRC", 
			x = "Recall", 
			y = "Precision"
		)

	# ── Combined patchwork plot ───────────────────────────────────────
	# Override the legend column count to `combined_legend_ncol` before
	# stitching so the shared legend at the bottom uses a layout suited
	# to the wider combined figure.
	roc_for_patch <- roc_plot + labs(title = NULL) +
		guides(
			color = guide_legend(ncol = combined_legend_ncol),
			alpha = guide_legend(ncol = combined_legend_ncol)
		)
	prc_for_patch <- prc_plot + labs(title = NULL) +
		guides(
			color = guide_legend(ncol = combined_legend_ncol),
			alpha = guide_legend(ncol = combined_legend_ncol)
		)

	roc_prc_plot <- (roc_for_patch + prc_for_patch +
		plot_layout(guides = "collect") &
		theme(legend.position = "bottom")) +
		plot_annotation(
			title = title,
			subtitle = subtitle,
			caption = caption,
			theme = theme(
				plot.title = element_text(size = (base_size * title_scale), face = "bold", hjust = 0.5, margin = margin(t = 10, b = 10)),
				plot.subtitle = element_text(size = (base_size * subtitle_scale), hjust = 0.5, margin = margin(b = 10)
				)
			)
		)

	if (!individual_plots) {
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
