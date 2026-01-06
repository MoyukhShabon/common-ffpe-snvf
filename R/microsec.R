#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# MicroSEC Wrapper Script (Single Sample Mode)
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
	library(argparser)
	library(MicroSEC)
	library(dplyr)
	library(readr)
	library(stringr)
	library(Rsamtools)
	library(BiocGenerics)
	library(Biostrings)
	library(GenomicRanges)
	library(IRanges)
	library(tictoc)
})


# ------------------------------------------------------------------------------
# Argument Parsing
# ------------------------------------------------------------------------------

p <- arg_parser("MicroSEC: Mutation filtering pipeline (Single Sample)")
p <- add_argument(p, "--sample_info", help = "Path to the sample information TSV file (1 row)", type = "character")
p <- add_argument(p, "--output_dir", help = "Directory to save final outputs", type = "character")
p <- add_argument(p, "--progress", help = "Show progress bar (Y/N)", default = "N", type = "character")

argv <- parse_args(p)

# Validate inputs
if (is.na(argv$sample_info) || !file.exists(argv$sample_info)) {
	stop("Error: Sample info file not found.")
}

dir.create(argv$output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------

#' Execute system command with error checking
run_cmd <- function(cmd) {
	exit_code <- system(cmd)
	if (exit_code != 0) {
		stop(paste("Command failed with exit code", exit_code, ":", cmd))
	}
}

#' Calculate download regions using GenomicRanges
get_target_regions <- function(df_mutation) {
	# Create GRanges object from mutations
	gr <- GRanges(
		seqnames = df_mutation$Chr,
		ranges = IRanges(start = df_mutation$Pos, end = df_mutation$Pos)
	)

	# Merge features within 400bp of each other
	gr_reduced <- reduce(gr, min.gapwidth = 400)

	# Add padding (200bp on both sides)
	start(gr_reduced) <- pmax(1, start(gr_reduced) - 200)
	end(gr_reduced) <- end(gr_reduced) + 200

	# Convert to data frame for BED format
	df_bed <- data.frame(
		chrom = seqnames(gr_reduced),
		chromStart = start(gr_reduced) - 1, # BED format uses 0-based start
		chromEnd = end(gr_reduced) # BED format uses 1-based end
	)

	return(df_bed)
}


#' Process BAM file: Sort, Index, Subset (Slim)
#' Saves files directly to inputs_dir/slim-bam/filename.bam
subset_bam <- function(bam_file, sample_name, df_mutation, ref_genome_path, inputs_dir) {

	bam_basename <- sub("\\.bam$|\\.cram$", "", basename(bam_file))
	bam_file_slim <- file.path(inputs_dir, paste0(bam_basename, ".SLIM.bam"))
	bam_file_tmp <- file.path(inputs_dir, paste0(bam_basename, ".tmp.bam"))
	bed_file <- file.path(inputs_dir, paste0(bam_basename, ".bed"))

	# Check if SLIM BAM already exists
	if (file.exists(bam_file_slim) && file.exists(paste0(bam_file_slim, ".bai"))) {
		message(paste("Slim BAM already exists:", bam_file_slim))
		return(bam_file_slim)
	}

	# 1. Handle Indexing / Sorting of original file if needed
	is_cram <- tools::file_ext(bam_file) == "cram"
	index_ext <- if (is_cram) ".crai" else ".bai"

	# Check for index file (either file.bam.bai or file.bai)
	if (!(file.exists(paste0(bam_file, index_ext)) || file.exists(gsub("\\.bam$|\\.cram$", index_ext, bam_file)))) {
		message(paste("Sorting/Indexing original file:", bam_file))

		sorted_path <- file.path(sorted_bam_dir, paste0(bam_basename, "_sort", if (is_cram) ".cram" else ".bam"))

		# Sort
		sort_cmd <- paste("samtools sort -@ 4",
						if (is_cram) "-O cram" else "",
						"-o", sorted_path, bam_file)
		run_cmd(sort_cmd)

		# Index
		run_cmd(paste("samtools index -@ 4", sorted_path))

		# Update bam_file to point to the new sorted file
		bam_file <- sorted_path
	}

	# 2. Create BED file
	message("Calculating target regions...")
	df_bed <- get_target_regions(df_mutation)
	# write_tsv is from readr, ensures tab separation
	write_tsv(df_bed, file = bed_file, col_names = FALSE, progress = FALSE)

	# 3. Subset (Slim) the BAM
	message(paste("Subsetting BAM to:", bam_file_slim))

	ref_opt <- if (is_cram && !is.null(ref_genome_path)) paste("-T", ref_genome_path) else ""

	# View regions to Temporary BAM
	cmd_view <- paste("samtools view -b -h --no-PG", ref_opt,
						"-L", bed_file,
						bam_file, ">", bam_file_tmp)
	run_cmd(cmd_view)

	# 4. Sort and Index the Slim BAM
	cmd_sort_slim <- paste("samtools sort -@ 4 -o", bam_file_slim, bam_file_tmp)
	run_cmd(cmd_sort_slim)

	cmd_index_slim <- paste("samtools index -@ 4", bam_file_slim)
	run_cmd(cmd_index_slim)

	# Cleanup
	if (file.exists(bam_file_tmp)) file.remove(bam_file_tmp)

	return(bam_file_slim)
}


# ------------------------------------------------------------------------------
# Main Execution
# ------------------------------------------------------------------------------

tic("Total Execution")

# Load Sample Info
sample_info <- read_tsv(argv$sample_info, show_col_types = FALSE)

if (nrow(sample_info) != 1) {
	stop(paste("Expected exactly 1 row in sample_info, found:", nrow(sample_info)))
}

# Extract Row Data
row <- sample_info[1, ]
s_name <- row$sample_name
mut_file <- row$mutation_file
bam_path <- row$bam_file
read_len <- as.integer(row$read_length)
adapt1   <- row$adapter_1
adapt2   <- if ("adapter_2" %in% colnames(row) && !is.na(row$adapter_2)) row$adapter_2 else adapt1
organism_id <- row$genome_build
ref_genome_path <- if ("ref_genome" %in% colnames(row)) row$ref_genome else NULL

message(paste0("Processing Sample: ", s_name))


# --- 1. Load Genome & Mutations ---
tic("Data Loading")

ref_genome <- fun_load_genome(organism_id)
chr_no <- fun_load_chr_no(organism_id)

df_mutation <- fun_load_mutation(
	mut_file,
	s_name,
	ref_genome,
	chr_no
)

toc(log = TRUE)

# --- 2. Prepare BAM ---
tic("BAM Processing")
# Determine inputs directory based on mutation file location
inputs_dir <- dirname(argv$sample_info)

final_bam_path <- subset_bam(
	bam_file = bam_path,
	sample_name = s_name,
	df_mutation = df_mutation,
	ref_genome_path = ref_genome_path,
	inputs_dir = inputs_dir
)

toc(log = TRUE)

# --- 3. MicroSEC Analysis ---
tic("Main Analysis")
print("Running Read Check...")
df_bam <- fun_load_bam(final_bam_path)

read_check_res <- fun_read_check(
	df_mutation = df_mutation,
	df_bam = df_bam,
	ref_genome = ref_genome,
	sample_name = s_name,
	read_length = read_len,
	adapter_1 = adapt1,
	adapter_2 = adapt2,
	short_homology_search_length = 4,
	min_homology_search = 40,
	progress_bar = argv$progress
)


msec <- read_check_res[[1]]
homology_search <- read_check_res[[2]]
mut_depth <- read_check_res[[3]]

print("Running Homology Search...")
homology_res <- fun_homology(
	msec,
	homology_search,
	min_homology_search = 40,
	ref_genome = ref_genome,
	chr_no = chr_no,
	progress_bar = argv$progress
)



print("Running Statistical Analysis...")
summary_res <- fun_summary(homology_res)


msec_final <- fun_analysis(
	msec = summary_res,
	mut_depth = mut_depth,
	short_homology_search_length = 4,
	min_homology_search = 40,
	threshold_p = 10^(-6),
	threshold_hairpin_ratio = 0.50,
	threshold_short_length = 0.75,
	threshold_distant_homology = 0.15,
	threshold_soft_clip_ratio = 0.50,
	threshold_low_quality_rate = 0.1,
	homopolymer_length = 15
)

toc(log = TRUE)

# --- 4. Save Results ---
out_file <- file.path(argv$output_dir, paste0(s_name, ".microsec.tsv"))
write.table(msec_final, out_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

message(paste("Finished:", s_name))

toc(log = TRUE)

# --- 5. Process Timing Logs ---
raw_log <- tic.log(format = FALSE)
timings_df <- do.call(rbind, lapply(raw_log, function(x) {
	data.frame(
		section = x$msg,
		seconds = x$toc - x$tic
	)
}))
timings_df$minutes <- round(timings_df$seconds / 60, 4)
timings_df$seconds <- round(timings_df$seconds, 4)

write.table(timings_df, file.path(argv$output_dir, paste0(s_name, ".runtime.tsv")), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
