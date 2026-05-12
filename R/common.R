library(vcfR)
library(dplyr)
library(tidyr)

source("strand-bias.R")


#' Calculate GATK Strand Bias (SB-GATK)
#'
#' @param ref_fwd (numeric) Reference forward reads
#' @param ref_rev (numeric) Reference reverse reads
#' @param alt_fwd (numeric) Alternate forward reads
#' @param alt_rev (numeric) Alternate reverse reads
#'
#' @return (numeric) SB-GATK score
calc_sb_gatk <- function(ref_fwd, ref_rev, alt_fwd, alt_rev) {
	tot_F <- ref_fwd + alt_fwd
	tot_R <- ref_rev + alt_rev
	tot_ref <- ref_fwd + ref_rev
	tot_all <- tot_F + tot_R
	
	val1 <- ((alt_fwd / tot_F) * (ref_rev / tot_R)) / (tot_ref / tot_all)
	val2 <- ((alt_rev / tot_R) * (ref_fwd / tot_F)) / (tot_ref / tot_all)
	
	pmax(val1, val2, na.rm = FALSE)
}

#' Calculate Guo Strand Bias (SB-GUO)
#'
#' @param ref_fwd (numeric) Reference forward reads
#' @param ref_rev (numeric) Reference reverse reads
#' @param alt_fwd (numeric) Alternate forward reads
#' @param alt_rev (numeric) Alternate reverse reads
#'
#' @return (numeric) SB-GUO score
calc_sb_guo <- function(ref_fwd, ref_rev, alt_fwd, alt_rev) {
	tot_F <- ref_fwd + alt_fwd
	tot_R <- ref_rev + alt_rev
	tot_alt <- alt_fwd + alt_rev
	tot_all <- tot_F + tot_R
	
	abs((alt_fwd / tot_F) - (alt_rev / tot_R)) / (tot_alt / tot_all)
}


#' Read a VCF into a data.frame
#'
#' @description
#' Read a VCF file using vcfR::read.vcfR and extract AD, F1R2, F2R1, and 
#' various Mutect2-specific INFO/FORMAT fields for a specific sample. 
#' Splits genotype and info fields into reference and alternate metrics.
#' Computes VAF, orientation-bias score (SOB), and converts Phred/Log-Odds 
#' scores to probabilities.
#'
#' @param vcf_path (string) Path to the VCF file.
#' @param sample_name (string) Name to assign to the sample/column in the output.
#'
#' @return (data.frame) A flattened dataframe with all relevant numeric features and probabilities.
#'
#' @examples
#' # df <- read_vcf("sample.vcf", "FFX_GZ_T_1h_1")
read_vcf <- function(vcf_path, sample_name){

	vcf0 <- read.vcfR(vcf_path, verbose = FALSE)
	
	# Extract basic fields, FORMAT fields, and INFO fields
	vcf1 <- data.frame(
		chrom = getCHROM(vcf0),
		pos = getPOS(vcf0),
		ref = getREF(vcf0),
		alt = getALT(vcf0),
		filter = getFILTER(vcf0),
		
		# FORMAT fields (Sample specific)
		ad = as.vector(extract.gt(vcf0, element = "AD")[, sample_name]),
		f1r1 = as.vector(extract.gt(vcf0, element = "F1R2")[, sample_name]),
		f2r1 = as.vector(extract.gt(vcf0, element = "F2R1")[, sample_name]),
		dp = as.vector(extract.gt(vcf0, element = "DP")[, sample_name]),
		fad = as.vector(extract.gt(vcf0, element = "FAD")[, sample_name]),
		sb = as.vector(extract.gt(vcf0, element = "SB")[, sample_name]),
		mutect_af = as.vector(extract.gt(vcf0, element = "AF")[, sample_name]),
		
		# INFO fields (Site specific, but usually represent the tumor in Mutect2 somatic calling)
		mfrl = extract.info(vcf0, element = "MFRL"),
		mbq = extract.info(vcf0, element = "MBQ"),
		mmq = extract.info(vcf0, element = "MMQ"),
		mpos = extract.info(vcf0, element = "MPOS"),
		tlod = extract.info(vcf0, element = "TLOD"),
		roq = extract.info(vcf0, element = "ROQ"),
		germq = extract.info(vcf0, element = "GERMQ")
	)

	# Split FORMAT fields into ref and alt components
	vcf1 <- separate(vcf1, ad, into = c("ad_ref", "ad_alt"), sep=",", extra="merge", convert = TRUE)
	vcf1 <- separate(vcf1, f1r1, into = c("f1r2_ref", "f1r2_alt"), sep=",", extra="merge", convert = TRUE)
	vcf1 <- separate(vcf1, f2r1, into = c("f2r1_ref", "f2r1_alt"), sep=",", extra="merge", convert = TRUE)
	vcf1 <- separate(vcf1, fad, into = c("fad_ref", "fad_alt"), sep=",", extra="merge", convert = TRUE)
	
	# Split Strand Bias (SB) into 4 components: ref_fwd, ref_rev, alt_fwd, alt_rev
	vcf1 <- separate(vcf1, sb, into = c("sb_ref_fwd", "sb_ref_rev", "sb_alt_fwd", "sb_alt_rev"), sep=",", extra="merge", convert = TRUE)

	# Split INFO fields that are reported as "Ref,Alt" arrays
	vcf1 <- separate(vcf1, mfrl, into = c("mfrl_ref", "mfrl_alt"), sep=",", extra="merge", convert = TRUE)
	vcf1 <- separate(vcf1, mbq, into = c("mbq_ref", "mbq_alt"), sep=",", extra="merge", convert = TRUE)
	vcf1 <- separate(vcf1, mmq, into = c("mmq_ref", "mmq_alt"), sep=",", extra="merge", convert = TRUE)

	# Split multiallelic sites to biallelic
	# Added the new alt-specific fields to ensure they duplicate correctly if a site has multiple ALTs
	vcf1 <- separate_rows(vcf1, alt, ad_alt, f1r2_alt, f2r1_alt, fad_alt, mfrl_alt, mbq_alt, mmq_alt, mpos, tlod, mutect_af, sep = ",")

	# Ensure all extracted columns are strictly numeric
	vcf1 <- mutate(
		vcf1,
		pos = as.numeric(pos),
		ad_ref = as.numeric(ad_ref),
		ad_alt = as.numeric(ad_alt),
		f1r2_ref = as.numeric(f1r2_ref),
		f1r2_alt = as.numeric(f1r2_alt),
		f2r1_ref = as.numeric(f2r1_ref),
		f2r1_alt = as.numeric(f2r1_alt),
		dp = as.numeric(dp),
		fad_ref = as.numeric(fad_ref),
		fad_alt = as.numeric(fad_alt),
		sb_ref_fwd = as.numeric(sb_ref_fwd),
		sb_ref_rev = as.numeric(sb_ref_rev),
		sb_alt_fwd = as.numeric(sb_alt_fwd),
		sb_alt_rev = as.numeric(sb_alt_rev),
		mutect_af = as.numeric(mutect_af),
		mfrl_ref = as.numeric(mfrl_ref),
		mfrl_alt = as.numeric(mfrl_alt),
		mbq_ref = as.numeric(mbq_ref),
		mbq_alt = as.numeric(mbq_alt),
		mmq_ref = as.numeric(mmq_ref),
		mmq_alt = as.numeric(mmq_alt),
		mpos = as.numeric(mpos),
		tlod = as.numeric(tlod),
		roq = as.numeric(roq),
		germq = as.numeric(germq)
	)

	# Calculate custom metrics
	vcf1 <- mutate(
		vcf1,
		vaf = ad_alt / (ad_ref + ad_alt),
		sob = abs((f1r2_alt - f2r1_alt) / (f1r2_alt + f2r1_alt)),
		sb_gatk = calc_sb_gatk(sb_ref_fwd, sb_ref_rev, sb_alt_fwd, sb_alt_rev),
		sb_guo = calc_sb_guo(sb_ref_fwd, sb_ref_rev, sb_alt_fwd, sb_alt_rev),
		# sb_mutect = post_strand_bias(sb_ref_fwd, sb_ref_rev, sb_alt_fwd, sb_alt_rev),
		fdeamc = case_when(
			ref == "C" & alt == "T" ~ f1r2_alt / (f1r2_alt + f2r1_alt),
			ref == "G" & alt == "A" ~ f2r1_alt / (f1r2_alt + f2r1_alt),
			TRUE ~ NA_real_ # Assigns NA to any other mutation types
		)
	)

	vcf1 <- mutate(
		vcf1,
		sb_mutect = mapply(
			post_strand_bias,
			sb_ref_fwd, sb_ref_rev, sb_alt_fwd, sb_alt_rev
		)
	)

	# Convert Phred-scaled qualities and Log-Odds to Probabilities
	vcf1 <- mutate(
		vcf1,
		# Phred conversions: P = 10^(-Q/10)
		mbq_ref_prob_err = 10^(-mbq_ref / 10),
		mbq_alt_prob_err = 10^(-mbq_alt / 10),
		mmq_ref_prob_err = 10^(-mmq_ref / 10),
		mmq_alt_prob_err = 10^(-mmq_alt / 10),
		roq_prob         = 10^(-roq / 10),
		germq_prob       = 10^(-germq / 10),
		
		# TLOD conversion: Log10-Odds to Probability: P = 10^TLOD / (1 + 10^TLOD)
		tlod_prob        = (10^tlod) / (1 + 10^tlod)
	)

	vcf1$sample_name <- sample_name

	return(vcf1)
}


#' Calculate True/False Positives/Negatives from prediction and ground truth
#' @param df [data.frame] with models prediction and ground truth columns
#' @param pred_col [string] prediction column name
#' @param truth_col [string] truth column name
#' @return [list] with numeric FP, TP, FN, and TN
calc_confusion_matrix <- function(df, pred_col = "pred", truth_col = "truth") {
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


#' Calculate evaluation metrics based on confusion matrix
#' @param df [data.frame] with models prediction and ground truth columns
#' @param pred_col [string] prediction column name
#' @param truth_col [string] truth column name
#' @return [list] with sensitivity, specificity, precision, and recall
calc_eval_metrics <- function(df, pred_col = "pred", truth_col = "truth") {

	confusion_matrix <- calc_confusion_matrix(df, pred_col, truth_col)

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