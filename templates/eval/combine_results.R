## This script is used to combine the results after evaluation of the ffpe filtering models
library(io)

## Dataset specific
eval_dirs = c(
		# "PRJEB8754/vcf_pass-orient-pos-sb_ad_filtered/roc-prc-auc/precrec",
		# "SRP044740/vcf_filtered_pass_orientation/roc-prc-auc/precrec",
		# "SRP065941/vcf_filtered_pass_orientation/roc-prc-auc/precrec",
		"PRJEB44073/vcf_filtered_pass_orientation/roc-prc-auc/precrec"
	)

for (eval_dir in eval_dirs){
	## List name of models that were evaluated. 
	models <- c("all-models", "mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm", "microsec")

	## Combine the AUC table for each sample and overall eval i.e all sample, liver samples, colon samples
	message("Combining AUCs")
	auc_paths <- list.files(eval_dir, pattern = "_auc_table.tsv", recursive = TRUE, full.names = TRUE)
	auc_paths <- auc_paths[!grepl("combined", auc_paths)]

	auc <- do.call(
		rbind,
		lapply(
			auc_paths, 
			read.delim
		)
	)

	auc <- auc[order(auc$sample_name, auc$model),]

	qwrite(auc, file.path(eval_dir, "combined_auc_table.tsv"))


	## List samples based on their path
	roc_path <- list.files(eval_dir, pattern = "_roc_coordinates.tsv", recursive = TRUE, full.names = TRUE)
	to_remove <- c(paste("", models, sep = "_"), '_roc_coordinates.tsv','_prc_coordinates.tsv')

	samples <- unique(sapply(
		roc_path, 
		function(i) {
			s <- basename(i)
			for (p in to_remove) {
				s <- gsub(p, "", s, fixed = TRUE)
			}
			s
		}, 
		USE.NAMES = FALSE
	))

	# Combine ROC and PRC results for all models within each sample
	message("Combining ROC and PRC results for the evaluated model within each sample:")
	for (sample in samples){

		message(sample)
		paths <- list.files(eval_dir, pattern = sample, recursive = TRUE, full.names = TRUE)

		## Combine ROC
		roc_paths <- paths[grepl("roc_coordinates.tsv", paths)]
		roc_paths <- roc_paths[!grepl("all-models_roc_coordinates.tsv", roc_paths)]
		
		roc <- do.call(
			rbind,
			lapply(roc_paths, read.delim)
		)

		dir <- dirname(roc_paths[1])
		qwrite(roc, file.path(dir, sprintf("%s_all-models_roc_coordinates.tsv", sample)))


		# Combine PRC
		prc_paths <- paths[grepl("prc_coordinates.tsv", paths)]
		prc_paths <- prc_paths[!grepl("all-models_prc_coordinates.tsv", prc_paths)]
		
		prc <- do.call(
			rbind,
			lapply(prc_paths, read.delim)
		)

		dir <- dirname(prc_paths[1])
		qwrite(prc, file.path(dir, sprintf("%s_all-models_prc_coordinates.tsv", sample)))

	}

	message("Done.")

}