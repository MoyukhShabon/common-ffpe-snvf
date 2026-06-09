#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(optparse)
	library(FIREVAT)
})

# ---- Default Mutect2 config added inside the package ----
default.config <- system.file("config",
                              "GATK4_Mutect2_config_Tumor_Normal_SB.json",
                              package = "FIREVAT")

# ---- Argument parsing ----
option_list <- list(
	make_option(c("-v", "--vcf"), 
				type = "character", default = NULL,
				help = "Path to input VCF file [required]", metavar = "FILE"),
	make_option(c("-c", "--config"), type = "character", 
				default = default.config,
				help = "Path to config file [default: packaged GATK4_Mutect2_config_Tumor_Normal_SB.json]",
				metavar = "FILE"),
	make_option(c("-o", "--outdir"), 
				type = "character", 
				default = NULL,
				help = "Output directory [required]", 
				metavar = "DIR"),
	make_option(c("-g", "--genome"), 
				type = "character", 
				default = "hg38",
				help = "Reference genome [default: %default]", metavar = "GENOME"),
	make_option(c("-n", "--num-cores"),
				type = "integer",
				dest = "num.cores",
				default = 1,
				help = "Number of cores [default: %default]", 
				metavar = "N"),
    make_option(c("--no-strandbias"),
				action = "store_false",
                default= TRUE,
				dest="strandbias",
                help="Don't perform strand bias analysis [default: FALSE]"),
    make_option(c("-p","--popsize"),
                type="integer",
                default= 200,
                dest="pop.size",
                help="Population size for GA algorithm [default\"%default\"]"),
    make_option(c("-i","--maxiter"),
                type="integer",
                default= 50,
                dest="max.iter",
                help="Maximum iteration for GA algorithm [default\"%default\"]"),
    make_option(c("-r","--run"),
                type="integer",
                default= 25,
                help="Limit for consecutive generations without changes [default\"%default\"]"),
    make_option("--pmutation",
                type="numeric",
                default= 0.25,
                dest="p.mutation",
                help="Mutation probability for GA algorithm [default\"%default\"]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---- Validate required arguments ----
required <- c("vcf", "outdir")
missing <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing) > 0) {
	stop("Missing required argument(s): ", paste0("--", missing, collapse = ", "),
			 call. = FALSE)
}

# config has a default, but make sure it actually resolved/exists
if (is.null(opt$config) || opt$config == "") {
	stop("No config supplied and packaged default could not be located. ",
			 "Pass one explicitly with --config.", call. = FALSE)
}

if (!file.exists(opt$vcf))    stop("VCF file not found: ", opt$vcf, call. = FALSE)
if (!file.exists(opt$config)) stop("Config file not found: ", opt$config, call. = FALSE)
if (!dir.exists(opt$outdir))  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ---- Run FIREVAT ----
message("Running FIREVAT")
message("  VCF:         ", opt$vcf)
message("  Config:      ", opt$config)
message("  Outdir:      ", opt$outdir)
message("  Genome:      ", opt$genome)
message("  Cores:       ", opt$num.cores)
message("  Strand bias: ", opt$strandbias)
message("  Pop size:    ", opt$pop.size)
message("  Max iter:    ", opt$max.iter)
message("  Run:         ", opt$run)
message("  P mutation:  ", opt$p.mutation)

res <- RunFIREVAT(
	vcf.file = opt$vcf,
	vcf.file.genome = opt$genome,
	config.file = opt$config,
	df.ref.mut.sigs = GetPCAWGMutSigs(),
	target.mut.sigs = GetPCAWGMutSigsNames(),
	sequencing.artifact.mut.sigs = PCAWG.All.Sequencing.Artifact.Signatures,
	output.dir = opt$outdir,
	objective.fn = Default.Obj.Fn,
	num.cores = opt$num.cores,
	ga.pop.size = opt$pop.size,
	ga.max.iter = opt$max.iter,
	ga.run = opt$run,
	ga.pmutation = opt$p.mutation,
	perform.strand.bias.analysis = opt$strandbias,
	ref.forward.strand.var = "TumorDPRefForward",
	ref.reverse.strand.var = "TumorDPRefReverse",
	alt.forward.strand.var = "TumorDPAltForward",
	alt.reverse.strand.var = "TumorDPAltReverse",
	annotate = FALSE,
	write.vcf = FALSE,
	report = FALSE,
	save.rdata = FALSE,
	save.tsv = TRUE,
	save.scores = TRUE,
	report.format = "html",
	verbose = TRUE
)

message("FIREVAT complete. Results written to: ", opt$outdir)
