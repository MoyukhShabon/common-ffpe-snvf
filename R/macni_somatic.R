log_sum_exp <- function(x) {
	x.max <- max(x);
	log(sum(exp(x - x.max))) + x.max
}

#' @param chom  chromosome
#' @param pos   position
#' @param ac    reference allele counts
#' @param bc    alternative allele counts
#' @param cnp   CNProfile object
call_somatic <- function(chrom, pos, ac, bc, cnp, hparams, scale=FALSE, log=FALSE) {
	aligned <- prepare_somatic(chrom, pos, cnp, scale=scale);
	somatic_model(ac, bc, aligned$alpha, aligned$beta, hparams=hparams, log=log)
}

find_gaps <- function(idx) {
	unaligned.idx <- is.na(idx);
	bounds <- diff(c(0, as.integer(unaligned.idx)));
	starts <- which(bounds > 0);
	ends <- which(bounds < 0) - 1;
	if (length(ends) < length(starts)) {
		ends <- c(ends, length(bounds));
	}
	cbind(start=starts, end=ends)
}

# Assign values to hyperparameter alpha and beta values to each query variant
# using the fitted segments in cnp.
# For query variants in gaps between segments, hyperparameter values may
# be assigned based on the nearest upstream and downstream loci
prepare_somatic <- function(chrom, pos, cnp, scale=FALSE) {
	if (is.character(chrom)) {
		chrom <- factor(chrom, levels=union(levels(cnp@chrom), unique(chrom)));
	}

	ranges <- IRanges(start = pos, end = pos);
	gr <- GRanges(seqnames = chrom, ranges = ranges, strand = NULL);

	# find segments that span each variant
	idx <- findOverlaps(gr, cnp@segments, select="first");
	gaps <- find_gaps(idx);

	# find index of chromosome starts and ends
	chrom <- as.integer(chrom);
	chrom.starts <- which(diff(c(0, as.integer(chrom))) > 0);
	chrom.ends <- c(chrom.starts[-1] - 1, length(chrom));

	# assign an upstream locus only if variant is not at the beginning of a chromsome
	upstreams <- ifelse(gaps[, 1] %in% chrom.starts, NA, gaps[, 1] - 1);
	# assign an dowstream locus only if variant is not at the end of a chromsome
	downstreams <- ifelse(gaps[, 2] %in% chrom.ends, NA, gaps[, 2] + 1);

	# fill in gaps with upstream match
	# fill in gaps with downstream match with alternative assignment
	idx2 <- rep(NA, length(idx));
	for (g in 1:nrow(gaps)) {
		gap <- gaps[g, 1] : gaps[g, 2];
		# assign to upstream locus
		if (!is.na(upstreams[g])) {
			idx[gap] <- idx[upstreams[g]];
		}
		# assign to downstream locus in alternative assignment
		if (!is.na(downstreams[g])) {
			idx2[gap] <- idx[downstreams[g]];
		}
	}

	cnp.alpha <- cnp@segments$alpha;
	cnp.beta <- cnp@segments$beta;

	if (scale) {
		ne <- cnp@segments$f_ne;
		scale.f <- ne / (cnp.alpha + cnp.beta);
		cnp.alpha <- cnp.alpha * scale.f;
		cnp.beta <- cnp.beta * scale.f;
	}
	
	# obtain hyperparameters for the prior on theta^{(0)}_j
	alpha <- cbind(
			cnp.alpha[idx],
			ifelse(is.na(idx2), NA, cnp.alpha[idx2])
	);
	beta <- cbind(
		cnp.beta[idx],
		ifelse(is.na(idx2), NA, cnp.beta[idx2])
	);


	list(idx = cbind(idx, idx2), alpha = alpha, beta = beta)
}


#' @param ac (vector(int)) reference allele supporting read count
#' @param bd (vector(int)) alternate allele suppporting read count
#' @param alpha  matrix of values of hyperameter for prior on germline alt
#'               allele frequency; each column represents one assignment
#' @param beta   matrix of values of hyperameter for prior on germline alt
#'               allele frequency; each column represents one assignment
#' @param log  whether to return the probability in log scale
somatic_model <- function(ac, bc, alpha, beta, hparams=NULL, log=TRUE) {

	if (length(ac) != length(bc)) {
		stop("ac and bc must have the same length");
	}

	if (nrow(alpha) != nrow(beta)) {
		stop("alpha and beta must have the same length");
	}

	if (length(ac) != nrow(alpha)) {
		stop("ac, bc, alpha, and beta must have the same length");
	}

	if (is.null(hparams)) {
		hparams <- list();
	}

	if (is.null(hparams$purity_alpha)) {
		hparams$purity_alpha <- 1;
	}
	if (is.null(hparams$purity_beta)) {
		hparams$purity_beta <- 1;
	}
	if (is.null(hparams$ccf_alpha)) {
		hparams$ccf_alpha <- 1;
	}
	if (is.null(hparams$ccf_beta)) {
		hparams$ccf_beta <- 1;
	}
	if (is.null(hparams$somatic_prob)) {
		hparams$somatic_prob <- 0.5;
	}
	if (is.null(hparams$epsilon)) {
		hparams$epsilon <- 1e-3;
	}
	if (is.null(hparams$V)) {
		hparams$V <- 6;
	}
	if (is.null(hparams$N)) {
		hparams$N <- 1e3;
	}

	if (
		hparams$purity_alpha == 1 &&
		hparams$purity_beta  == 1 &&
		hparams$ccf_alpha    == 1 &&
		hparams$ccf_beta     == 1
	) {
		# simplified model with where \theta^{(1)}_j \sim Beta(1, 1)
		unif_theta1 <- TRUE;
	} else {
		unif_theta1 <- FALSE;
	}

	# number of germline classes
	K <- 1 + ncol(alpha)*2;

	# log prior on genotype
	lprior <- numeric(K + 1);
	# 1: somatic
	lprior[1] <- log(hparams$somatic_prob);
	# 2: germline homozygous
	# 3: germline
	# 4: germline with flipping
	# 3: germline with alternative hyperparameter
	# 4: germline with flipping alternative hyperparameter
	lprior[1 + (1:K)] <- log((1 - hparams$somatic_prob)/K);

	J <- length(ac);
	tc <- ac + bc;

	# log likelihood

	llike <- matrix(-Inf, nrow=J, ncol=K+1);

	# somatic class
	if (unif_theta1) {

		llike[, 1] <- -log(tc);

	} else {

		rho <- with(hparams, rbeta(N, purity_alpha, purity_beta));
		phi <- with(hparams, rbeta(N, ccf_alpha, ccf_beta));
		f.domain <- with(hparams, unique(unlist(lapply(1:V, function(m) (1:m)/m))));
		f <- sample(f.domain, replace=TRUE, hparams$N);
		theta <- get_theta(f, phi, rho, hparams$epsilon);

		# compute marginal likelihood by sampling
		llike[, 1] <- unlist(mapply(
			function(x, n) {
				-log(hparams$N) + log_sum_exp(dbinom(x, n, theta, log=TRUE))
			},
			bc, tc,
			SIMPLIFY=FALSE
		));

	}

	# germline homozygous class
	homo.beta <- (1 - hparams$epsilon) * tc;
	homo.alpha <- hparams$epsilon * tc;
	llike[, 2] <- -lbeta(homo.beta, homo.alpha) + lchoose(tc, bc) + 
		lbeta(homo.beta + bc, homo.alpha + ac);

	# calculate the remaining germline likelihoods
	# germline heterozygous class (non-flipped and flipped)
	for (k in 1:ncol(alpha)) {
		beta.k <- beta[, k];
		alpha.k <- alpha[, k];
		valid <- which(!is.na(beta.k) & !is.na(alpha.k));
		if (length(valid) > 0) {	
			beta.k <- beta.k[valid];
			alpha.k <- alpha.k[valid];
			tc2 <- tc[valid];
			bc2 <- bc[valid];
			ac2 <- ac[valid];
			llike[valid, 1 + 2*k] <- -lbeta(beta.k, alpha.k) + 
				lchoose(tc2, bc2) + 
				lbeta(beta.k + bc2, alpha.k + ac2);
			llike[valid, 2 + 2*k] <- -lbeta(beta.k, alpha.k) + 
				lchoose(tc2, ac2) + 
				lbeta(beta.k + ac2, alpha.k + bc2);
		}
	}

	# log joint
	lp <- t(lprior + t(llike));
	lpost <- numeric(J);

	# calculuate log p(genotype = "somatic" | data)
	for (j in 1:J) {
		lpost[j] <- lp[j, 1] - log_sum_exp(lp[j, ]);
	}

	if (log) {
		lpost
	} else {
		exp(lpost)
	}
}

get_theta <- function(f, phi, rho, e) {
	f2 <- f * phi * rho;
	(1 - f2) * e  +  f2 * (1 - e)
}

#' Extract Allelic Depth (count of reads supporting ref and alt allele)
#' @param vcf (data.frame) VCF read using qread from io package
#' @param sample_name (string) name of tumor sample name in the VCF
#' @param annotate_vaf (boolean) calcualtes variant allele frequency using AD
#' @return data.frame with referene and alternate read count extracted to column level
#' columns: [chrom, pos, ref, alt, ref_ad, alt_ad]
extract_ad <- function(vcf, sample_name, annotate_vaf=FALSE) {
	ad.list <- lapply(1:nrow(vcf),
		function(i) {
			ad <- vcf[[sample_name]][[i]]$AD;
			if (length(ad) > 2) {
				# first allele is reference; rest are alternative
				# choose alternative count that is the maximum
				c(ad[1], max(ad[-1]))
			} else {
				ad
			}
		}
	)

	# some allelic counts have 3 elements
	#idx <- which(unlist(lapply(ad.list, function(x) length(x) > 2)));

	ad <- matrix(unlist(ad.list), ncol = 2, byrow=TRUE);
	d <- cbind(vcf[, c("chrom", "pos", "ref", "alt")], ad_ref=ad[,1], ad_alt=ad[,2]);
	# d$chrom <- sub("chr", "", d$chrom, fixed=TRUE);
	if (annotate_vaf){
		d$vaf <- d$ad_alt / (d$ad_ref + d$ad_alt)
	}
	
	return(d)
}

#' Extract Allelic Depth (count of reads supporting ref and alt allele)
#' @param vcf (data.frame) VCF read using qread from io package.
#' @param sample_name (string) name of tumor sample name in the VCF.
#' @param tresh (float) posterior probability threshold above which 
#' the varaint is to be considerd somatic.
#' @return (data.frame) variants with macni posterior probability annotated.
#' columns:  [chrom, pos, ref, alt, ref_ad, alt_ad, "macni_pp", "is_somatic"]
run_macni <- function(vcf, sample_name, thresh=0.9, alpha.ghomo = 100, alpha.ghet = 100, beta.ghomo = 0, beta.ghet = 100) {

	data <- extract_ad(vcf, sample_name, annotate_vaf=TRUE)

	alpha <- data.frame(
		ghomo = rep(alpha.ghomo, times = nrow(data)), 
		ghet = rep(alpha.ghet, times = nrow(data))
	)

	beta <- data.frame(
		ghomo = rep(beta.ghomo, times = nrow(data)), 
		ghet = rep(beta.ghet, times = nrow(data))
	)

	data$macni_pp <- somatic_model(
		ac = data$ad_ref,
		bc = data$ad_alt,
		alpha = alpha,
		beta = beta,
		log = FALSE
	)

	data$is_somatic <- ifelse(data$macni_pp > thresh & !(is.nan(data$macni_pp) & data$vaf == 1), TRUE, FALSE)
	# data$chrom <- paste0("chr", data$chrom)

	return(data)
}
