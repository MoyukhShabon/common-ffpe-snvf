# Strand bias model
# https://github.com/broadinstitute/gatk/blob/master/docs/mutect/mutect.pdf

library(matrixStats)

# log pmdf for beta-binomial distribution
ldbetabinom <- function(x, n, a, b) {
	lchoose(n, x) + lbeta(x + a, n - x + b) - lbeta(a, b);
}

# posterior probability of strand bias
# notations are modified from pdf
# z \in {0, 1, 2} s.t.
#		no artifact: z = 0
#		artifact on forward strand: z = 1
#   artifact on reverse strand: z = 2
# eps: strand-biased error rate
# theta: non-strand-baised error rate
# f: true allele fraction
post_strand_bias <- function(n_ref_f, n_ref_r, n_alt_f, n_alt_r,
	p_artifact = 0.5,
	a_eps=1, b_eps=100,
	# low read error rate on Illumina (maybe higher for other platforms)
	a_theta=1, b_theta=1000,
	a_f=1, b_f=1,
	log=FALSE
) {

	n_f <- n_ref_f + n_alt_f;
	n_r <- n_ref_r + n_alt_r;
	n <- n_f + n_r;
	x_f <- n_alt_f;
	x_r <- n_alt_r;
	x <- x_f + x_r;

	# log likelihoods

	llike <- c(
		# no artifact
		lchoose(n_f, x_f) + lchoose(n_r, x_r) - lchoose(n, x) + 
			ldbetabinom(x, n, a_f, b_f),
		# artifact on forward strand
		ldbetabinom(x_f, n_f, a_eps, b_eps) + 
			ldbetabinom(x_r, n_r, a_theta, b_theta),
		# artifact on reverse strand
		ldbetabinom(x_f, n_f, a_theta, b_theta) + 
			ldbetabinom(x_r, n_r, a_eps, b_eps)
	);

	llike + log(1 - p_artifact)

	lprior <- log(c(1 - p_artifact, p_artifact/2, p_artifact/2));

	lp <- llike + lprior;
	lp <- lp - logSumExp(lp);
	p_sb <- 1 - exp(lp[1]);

	if (log) log(p_sb) else p_sb
}

# Remarks  MuTect2 optimizes pi = p_artifact and posterior probability of z
#          using the EM algorithm. This doesn't make sense because p_artifact
#          is a parameter on the prior of z.


# testing

# no strand bias
# post_strand_bias(100, 100, 10, 10)
# post_strand_bias(90, 100, 90, 100)
# post_strand_bias(90, 100, 9, 10)

# # strand bias on WGS
# post_strand_bias(13, 31, 4, 0)
# post_strand_bias(8, 33, 3, 0)

# # strand bias on WES
# post_strand_bias(86, 108, 0, 0)
# post_strand_bias(19, 37, 1, 0)
# post_strand_bias(824, 475, 9, 36)