`adk.test` <-
function (...) 
{
#################################################################################
# This function "adk.test" tests whether k samples (k>1) come from a common
# continuous distribution, using the nonparametric (rank) test described in
# Scholz F.W. and Stephens M.A. (1987), K-sample Anderson-Darling Tests,
# Journal of the American Statistical Association, Vol 82, No. 399, pp. 918-924. 
# This test is consistent against all alternatives. 
# There is an adjustment for ties, and for a moderate number of tied observations
# the P-value calculation should be reasonable.
# 
# The inputs can either be a sequence or a list of k (>1) sample vectors.
#
# An example:
# x1 <- c(1, 3, 2, 5, 7), x2 <- c(2, 8, 1, 6, 9, 4), 
# and x3 <- c(12,  5,  7,  9, 11)
# adk.test(x1,x2,x3) # or adk.test(list(x1,x2,x3)) produces the output below.
#################################################################################
# Anderson-Darling k-sample test.
#
# Number of samples:  3
# Sample sizes: 5 6 5
# Total number of values: 16
# Number of unique values: 11
#
# Mean of Anderson Darling Criterion: 2
# Standard deviation of Anderson Darling Criterion: 0.92837
#
# T = (Anderson Darling Criterion - mean)/sigma
#
# Null Hypothesis: All samples come from a common population.
#
#                         T P-value extrapolation
# not adj. for ties 1.41756 0.08956             0
# adj. for ties     1.62856 0.07084             0
#################################################################################
# In order to get the output list, call out.adk <- adk.test(x1,x2,x3)
# then out.adk is of class adk and has components 
# > names(out.adk)
# [1] "k"       "ns"      "n"       "n.ties"  "sig"    "adk"     "warning"
#
# where
# k = number of samples being compared
# ns = vector of the k sample sizes
# n = total sample size = n_1+...+n_k
# n.ties = number of ties in the combined set of all n observations
# sig = standard deviation of the AD statistic
# adk = 2 x 3 matrix containing t.obs, P-value, extrapolation, adjusting for
#       ties and not adjusting for ties. extrapolation is TRUE when 
#       the P-value was extrapolated. Here t.obs is the observed value of T
#       and P-value = P(T > t.obs).
# warning = logical indicator, warning = TRUE indicates that at least  
# one of the sample sizes is < 5.
#
# Fritz Scholz, January 2008
#
#################################################################################
    if (nargs() == 1 & is.list(list(...)[[1]])) {
        samples <- list(...)[[1]]
    }
    else {
        samples <- list(...)
    }
    k <- length(samples)
    if (k < 2) 
        stop("Must have at least two samples.")
    ns <- sapply(samples, length)
    if (any(ns == 0)) 
        stop("One or more samples have no observations.")
    x <- NULL
    for (i in 1:k) x <- c(x, samples[[i]])
    n <- length(x)
    Z.star <- sort(unique(x))
    L <- length(Z.star)
    AkN2 <- 0
    AakN2 <- 0
    l.vec <- NULL
    for (j in 1:L) {
        fij <- NULL
        for (i in 1:k) {
            fij <- c(fij, sum(samples[[i]] == Z.star[j]))
        }
        l.vec <- c(l.vec, sum(fij))
    }
    for (i in 1:k) {
        Mij <- 0
        Maij <- 0
        inner.sum <- 0
        inner.sum.a <- 0
        for (j in 1:L) {
            fij <- sum(samples[[i]] == Z.star[j])
            Mij <- Mij + fij
            Maij <- Mij - fij/2
            Bj <- sum(l.vec[1:j])
            Baj <- Bj - l.vec[j]/2
            if (j < L) {
                inner.sum <- inner.sum + l.vec[j] * (n * Mij - 
                  ns[i] * Bj)^2/(Bj * (n - Bj))
            }
            inner.sum.a <- inner.sum.a + l.vec[j] * (n * Maij - 
                ns[i] * Baj)^2/(Baj * (n - Baj) - n * l.vec[j]/4)
        }
        AkN2 <- AkN2 + inner.sum/ns[i]
        AakN2 <- AakN2 + inner.sum.a/ns[i]
    }
    AkN2 <- AkN2/n
    AakN2 <- (n - 1) * AakN2/n^2
    coef.d <- 0
    coef.c <- 0
    coef.b <- 0
    coef.a <- 0
    H <- sum(1/ns)
    h <- sum(1/(1:(n - 1)))
    g <- 0
    for (i in 1:(n - 2)) {
        g <- g + (1/(n - i)) * sum(1/((i + 1):(n - 1)))
    }
    coef.a <- (4 * g - 6) * (k - 1) + (10 - 6 * g) * H
    coef.b <- (2 * g - 4) * k^2 + 8 * h * k + (2 * g - 14 * h - 
        4) * H - 8 * h + 4 * g - 6
    coef.c <- (6 * h + 2 * g - 2) * k^2 + (4 * h - 4 * g + 6) * 
        k + (2 * h - 6) * H + 4 * h
    coef.d <- (2 * h + 6) * k^2 - 4 * h * k
    sig2 <- (coef.a * n^3 + coef.b * n^2 + coef.c * n + coef.d)/((n - 
        1) * (n - 2) * (n - 3))
    sig <- sqrt(sig2)
    TkN <- (AkN2 - (k - 1))/sig
    TakN <- (AakN2 - (k - 1))/sig
    pvalTkN <- adk.pval(TkN, k - 1)
    pvalTakN <- adk.pval(TakN, k - 1)
    warning <- min(ns) < 5
    ad.mat <- matrix(c(signif(TkN, 7), round(pvalTkN[[1]][1], 
        7), pvalTkN[[2]][1], signif(TakN, 7), round(pvalTakN[[1]][1], 
        7), pvalTakN[[2]][1]), byrow = T, ncol = 3)
    dimnames(ad.mat) <- list(c("not adj. for ties", "adj. for ties"), 
        c("t.obs", "P-value", "extrapolation"))
    object <- list(k = k, ns = ns, n = n, n.ties = n - L, sig = round(sig, 
        5), adk = round(ad.mat, 5), warning = warning)
    class(object) <- "adk"
    object
}
