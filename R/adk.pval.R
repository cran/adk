adk.pval <-
function (tx,m) 
{
# This function "adk.pval" evaluates the p-value of the observed value 
# tx of the T_m statistic in "K-Sample Anderson-Darling Tests" by F.W. Scholz 
# and M.A. Stephens (1987), Journal of the American Statistical Association, 
# Vol 82, No. 399, pp 918-924. Thus this p-value is P(T_m >= tx).
#
# Input: tx = observed value of T_m, tx > 0.
#         m = the index of T_m, m >= 1.
#
# Output: a list with components
#         p0 = p-value of tx, i.e., p0 = P(T_m >= tx)
#         extrap = a logical indicator
#                    extrap = TRUE means that linear extrapolation took place
#                    extrap = FALSE means that quadratic interpolation was used.
#
# Computational Details:
#
# This function first interpolates the upper T_m quantiles as given in Table 1
# of the above reference to the given value of m by fitting a quadratic 
# in 1/sqrt(m) to the quantiles as tabulated for the upper quantile
# levels .25, .10, .05, .025, .01.
#
# Next a quadratic in the interpolated quantiles (for m) is fitted to the  
# log-odds of the upper probability levels defining these quantiles
# and the fitted log-odds value at tx is converted back to the calculated upper 
# probability value, i.e., the p-value. p-values outside the tabulated range
# [.01,.25] are obtained by linear extrapolation of the fitted quadratic.
#  
# Fritz Scholz, May 2011
#=========================================================================================
table1.adk <- cbind(c(1, 2, 3, 4, 6, 8, 10, Inf), c(0.326, 
        0.449, 0.498, 0.525, 0.557, 0.576, 0.59, 0.674), c(1.225, 
        1.309, 1.324, 1.329, 1.332, 1.33, 1.329, 1.282), c(1.96, 
        1.945, 1.915, 1.894, 1.859, 1.839, 1.823, 1.645), c(2.719, 
        2.576, 2.493, 2.438, 2.365, 2.318, 2.284, 1.96), c(3.752, 
        3.414, 3.246, 3.139, 3.005, 2.92, 2.862, 2.326))
    extrap <- FALSE
    mt <- table1.adk[, 1]
    sqm1 <- 1/sqrt(mt)
    sqm2 <- sqm1^2
    tm <- NULL
    p <- c(0.25, 0.1, 0.05, 0.025, 0.01)
    lp <- log(p/(1 - p))
    for (i in 1:5) {
        out <- lsfit(cbind(sqm1, sqm2), table1.adk[, i + 1])
        x <- 1/sqrt(m)
        coef <- out$coef
        y <- coef[1] + coef[2] * x + coef[3] * x^2
        tm <- c(tm, y)
    }
    out <- lsfit(cbind(tm, tm^2), lp)
    coef <- out$coef
    if (tx <= max(tm) & tx >= min(tm)) {
        lp0 <- coef[1] + coef[2] * tx + coef[3] * tx^2
    }
    if (tx > max(tm)) {
        extrap <- TRUE
        lp0 <- min(lp) + (tx - max(tm)) * (coef[2] + 2 * coef[3] * 
            max(tm))
    }
    if (tx < min(tm)) {
        extrap <- TRUE
        lp0 <- max(lp) + (tx - min(tm)) * (coef[2] + 2 * coef[3] * 
            min(tm))
    }
    p0 <- exp(lp0)/(1 + exp(lp0))
    names(p0) <- NULL
    list(p0 = p0, extrap = extrap)
}

