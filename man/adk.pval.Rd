\name{adk.pval}
\alias{adk.pval}
\title{
Function for computing p-values of Anderson-Darling 
test statistics 
}
\description{
This function "adk.pval" evaluates the p-value of the 
observed value tx of T_m, the standardized Anderson-Darling
statistic. It is used by \code{adk.test} and 
\code{adk.combined.test}.
}
\usage{
adk.pval(tx,m)
}
\arguments{
\item{tx}{threshold for which the right tail
probability of the standardized Anderson-Darling 
statistics T_m is to be calculated
}
\item{m}{
index of the standardized T_m statistic (see reference)
}
}
\details{
This function first interpolates the upper 
T_m quantiles as given in Table 1
(see reference below) to the given value of m by fitting a quadratic in 1/sqrt(m) to the quantiles as tabulated for the upper quantile levels .25, .10, .05, .025, .01.

Next a quadratic in the interpolated quantiles (for m) 
is fitted to the  
log-odds of the upper probability 
levels defining these quantiles
and the fitted log-odds value at tx is converted back to the calculated upper 
probability value, i.e., the p-value. p-values outside the tabulated range
[.01,.25] are obtained by linear extrapolation of the fitted quadratic.
}
\value{
A list with components 
 \item{p0}{p-value of tx, i.e., p0=P(T_m>= tx)}
 \item{extrap}{a logical indicator, \code{extrap = TRUE}
indicates that linear extrapolation took place,
otherwise quadratic interpolations was used.}
}
\references{
Scholz, F. W. and Stephens, M. A. (1987), K-sample Anderson-Darling Tests, \emph{Journal of the American Statistical Association}, 
\bold{Vol 82, No. 399}, 918--924. 
}
\author{
Fritz Scholz
}
\seealso{
\code{\link{adk.test}}  and 
\code{\link{adk.combined.test}}
}
\examples{
## compute p-value when tx = 3 and m = 6
adk.pval(3,6)
}

