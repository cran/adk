\name{adk.test}
\alias{adk.test}
\title{
Anderson-Darling K-Sample Test
}
\description{
The Anderson-Darling k-sample test may be used to test the hypothesis that k samples of various sizes ( > 4 ) come from one common continuous distribution. 
It is a rank test and it is consistent against all alternatives, a property not shared by the Kruskal-Wallis k-sample rank test. Also provided is a version 
that adjusts for a moderate number of ties (due to rounding). 

NA values are removed and the user is alerted with the total NA count.
It is up to the user to judge whether the removal of NA's is appropriate.
}
\usage{
adk.test(\dots)
}
\arguments{
\item{\dots}{Either several sample vectors of respective sizes n.1, ... , n.k, 
with n.i > 4 recommended,

or a list of such sample vectors 
}
}
\details{
See the given reference for details on the Anderson-Darling k-sample 
criterion AD and its modification in case of ties. The standardized value 
of AD, i.e., T = (AD - mu)/sig, is used as test statistic. Here mu = k-1 
and sig are the mean and standard deviation of AD. 
The P-value = P( T >= t.obs )
corresponding to an observed t.obs of T is computed by quadratic interpolation
w.r.t. 1/sqrt(mu) and by quadratic interpolation w.r.t. log(p/(1-p)), where p 
is the tail probability corresponding to the quantiles given in Table 1 of the
cited reference. Both interpolations are reasonably accurate. For p beyond the
range [.01,.25] of Table 1 linear exptrapolation is used w.r.t. the log(p/(1-p))
fit. Such extrapolation affects the accuracy of the P-value calculation to some
extent but this should not strongly affect any decisions regarding the tested
hypothesis. 
}
\value{
A list of class adk with components 
 \item{k}{number of samples being compared}
 \item{ns}{	vector of the k sample sizes c(n.1, ...,n.k)}
 \item{n}{total sample size = n.1 + ... + n.k}
 \item{n.ties}{number of ties in the combined set of all n observations}
 \item{sig}{standard deviation of the AD statistic}
 \item{adk}{2 x 3 matrix containing t.obs, P-value, extrapolation, not adjusting for ties and adjusting for ties. extrapolation = 1 when the P-value was extrapolated.}
 \item{warning}{logical variable, warning = TRUE if n.i < 5 for at least one of the samples, otherwise warning = FALSE .}
}
\references{
Scholz, F. W. and Stephens, M. A. (1987), K-sample Anderson-Darling Tests, \emph{Journal of the American Statistical Association}, 
\bold{Vol 82, No. 399}, 918--924. 
}
\author{
Fritz Scholz
}
\seealso{
\code{\link{kruskal.test}} as a nonparametric alternative to \code{adk.test} and 
\code{\link{adk.combined.test}} for combining several such tests for different and independent groups of samples
}
\examples{
## Create input list of 3 sample vectors.
x <- list(c(1,3,2,5,7),c(2,8,1,6,9,4),c(12,5,7,9,11))
out <- adk.test(x) 
# or out <- adk.test(c(1,3,2,5,7),c(2,8,1,6,9,4), c(12,5,7,9,11))
## Examine the component names of out
names(out)

## Examine the matrix adk of out.
out$adk

## Fully print formatted object out of class adk.
out
}

