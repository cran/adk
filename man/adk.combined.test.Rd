\name{adk.combined.test}
\alias{adk.combined.test}
\title{
Combined Anderson-Darling K-Sample Tests
}
\description{
This function combines several independent Anderson-Darling k-sample tests
into one overall test of the hypothesis that the independent samples within 
each group come from a common unspecified distribution, while the common
distributions may vary from group to group. The k for each group of k
independent samples may also change from group to group. All samples from
all groups are independent and are assumed to come from continuous
distributions. Also provided is a version that adjusts for a moderate number 
of ties (due to rounding). 

NA values are removed and the user is alerted with the total NA count.
It is up to the user to judge whether the removal of NA's is appropriate.
}
\usage{
adk.combined.test(\dots)
}
\arguments{
  \item{\dots}{
Either several lists, say L.1, ... , L.M, where list L.i contains K.i sample vectors of respective sizes n.i[1], ..., n.i[K.i], where n.i[j] > 4 is recommended

or a list of such lists. 
}
}
\details{
If AD.i is the Anderson-Darling criterion for the i-th group of K.i samples, 
its standardized test statistic is T.i = (AD.i - mu.i)/sig.i, with mu.i and
sig.i representing mean and standard deviation of AD.i. This test statistic 
is used to test the hypothesis that the samples in the i-th group all come 
from the same but unspecified continuous distribution function F(x).

The combined Anderson-Darling criterion is AD.c=AD.1 + ... + AD.M and 
T.combined = (AD.c - mu.c)/sig.c is the standardized test statistic, 
where mu.c and sig.c represent the mean and standard deviation of AD.c. 
This test statistic is used to simultaneously test whether the samples 
in each group come from the same continuous distribution function F(x). 
However, the unspecified common distribution function F(x) may change 
from group to group. 
}
\value{
A list with components 
\item{M}{number of groups of samples being compared}
\item{n.samples}{list of M vectors, each vector giving the sample sizes for 
each group of samples being compared}
\item{nt}{vector of length M of total sample sizes involved in each of the 
M comparisons of K.i samples, respectively}
\item{n.ties}{vector giving the number of ties in each the M
comparison groups}
\item{adk.i}{(2*M) * 3 matrix containing the T.kN, P-value, 
and extrapolation for the M individual Anderson-Darling tests, 
not adjusted for ties and adjusted for ties}
\item{mu}{vector of means of the M AD statistics}
\item{sig}{vector of standard deviations of the M AD statistics}
\item{adk.c}{2*3 matrix containing T.combined, P-value, and extrapolation 
for the combined test not adjusted for ties and adjusted for ties}
\item{mu.c}{mean of the combined AD statistic}
\item{sig.c}{standard deviation of the combined AD statistic}
\item{warning}{logical indicator, warning = TRUE when at least one of 
the sample sizes is < 5.}
}
\references{
Scholz, F. W. and Stephens, M. A. (1987), K-sample Anderson-Darling Tests, 
\emph{Journal of the American Statistical Association}, 
\bold{Vol 82, No. 399}, 918--924. 
}
\author{
Fritz Scholz 
}
\note{
This test is useful in analyzing treatment effects in randomized 
(incomplete) block experiments and in examining performance 
equivalence of several laboratories when presented with different 
test materials for comparison.
}
\seealso{
\code{\link{adk.test}}
}
\examples{
## Create two lists of sample vectors.
x1 <- list( c(1, 3, 2, 5, 7), c(2, 8, 1, 6, 9, 4), c(12, 5, 7, 9, 11) )
x2 <- list( c(51, 43, 31, 53, 21, 75), c(23, 45, 61, 17, 60) )

## Run adk.combined.test.
adk.combined.out <- adk.combined.test(x1,x2)  
# or out <- adk.combined.test(list(x1,x2))

## Examine list objects in adk.combined.out.
names(adk.combined.out)

## Extract matrix components adk.i and adk.c
adk.combined.out$adk.i
adk.combined.out$adk.c


## Fully formatted output 
adk.combined.out
}

