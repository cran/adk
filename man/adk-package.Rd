\name{adk-package}
\alias{adk-package}
\docType{package}
\title{
The Package adk Contains a K-Sample Anderson-Darling Test and its Combinations
}
\description{
The K-sample Anderson-Darling test \code{adk.test} is used to test the hypothesis that K samples of various sizes come from a common continuous distribution that is otherwise unspecified. It is a rank test and it is consistent against all alternatives. The combined version of the test 
\code{adk.combined.test} is used to test several such hypotheses at the same time and the common distribution may vary from hypothesis to hypothesis.
}
\details{
\tabular{ll}{
Package: \tab adk\cr
Type: \tab Package\cr
Version: \tab 1.0-1\cr
Date: \tab 2011-05-21\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
}
\author{
Fritz Scholz

Maintainer: \email{fscholz@uw.edu}
}
\references{
Scholz, F. W. and Stephens, M. A. (1987), K-sample Anderson-Darling Tests, \emph{Journal of the American Statistical Association}, 
\bold{Vol 82, No. 399}, 918--924. 
}

\examples{
## Example using adk.test
x <- list(c(1,3,2,5,7),c(2,8,1,6,9,4),c(12,5,7,9,11))
out <- adk.test(x) # or out <- adk.test(c(1,3,2,5,7),c(2,8,1,6,9,4),c(12,5,7,9,11))

## Example using adk.combined.test
x1 <- list( c(1, 3, 2, 5, 7), c(2, 8, 1, 6, 9, 4), c(12, 5, 7, 9, 11) )
x2 <- list( c(51, 43, 31, 53, 21, 75), c(23, 45, 61, 17, 60) )
adk.combined.out <- adk.combined.test(x1,x2)  # or out <- adk.combined.test(list(x1,x2))
}
