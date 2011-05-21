adk.combined.test <-
function (...) 
{
# This function adk.combined.test combines several Anderson-Darling 
# K-sample test statistics AD.m, m = 1,...,M, into one overall test 
# statistic AD.combined as suggested in Section 8 of Scholz F.W. and 
# Stephens M.A. (1987), K-sample Anderson-Darling Tests,
# Journal of the American Statistical Association, Vol 82, No. 399, 
# pp. 918-924.
# See also the documentation of adk.test for the comparison of a single 
# set of K samples.
# Here each application of the Anderson-Darling K-sample test can be 
# based on a different K > 1. This combined version tests the hypothesis 
# that all the hypotheses underlying the individual K-sample tests are 
# true simultaneously.
# The individual K-sample test hypothesis is that all samples from 
# the m-th group come from a common population. However, that population 
# may be different from one individual K-sample situation to the next. 
# Such a combined test is useful in 
#
# 1) examining intra-laboratory measurement equivalence, when samples from 
#    the same material or batch are compared for M laboratories and such
#    comparisons are made for samples from several different materials or
#    batches and one assessment over all materials/batches is desired.
#
# 2) analyzing treatment effects in randomized complete or incomplete 
#    block designs.
#
# When there are NA's among the sample values they are removed,
# with a warning message indicating the number of NA's.
# It is up to the user to judge whether such removals make sense.
#
# Input: Either several lists, say L.1,...,L.M, where list L.i contains 
#        K.i sample vectors of respective sizes n.i[1], ..., n.i[K.i] 
#        (n.i[j] > 4 is recommended)
#
#        or a single list of such lists.
#        
#
# An example:
# x1 <- c(1, 3, 2, 5, 7), x2 <- c(2, 8, 1, 6, 9, 4), and 
# x3 <- c(12,  5,  7,  9, 11)
# y1 <- c(51, 43, 31, 53, 21, 75) and y2 <- c(23, 45, 61, 17, 60)
# then adk.combined.test(list(x1,x2,x3),list(y1,y2)) 
# or equivalently adk.combined.test(list(list(x1,x2,x3),list(y1,y2)))
# produces the outcome below.
##########################################################################
#Combination of Anderson-Darling K-Sample Tests.
#
# Number of data sets = 2 
#
# Sample sizes within each data set:
# Data set 1 :  5 6 5
# Data set 2 :  6 5
# Total sample size per data set: 16 11 
# Number of unique values per data set: 11 11 
#
# AD.i = Anderson-Darling Criterion for i-th data set
# Means: 2 1 
# Standard deviations: 0.92837 0.64816 
#
# T.i = (AD.i - mean.i)/sigma.i
#
#
# Null Hypothesis:
# All samples within a data set come from a common distribution.
# The common distribution may change between data sets.
#
#                        T.i P-value extrapolation
# not adj. for ties  1.41756 0.08956             0
# adj. for ties      1.62856 0.07084             0
#                                                
# not adj. for ties -0.96786 0.61831             1
# adj. for ties     -1.02925 0.63586             1
#                                                
#
# Combined Anderson-Darling Criterion: AD.combined = AD.1+AD.2 
# Mean = 3    Standard deviation = 1.13225 
#
# T.combined = (AD.combined - mean)/sigma
#
#                  T.combined P-value extrapolation
# not adj. for ties    0.60825 0.22302             0
# adj. for ties        0.74612 0.19293             0
#
###############################################################################
# For out.adk.combined <- adk.combined.test(list(x1,x2,x3),list(y1,y2))
# or out.adk.combined <- adk.combined.test(list(list(x1,x2,x3),list(y1,y2)))
# we get the object out.adk.combined of class adk with the following components
# > names(out)
# [1] "M"         "n.samples" "nt"        "n.ties"    "adk.i"     "mu"       
# [7] "sig"       "adk.c"     "mu.c"      "sig.c"     "warning"  
# where 
# M = number of sets of samples being compared
# n.samples = is a list of the vectors giving the sample sizes for each 
#             set of samples being compared
# nt = vector of total sample sizes involved in each of the M comparisons
# n.ties = vector of lenth M giving the number of ties in each comparison group
# adk.i = (2*M) * 3 matrix containing the ti.obs, P-value, and 
#            extrapolation for the M individual Anderson-Darling tests, 
#            not adjusted for ties and adjusted for ties. 
#            Here ti.obs is the observed value of T.i and
#            the corresponding P-value = P(T.i > ti.obs). 
#            Further, extrapolation = TRUE when the P-value is linearly
#            extrapolated outside of [.01,.25].
# mu = vector of means of the M AD statistics
# sig = vector of standard deviations of the M AD statistics
# adk.c = 2 * 3 matrix containing tc.obs, P-value, and extrapolation 
#            for the combined test, not adjusted for ties and adjusted for ties.
#            Here tc.obs is the observed value of Tc and 
#            P-value = P(Tc > tc.obs).
#            Further, extrapolation = TRUE when the P-value is 
#            linearly extrapolated outside of [.01,.25].
# mu.c = mean of the combined AD statistic
# sig.c = standard deviation of the combined AD statistic
# warning = logical indicator, warning = TRUE when at least one of 
#           the sample sizes is < 5.
#
# Fritz Scholz, May 2011
#####################################################################################
    na.remove <- function(x){
#
# This function removes NAs from a list and counts the total 
# number of NAs in na.total.
# Returned is a list with the cleaned list x.new and with 
# the count na.total of NAs.
#
	na.status <- sapply(x,is.na)
	k <- length(x)
	x.new <- list()
	na.total <- 0
	for( i in 1:k ){
		x.new[[i]] <- x[[i]][!na.status[[i]]]
		na.total <- na.total + sum(na.status[[i]])
	}
	list(x.new=x.new,na.total=na.total)
	}
   	if(nargs() == 1 & is.list(list(...)[[1]])) {
        data.sets <- list(...)[[1]]
   	}
   	else {
        data.sets <- list(...)
   	}
   	n.sizes <- NULL
   	M <- length(data.sets)
   	if(M < 2) 
        stop("To combine test results you must have at least two data sets.")
   	n.data <- sapply(data.sets, length)
   	if(any(n.data <= 1)) 
        stop("One or more of the data sets consists of less than 2 samples.")
   	n.samples <- list()
  	na.total <- 0
   	for(i in 1:M){
  	  	out <- na.remove(data.sets[i])
    	na.total <- na.total + out$na.total
        data.sets[i] <- out$x.new
		n.sample <- sapply(data.sets[[i]], length)
      	n.sizes <- c(n.sizes, n.sample)
      	if(any(n.sample==0))
      		stop(paste("one or more samples in data set", i,
                       " has no observations"))
      	n.samples[[i]] <- n.sample
   	}
	na.t <- na.total
   	if( na.t > 1) cat(paste("\n",na.t," NAs were removed!\n\n"))
   	if( na.t == 1) cat(paste("\n",na.t," NA was removed!\n\n"))
   	warning <- min(n.sizes) < 5
   	AD <- 0
   	sig <- NULL
   	n.ties <- NULL
   	nt <- NULL
   	mu <- NULL
   	adk.i <- NULL
   	mu.c <- 0
   	for(i in 1:M){
    	  out <- adk.test(data.sets[[i]])
    	  out.adk=out$adk
    	  dimnames(out.adk) <- list(c("not adj. for ties","adj. for ties"),
                               c("ti.obs","P-value","extrapolation"))
    	  adk.i <- rbind(adk.i, out.adk)
    	  sig.i <- out$sig
    	  mu <- c(mu, length(data.sets[[i]])-1)
    	  AD.i <- out$adk[,1]*sig.i + length(data.sets[[i]]) - 1
    	  sig <- c(sig, sig.i)
    	  AD <- AD+AD.i
    	  mu.c <- mu.c + length(data.sets[[i]]) - 1
    	  n.ties <- c(n.ties, out$n.ties)
    	  nt <- c(nt, sum(out$ns))
   	}
   	sig.c <- sqrt(sum(sig^2))
   	tc.obs <- (AD - mu.c)/sig.c
   	adk.pval1 <- adk.pval(tc.obs[1], mu.c)
   	adk.pval2 <- adk.pval(tc.obs[2], mu.c)
	adk.c <- matrix(c(signif(tc.obs[1],7),
           round(adk.pval1[[1]],7), adk.pval1[[2]],
           signif(tc.obs[2],7), round(adk.pval2[[1]],7),
           adk.pval2[[2]]), byrow=TRUE, ncol=3)
	dimnames(adk.c) <- list(c("not adj. for ties", "adj. for ties"),
                        c("tc.obs", "P-value", "extrapolation"))
	object <- list(M=M, n.samples=n.samples, nt=nt, n.ties=n.ties, adk.i=adk.i,
           mu=mu, sig=sig, adk.c = round(adk.c,5), mu.c=mu.c,
           sig.c=round(sig.c,5), warning=warning)
	class(object) <-  "adk"
	object
}

