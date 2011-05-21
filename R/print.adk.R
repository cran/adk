print.adk <-
function (x,...) 
{
######################################################
#
# This is a print function for objects of class adk,
# as they are produced by adk.test and adk.combined.test.
#
# Fritz Scholz, May 2011
#
#######################################################
    if(names(x)[1]=="k"){# checking whether the object x came from adk.test
    cat("Anderson-Darling k-sample test.\n")
    cat(paste("\nNumber of samples: ", x$k))
    cat("\nSample sizes: ")
    cat(x$ns)
    cat(paste("\nTotal number of values:", x$n))
    cat(paste("\nNumber of unique values:", x$n-x$n.ties))
    cat(paste("\n\nMean of Anderson-Darling Criterion:", 
        x$k-1))
    cat(paste("\nStandard deviation of Anderson-Darling Criterion:", 
        x$sig))
    cat("\n\nT.AD = (Anderson-Darling Criterion - mean)/sigma")
    cat("\n\nNull Hypothesis: All samples come from a common population.\n\n")
    print(x$adk)
    if (x$warning) {
        cat("\n\nWarning: At least one sample size is less than 5.\n")
        cat("   p-values may not be very accurate.\n")
    }
    invisible(x)
    }
  if(names(x)[1]=="M"){# checking whether the object x came from adk.combined.test
    cat("Combination of Anderson-Darling K-Sample Tests.\n")
    cat(paste("\nNumber of data sets =", x$M,"\n"))
    cat("\nSample sizes within each data set:\n")
    ns <- NULL
    k <- length(x$n.samples)
    d.sets <- paste("Data set",1:k)
    for(i in 1:k){
       cat(d.sets[i],": ",x$n.samples[[i]])
       cat("\n")
    }
    if(k>3) AD.name=paste("AD.1","...",paste("AD.",k,sep=""),sep="+")
    if(k==2)AD.name=paste("AD.1+AD.2")
    if(k==3)AD.name=paste("AD.1+AD.2+AD.3")
    cat("Total sample size per data set: ")
    cat(x$nt,"\n")
    cat("Number of unique values per data set: ")
    cat(x$nt-x$n.ties,"\n")
    cat("\nAD.i = Anderson-Darling Criterion for i-th data set\n")
    cat("Means:",x$mu,"\n")
    cat("Standard deviations:", x$sig,"\n")
    cat("\nT.i = (AD.i - mean.i)/sigma.i\n")
    cat("\nNull Hypothesis:\nAll samples within a data set come from a common distribution.\n")
    cat("The common distribution may change between data sets.\n")
    cat("\n")
    adk.mat=NULL
    nx <- nrow(x$adk.i)/2
    for(i in 1:nx){
    adk.mat <- rbind(adk.mat,x$adk.i[c(2*i-1,2*i),],c(NA,NA,NA))
    }
    print(adk.mat,na.print=" ")
    cat("Combined Anderson-Darling Criterion: AD.combined =",AD.name,"\n")
    cat("Mean =",x$mu.c,"   Standard deviation =",round(x$sig.c,5),"\n")
    cat("\nT.combined = (AD.combined - mean)/sigma\n")
    cat("\n")
    print(x$adk.c)
    if (x$warning) {
        cat("\n\nWarning: At least one sample size is less than 5.\n         p-values may not be very accurate.\n")
    }
    invisible(x)
  }
}

