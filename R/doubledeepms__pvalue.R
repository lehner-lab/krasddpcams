doubledeepms__pvalue <- function(
  av, 
  se, 
  degreesFreedom = 5, 
  mu = 0,
  testType = "ztest"){
  
  #Perform z test
  if(testType=="ztest"){
    zscore <- (av - mu)/se
    pval <- 2*pnorm(abs(zscore), lower.tail = FALSE)
    return(pval)
  }
  
  #Perform t test
  if(testType=="ttest"){
    tstat <- (av - mu)/se
    pval <- 2*pt(abs(tstat), degreesFreedom, lower = FALSE)
    return(pval)
  }
  
}