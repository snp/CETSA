fitSigmoid <- function(data, startPars=c("Tm"=50, "Pl"=0, "b" = 0.05),maxAttempts=100){
  fctStr <- "(1 - Pl) * 1 / (1+exp((x-Tm)/b/x)) + Pl"
  fitFct <- as.formula(paste("y ~", fctStr))
  xVec <- data[[1]]
  yVec <- data[[2]]
  varyPars <- 0

  attempts <- 0
  repeatLoop <- TRUE

  ## Check if number of non-missing values is sufficient
  ## (NLS can only handle data with at least three non-missing values)
  validValues <- !is.na(yVec)
  if (sum(validValues) <=2){
    m <- NA
    class(m) <- "try-error"
  } else{
    ## Perform fit
    while(repeatLoop & attempts < maxAttempts){
      parTmp <- startPars * (1 + varyPars*(runif(1, -0.2, 0.2)))
      m <- try(nls(formula=fitFct, start=parTmp, data=list(x=xVec, y=yVec), na.action=na.exclude,
                   algorithm="port", lower=c(30.0, 0, 1e-8), upper=c(90.0, 1, 1)),
               silent=TRUE)
      attempts <- attempts + 1
      varyPars <- 1
      if (class(m)!="try-error") {
        repeatLoop <- FALSE
      }
    }
  }
  return(m)
}
