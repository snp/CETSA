rSquared <- function(model, y) {
  ## Determine R2 of a fitted model
  if (class(model)!="try-error"){
    ssTot <- sum((y - mean(y, na.rm=TRUE))^2, na.rm=TRUE)
    ssRes <- sum( (y - predict(model))^2 , na.rm=TRUE)
    r2 <- 1 - ssRes/ssTot
  } else {
    r2 <- NA_real_
  }
  return(r2)
}
