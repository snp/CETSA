fitPeptide <- function(pepdata, startPars=c("Tm"=50, "Pl"=0, "b" = 0.05)){
  pepdata %>%
    group_by(Sample) %>%
    do(
      model=fitSigmoid(.[,c("Temperature","Value")],startPars),
      yVec = .$Value
    ) %>%
    filter(class(model)=='nls') %>%
    group_by(Sample) %>%
    summarize(
      sigma=sigma(model[[1]]),
      model=model,
      rSquared = rSquared(model[[1]], yVec[[1]])) %>%
    arrange(sigma)
}
