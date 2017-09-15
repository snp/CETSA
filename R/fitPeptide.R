fitPeptide <- function(pepdata, startPars=c("Tm"=50, "Pl"=0, "b" = 0.05)){
  pepdata %>%
    group_by(Sample) %>%
    do(
      model=fitSigmoid(.[,c("Temperature","Value")],startPars),
      yVec = .$Value
    ) %>%
    filter(class(model)=='nls') %>%
    rowwise() %>%
    summarise(
      Sample=Sample,
      sigma=sigma(model),
      rSquared = rSquared(model, yVec),
      model=model) %>%
    arrange(sigma)
}
