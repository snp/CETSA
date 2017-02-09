filename <- "../TPPQC/NTUB1P_MTX/proteinGroups_M28.txt"
analyzeTR <- function(filename, vehicle=c('C1'), treatment=c('M2')){
  data <- importTR_MQ(filename)
  #data <- importTR_MQ("../TPPQC/NTUB1P_MTX/peptides_M28.txt", idVar = "Sequence", qPrefix="Reporter intensity corrected")
  normdata <- normalizeTR(data, limits=list("low"=c(42,0.8), "high"=c(62,0.2)))
  fitted <- fitPeptides(normdata, plotCurves = FALSE)
  fitted %>%
    gather(Measure, Value, estimate:p.value) %>%
    unite(temp, term, Measure) %>%
    spread(temp,Value) -> gathered
  gathered %>%
    group_by(id) %>%
    do({
      pid <- unique(.$id)
      pdata <- .
      vdata <- pdata %>% filter(Sample %in% vehicle)
      tdata <- pdata %>% filter(Sample %in% treatment)
      ret <- data.frame()
      if((nrow(vdata)>0) & (nrow(tdata)>0))
        ret <- data.frame(
          id = pid,
          Tm_vehicle = mean(vdata$Tm_estimate),
          Tm_treatment = mean(tdata$Tm_estimate),
          Tm_diff = mean(tdata$Tm_estimate)-mean(vdata$Tm_estimate),
          Tm_se = max(c(vdata$Tm_std.error, tdata$Tm_std.error)),
          Pl_vehicle = mean(vdata$Pl_estimate),
          Pl_treatment = mean(tdata$Pl_estimate),
          Pl_diff = mean(tdata$Pl_estimate)-mean(vdata$Pl_estimate),
          Pl_se = max(c(vdata$Pl_std.error, tdata$Pl_std.error)),
          b_vehicle = mean(vdata$b_estimate),
          b_treatment = mean(tdata$b_estimate),
          b_se = max(c(vdata$b_std.error, tdata$b_std.error))
        )
      ret
    }) -> result
  result %>%
    filter(Tm_se < Tm_vehicle) %>%
    filter(Tm_se < Tm_treatment) %>%
    filter(b_se < 2*max(b_vehicle,b_treatment)) %>%
    filter(Pl_se < 0.5) %>%
    write_csv("../TPPQC/NTUB1P_MTX/result.csv")

}

