summarizeTR <- function(fitResult,
                        vehicle = c('C1'),
                        treatment = c('M2')){
  fitResult %>%
    gather(Measure, Value, estimate:p.value) %>%
    unite(temp, term, Measure) %>%
    spread(temp, Value) -> gathered
  gathered %>%
    gather(Measure, Value, -id, -Sample) %>%
    unite(temp, Sample, Measure) %>%
    spread(temp, Value) -> regathered
  gathered %>%
    group_by(id) %>%
    do({
      pid <- unique(.$id)
      pdata <- .
      vdata <- pdata %>% filter(Sample %in% vehicle)
      tdata <- pdata %>% filter(Sample %in% treatment)
      ret <- data.frame()
      if ((nrow(vdata) > 0) & (nrow(tdata) > 0)){
        df_ <- 9 * sum(pdata$Tm_std.error**2, na.rm = T)**2 / sum(pdata$Tm_std.error**4, na.rm = T)
        se_ <- sqrt(sum(pdata$Tm_std.error**2))
        t_ <- abs(mean(tdata$Tm_estimate) - mean(vdata$Tm_estimate))/se_
        p_ <- pt(t_, df=df_,lower.tail=F)

        pv_ = 1
        pt_ = 1
        pn_ = 1
        if(nrow(vdata)==2){
          se_ <- sqrt(mean(vdata$Tm_std.error**2))
          df_ <- 9 * sum(vdata$Tm_std.error**2, na.rm = T)**2 / sum(vdata$Tm_std.error**4, na.rm = T)
          t_ <- abs(vdata$Tm_estimate[2] - vdata$Tm_estimate[1])/se_
          pv_ <- pt(t_, df=df_, lower.tail=F)
        }
        if(nrow(tdata)==2){
          se_ <- sqrt(mean(tdata$Tm_std.error**2))
          df_ <- 9 * sum(tdata$Tm_std.error**2, na.rm = T)**2 / sum(tdata$Tm_std.error**4, na.rm = T)
          t_ <- abs(tdata$Tm_estimate[2] - tdata$Tm_estimate[1])/se_
          pt_ <- pt(t_, df=df_, lower.tail=F)
        }
        if(nrow(vdata)>1 & nrow(tdata)>1){
          tt <- try(t.test(vdata$Tm_estimate, tdata$Tm_estimate), silent=T)
          if (class(tt)!="try-error")
            pn_ <- tt$p.value
        }


        ret <- data.frame(
          id = pid,
          sigma_vehicle = mean(vdata$sigma),
          sigma_treatment = mean(tdata$sigma),
          N_vehicle = nrow(vdata),
          N_treatment = nrow(tdata),
          Tm_vehicle = mean(vdata$Tm_estimate),
          Tm_vehicle_diff = ifelse(nrow(vdata)>1, vdata$Tm_estimate[-1] - vdata$Tm_estimate[1], NA),
          Tm_vehicle_se = min(vdata$Tm_std.error),
          Tm_vehicle_pval = pv_,
          Tm_treatment = mean(tdata$Tm_estimate),
          Tm_treatment_diff = ifelse(nrow(tdata)>1,tdata$Tm_estimate[-1] - tdata$Tm_estimate[1],NA),
          Tm_treatment_se = min(tdata$Tm_std.error),
          Tm_treatment_pval = pt_,
          Tm_diff = mean(tdata$Tm_estimate) - mean(vdata$Tm_estimate),
          Tm_pval = p_,
          Tm_pval_naive = pn_,
          Tm_se = max(c(
            vdata$Tm_std.error, tdata$Tm_std.error
          )),
          Pl_vehicle = mean(vdata$Pl_estimate),
          Pl_treatment = mean(tdata$Pl_estimate),
          Pl_diff = mean(tdata$Pl_estimate) - mean(vdata$Pl_estimate),
          Pl_se = max(c(
            vdata$Pl_std.error, tdata$Pl_std.error
          )),
          b_vehicle = mean(vdata$b_estimate),
          b_treatment = mean(tdata$b_estimate),
          b_se = max(c(
            vdata$b_std.error, tdata$b_std.error
          ))
        )
      }
      ret
    }) -> result
  result[,"Tm_padj"] <- p.adjust(result$Tm_pval, method='BH')
  result[,"Tm_padj_naive"] <- p.adjust(result$Tm_pval_naive, method='BH')
  result %>%
    full_join(regathered, by="id")
}
