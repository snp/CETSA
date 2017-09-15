fitPeptides <-
  function(data,
           startPars = c("Tm" = 50, "Pl" = 0, "b" = 0.05),
           plotCurves = FALSE,
           resultPath = '.', vehicle=c(), treatment=c()) {
    require(multidplyr)
    cl <- create_cluster()
    cluster_library(cl, "CETSA")
    cluster_library(cl, "ggplot2")
    cluster_library(cl, "gridExtra")
    cluster_library(cl, "broom")

    cluster_copy(cl, startPars)
    cluster_copy(cl, plotCurves)
    cluster_copy(cl, resultPath)
    cluster_copy(cl, fitPeptide)
    cluster_copy(cl, vehicle)
    cluster_copy(cl, treatment)
    set_default_cluster(cl)
    if (!dir.exists(file.path(resultPath, "plots")))
      dir.create(file.path(resultPath, "plots"), recursive = TRUE)

    data %>%
      # group_by(id) %>%
      partition(id) %>%
      do({
        pepdata <- .
        pid <- unique(pepdata$id)
        models <- fitPeptide(pepdata, startPars)
        result <- data_frame()
        if (nrow(models) > 1) {
          models %>%
            group_by(Sample) %>%
            do({
              m = .$model[[1]]
              res <- try(tidy(m), silent = TRUE)
              if (class(res) == 'try-error')
                res <- data.frame()
              else{
                res[, 'sigma'] = .$sigma
                res[, 'rSquared'] = .$rSquared
              }
              res
            }) %>%
            ungroup() %>%
            mutate(id = pid) -> result
          if (plotCurves & (nrow(result) > 1)) {
            temps <- unique(pepdata$Temperature)
            xtemps <- seq(min(temps), max(temps), length.out = 100)
            pepdata_m =  models %>%
              group_by(Sample) %>%
              do(data.frame(
                Sample = .$Sample,
                Temperature = xtemps,
                prValue = predict(.$model[[1]], list(x = xtemps))
              )) %>%
              full_join(pepdata, by = c("Sample", "Temperature")) %>%
              arrange(Sample, Temperature) %>%
              mutate(Type=Sample)

            if(length(treatment)>0){
              pepdata_m <- pepdata_m %>%
                mutate(Type = ifelse(Sample %in% treatment, "Treatment","Vehicle"))
            }
            gg1 <-
              pepdata_m %>% ggplot(aes(x = Temperature, color = Type, group=Sample))  +
                              geom_line(aes(y=prValue)) +
                              geom_point(aes(y = Value)) +
                              labs(x = "Temperature",
                                   y = "Value",
                                   title = .$id) +
                              theme_minimal()
            gg2 <- tableGrob(
                result %>% mutate(
                  V = sprintf("%.2f", estimate),
                  se = sprintf("%.2f", std.error),
                  R2 = sprintf("%.2f", rSquared)
                )  %>% filter(term=='Tm')%>% select(Sample, V, se) %>% unite(Val, V, se, sep="±") %>% spread(Sample, Val), theme=ttheme_minimal(base_size=6), rows="")
            m1 <- arrangeGrob(gg1, gg2,  heights=c(13,1))
            ggsave(file.path(resultPath, "plots", sprintf("fit_%s.pdf", .$id)), m1, device =
                     'pdf', width=6, height=4)
          }
        }
        result
      }) -> fitted
    message("fitted")
    fitted %>% collect()
  }
