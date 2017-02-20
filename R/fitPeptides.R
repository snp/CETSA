fitPeptides <-
  function(data,
           startPars = c("Tm" = 50, "Pl" = 0, "b" = 0.05),
           plotCurves = FALSE,
           resultPath = '.') {
    # cl <- create_cluster()
    # cluster_library(cl, "CETSA")
    # cluster_copy(cl, startPars)
    # cluster_copy(cl, plotCurves)
    # cluster_copy(cl, resultPath)

    data %>%
      group_by(id) %>%
      # partition(id, cluster=cl) %>%
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
              arrange(Sample, Temperature)
            if (!dir.exists(file.path(resultPath, "plots")))
              dir.create(file.path(resultPath, "plots"))
            gg1 <-
              pepdata_m %>% ggplot(aes(x = Temperature, color = Sample))  +
                              geom_line(aes(y=prValue)) +
                              geom_point(aes(y = Value)) +
                              labs(x = "Temperature",
                                   y = "Value",
                                   title = .$id) +
                              theme_minimal()
            gg1 <- gg1 + annotation_custom(
              tableGrob(
                result %>% mutate(
                  V = sprintf("%.2f", estimate),
                  se = sprintf("%.2f", std.error),
                  R2 = sprintf("%.2f", rSquared)
                ) %>% select(term, V, se),
                rows = result$Sample
              ),
              xmin = 45,
              xmax = 80,
              ymin = -0.5,
              ymax = 2
            )
            ggsave(file.path(resultPath, "plots", sprintf("fit_%s.pdf", .$id)), gg1, device =
                     'pdf')
          }
        }
        result
      }) -> fitted
    message("fitted")
    fitted %>% glimpse()
    message("collected")
    fitted %>% collect() %>% glimpse()
    fitted
  }
