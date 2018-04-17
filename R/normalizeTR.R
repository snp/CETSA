normalizeTR <-
  function(data, limits = list("low" = c(42, 0.8), "high" = c(60, 0.2))) {
    lowPass <- data %>%
      filter(Temperature <= limits$low[1]) %>%
      group_by(Sample, id) %>%
      summarise(minVal = min(Value))
    numSamples <- length(unique(data$Sample))
    highPass <- data %>%
      filter(Temperature >= limits$high[1]) %>%
      group_by(Sample, id) %>%
      summarise(maxVal = max(Value))

    lowPass %>%
      full_join(highPass, by = c("Sample", "id")) %>%
      mutate(good = (minVal > limits$low[2]) &
               (maxVal < limits$high[2])) %>%
      group_by(id) %>%
      summarize(score = sum(good)) %>%
      filter(score == numSamples) -> good_proteins

    message("Proteins for normalization: ", nrow(good_proteins))

    models <- data %>%
      filter(id %in% good_proteins$id) %>%
      group_by(Sample, Temperature) %>%
      summarize(total = mean(Value, na.rm=T)) %>%
      do(model = fitSigmoid(.[, 2:3]),
         yVec = .$total) %>%
      filter(class(model) == 'nls') %>%
#      group_by(Sample) %>%
      ungroup() %>%
      rowwise() %>%
      summarise(
        sigma = sigma(model),
        Rsq = rSquared(model, yVec),
        model = list(model)
      ) %>%
      arrange(desc(Rsq))

    model_best <- models$model[[1]]
    message(sprintf("Normalization curve R-squared: %.2f", models$Rsq[[1]]))

    temps <- unique(data$Temperature)
    profile <-
      data.frame(Temperature = temps,
                 Predicted = predict(model_best, temps))

    data %>%
      filter(id %in% good_proteins$id) %>%
      group_by(Sample, Temperature) %>%
      summarize(total = mean(Value)) %>%
      ungroup() %>%
      full_join(profile, by="Temperature") %>%
      mutate(norm_coef = Predicted / total) %>%
      select(Sample, Temperature, norm_coef) %>%
      full_join(data, by=c("Sample","Temperature")) %>%
      mutate(Value = Value * norm_coef) %>%
      mutate(normProtein = id %in% good_proteins$id)
  }
