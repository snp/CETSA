normalizeTR <-
  function(data, limits = list("low" = c(42, 0.8), "high" = c(60, 0.2))) {
    maxN = nrow(data %>% distinct(Sample,Temperature))
    data %>% group_by(id) %>% summarize(n=sum(!is.na(Value))) %>% mutate(all_points = (n>=maxN)) -> npoints
    lowPass <- data %>%
      filter(Temperature <= limits$low[1]) %>%
      group_by( id) %>%
      summarise(minVal = min(Value, na.rm=T))
    numSamples <- length(unique(data$Sample))
    highPass <- data %>%
      filter(Temperature >= limits$high[1]) %>%
      group_by(id) %>%
      summarise(maxVal = max(Value, na.rm=T))

    lowPass %>%
      full_join(highPass, by = c( "id")) %>%
      mutate(good = (minVal > limits$low[2]) &
               (maxVal < limits$high[2])) %>%
      left_join(npoints, by='id') %>%
      filter(all_points & good) -> good_proteins

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
    print(summary(models$model[[1]]))
    model_best <- models$model[[1]]
    message(sprintf("Normalization curve R-squared: %.2f", models$Rsq[[1]]))

    temps <- unique(data$Temperature)
    profile <-
      data.frame(Temperature = temps,
                 Predicted = predict(model_best, temps))

    data %>%
      filter(id %in% good_proteins$id) %>%
      group_by(Sample, Temperature) %>%
      summarize(total = mean(Value, na.rm=T)) %>%
      ungroup() %>%
      full_join(profile, by="Temperature") %>%
      mutate(norm_coef = Predicted / total) %>%
      select(Sample, Temperature, norm_coef) %>%
      full_join(data, by=c("Sample","Temperature")) %>%
      mutate(Value = Value * norm_coef) %>%
      mutate(normProtein = id %in% good_proteins$id)
  }

normalizeTR_notscaled <-
  function(data, limits = list("low" = c(42, 0.8), "high" = c(60, 0.2)), lowPl=42) {
    maxN = nrow(data %>% distinct(Sample,Temperature))
    data %>% group_by(id) %>% summarize(n=sum(!is.na(Value))) %>% mutate(all_points = (n>=maxN)) -> npoints
    lowT <- data %>%
      filter(Temperature<lowPl) %>%
      group_by(Sample, id) %>%
      summarize(lowT = mean(Value_)) %>% ungroup()
    lowPass <- data %>%
      select(-one_of("lowT")) %>%
      left_join(lowT) %>%
      mutate(Value__=Value_/lowT) %>%
      filter(Temperature <= limits$low[1]) %>%
      group_by(id) %>%
      summarise(minVal = min(Value__, na.rm=T)) %>% ungroup()
    numSamples <- length(unique(data$Sample))
    highPass <- data %>%
      select(-one_of("lowT")) %>%
      left_join(lowT) %>%
      mutate(Value__=Value_/lowT) %>%
      filter(Temperature >= limits$high[1]) %>%
      group_by(id) %>%
      summarise(maxVal = max(Value__, na.rm=T)) %>% ungroup()

    lowPass %>%
      full_join(highPass, by = c("id")) %>%
      mutate(good = (minVal > limits$low[2]) &
               (maxVal < limits$high[2])) %>%
      left_join(npoints, by='id') %>%
      filter(all_points & good) -> good_proteins

    message("Proteins for normalization: ", nrow(good_proteins))

    models <- data %>%
      filter(id %in% good_proteins$id) %>%
      group_by(Sample, Temperature) %>%
      summarize(total = mean(Value_, na.rm=T)) %>%
      do(model = fitSigmoid_notscaled(.[, 2:3]),
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
    print(summary(models$model[[1]]))
    model_best <- models$model[[1]]
    message(sprintf("Normalization curve R-squared: %.2f", models$Rsq[[1]]))

    temps <- unique(data$Temperature)
    profile <-
      data.frame(Temperature = temps,
                 Predicted = predict(model_best, temps))

    data_norm <- data %>%
      filter(id %in% good_proteins$id) %>%
      group_by(Sample, Temperature) %>%
      summarize(total = mean(Value_, na.rm=T)) %>%
      ungroup() %>%
      full_join(profile, by="Temperature") %>%
      mutate(norm_coef = Predicted / total) %>%
      select(Sample, Temperature, norm_coef) %>%
      full_join(data, by=c("Sample","Temperature")) %>%
      mutate(Value_ = Value_ * norm_coef) %>%
      mutate(normProtein = id %in% good_proteins$id)
    lowT <- data_norm %>%
      filter(Temperature<lowPl) %>%
      group_by(Sample, id) %>%
      summarize(lowT = mean(Value_)) %>% ungroup()
    data_norm <- data_norm %>%
      select(-one_of("lowT")) %>%
      full_join(lowT, by=c("Sample","id")) %>%
      filter(lowT > 0) %>%
      mutate(Value=Value_/lowT)
    data_norm
  }
