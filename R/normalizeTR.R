normalizeTR <- function(data, limits=list("low"=c(42,0.8), "high"=c(60,0.2))){
  lowPass <- data %>%
    filter(Temperature <= limits$low[1]) %>%
    group_by(Sample, id) %>%
    summarise(minVal = min(Value))
#    filter(minVal > limits$low[2])

  highPass <- data %>%
    filter(Temperature >= limits$high[1]) %>%
    group_by(Sample, id) %>%
    summarise(maxVal = max(Value))
#    filter(maxVal < limits$high[2])

  lowPass %>%
    full_join(highPass) %>%
    mutate(good=(minVal > limits$low[2]) & (maxVal < limits$high[2])) %>%
    group_by(id) %>%
    summarize(score=sum(good)) %>%
    filter(score==2) -> good_proteins


  models <- data %>%
    full_join(inner_join(lowPass, highPass)) %>%
    filter(minVal > limits$low[2], maxVal < limits$high[2]) %>%
    filter(id %in% good_proteins$id) %>%
    group_by(Sample, Temperature) %>%
    summarize(total=mean(Value)) %>%
    do(
      model=fitSigmoid(.[,2:3])
    ) %>%
    filter(class(model)=='nls') %>%
    group_by(Sample) %>%
    summarize(sigma=sigma(model[[1]]), model=model) %>%
    arrange(sigma)

  model_best <- models$model[[1]]
  temps <- unique(data$Temperature)
  profile <- data.frame(Temperature = temps, Predicted=predict(model_best,temps))

  data %>%
    full_join(inner_join(lowPass, highPass)) %>%
    filter(minVal > limits$low[2], maxVal < limits$high[2]) %>%
    group_by(Sample, Temperature) %>%
    summarize(total=mean(Value)) %>%
    ungroup() %>%
    full_join(profile) %>%
    mutate(norm_coef = Predicted/total) %>%
    select(Sample, Temperature, norm_coef) %>%
    full_join(data) %>%
    mutate(Value=Value*norm_coef)
}
