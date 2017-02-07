importTR_MQ <- function(proteinGroups="proteinGroups.txt",
                        idVar="Majority protein IDs",
                        controls,
                        treatments,
                        qPrefix="Reporter intensity corrected",
                        temperatures=c(37,41,44,47,50,53,57,61,64,67)){

  get_temp <- function(x, t=temperatures){sapply(x, function(xx){t[xx]})}
  # Read datafile
  data <- read_tsv(proteinGroups) %>%
    gather(Sample, Value, matches(sprintf("^%s (\\d) (.+)$", qPrefix))) %>%
    mutate(
      Channel = as.numeric(sub(sprintf("^%s (\\d) (.+)$", qPrefix), "\\1", Sample)),
      Sample = sub(sprintf("^%s (\\d) (.+)$", qPrefix), "\\2", Sample)) %>%
    mutate(Temperature=temperatures[Channel+1])
  # Calculate ratios from reporter intensities
  data %>%
    filter(Temperature<42) %>%
    group_by_("Sample", sprintf("`%s`", idVar)) %>%
    summarize(lowT = mean(Value)) -> lowT
  data <- data %>%
    full_join(lowT) %>%
    filter(lowT > 0) %>%
    mutate(Value=Value/lowT)
  # data %>%
  #   group_by(Sample, Temperature) %>%
  #   summarize(total=median(Value)) %>%
  #   ggplot(aes(x=factor(Temperature),y=total, fill=Sample)) + geom_bar(stat='identity', position='dodge')
  return(data %>% select(-contains("Reporter")) %>% select(-contains("peptides")))
}
