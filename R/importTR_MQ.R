importTR_MQ <- function(proteinGroups="proteinGroups.txt",
                        idVar="Majority protein IDs",
                        qPrefix="Reporter intensity corrected",
                        temperatures=c(37,41,44,47,50,53,57,61,64,67)){

  # Read datafile
  data <- read_tsv(proteinGroups) %>%
    mutate_("id" = sprintf("`%s`", idVar)) %>%
    mutate(id = sub(";.*","", id)) %>%
    mutate(id = sub("[^a-zA-Z0-9-]*","", id)) %>%
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
  return(data %>% select(-contains("Reporter")) %>% select(-contains("peptides")))
}
