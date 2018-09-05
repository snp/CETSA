importTR_MQ <- function(proteinGroups="proteinGroups.txt",
                        idVar="Majority protein IDs",
                        qPrefix="Reporter intensity corrected",
                        temperatures=c(37,41,44,47,50,53,56,59,63,67),
                        lowPl = 44){

  # Read datafile
  data <- read_tsv(proteinGroups) %>%
    mutate_("id" = sprintf("`%s`", idVar)) %>%
    mutate(id = sub(";.*","", id)) %>%
    mutate(id = gsub("[^a-zA-Z0-9-]+","_", id)) %>%
    gather(Sample, Value, matches(sprintf("^%s (\\d) (.+)$", qPrefix))) %>%
    mutate(
      Channel = as.numeric(sub(sprintf("^%s (\\d) (.+)$", qPrefix), "\\1", Sample)),
      Sample = sub(sprintf("^%s (\\d) (.+)$", qPrefix), "\\2", Sample)) %>%
    mutate(Temperature=temperatures[Channel+1])
  # Calculate ratios from reporter intensities
  data %>%
    filter(Temperature<lowPl) %>%
    group_by(Sample, id) %>%
    summarize(lowT = mean(Value)) -> lowT
  data <- data %>%
      full_join(lowT, by=c("Sample","id")) %>%
      filter(lowT > 0) %>%
      mutate(Value_ = Value,
             Value=Value_/lowT)
  return(data)
}

importTR_PD <- function(filename="proteinGroups.txt",
                        idVar="Accession",
                        qPrefix="Abundance:",
                        resultColumns=c("Description"),
                        temperatures=c(37,41,44,47,50,53,56,59,63,67),
                        lowPl = 44){

  # Read datafile
  if (!file.exists( filename )){
    warning(sprintf("Can't find input file: [%s]", filename))
    return(data.frame())
  }
  TMT_channels_ <- c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131")
  message("Reading from ", filename)
  data <- read_tsv(filename) %>%
    mutate_("id" = sprintf("`%s`", idVar)) %>%
    mutate(id = sub(";.*","", id)) %>%
    mutate(id = gsub("[^a-zA-Z0-9-]+","_", id))
  if("Master" %in% names(data)) {
    data <- data %>% filter(Master == 'IsMasterProtein')
  }
  if("Contaminant" %in% names(data)){
    data <- data %>% filter(Contaminant == 'False')
  }
  if("Description" %in% names(data)) {
    data <- data %>%
      mutate(Gene = sub(".*GN=([A-Z0-9a-z]+)","\\1", Description)) %>%
      mutate(Gene = sub("\\s+.*$","", Gene))

  }
  data <- data %>%
    select(id, one_of(c("Description", "Gene", resultColumns)), matches(sprintf("^%s? (F\\d+)[\\:,]? (\\d+[NC]?).+?", qPrefix))) %>%
    gather(Sample, Value, matches(sprintf("^%s? (F\\d+)[\\:,]? (\\d+[NC]?).+?", qPrefix))) %>%
    # tidyr::gather(Sample, Value, matches(sprintf("^%s? (F\\d+)[\\:,]? (\\d+[NC]?).+?", qPrefix))) %>%
    filter(!is.na(Value)) %>%
    tidyr::extract(
      Sample,
      c("File", "Channel"),
      sprintf("^%s? (F\\d+)[\\:,]? (\\d+[NC]?).+?$", qPrefix)) %>%
    mutate(Channel = as.numeric(factor(Channel, levels=TMT_channels_))) %>%
    mutate(Sample = File) %>%
    mutate(Temperature=temperatures[Channel])
  # Calculate ratios from reporter intensities
  data %>%
    filter(Temperature<lowPl) %>%
    group_by(Sample, id) %>%
    summarize(lowT = mean(Value)) %>% ungroup() -> lowT
  data <- data %>%
    left_join(lowT, by=c("Sample","id")) %>%
    filter(lowT > 0) %>%
    mutate(Value_ = Value,
           Value=Value_/lowT)
  return(data)
}
