#' @title Analyse CETSA experiment
#'
#' @description Reads MaxQuant output, normalizes data, fits melting curves and
#'     exports results
#'
#' @param filename Path to proteinGroups.txt or peptides.txt file
#' @param vehicle Array of names of vehicle experiments
#' @param treatment Array of names of treatment experiments
#' @param resultPath Path where the results will be exported
#' @param plotCurves If TRUE pdf plots for each protein will be created
#' @param idVar Column name to use for protein id
#' @param qPrefix Prefix for column with reporter intensities
#' @param temperatures array of temperatures for the respective TMT channels
#'
#' @export
#'
#' @examples
analyzeTR_pep <- function(filename = "peptides.txt",
                      vehicle = c('C1'),
                      treatment = c('M2'),
                      resultPath = "CETSA_result",
                      plotCurves = TRUE,
                      idVar = "Sequence",
                      protVar = "Leading razor protein",
                      qPrefix = "Reporter intensity corrected",
                      temperatures = c(37, 41, 44, 47, 50, 53, 57, 61, 64, 67),
                      resultColumns = c("Protein names", "Gene names"),
                      normLimits = list("low" = c(42, 0.7), "high" = c(60, 0.4))) {
  if (!dir.exists(resultPath))
    dir.create(resultPath, recursive = TRUE)

  message("Reading data from ", filename)
  data <- importTR_MQ(
    proteinGroups = filename,
    idVar = idVar,
    qPrefix = qPrefix,
    temperatures = temperatures
  ) %>% filter(Sample %in% c(vehicle, treatment)) %>%
    mutate_("pid" = sprintf("`%s`", protVar))
  save(data, file = file.path(resultPath, "data.RData"))
  #data <- importTR_MQ("../TPPQC/NTUB1P_MTX/peptides_M28.txt", idVar = "Sequence", qPrefix="Reporter intensity corrected")

  message("Normalizing data")
  normdata <-
    normalizeTR(data, limits = normLimits)
  save(normdata, file = file.path(resultPath, "normdata.RData"))

  message("Grouping peptides")
  normdata %>%
    #  group_by(pid) %>%
    group_by(pid, Sample) %>%
    do({
      res <- data.frame()
      pdata <- .
      if(nrow(pdata)>1){
        ppid = pdata$pid[1]
        s = pdata$Sample[1]
        x_ <- unique(pdata$Temperature)
        y_ <- pdata %>% select(id, Temperature, Value) %>% spread(Temperature, Value) %>% select(-id)
        x__ <- c()
        y__ <- data.frame(row.names = rownames(y_))
        for(fc in names(y_)){
          col_ <- y_[,fc]
          col_[col_<1e-4] <- NA
          col__ <- zoo::na.locf(col_)
          if(sum(is.na(col_))==0){
            y__[,fc] <- col__
            x__ <- c(x__, as.numeric(fc))
          }
          
        }
        y_ <-as.matrix(y__)
        
        
        y_ <- 2**y_
        farms.res <- try(generateExprVal.method.farms(y_,weighted.mean = T))

        if(class(farms.res) != "try-error")
          res <- data.frame(Sample=s, pid=ppid, Temperature=x__, Value=farms.res$exprs)
      }
      res
    }) -> protdata
  if(length(resultColumns) > 0){
    protdata %>%
      left_join(
        normdata %>%select(pid,one_of(resultColumns)),
        by="pid"
      )
  }
  save(protdata, file = file.path(resultPath, "proteins.RData"))
  protdata <- protdata %>%ungroup() %>%  mutate(id=pid)
  # normdata <- protdata

  message("Fitting individual proteins")
  fitted <-
    fitPeptides(protdata, plotCurves = plotCurves, resultPath = resultPath, vehicle=vehicle, treatment=treatment)
  save(fitted, file = file.path(resultPath, "fitted.RData"))

  message("Summarizing results")
  result <- summarizeTR(fitted, vehicle=vehicle, treatment=treatment) %>% ungroup()

  if(length(resultColumns) > 0){
    result <- result %>%
      mutate(pid=id) %>%
      left_join(data %>% select(pid, one_of(resultColumns)) %>% distinct(), by='pid')
  }
  message("Saving results to ", resultPath)
  result %>% write_csv(file.path(resultPath, "result_all.csv"))
  save(result, file = file.path(resultPath, "result.RData"))

  result %>%
    filter(N_vehicle == length(vehicle)) %>%
    filter(N_treatment == length(treatment)) %>%
    filter(Tm_se < Tm_vehicle) %>%
    filter(Tm_se < Tm_treatment) %>%
    filter(b_se < 2 * max(b_vehicle, b_treatment)) %>%
    filter(Pl_se < 0.5) %>%
    #    ggplot(aes(Tm_diff)) + geom_histogram(bins=200)
    write_csv(file.path(resultPath, "result_filtered.csv"))
  result %>%
    filter(Tm_se < Tm_vehicle) %>%
    filter(Tm_se < Tm_treatment) %>%
    filter((sigma_treatment + sigma_vehicle) < 0.9) %>%
    ggplot(aes(Tm_diff)) + geom_histogram(binwidth = 0.2)
  ggsave(file.path(resultPath, "Tm_diff_histogram.pdf"), device = 'pdf')
  render(system.file("Rmd/Report_pep.Rmd", package =getPackageName()), envir= sys.frame(sys.nframe()), output_file=file.path(resultPath,"CETSA_report.html"))
  message("Done!")
}
