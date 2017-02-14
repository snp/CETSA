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
analyzeTR <- function(filename = "proteinGroups.txt",
                      vehicle = c('C1'),
                      treatment = c('M2'),
                      resultPath = "CETSA_result",
                      plotCurves = TRUE,
                      idVar = "Majority protein IDs",
                      qPrefix = "Reporter intensity corrected",
                      temperatures = c(37, 41, 44, 47, 50, 53, 57, 61, 64, 67),
                      resultColumns = c("Protein names", "Gene names", "Unique peptides")) {
  if (!dir.exists(resultPath))
    dir.create(resultPath, recursive = TRUE)

  message("Reading data from ", filename)
  data <- importTR_MQ(
    proteinGroups = filename,
    idVar = idVar,
    qPrefix = qPrefix,
    temperatures = temperatures
  )

  save(data, file = file.path(resultPath, "data.RData"))
  #data <- importTR_MQ("../TPPQC/NTUB1P_MTX/peptides_M28.txt", idVar = "Sequence", qPrefix="Reporter intensity corrected")

  message("Normalizing data")
  normdata <-
    normalizeTR(data, limits = list("low" = c(42, 0.8), "high" = c(60, 0.2)))
  save(normdata, file = file.path(resultPath, "normdata.RData"))

  message("Fitting individual proteins")
  fitted <-
    fitPeptides(normdata, plotCurves = plotCurves, resultPath = resultPath)
  save(fitted, file = file.path(resultPath, "fitted.RData"))

  message("Summarizing results")
  result <- summarizeTR(fitted, vehicle=vehicle, treatment=treatment)

  result <- data[,c("id", resultColumns)] %>%
    distinct(id, .keep_all = TRUE) %>%
    full_join(result, by='id')

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
  render(system.file("Rmd/Report.Rmd", package =getPackageName()), envir= parent.frame(), output_file=file.path(resultPath,"CETSA_report.html"))
  message("Done!")
}
