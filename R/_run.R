# can be run from either laptop or Rstudio on the Hutch server, but make_tables.R is too slow on laptop

# set working directory to the root dir of the repo
{
if (Sys.info()[["sysname"]] == "Windows") setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/.."))

dir.create("reports", showWarnings = FALSE, recursive=TRUE)
dir.create("output", showWarnings = FALSE, recursive=TRUE)

source("R/data.R") # writes a Rdata file with the data
}


# density curves
source("R/compare_intact_ulcerative.R")
source("R/compare_markers.R")
#source("R/compare_biopsies.R")
rmarkdown::render("report_32bit_v2_density_curves.Rmd", output_file = sprintf("reports/report_32bit_v2_density_curves_%s.pdf", format(Sys.Date(), "%Y%m%d")))


source("R/make_counts_tables.R") 
rmarkdown::render("report_32bit_v2_density_tables.Rmd", output_file = sprintf("reports/report_32bit_v2_density_tables_%s.pdf", format(Sys.Date(), "%Y%m%d")))


# correlation study
rmarkdown::render("report_32bit_v2_cor.Rmd", output_file = sprintf("reports/report_32bit_v2_cor_%s.pdf", format(Sys.Date(), "%Y%m%d")))


# clinical 
source("R/plot_shedding_lesion_history.R")


# for manus
rmarkdown::render("count_size_cor_epi_lesion.Rmd", output_file = sprintf("reports/count_size_cor_epi_lesion_%s.pdf", format(Sys.Date(), "%Y%m%d")))
