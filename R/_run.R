# can be run from either laptop or Rstudio on the Hutch server, but make_tables.R is too slow on laptop

# set working directory to the root dir of the repo
if (Sys.info()[["sysname"]] == "Windows") {
  setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/.."))
} 

source("R/data.R") # writes a Rdata file with the data

source("R/make_counts_tables.R") 

source("R/compare_intact_ulcerative.R")

source("R/compare_markers.R")

#source("R/compare_biopsies.R")

# correlates study is done in the Rmd

rmarkdown::render(
  "report_32bit_v2.Rmd",
  output_file = sprintf(
    "report_32bit_v2_%s.pdf",
    format(Sys.Date(), "%Y%m%d")
  )
)

# clinical figure is drawn with plot_shedding_lesion_history.R