# can be run from either laptop or Rstudio on the Hutch server, but make_tables.R is too slow on laptop

if (Sys.info()[["sysname"]] == "Windows") {
  # set working directory to the root dir of the repo
  setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/.."))
} 

# this write a Rdata file with the data
source("R/data.R")

source("R/make_counts_tables.R") 

source("R/compare_intact_ulcerative.R")

source("R/compare_markers.R")

#source("R/compare_biopsies.R")

rmarkdown::render("report_32bit_v2.Rmd")

