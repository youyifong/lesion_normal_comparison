# can be run from either laptop or Rstudio on the Hutch server
# it runs faster on the Hutch server
# check if os is windows or linux
if (Sys.info()[["sysname"]] == "Windows") {
  # set working directory to the parent of the location of this file
  setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/.."))
} else {
}


source("code/data.R")

source("code/compare_biopsies.R")
source("code/compare_intact_ulcerative.R")
source("code/compare_markers.R")
source("code/make_tables.R") # dist.bt.peaks is computed in compare_markers.R
