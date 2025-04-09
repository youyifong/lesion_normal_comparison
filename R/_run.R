# can be run from either laptop or Rstudio on the Hutch server
# make_tables.R runs too slow on laptop and better be run on the Hutch server

# check if os is windows or linux and set workig dir accordingly
if (Sys.info()[["sysname"]] == "Windows") {
  # set working directory to the parent of the location of this file
  setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/.."))
} else {
  # automatically set working directory to the parent of the location of this file
}


source("R/data.R")

source("R/make_tables.R") 

#source("R/compare_biopsies.R")

source("R/compare_intact_ulcerative.R")

source("R/compare_markers.R")
