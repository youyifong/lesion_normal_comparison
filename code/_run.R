# table 1 500um and 1000 um from DEJ in dermis
# add a third group with normalized density which is the average of area under density curve

# run from Rstudio on the Hutch server

source("code/data.R")

source("code/compare_biopsies.R")
source("code/compare_intact_ulcerative.R")
source("code/compare_markers.R")
source("code/make_tables.R") # depends on output from compare_markers.R
