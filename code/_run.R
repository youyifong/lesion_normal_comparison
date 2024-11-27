source("code/data.R")


source("code/compare_biopsies.R")
source("code/compare_intact_ulcerative.R")
source("code/compare_markers.R")
source("code/make_tables.R") # depends on output from compare_markers.R







################################################################################
# compare normal vs lesion tissues, one png per marker

# compare.regions=function(dat, legend, ...) {
#   hist(dat[[dist]][dat$lesion], border="red", col=NULL, xlab=dist %.% " (um)", ylab="Number of cells", breaks=20
#        , ylim=range(c(hist(dat[[dist]][dat$lesion], breaks=20, plot=F)$counts, hist(dat[[dist]][dat$normal], breaks=20, plot=F)$counts))
#        , ...)
#   hist(dat[[dist]][dat$normal], border="blue", col=NULL, breaks=20, add=T)
#   mylegend(x=3, legend=legend, lty=1, col=c("blue","red"), y.intersp=1.5, text.width=700)
#   
#   plot(density(dat[[dist]][dat$lesion]), col="red", xlab=dist %.% " (um)", ylab="Density"
#        , ylim=range(c(density(dat[[dist]][dat$lesion])$y, density(dat[[dist]][dat$normal])$y))
#        , ...)
#   lines(density(dat[[dist]][dat$in_Les_normal_side]), col="blue")
#   mylegend(x=3, legend=legend, lty=1, col=c("blue","red"), y.intersp=1, text.width=700)
# }


# for (s in c("DD","CF")) {
#   if (s=="DD") markers=c(        "CD3", "CD4", "CD8a",          "CD15", "CD45", "CD68", "CD69", "CD123", "CD279")
#   if (s=="CF") markers=c("CD1a", "CD3", "CD4", "CD8a", "CD11c", "CD15", "CD45",                 "CD123", "CD279")
#   for(m in markers) {
#     #s="CF"; m="CD1a"
#     myfigure(mfcol=c(1,2))
#     compare.regions(get.dat(s%.%"_Les_"%.%m%.%"_masks_d2b_regmem.csv", s), xlim=xlim, main=s%.%" Les, "%.%m, legend=region.labels)    
#     mydev.off(file="figures/"%.%s%.%"_dist2skin_compare_regions_"%.%m, ext="png")
#   }
# }