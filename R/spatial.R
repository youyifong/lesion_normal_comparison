# dist="dist2DEJ"
dist="dist2skin"

library(kyotil)

get.dat=function(filename, sample.label, dist) {
    dat=read.csv("data/"%.%filename)
    
    if (sample.label=="DD") {
      dat$normal.epi =    dat$in_DD_Les_epi_region_1_normal=="True" 
      dat$normal.dermis = dat$in_DD_Les_dermis_region_1_normal=="True" 
      dat$lesion.epi =    dat$in_DD_Les_epi_region_2_lesion=="True" 
      dat$lesion.dermis = dat$in_DD_Les_dermis_region_2_lesion=="True" 

      dat$epi = dat$normal.epi | dat$lesion.epi
      dat$dermis = dat$normal.dermis | dat$lesion.dermis
      dat$normal = dat$normal.epi | dat$normal.dermis
      dat$lesion = dat$lesion.epi | dat$lesion.dermis
      
      # each pixel is 1/2 micron, negative dist in epi
      dat$dist2skin = dat$dist2DD_Les_skin_boundary/2
      dat$dist2DEJ  = dat$dist2DD_Les_DEJviaDAPI_v2/2
      
      if (dist=="dist2DEJ") {
        dat$dist2skin=dat$dist2skin * ifelse(dat$epi, -1, 1)
        dat$dist2DEJ =dat$dist2DEJ  * ifelse(dat$epi, -1, 1)
      }
      
    } else if (sample.label=="CF") {
      dat$normal.epi =    dat$in_CF_Les_epi_region_1_normal=="True" | dat$in_CF_Les_epi_region_3_normal=="True" 
      dat$normal.dermis = dat$in_CF_Les_dermis_region_1_normal=="True" | dat$in_CF_Les_dermis_region_3_normal=="True" 
      dat$lesion.epi =    dat$in_CF_Les_epi_region_2_lesion=="True" 
      dat$lesion.dermis = dat$in_CF_Les_dermis_region_2_lesion=="True" 
      
      dat$epi = dat$normal.epi | dat$lesion.epi
      dat$dermis = dat$normal.dermis | dat$lesion.dermis
      dat$normal = dat$normal.epi | dat$normal.dermis
      dat$lesion = dat$lesion.epi | dat$lesion.dermis
      
      # each pixel is 1/2 micron, negative dist in epi
      dat$dist2skin = dat$dist2CF_Les_skin_boundary/2 * ifelse(dat$epi, -1, 1)
      dat$dist2DEJ  = dat$dist2CF_Les_DEJ_v2/2 * ifelse(dat$epi, -1, 1)
    }
    dat
}

compare.regions=function(dat, legend, ...) {
    hist(dat[[dist]][dat$lesion], border="red", col=NULL, xlab=dist %.% " (um)", ylab="Number of cells", breaks=20
        , ylim=range(c(hist(dat[[dist]][dat$lesion], breaks=20, plot=F)$counts, hist(dat[[dist]][dat$normal], breaks=20, plot=F)$counts))
        , ...)
    hist(dat[[dist]][dat$normal], border="blue", col=NULL, breaks=20, add=T)
    mylegend(x=3, legend=legend, lty=1, col=c("blue","red"), y.intersp=1.5, text.width=700)
    
    plot(density(dat[[dist]][dat$lesion]), col="red", xlab=dist %.% " (um)", ylab="Density"
        , ylim=range(c(density(dat[[dist]][dat$lesion])$y, density(dat[[dist]][dat$normal])$y))
        , ...)
    lines(density(dat[[dist]][dat$in_Les_normal_side]), col="blue")
    mylegend(x=3, legend=legend, lty=1, col=c("blue","red"), y.intersp=1, text.width=700)
}


compare.markers=function(dist.ls, legend, ...) {
    cnts.all = sapply(dist.ls, function (dist) hist(dist, breaks=20, plot=F)$counts)
    ylim = range(cnts.all)    
    
    hist(dist.ls[[1]], border=1, xlab=dist %.% " (um)", ylab="Number of cells", breaks=20, col=NULL, ylim=ylim, ...)
    for (i in 2:length(dist.ls)) hist(dist.ls[[i]], border=i, breaks=20, col=NULL, add=T)
    mylegend(x=3, legend=legend, lty=1, col=1:length(dist.ls), y.intersp=1, text.width=700)
    
    y.all = sapply(dist.ls, function (dist) density(dist)$y)
    ylim = range(y.all)    
    
    plot(density(dist.ls[[1]]), col=1, xlab=dist %.% " (um)", ylab="Density", ylim=ylim, ...)
    for (i in 2:length(dist.ls)) lines(density(dist.ls[[i]]), col=i)
    mylegend(x=3, legend=legend, lty=1, col=1:i, y.intersp=1, text.width=700)
}


region.labels=c(normal="intact epi", lesion="ulcerative")
region.labels.title=c(normal="Intact Epidermis", lesion="Ulcerative")

red=rgb(red=1, blue=0, green=0, alpha=1)
blue=rgb(red=0, blue=1, green=0, alpha=1)

xlim=c(-1000,2000)



#### compare markers ####

## DD
legend=c("CD3","CD15","CD123","CD11c","CD68")

dat.cd3=get.dat("DD_Les_CD3_masks_d2b_regmem.csv", "DD", dist)
dat.cd15=get.dat("DD_Les_CD15_masks_d2b_regmem.csv", "DD", dist)
dat.cd123=get.dat("DD_Les_CD123_masks_d2b_regmem.csv", "DD", dist)
dat.cd11c=get.dat("DD_Les_CD11c_masks_d2b_regmem.csv", "DD", dist)
dat.cd68=get.dat("DD_Les_CD68_masks_d2b_regmem.csv", "DD", dist)


dist.ls.normal = list(subset(dat.cd3, normal)[[dist]], 
                      subset(dat.cd15, normal)[[dist]], 
                      subset(dat.cd123, normal)[[dist]], 
                      subset(dat.cd11c, normal)[[dist]],
                      subset(dat.cd68, normal)[[dist]])
dist.ls.lesion = list(subset(dat.cd3, lesion)[[dist]], 
                      subset(dat.cd15, lesion)[[dist]], 
                      subset(dat.cd123, lesion)[[dist]], 
                      subset(dat.cd11c, lesion)[[dist]],
                      subset(dat.cd68, lesion)[[dist]])

mypng(mfcol=c(2,2), file="figures/DD_"%.%dist%.%"_compare_markers")
    compare.markers(dist.ls.normal, xlim=xlim, main="DD Les, "%.%region.labels.title["normal"], legend=legend)    
    compare.markers(dist.ls.lesion, xlim=xlim, main="DD Les, "%.%region.labels.title["lesion"], legend=legend)    
dev.off()


## CF
legend=c("CD3","CD15","CD123","CD11c") # CD68 did not pass QC

dat.cd3=get.dat("CF_Les_CD3_masks_d2b_regmem.csv", "CF", dist)
dat.cd15=get.dat("CF_Les_CD15_masks_d2b_regmem.csv", "CF", dist)
dat.cd123=get.dat("CF_Les_CD123_masks_d2b_regmem.csv", "CF", dist)
dat.cd11c=get.dat("CF_Les_CD11c_masks_d2b_regmem.csv", "CF", dist)

dist.ls.normal = list(subset(dat.cd3, normal)[[dist]], 
                      subset(dat.cd15, normal)[[dist]], 
                      subset(dat.cd123, normal)[[dist]], 
                      subset(dat.cd11c, normal)[[dist]])
dist.ls.lesion = list(subset(dat.cd3, lesion)[[dist]], 
                      subset(dat.cd15, lesion)[[dist]], 
                      subset(dat.cd123, lesion)[[dist]], 
                      subset(dat.cd11c, lesion)[[dist]])
    
mypng(mfcol=c(2,2), file="figures/CF_"%.%dist%.%"_compare_markers")
    compare.markers(dist.ls.normal, xlim=xlim, main="CF Les, "%.%region.labels.title["normal"], legend=legend)    
    compare.markers(dist.ls.lesion, xlim=xlim, main="CF Les, "%.%region.labels.title["lesion"], legend=legend)    
dev.off()


# 
# 
# ### compare normal vs lesion tissues ####
# # generate one png per marker
# 
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


