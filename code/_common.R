library(kyotil)

region.labels=c(normal="intact", lesion="ulcerative")
region.labels.title=c(normal="Intact", lesion="Ulcerative")

red=rgb(red=1, blue=0, green=0, alpha=1)
blue=rgb(red=0, blue=1, green=0, alpha=1)

sample.labels=c("DD","CF","AS","CJ","JM","LG","LS"); names(sample.labels)=sample.labels
dist.labels=c("dist2skin", "dist2DEJ"); names(dist.labels)=dist.labels
tissue.types=c("lesion","normal"); names(tissue.types)=tissue.types

markers=c("CD3","CD15","CD123","CD11c"); names(markers)=markers
# CD68 did not pass QC for CF
markers.ls=list(DD=c(markers,CD68="CD68"), CF=markers)

xlim.ls=list(dist2DEJ=c(-750,1750), dist2skin=c(0,2500))

# computed from analysis.py
tissue.widths=list(
  lesion=c(DD=5723/2, CF=2404/2),
  normal=c(DD= 502/2, CF=3062/2)
)

# load data
get.dat=function(filename, sample.label, dist) {
  dat=read.csv("data/"%.%filename)
  
  if (sample.label=="DD") {
    dat$normal.epi =    dat$in_DD_Les_epi_region_1_normal=="True" 
    dat$normal.dermis = dat$in_DD_Les_dermis_region_1_normal=="True" 
    dat$lesion.epi =    dat$in_DD_Les_epi_region_2_lesion_alternate=="True" 
    dat$lesion.dermis = dat$in_DD_Les_dermis_region_2_lesion_alternate=="True" 
    
    # each pixel is 1/2 micron, negative dist in epi
    dat$dist2skin = dat$dist2DD_Les_skin_boundary/2
    dat$dist2DEJ  = dat$dist2DD_Les_DEJviaDAPI.WGA_alternate/2
    
  } else if (sample.label=="CF") {
    dat$normal.epi =    dat$in_CF_Les_epi_region_1_normal=="True" | dat$in_CF_Les_epi_region_3_normal=="True" 
    dat$normal.dermis = dat$in_CF_Les_dermis_region_1_normal=="True" | dat$in_CF_Les_dermis_region_3_normal=="True" 
    dat$lesion.epi =    dat$in_CF_Les_epi_region_2_lesion=="True" 
    dat$lesion.dermis = dat$in_CF_Les_dermis_region_2_lesion=="True" 
    
    # each pixel is 1/2 micron, negative dist in epi
    dat$dist2skin = dat$dist2CF_Les_skin_boundary/2
    dat$dist2DEJ  = dat$dist2CF_Les_DEJ_v2/2
    
  } else {
    stop("tissue not recognized: "%.%sample.label)
  }
  
  dat$epi = dat$normal.epi | dat$lesion.epi
  dat$dermis = dat$normal.dermis | dat$lesion.dermis
  dat$normal = dat$normal.epi | dat$normal.dermis
  dat$lesion = dat$lesion.epi | dat$lesion.dermis
  
  if (dist=="dist2DEJ") {
    dat$dist = dat$dist2DEJ  * ifelse(dat$epi, -1, 1)
  } else {
    dat$dist = dat$dist2skin 
  }
  dat
}


alldat = lapply (sample.labels, function(s) 
  lapply (dist.labels, function(dist) 
    lapply (markers.ls[[s]], function (l) 
      get.dat(s%.%"_Les_"%.%l%.%"_masks_d2b_regmem.csv", s, dist)
    )))