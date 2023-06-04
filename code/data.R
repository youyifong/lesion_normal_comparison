library(kyotil)
library(glue)

data.dir = "/fh/fast/fong_y/images_from_K/ROI_lesion_vs_normal_TSP"

region.labels=c(normal="intact", lesion="ulcerative")
region.labels.title=c(normal="Intact", lesion="Ulcerative")

dist.labels=c("dist2skin", "dist2DEJ"); names(dist.labels)=dist.labels

tissue.types=c("ulcerative","intact"); names(tissue.types)=tissue.types

biopsies=c("DD","CF","JM","AS","CJ","LG","LS"); names(biopsies)=biopsies
biopsies.1=c("DD","CF","JM"); names(biopsies.1)=biopsies.1

xlim.ls=list(dist2DEJ=c(-750,1750), dist2skin=c(0,2500))
ylim=c(0,175)

markers.1=c("CD3","CD15","CD123","CD11c"); names(markers.1)=markers.1
# CD68 did not pass QC for CF

markers.ls = lapply(biopsies, function(s) {
  files <- list.files(path = data.dir)
  tmp = list.files(data.dir%.%"/"%.%files[startsWith(files, s%.%"_Les_")], pattern=".*_Les_CD.*.masks.png")
  out = sapply(tmp, function(n) sub("_masks.png", "", sub(s%.%"_Les_","",n)))
  names(out)=out
  out
})
markers.ls

red=rgb(red=1, blue=0, green=0, alpha=1)
blue=rgb(red=0, blue=1, green=0, alpha=1)


tissue.widths = lapply(biopsies, function(s) {
# s="LS"
  dat=read.csv(glue("{data.dir}/{s}_Les_stitched_alignedtoCD45_cropped/{s}_Les_region_area_length.csv"))
  
  out = if (s=="DD") {
    c(ulcerative=dat$size[dat$name=='DD_Les_length_2_lesion'], 
      intact=dat$size[dat$name=='DD_Les_length_1_normal'])
    
  } else if (s=="CF") {
    c(ulcerative=dat$size[dat$name=='CF_Les_length_2_lesion'], 
      intact=dat$size[dat$name=='CF_Les_length_1_normal'] + dat$size[dat$name=='CF_Les_length_3_normal'])
    
  } else if (s=="JM") {
    c(ulcerative=dat$size[dat$name=='JM_Les_length_2_lesion'] + dat$size[dat$name=='JM_Les_length_4_lesion'], 
      intact=dat$size[dat$name=='JM_Les_length_1_normal'] + dat$size[dat$name=='JM_Les_length_3_normal'])
    
  } else if (s=="AS") {
    c(ulcerative=NA, 
      intact=1000)
    
  } else if (s=="CJ") {
    c(ulcerative=NA,
      intact=dat$size[dat$name=='CJ_Les_length_1_normal'])
    
  } else if (s=="LG") {
    c(ulcerative=NA, 
      intact=1000)
    
  } else if (s=="LS") {
    c(ulcerative=NA, 
      intact=1000)
  } 
  
  # each pixel is 1/2 micron
  out/2
})


# load mask-level data
get.dat=function(s, marker) {
# s="JM"; dist="dist2skin"; marker="CD3"
  dat=read.csv(glue("{data.dir}/{s}_Les_stitched_alignedtoCD45_cropped/{s}_Les_{marker}_masks_d2b_regmem.csv"))
  
  if (s=="DD") {
    dat$normal.epi =    dat$in_DD_Les_epi_region_1_normal=="True" 
    dat$normal.dermis = dat$in_DD_Les_dermis_region_1_normal=="True" 
    dat$lesion.epi =    dat$in_DD_Les_epi_region_2_lesion_alternate=="True" 
    dat$lesion.dermis = dat$in_DD_Les_dermis_region_2_lesion_alternate=="True" 
    
    dat$dist2skin = dat$dist2DD_Les_skin_boundary/2
    dat$dist2DEJ  = dat$dist2DD_Les_DEJviaDAPI.WGA_alternate/2
    
  } else if (s=="CF") {
    dat$normal.epi =    dat$in_CF_Les_epi_region_1_normal=="True" | dat$in_CF_Les_epi_region_3_normal=="True" 
    dat$normal.dermis = dat$in_CF_Les_dermis_region_1_normal=="True" | dat$in_CF_Les_dermis_region_3_normal=="True" 
    dat$lesion.epi =    dat$in_CF_Les_epi_region_2_lesion=="True" 
    dat$lesion.dermis = dat$in_CF_Les_dermis_region_2_lesion=="True" 

    dat$dist2skin = dat$dist2CF_Les_skin_boundary/2
    dat$dist2DEJ  = dat$dist2CF_Les_DEJ_v2/2
    
  } else if (s=="JM") {
    dat$normal.epi =    dat$in_JM_Les_epi_region_1_normal   =="True" | dat$in_JM_Les_epi_region_3_normal      =="True"
    dat$normal.dermis = dat$in_JM_Les_dermis_region_1_normal=="True" | dat$in_JM_Les_dermis_region_3_normal   =="True"
    dat$lesion.epi =    dat$in_JM_Les_epi_region_2_lesion   =="True" | dat$in_JM_Les_epi_region_4_lesion      =="True"
    dat$lesion.dermis = dat$in_JM_Les_dermis_region_2_lesion=="True" | dat$in_JM_Les_dermis_region_4_lesion   =="True"

    dat$dist2skin = dat$dist2JM_Les_skin_boundary              
    dat$dist2DEJ  = dat$dist2JM_Les_DEJ_v3                     
    
  } else if (s=="AS") {
    dat$normal.epi =    dat$in_AS_Les_epi_region_1_normal   =="True" 
    dat$normal.dermis = dat$in_AS_Les_dermis_region_1_normal=="True" 
    dat$lesion.epi =    FALSE
    dat$lesion.dermis = FALSE

    dat$dist2skin = dat$dist2AS_Les_skin_boundary_v2    
    dat$dist2DEJ  = dat$dist2AS_Les_DEJ_v2              
    
  } else if (s=="CJ") {
    dat$normal.epi =    dat$in_CJ_Les_epi_region_1_normal_v3   =="True" 
    dat$normal.dermis = dat$in_CJ_Les_dermis_region_1_normal_v3=="True"
    dat$lesion.epi =    FALSE
    dat$lesion.dermis = FALSE

    dat$dist2skin = dat$dist2CJ_Les_skin_boundary_v3       
    dat$dist2DEJ  = dat$dist2CJ_Les_DEJ_v3                 
    
  } else if (s=="LG") {
    dat$normal.epi =    dat$in_LG_Les_epi_region_1_normal   =="True" 
    dat$normal.dermis = dat$in_LG_Les_dermis_region_1_normal=="True" & dat$in_LG_Les_dermis_region_1_void  !="True" 
    dat$lesion.epi =    FALSE
    dat$lesion.dermis = FALSE

    dat$dist2skin = dat$dist2LG_Les_skin_boundary       
    dat$dist2DEJ  = dat$dist2LG_Les_DEJ                 
    
  } else if (s=="LS") {
    dat$normal.epi =    dat$in_LS_Les_epi_region_1_normal   =="True" 
    dat$normal.dermis = dat$in_LS_Les_dermis_region_1_normal=="True" 
    dat$lesion.epi =    FALSE
    dat$lesion.dermis = FALSE

    dat$dist2skin = dat$dist2LS_Les_skin_boundary       
    dat$dist2DEJ  = dat$dist2LS_Les_DEJ_v2              
    
  } else {
    stop("tissue not recognized: "%.%s)
  }
  
  dat$epi = dat$normal.epi | dat$lesion.epi
  dat$dermis = dat$normal.dermis | dat$lesion.dermis

  dat$region = NA
  dat$region[dat$normal.epi | dat$normal.dermis] = "intact"
  dat$region[dat$lesion.epi | dat$lesion.dermis] = "ulcerative"

  # each pixel is 1/2 micron, negative dist in epi
  dat$dist2DEJ = dat$dist2DEJ * ifelse(dat$epi, -1, 1)/2
  dat$dist2skin = dat$dist2skin/2
  dat
}


alldat = 
lapply (biopsies, function(s) {
lapply (markers.ls[[s]], function (marker) {
  # myprint(s, dist, marker)
  get.dat(s, marker)
})
})
