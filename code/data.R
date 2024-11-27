library(kyotil)
library(glue)

data.dir = "../Chipcytometry_images_from8bit"

{
region.labels=c(normal="intact", lesion="ulcerative")
region.labels.title=c(normal="Intact", lesion="Ulcerative")

dist.labels=c("dist2skin", "dist2DEJ"); names(dist.labels)=dist.labels

tissue.types=c("ulcerative","intact"); names(tissue.types)=tissue.types

biopsies=c("DD_Les","CF_Les","JM_Les","AS_Les","CJ_Les","LG_Les","LS_Les"); names(biopsies)=biopsies
biopsies.1=c("DD_Les","CF_Les","JM_Les"); names(biopsies.1)=biopsies.1

xlim.ls=list(dist2DEJ=c(-750,1750), dist2skin=c(0,2500))
ylim=c(0,175)

markers.1=c("CD3","CD15","CD123","CD11c"); names(markers.1)=markers.1
# CD68 did not pass QC for CF

markers.ls = lapply(biopsies, function(s) {
  files <- list.files(path = data.dir)
  tmp = list.files(data.dir%.%"/"%.%files[startsWith(files, s)]%.%"/ZKWstitched_alignedtoCD45_cropped", pattern=".*_CD.*.masks.png")
  out = sapply(tmp, function(n) sub("_masks.png", "", sub(s%.%"_","",n)))
  names(out)=out
  out
})
markers.ls

red=rgb(red=1, blue=0, green=0, alpha=1)
blue=rgb(red=0, blue=1, green=0, alpha=1)
}


tissue_size_by_depth_and_region=read.csv(paste0(data.dir, "/tissue_size_by_depth_and_region.csv"))



# tissue.widths is not used for normalzing cell density because we are now using depth-level biopsy size data.
# it is used to give a rough indication of the size of the tissue in the biopsy
tissue.widths = lapply(biopsies, function(s) {
# s="DD"
  dat=read.csv(glue("{data.dir}/{s}/ZKWstitched_alignedtoCD45_cropped/{s}_region_area_length.csv"))

  out = if (s=="DD_Les") {
    c(ulcerative=dat$size[dat$name=='DD_Les_length_2_lesion'],
      intact=dat$size[dat$name=='DD_Les_length_1_normal'])

  } else if (s=="CF_Les") {
    c(ulcerative=dat$size[dat$name=='CF_Les_length_2_lesion'],
      intact=dat$size[dat$name=='CF_Les_length_1_normal'] + dat$size[dat$name=='CF_Les_length_3_normal'])

  } else if (s=="JM_Les") {
    c(ulcerative=dat$size[dat$name=='JM_Les_length_4_lesion'],# + dat$size[dat$name=='JM_Les_length_2_lesion']
      intact=dat$size[dat$name=='JM_Les_length_1_normal'] + dat$size[dat$name=='JM_Les_length_3_normal'])

  } else if (s=="AS_Les") {
    c(ulcerative=NA,
      intact=dat$size[dat$name=='AS_Les_length_1_normal'])

  } else if (s=="CJ_Les") {
    c(ulcerative=NA,
      intact=dat$size[dat$name=='CJ_Les_length_1_normal'])

  } else if (s=="LG_Les") {
    c(ulcerative=NA,
      intact=dat$size[dat$name=='LG_Les_length_1_normal'])

  } else if (s=="LS_Les") {
    c(ulcerative=NA,
      intact=dat$size[dat$name=='LS_Les_length_1_normal'])
  }

  # each pixel is 1/2 micron
  out/2
})


# load mask-level data
get.dat=function(s, marker) {
# s="JM"; dist="dist2skin"; marker="CD3"
  dat=read.csv(glue("{data.dir}/{s}/ZKWstitched_alignedtoCD45_cropped/{s}_{marker}_masks_d2b_regmem.csv"))
  
  if (s=="DD_Les") {
    dat$normal.epi =    dat$in_DD_Les_epi_region_1_normal=="True" 
    dat$normal.dermis = dat$in_DD_Les_dermis_region_1_normal=="True" 
    dat$lesion.epi =    dat$in_DD_Les_epi_region_2_lesion_alternate=="True" 
    dat$lesion.dermis = dat$in_DD_Les_dermis_region_2_lesion_alternate=="True" 
    
    dat$dist2skin = dat$dist2DD_Les_skin_boundary
    dat$dist2DEJ  = dat$dist2DD_Les_DEJviaDAPI.WGA_alternate
    
  } else if (s=="CF_Les") {
    dat$normal.epi =    dat$in_CF_Les_epi_region_1_normal=="True" | dat$in_CF_Les_epi_region_3_normal=="True" 
    dat$normal.dermis = dat$in_CF_Les_dermis_region_1_normal=="True" | dat$in_CF_Les_dermis_region_3_normal=="True" 
    dat$lesion.epi =    dat$in_CF_Les_epi_region_2_lesion=="True" 
    dat$lesion.dermis = dat$in_CF_Les_dermis_region_2_lesion=="True" 

    dat$dist2skin = dat$dist2CF_Les_skin_boundary
    dat$dist2DEJ  = dat$dist2CF_Les_DEJ_v2
    
  } else if (s=="JM_Les") {
    dat$normal.epi =    dat$in_JM_Les_epi_region_1_normal   =="True" | dat$in_JM_Les_epi_region_3_normal      =="True"
    dat$normal.dermis = dat$in_JM_Les_dermis_region_1_normal=="True" | dat$in_JM_Les_dermis_region_3_normal   =="True"
    dat$lesion.epi =    dat$in_JM_Les_epi_region_4_lesion      =="True" # | dat$in_JM_Les_epi_region_2_lesion   =="True" 
    dat$lesion.dermis = dat$in_JM_Les_dermis_region_4_lesion   =="True" # | dat$in_JM_Les_dermis_region_2_lesion=="True"

    dat$dist2skin = dat$dist2JM_Les_skin_boundary              
    dat$dist2DEJ  = dat$dist2JM_Les_DEJ_v3                     
    
  } else if (s=="AS_Les") {
    dat$normal.epi =    dat$in_AS_Les_epi_region_1_normal   =="True" 
    dat$normal.dermis = dat$in_AS_Les_dermis_region_1_normal=="True" 
    dat$lesion.epi =    FALSE
    dat$lesion.dermis = FALSE

    dat$dist2skin = dat$dist2AS_Les_skin_boundary_v2    
    dat$dist2DEJ  = dat$dist2AS_Les_DEJ_v2              
    
  } else if (s=="CJ_Les") {
    dat$normal.epi =    dat$in_CJ_Les_epi_region_1_normal_v3   =="True" 
    dat$normal.dermis = dat$in_CJ_Les_dermis_region_1_normal_v3=="True"
    dat$lesion.epi =    FALSE
    dat$lesion.dermis = FALSE

    dat$dist2skin = dat$dist2CJ_Les_skin_boundary_v3       
    dat$dist2DEJ  = dat$dist2CJ_Les_DEJ_v3                 
    
  } else if (s=="LG_Les") {
    dat$normal.epi =    dat$in_LG_Les_epi_region_1_normal   =="True" 
    dat$normal.dermis = dat$in_LG_Les_dermis_region_1_normal=="True" & dat$in_LG_Les_dermis_region_1_void  !="True" 
    dat$lesion.epi =    FALSE
    dat$lesion.dermis = FALSE

    dat$dist2skin = dat$dist2LG_Les_skin_boundary       
    dat$dist2DEJ  = dat$dist2LG_Les_DEJ                 
    
  } else if (s=="LS_Les") {
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
