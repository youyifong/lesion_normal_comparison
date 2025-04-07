library(kyotil)
library(glue)

{
data.dir = "../Chipcytometry_images_from32bit"
tissue_size_by_depth_and_region=read.csv(paste0(data.dir, "/tissue_size_by_depth_and_region.csv"))

# list unique biopsies
biopsies = unique(tissue_size_by_depth_and_region$biopsy)
biopsies = sort(biopsies)
names(biopsies)=biopsies

# a select subset
biopsies.1=c("DD/DD_Les","CF_Les","JM_Les")
; names(biopsies.1)=biopsies.1
}


{
region.labels=c(normal="intact", lesion="ulcerative")
region.labels.title=c(normal="Intact", lesion="Ulcerative")

dist.labels=c("dist2skin", "dist2DEJ"); names(dist.labels)=dist.labels

tissue.types=c("ulcerative","intact"); names(tissue.types)=tissue.types

xlim.ls=list(dist2DEJ=c(-750,1750), dist2skin=c(0,2500))
ylim=list(dist2DEJ=c(0,100), dist2skin=c(0,100))

markers.1=c("CD3","CD15","CD123","CD11c"); names(markers.1)=markers.1
# CD68 did not pass QC for CF

markers.ls = lapply(biopsies, function(s) {
  # s is like AS/AS_Les1, extract AS_Les1 from s
  s1 = sub(".*?/(.*)", "\\1", s)
  tmp = list.files(paste0(data.dir, "/", s, "/TSP_stitched"), pattern=".*_CD.*_m_c7.png")
  out = sapply(tmp, function(n) sub("_m_.*.png", "", sub(s1%.%"_","",n)))
  names(out)=out
  out
})
markers.ls

# list region ROI files
roi.ls = lapply(biopsies, function(s) {
  list.files(paste0(data.dir, "/", s, "/TSP_stitched/ROI_regions"), pattern=".*_region_.*.roi")
})
roi.ls

red=rgb(red=1, blue=0, green=0, alpha=1)
blue=rgb(red=0, blue=1, green=0, alpha=1)
}


# load mask-level data
get.dat=function(s, marker) {
# s="LG/LG2_Les"; dist="dist2skin"; marker="CD3"
  
  s1 = sub(".*?/(.*)", "\\1", s)
  dat=read.csv(glue("{data.dir}/{s}/TSP_stitched/{s1}_{marker}_m_c7_d2b_regmem.csv"))
  
  # find the column with epi_region and normal in the name
  normal_epi_col = grep("epi_region", grep("normal", colnames(dat), value=TRUE), value=TRUE)
  # find the column with dermis_region and normal in the name
  normal_dermis_col = grep("dermis_region", grep("normal", colnames(dat), value=TRUE), value=TRUE)
  # find the column with epi_region and lesion in the name
  lesion_epi_col = grep("epi_region", grep("lesion", colnames(dat), value=TRUE), value=TRUE)
  # find the column with dermis_region and lesion in the name
  lesion_dermis_col = grep("dermis_region", grep("lesion", colnames(dat), value=TRUE), value=TRUE)
  
  # define column normal.epi as any of the normal_epi_col columns equals "True"
  dat$normal.epi = rowSums(dat[, normal_epi_col, drop=F] == "True", na.rm=TRUE) > 0
  # define column normal.dermis as any of the normal_dermis_col columns equals "True"
  dat$normal.dermis = rowSums(dat[, normal_dermis_col, drop=F] == "True", na.rm=TRUE) > 0
  # define column lesion.epi as any of the lesion_epi_col columns equals "True"
  dat$lesion.epi = rowSums(dat[, lesion_epi_col, drop=F] == "True", na.rm=TRUE) > 0
  # define column lesion.dermis as any of the lesion_dermis_col columns equals "True"
  dat$lesion.dermis = rowSums(dat[, lesion_dermis_col, drop=F] == "True", na.rm=TRUE) > 0
  
  # find the column with skin_boundary in the name
  skin_boundary_col = grep("skin_boundary", colnames(dat), value=TRUE)
  # assert there is only one column
  if (length(skin_boundary_col) != 1) {
    stop("more than one skin boundary column found")
  }
  # find the column with DEJ in the name
  DEJ_col = grep("DEJ", colnames(dat), value=TRUE)
  # assert there is only one column
  if (length(DEJ_col) != 1) {
    stop("more than one DEJ column found")
  }
  # define column dist2skin as the skin_boundary_col
  dat$dist2skin = dat[, skin_boundary_col]
  # define column dist2DEJ as the DEJ_col
  dat$dist2DEJ = dat[, DEJ_col]
  
  
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


# skip a biopsy if markers.1 are not all present
biopsies.tmp = sapply (biopsies, function(s) {
  markers=markers.ls[[s]]
  if (any(!markers.1 %in% markers)) {
    return("")
  } else {
    return(s)
  }
})
biopsies=biopsies.tmp[biopsies.tmp!=""]

alldat = 
lapply (biopsies, function(s) {
lapply (markers.ls[[s]], function (marker) {
  # myprint(s, dist, marker)
  get.dat(s, marker)
})
})
