library(kyotil)
library(glue)

{
data.dir = "../Chipcytometry_images_from32bit"

tissue_size_by_depth_and_region=read.csv(paste0(data.dir, "/tissue_size_by_depth_and_region.csv"))

all.slides = sort(unique(tissue_size_by_depth_and_region$biopsy)); names(all.slides)=all.slides

markers.ls = lapply(all.slides, function(s) {
  # s is like AS/AS_Les1, extract AS_Les1 from s
  s1 = sub(".*?/(.*)", "\\1", s)
  tmp = list.files(paste0(data.dir, "/", s, "/TSP_stitched"), pattern=".*_CD.*_m_c7.png")
  out = sapply(tmp, function(n) sub("_m_.*.png", "", sub(s1%.%"_","",n)))
  names(out)=out
  out
})

all.slides = all.slides[sapply(markers.ls, function(x) length(x)>0)]

# subject initials
initials = unique(sort(sapply(all.slides, function(s) sub("/.*", "", s))))

# choose one slide for each biopsie
slides = c("AS/AS1_Les", "CF/CF1_Les", "CJ/CJ1_Les", "DD/DD1_Les",  "EB/EB_Les", "JK/JK1_Les", 
          "JM/JM1_Les",  "JT/JT_Les",  "LD/LD_Les", "LG/LG2_Les",  "LM/LM_Les", 
          "LS/LS2_Les", "SK/SK2_Les")
# use biopsies as an alias for slides
biopsies = slides

key.markers=c("CD11c", "CD123", "CD3", "CD15", "CD56", "CD68")

# make sure all slides have all key markers
all(sapply(slides, function(b) length(intersect(markers.ls[[b]], key.markers)))==6)

# define a function that returns c("DD1_Les") from c("DD/DD1_Les")
get.slide.name=function(ss)  sapply(ss, function (s) sub(".*?/(.*)", "\\1", s))
# define a function that returns c("DD") from c("DD/DD1_Les")
get.initial=function(ss)  sapply(ss, function (s) sub("/.*", "", s))

# read severity data
severity = read.csv("R/Lesion_HSV_rates_04082025.csv")
# add a column for the initial, which are the letter in front of the arabic number
severity$initial = sapply(severity$PID, function(x) sub("[0-9]+$", "", x))

}




{
region.labels=c(normal="intact", lesion="ulcerative")
region.labels.title=c(normal="Intact", lesion="Ulcerative")

dist.labels=c("dist2skin", "dist2DEJ"); names(dist.labels)=dist.labels

tissue.types=c("intact","ulcerative"); names(tissue.types)=tissue.types

xlim.ls=list(dist2DEJ=c(-750,1750), dist2skin=c(0,2500))
ylim=list(dist2DEJ=c(0,100), dist2skin=c(0,100))


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
# s="LG/LG2_Les"; dist="dist2skin"; marker="CD3
  
  # use CF5 CD68 instead of CF1 CD68 due to quality
  if (s == "CF/CF1_Les" && marker == "CD68") {
    s = "CF/CF5_Les"
  }
  
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


alldat = 
lapply (all.slides, function(s) {
lapply (markers.ls[[s]], function (marker) {
  # myprint(s, dist, marker)
  get.dat(s, marker)
})
})
