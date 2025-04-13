library(kyotil)
library(glue)

{
data.dir = "../Chipcytometry_images_from32bit"

tissue_size_by_depth_and_region=read.csv(paste0(data.dir, "/tissue_size_by_depth_and_region.csv"))

# define a function that returns c("DD1_Les") from c("DD/DD1_Les")
get.slide.name=function(ss)  sapply(ss, function (s) sub(".*?/(.*)", "\\1", s))

# define a function that returns c("DD") from c("DD1_Les")
get.initial=function(ss)  sapply(ss, function (s) sub("^([^0-9_]+).*", "\\1", s) )

all.slides = sort(unique(tissue_size_by_depth_and_region$biopsy))
all.slides=get.slide.name(all.slides)
names(all.slides)=all.slides

markers.ls = lapply(all.slides, function(s) {
  # s = 'JT_Les'
  tmp = list.files(paste0(data.dir, "/", get.initial(s), "/", s, "/TSP_stitched"), pattern=".*_CD.*_m_c7.png")
  tmp = sort(tmp)
  out = sapply(tmp, function(n) sub("_m_.*.png", "", sub(s%.%"_","",n)))
  names(out)=out
  out
})

all.slides = all.slides[sapply(markers.ls, function(x) length(x)>0)]

# choose one slide for each biopsie
select.slides = c("AS1_Les", "CF1_Les", "CJ1_Les", "DD1_Les",  "EB_Les", "JK1_Les", 
          "JM1_Les", "JT_Les",  "LD_Les", "LG2_Les",  "LM_Les", 
          "LS2_Les", "SK2_Les")
names(select.slides)=get.slide.name(select.slides)

# average reg.densities of all slides from a single biopsy
avg_slides = glue("{biopsies}_Les_avg")

biopsies=sapply(select.slides, function(x) get.initial(x)); names(biopsies)=biopsies

slides.by.biopsy = lapply(biopsies, function(b) all.slides[startsWith(all.slides,b)])

key.markers=c("CD3", "CD11c", "CD123", "CD15", "CD56", "CD68")

# make sure all select.slides have all key markers
all(sapply(select.slides, function(b) length(intersect(markers.ls[[b]], key.markers)))==6)

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


red=rgb(red=1, blue=0, green=0, alpha=1)
blue=rgb(red=0, blue=1, green=0, alpha=1)
}


# load mask-level data
get.dat=function(s, marker) {
# s="LG2_Les"; dist="dist2skin"; marker="CD3
  
  # use CF5 CD68 instead of CF1 CD68 due to quality
  if (s == "CF1_Les" && marker == "CD68") {
    s = "CF5_Les"
  }
  
  dat=read.csv(glue("{data.dir}/{get.initial(s)}/{s}/TSP_stitched/{s}_{marker}_m_c7_d2b_regmem.csv"))
  
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
  get.dat(s, marker)
})
})




################################################################################
## tabulate cell counts in epi/dermis x normal/lesion

tissue_size = read.csv(glue("{data.dir}/tissue_size_by_depth_and_region.csv"))

tissue_size$biopsy = get.slide.name(tissue_size$biopsy)

# count cells to different depths below DEJ
count.tables = lapply (all.slides, function(s) {
  # s="DD1_Les"
  
  markers=markers.ls[[s]]; names(markers)=markers
  
  depths=c(500, 1000, Inf)
  names(depths)=paste0("DEJ", depths)
  lapply (depths, function(max_depth) { # um
    
    myprint(s, max_depth)
    
    counts = sapply(markers, function (a) 
      c(EI=nrow(subset(alldat[[s]][[a]], dist2DEJ<max_depth & epi &    region=="intact")),
        DI=nrow(subset(alldat[[s]][[a]], dist2DEJ<max_depth & dermis & region=="intact")),
        EU=nrow(subset(alldat[[s]][[a]], dist2DEJ<max_depth & epi &    region=="ulcerative")),
        DU=nrow(subset(alldat[[s]][[a]], dist2DEJ<max_depth & dermis & region=="ulcerative")))
    )
    
    
    # compute sizes of epi/dermis, normal/lesion
    dat = tissue_size[tissue_size$biopsy==s & tissue_size$dist_type=="dist2DEJ",] # max_depth is relative to DEJ
    dat = dat[dat$depth<max_depth,] # < b/c depth is the beginning of the bin
    sizes = with(dat, c(EI=sum(area[depth<0 & region=="intact"]),
                        DI=sum(area[depth>0 & region=="intact"]),
                        EU=sum(area[depth<0 & region=="ulcerative"]),
                        DU=sum(area[depth>0 & region=="ulcerative"])))
    # convert size from 100um x 100um to 1mm x 1mm
    sizes = sizes/100
    
    # density
    counts2=counts/sizes
    mode(counts2) <- "integer"
    
    ## depth-specific density
    
    # first, compute sizes of epi/dermis, normal/lesion
    dat = tissue_size[tissue_size$biopsy==s & tissue_size$dist_type=="dist2DEJ",] # max_depth is relative to DEJ
    dat = dat[dat$depth<max_depth,] # use < instead of <= b/c depth is the beginning of the bin
    dat$area = dat$area/100    # convert size from 100um x 100um to 1mm x 1mm
    size_epi_normal = dat[dat$depth<0 & dat$region=="intact",]
    size_epi_lesion = dat[dat$depth<0 & dat$region=="ulcerative",]
    size_dermis_normal = dat[dat$depth>0 & dat$region=="intact",]
    size_dermis_lesion = dat[dat$depth>0 & dat$region=="ulcerative",]
    
    counts3 = sapply(markers, function (a) {
      # epi normal
      epi_normal = subset(alldat[[s]][[a]], dist2DEJ<max_depth & epi &    region=="intact")
      # add a depth column computed as the floor of dist2DEJ/100
      epi_normal$depth = floor(epi_normal$dist2DEJ/100)*100
      count_by_depth = mytable(epi_normal$depth)
      # add count_by_depth column to size_epi_normal by match the depth in size_epi_normal and count_by_depth
      size_epi_normal$count_by_depth = count_by_depth[size_epi_normal$depth%.%""]
      epi_normal_density = sum(size_epi_normal$count_by_depth/size_epi_normal$area, na.rm=TRUE)
      
      # dermis normal
      dermis_normal = subset(alldat[[s]][[a]], dist2DEJ<max_depth & dermis & region=="intact")
      # add a depth column computed as the floor of dist2DEJ/100
      dermis_normal$depth = floor(dermis_normal$dist2DEJ/100)*100
      count_by_depth = mytable(dermis_normal$depth)
      # add count_by_depth column to size_dermis_normal by match the depth in size_dermis_normal and count_by_depth
      size_dermis_normal$count_by_depth = count_by_depth[size_dermis_normal$depth%.%""]
      dermis_normal_density = sum(size_dermis_normal$count_by_depth/size_dermis_normal$area, na.rm=TRUE)
      
      # epi lesion
      epi_lesion = subset(alldat[[s]][[a]], dist2DEJ<max_depth & epi &    region=="ulcerative")
      # add a depth column computed as the floor of dist2DEJ/100
      epi_lesion$depth = floor(epi_lesion$dist2DEJ/100)*100
      count_by_depth = mytable(epi_lesion$depth)
      # add count_by_depth column to size_epi_lesion by match the depth in size_epi_lesion and count_by_depth
      size_epi_lesion$count_by_depth = count_by_depth[size_epi_lesion$depth%.%""]
      epi_lesion_density = sum(size_epi_lesion$count_by_depth/size_epi_lesion$area, na.rm=TRUE)
      
      # dermis lesion
      dermis_lesion = subset(alldat[[s]][[a]], dist2DEJ<max_depth & dermis & region=="ulcerative")
      # add a depth column computed as the floor of dist2DEJ/100
      dermis_lesion$depth = floor(dermis_lesion$dist2DEJ/100)*100
      count_by_depth = mytable(dermis_lesion$depth)
      # add count_by_depth column to size_dermis_lesion by match the depth in size_dermis_lesion and count_by_depth
      size_dermis_lesion$count_by_depth = count_by_depth[size_dermis_lesion$depth%.%""]
      dermis_lesion_density = sum(size_dermis_lesion$count_by_depth/size_dermis_lesion$area, na.rm=TRUE)
      
      # return the counts
      c(EI=epi_normal_density, DI=dermis_normal_density, EU=epi_lesion_density, DU=dermis_lesion_density)
    })
    mode(counts3) <- "integer"
    
    
    list(sizes=sizes, counts=counts, densities=counts2, reg.densities=counts3)
    
  })
}) 


# count cells to different depths below skin
count.tables.2 = lapply (all.slides, function(s) {
  # s="DD1_Les" 
  
  markers=markers.ls[[s]]; names(markers)=markers
  
  depths=c(1000, Inf)
  names(depths)=paste0("Skin", depths)
  lapply (depths, function(max_depth) { # um
    
    myprint(s, max_depth)
    
    counts = sapply(markers, function (a) 
      c(EI=nrow(subset(alldat[[s]][[a]], dist2skin<max_depth & epi &    region=="intact")),
        DI=nrow(subset(alldat[[s]][[a]], dist2skin<max_depth & dermis & region=="intact")),
        EU=nrow(subset(alldat[[s]][[a]], dist2skin<max_depth & epi &    region=="ulcerative")),
        DU=nrow(subset(alldat[[s]][[a]], dist2skin<max_depth & dermis & region=="ulcerative")))
    )
    
    
    # compute sizes of epi/dermis, normal/lesion
    dat = tissue_size[tissue_size$biopsy==s & tissue_size$dist_type=="dist2skin",] # max_depth is relative to skin
    dat = dat[dat$depth<max_depth,] # < b/c depth is the beginning of the bin
    sizes = with(dat, c(I=sum(area[region=="intact"]),
                        U=sum(area[region=="ulcerative"])))
    # convert size from 100um x 100um to 1mm x 1mm
    sizes = sizes/100
    # rep to match the dimension of counts: I I U U
    sizes = rep(sizes, each = 2)
    
    # density
    counts2=counts/sizes
    mode(counts2) <- "integer"
    
    ## depth-specific density
    
    # first, compute sizes of normal/lesion
    dat = tissue_size[tissue_size$biopsy==s & tissue_size$dist_type=="dist2skin",] # max_depth is relative to skin
    dat = dat[dat$depth<max_depth,] # use < instead of <= b/c depth is the beginning of the bin
    dat$area = dat$area/100    # convert size from 100um x 100um to 1mm x 1mm
    size_normal = dat[dat$region=="intact",]
    size_lesion = dat[dat$region=="ulcerative",]
    size_normal = dat[dat$region=="intact",]
    size_lesion = dat[dat$region=="ulcerative",]
    
    counts3 = sapply(markers, function (a) {
      # normal
      normal = subset(alldat[[s]][[a]], dist2skin<max_depth & region=="intact")
      # add a depth column computed as the floor of dist2skin/100
      normal$depth = floor(normal$dist2skin/100)*100
      count_by_depth = mytable(normal$depth)
      # add count_by_depth column to size_normal by match the depth in size_normal and count_by_depth
      size_normal$count_by_depth = count_by_depth[size_normal$depth%.%""]
      normal_density = sum(size_normal$count_by_depth/size_normal$area, na.rm=TRUE)
      
      # epi lesion
      lesion = subset(alldat[[s]][[a]], dist2skin<max_depth & region=="ulcerative")
      # add a depth column computed as the floor of dist2skin/100
      lesion$depth = floor(lesion$dist2skin/100)*100
      count_by_depth = mytable(lesion$depth)
      # add count_by_depth column to size_lesion by match the depth in size_lesion and count_by_depth
      size_lesion$count_by_depth = count_by_depth[size_lesion$depth%.%""]
      lesion_density = sum(size_lesion$count_by_depth/size_lesion$area, na.rm=TRUE)
      
      # return the counts
      c(I=normal_density, U=lesion_density)
    })
    mode(counts3) <- "integer"
    
    
    list(sizes=sizes, counts=counts, densities=counts2, reg.densities=counts3)
    
  })
}) 


# create avg.slides that average regularized densities across slides for a biopsy
out = lapply(biopsies, function(b){
  
  depths=c(500, 1000, Inf)
  depths=paste0("DEJ", depths)
  names(depths)=depths
  out=lapply (depths, function(max_depth) { # um
    tab=sapply(slides.by.biopsy[[b]],simplify="array", function(s){
      count.tables[[s]][[max_depth]]$reg.densities[,"CD3",drop=F]
    })
    out=apply(tab,c(1,2), mean)
    list(reg.densities=out, sizes=NA, counts=NA, densities=NA)
  })
  
})
names(out) = glue("{biopsies}_Les_avg")
count.tables = c(count.tables, out)

# repeat for count.tables.2
out = lapply(biopsies, function(b){
  
  depths=c(1000, Inf)
  depths=paste0("Skin", depths)
  names(depths)=depths
  out=lapply (depths, function(max_depth) { # um
    tab=sapply(slides.by.biopsy[[b]],simplify="array", function(s){
      count.tables.2[[s]][[max_depth]]$reg.densities[,"CD3",drop=F]
    })
    out=apply(tab,c(1,2), mean)
    list(reg.densities=out, sizes=NA, counts=NA, densities=NA)
  })
  
})
names(out) = glue("{biopsies}_Les_avg")
count.tables.2 = c(count.tables.2, out)


# save the objects except alldat to a file
# alldat is large and only used for plotting histograms
save(list = setdiff(ls(), "alldat"), file = "R/forRmd.RData")
