library(kyotil)
library(glue)

{ # slide level data
data.dir = "../Chipcytometry_images_from32bit"
key.markers=c("CD3", "CD11c", "CD123", "CD15", "CD56", "CD68")
cell_types=c("CD4T", "CD8T", "Treg")
get.slide.name=function(ss)  sapply(ss, function (s) sub(".*?/(.*)", "\\1", s)) # "DD/DD1_Les" -> DD1_Les
get.initial=function(ss)  sapply(ss, function (s) sub("^([^0-9_]+).*", "\\1", s) ) # "DD1_Les" -> "DD"

# in tissue_size, the depth column is the starting point of the bin.
tissue_size=read.csv(paste0(data.dir, "/tissue_size_by_depth_and_region.csv"))
tissue_size$biopsy = get.slide.name(tissue_size$biopsy)
all.slides = sort(unique(tissue_size$biopsy)); names(all.slides)=all.slides

markers.ls = lapply(all.slides, function(s) {  # s = 'JT_Les'
  tmp = sort(list.files(paste0(data.dir, "/", get.initial(s), "/", s, "/TSP_stitched"), pattern=".*_d2b_regmem.csv"))
  out = sapply(tmp, function(n) sub("_m.*.csv", "", sub(s%.%"_","",n)))
  # move CD4T, CD8T, Treg to the end
  out = c(setdiff(out, cell_types), setdiff(out, setdiff(out, cell_types)))
  names(out)=out; out
})

# remove slides with no markers
all.slides = all.slides[sapply(markers.ls, function(x) length(x)>0)]
# choose one slide for each biopsy
select.slides = c("AS1_Les", "CF1_Les", "CJ1_Les", "DD1_Les",  "EB_Les", "JK1_Les", "JM1_Les", "JT_Les",  "LD_Les", "LG2_Les",  "LM_Les", "LS2_Les", "SK2_Les"); names(select.slides)=get.slide.name(select.slides)
stopifnot(all(sapply(select.slides, function(b) length(intersect(markers.ls[[b]], key.markers)))==6))# make sure all select.slides have all key markers
biopsies=sapply(select.slides, function(x) get.initial(x)); names(biopsies)=biopsies
avg_slides = glue("{biopsies}_Les_avg")# average norm.counts of all slides from a single biopsy
slides.by.biopsy = lapply(biopsies, function(b) all.slides[startsWith(all.slides,b)])

severity = read.csv("R/Lesion_HSV_rates_04082025.csv")
severity$initial = sapply(severity$PID, function(x) sub("[0-9]+$", "", x))# add a column for the initial, which are the letter in front of the arabic number

region.labels=c(normal="intact", lesion="ulcerative")
region.labels.title=c(normal="Intact", lesion="Ulcerative")
dist.labels=c("dist2skin", "dist2DEJ"); names(dist.labels)=dist.labels
tissue.types=c("intact","ulcerative"); names(tissue.types)=tissue.types
xlim.ls=list(dist2DEJ=c(-750,1750), dist2skin=c(0,2500))
ylim=list(dist2DEJ=c(0,100), dist2skin=c(0,100))
red=rgb(red=1, blue=0, green=0, alpha=1)
blue=rgb(red=0, blue=1, green=0, alpha=1)
}



################################################################################
# cell-level data

get.dat=function(s, marker) {
# s="LG2_Les"; dist="dist2skin"; marker="CD3
    # use CF5 CD68 instead of CF1 CD68 due to quality
  if (s == "CF1_Les" && marker == "CD68") s = "CF5_Les"
  
  fn = glue("{data.dir}/{get.initial(s)}/{s}/TSP_stitched/{s}_{marker}_m_c7_d2b_regmem.csv")
  if (!file.exists(fn)) {
    fn = glue("{data.dir}/{get.initial(s)}/{s}/TSP_stitched/{s}_{marker}_masks_d2b_regmem.csv")
    if (!file.exists(fn)) return (NULL)
  }
  dat=read.csv(fn)
  
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


alldat =  lapply (all.slides, function(s) {
          lapply (markers.ls[[s]], function (marker) {
            get.dat(s, marker)
          })
          })


################################################################################
## count cells in epi/dermis x normal/lesion
## count.tables stratifies by distance from DEJ
## count.tables.2 stratifies by distance from skin

count.tables = lapply (all.slides, function(s) {
  # s="LG3_Les"; max_depth=500
  
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
                        DI=sum(area[depth>=0 & region=="intact"]),
                        EU=sum(area[depth<0 & region=="ulcerative"]),
                        DU=sum(area[depth>=0 & region=="ulcerative"])))
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
    size_dermis_normal = dat[dat$depth>=0 & dat$region=="intact",]
    size_dermis_lesion = dat[dat$depth>=0 & dat$region=="ulcerative",]
    
    counts3 = NULL
    density.by.depth = list()
    for (a in markers) {
      
      # epi normal
      
      epi_normal = subset(alldat[[s]][[a]], dist2DEJ<max_depth & epi &    region=="intact")
      # add a depth column computed as the floor of dist2DEJ/100
      epi_normal$depth = floor(epi_normal$dist2DEJ/100)*100
      count_by_depth = mytable(epi_normal$depth)
      # add count_by_depth column to size_epi_normal by matching the depth in size_epi_normal and count_by_depth
      size_epi_normal$count_by_depth = count_by_depth[size_epi_normal$depth%.%""]
      # set NA in size_epi_normal$count_by_depth to 0 because there are no cells at the depth
      size_epi_normal$count_by_depth[is.na(size_epi_normal$count_by_depth)] = 0
      # unname is needed otherwise rownames of counts becomes weird
      # set to NA when there is no lesion or normal
      epi_normal_density = unname(ifelse (sizes["EI"]==0, NA, sum(size_epi_normal$count_by_depth/size_epi_normal$area, na.rm=TRUE)))
      
      # dermis normal
      
      dermis_normal = subset(alldat[[s]][[a]], dist2DEJ<max_depth & dermis & region=="intact")
      # add a depth column computed as the floor of dist2DEJ/100
      dermis_normal$depth = floor(dermis_normal$dist2DEJ/100)*100
      count_by_depth = mytable(dermis_normal$depth)
      # add count_by_depth column to size_dermis_normal by match the depth in size_dermis_normal and count_by_depth
      size_dermis_normal$count_by_depth = count_by_depth[size_dermis_normal$depth%.%""]
      size_dermis_normal$count_by_depth[is.na(size_dermis_normal$count_by_depth)]=0
      dermis_normal_density = unname(ifelse (sizes["DI"]==0, NA, sum(size_dermis_normal$count_by_depth/size_dermis_normal$area, na.rm=TRUE)))
      
      # epi lesion
      
      epi_lesion = subset(alldat[[s]][[a]], dist2DEJ<max_depth & epi &    region=="ulcerative")
      # add a depth column computed as the floor of dist2DEJ/100
      epi_lesion$depth = floor(epi_lesion$dist2DEJ/100)*100
      count_by_depth = mytable(epi_lesion$depth)
      # add count_by_depth column to size_epi_lesion by match the depth in size_epi_lesion and count_by_depth
      size_epi_lesion$count_by_depth = count_by_depth[size_epi_lesion$depth%.%""]
      size_epi_lesion$count_by_depth[is.na(size_epi_lesion$count_by_depth)] = 0
      epi_lesion_density = unname(ifelse (sizes["EU"]==0, NA, sum(size_epi_lesion$count_by_depth/size_epi_lesion$area, na.rm=TRUE)))
      
      # dermis lesion
      
      dermis_lesion = subset(alldat[[s]][[a]], dist2DEJ<max_depth & dermis & region=="ulcerative")
      # add a depth column computed as the floor of dist2DEJ/100
      dermis_lesion$depth = floor(dermis_lesion$dist2DEJ/100)*100
      count_by_depth = mytable(dermis_lesion$depth)
      # add count_by_depth column to size_dermis_lesion by match the depth in size_dermis_lesion and count_by_depth
      size_dermis_lesion$count_by_depth = count_by_depth[size_dermis_lesion$depth%.%""]
      size_dermis_lesion$count_by_depth[is.na(size_dermis_lesion$count_by_depth)] = 0
      dermis_lesion_density = unname(ifelse (sizes["DU"]==0, NA, sum(size_dermis_lesion$count_by_depth/size_dermis_lesion$area, na.rm=TRUE)))
      
      # return the counts
      counts3 = cbind(counts3, c(EI=epi_normal_density, DI=dermis_normal_density, EU=epi_lesion_density, DU=dermis_lesion_density))
    
      density.by.depth[[a]] = list(
        EI = cbind(density=size_epi_normal$count_by_depth/size_epi_normal$area, depth=size_epi_normal$depth), 
        DI = cbind(density=size_dermis_normal$count_by_depth/size_dermis_normal$area, depth=size_dermis_normal$depth),
        EU = cbind(density=size_epi_lesion$count_by_depth/size_epi_lesion$area, depth=size_epi_lesion$depth),
        DU = cbind(density=size_dermis_lesion$count_by_depth/size_dermis_lesion$area, depth=size_dermis_lesion$depth)
      )
    }
    mode(counts3) <- "integer"; colnames(counts3)=markers
    
    # add ratio to counts3
    density.ratio.by.depth=list()
    # add the same number of columns to counts and counts2
    counts=cbind(counts, NA, NA, NA, NA, NA, NA)
    counts2=cbind(counts2, NA, NA, NA, NA, NA, NA)
      
    # make density.ratio.by.depth
    for(i in 1:6) {
      if (i==1) {a="Treg"; b="CD3"}
      if (i==2) {a="Treg"; b="CD4T"}
      if (i==3) {a="Treg"; b="CD8T"}
      if (i==4) {a="CD4T"; b="CD3"}
      if (i==5) {a="CD8T"; b="CD3"}
      if (i==6) {a="CD4T"; b="CD8T"}
      
      if(a %in% colnames(counts3) & b %in% colnames(counts3)) {
        counts3 = cbind(counts3, counts3[,a]/counts3[,b])
        colnames(counts3)[ncol(counts3)]=glue("{a}/{b}")
        
        tmp = list(
          EI=cbind(perc = 100*density.by.depth[[a]]$EI[,"density"]/density.by.depth[[b]]$EI[,"density"], depth = density.by.depth[[b]]$EI[,"depth"]),
          DI=cbind(perc = 100*density.by.depth[[a]]$DI[,"density"]/density.by.depth[[b]]$DI[,"density"], depth = density.by.depth[[b]]$DI[,"depth"]),
          EU=cbind(perc = 100*density.by.depth[[a]]$EU[,"density"]/density.by.depth[[b]]$EU[,"density"], depth = density.by.depth[[b]]$EU[,"depth"]),
          DU=cbind(perc = 100*density.by.depth[[a]]$DU[,"density"]/density.by.depth[[b]]$DU[,"density"], depth = density.by.depth[[b]]$DU[,"depth"])
        )
        density.ratio.by.depth[[glue("{a}/{b}")]]$I = rbind(EI=tmp$EI, DI=tmp$DI)
        density.ratio.by.depth[[glue("{a}/{b}")]]$U = rbind(EU=tmp$EU, DU=tmp$DU)
        # adjust depth by 50 so that the depth is in the middle of the bin
        density.ratio.by.depth[[glue("{a}/{b}")]]$I[,"depth"] = density.ratio.by.depth[[glue("{a}/{b}")]]$I[,"depth"] + 50
        density.ratio.by.depth[[glue("{a}/{b}")]]$U[,"depth"] = density.ratio.by.depth[[glue("{a}/{b}")]]$U[,"depth"] + 50
        
      } else {
        counts3 = cbind(counts3, NA)
        colnames(counts3)[ncol(counts3)]=glue("{a}/{b}")
        
      }
    }

    list(sizes=sizes, counts=counts, densities=counts2, norm.counts=counts3, density.ratio.by.depth = density.ratio.by.depth)
    
  })
}) 


# count cells to different depths below skin
count.tables.2 = lapply (all.slides, function(s) {
  # s="LG3_Les" 
  
  markers=markers.ls[[s]]; names(markers)=markers
  
  depths=c(1000, 500, Inf)
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
    
    counts3 = NULL
    density.by.depth = list()
    
    for (a in markers) {
      
      # normal
      
      normal = subset(alldat[[s]][[a]], dist2skin<max_depth & region=="intact")
      # add a depth column computed as the floor of dist2skin/100
      normal$depth = floor(normal$dist2skin/100)*100
      count_by_depth = mytable(normal$depth)
      # add count_by_depth column to size_normal by match the depth in size_normal and count_by_depth
      size_normal$count_by_depth = count_by_depth[size_normal$depth%.%""]
      size_normal$count_by_depth[is.na(size_normal$count_by_depth)] = 0
      normal_density = unname(ifelse (sizes["I"]==0, NA, sum(size_normal$count_by_depth/size_normal$area, na.rm=TRUE)))
      
      # epi lesion
      
      lesion = subset(alldat[[s]][[a]], dist2skin<max_depth & region=="ulcerative")
      # add a depth column computed as the floor of dist2skin/100
      lesion$depth = floor(lesion$dist2skin/100)*100
      count_by_depth = mytable(lesion$depth)
      # add count_by_depth column to size_lesion by match the depth in size_lesion and count_by_depth
      size_lesion$count_by_depth = count_by_depth[size_lesion$depth%.%""]
      size_lesion$count_by_depth[is.na(size_lesion$count_by_depth)] = 0
      lesion_density = unname(ifelse (sizes["U"]==0, NA, sum(size_lesion$count_by_depth/size_lesion$area, na.rm=TRUE)))
      
      # return the counts
      counts3 = cbind(counts3, c(I=normal_density, U=lesion_density))
      
      density.by.depth[[a]] = list(
        I = cbind(density=size_normal$count_by_depth/size_normal$area, depth=size_normal$depth), 
        U = cbind(density=size_lesion$count_by_depth/size_lesion$area, depth=size_lesion$depth)
      )
    }
    mode(counts3) <- "integer"; colnames(counts3)=markers
    
    # add the same number of columns to counts and counts2
    counts=cbind(counts, NA, NA, NA, NA, NA, NA)
    counts2=cbind(counts2, NA, NA, NA, NA, NA, NA)

    # add ratio to counts3
    density.ratio.by.depth=list()
    
    # make density.ratio.by.depth
    for(i in 1:6) {
      if (i==1) {a="Treg"; b="CD3"}
      if (i==2) {a="Treg"; b="CD4T"}
      if (i==3) {a="Treg"; b="CD8T"}
      if (i==4) {a="CD4T"; b="CD3"}
      if (i==5) {a="CD8T"; b="CD3"}
      if (i==6) {a="CD4T"; b="CD8T"}
      
      if(a %in% colnames(counts3) & b %in% colnames(counts3)) {
        counts3 = cbind(counts3, counts3[,a]/counts3[,b])
        colnames(counts3)[ncol(counts3)]=glue("{a}/{b}")
        
        tmp = list(
          I=cbind(perc = 100*density.by.depth[[a]]$I[,"density"]/density.by.depth[[b]]$I[,"density"], depth = density.by.depth[[b]]$I[,"depth"]),
          U=cbind(perc = 100*density.by.depth[[a]]$U[,"density"]/density.by.depth[[b]]$U[,"density"], depth = density.by.depth[[b]]$U[,"depth"])
        )
        density.ratio.by.depth[[glue("{a}/{b}")]]$I = tmp$I
        density.ratio.by.depth[[glue("{a}/{b}")]]$U = tmp$U
        # add 50 to depth because the first depth is 0
        density.ratio.by.depth[[glue("{a}/{b}")]]$I[,"depth"] = density.ratio.by.depth[[glue("{a}/{b}")]]$I[,"depth"] + 50
        density.ratio.by.depth[[glue("{a}/{b}")]]$U[,"depth"] = density.ratio.by.depth[[glue("{a}/{b}")]]$U[,"depth"] + 50
        
      } else {
        counts3 = cbind(counts3, NA)
        colnames(counts3)[ncol(counts3)]=glue("{a}/{b}")
      
      }
      
      
    }
    
    list(sizes=sizes, counts=counts, densities=counts2, norm.counts=counts3, density.ratio.by.depth = density.ratio.by.depth)
    
  })
}) 


# create avg.slides that average normalized densities across slides for a biopsy
# once for adding to count.tables and once for adding to count.tables.2
for (iter in 1:2) {
  out = lapply(biopsies, function(b){
    if (iter==1){
      depths=c(500, 1000, Inf); depths=paste0("DEJ", depths)
    } else {
      depths=c(500, 1000, Inf); depths=paste0("Skin", depths)
    }
    names(depths)=depths
    out=lapply (depths, function(max_depth) { 
      # iter=1; b=biopsies[1]; max_depth="DEJ1000"
      
      # find the union of markers for slides.by.biopsy[[b]]
      markers=c(); for (s in slides.by.biopsy[[b]]) markers=c(markers, markers.ls[[s]])
      markers=sort(unique(markers))
      markers=markers[!endsWith(markers, "_raw")]      # remove markers ending in raw
      names(markers)=markers
      
      tab=sapply(slides.by.biopsy[[b]],simplify="array", function(s){
        # myprint(b, max_depth, s)
        norm.counts=(if (iter==1) count.tables else count.tables.2)[[s]][[max_depth]]$norm.counts
        # add a column of NAs to norm.counts for each marker in markers that is not in colnames(norm.counts)
        missing_markers <- setdiff(markers, colnames(norm.counts))
        if (length(missing_markers) > 0) {
          na_matrix <- matrix(NA, nrow = nrow(norm.counts), ncol = length(missing_markers))
          colnames(na_matrix) <- missing_markers
          norm.counts <- cbind(norm.counts, na_matrix)
        }
        norm.counts[,markers,drop=F]
      })
      out=apply(tab,c(1,2), mean, na.rm=T)
      
      # add ratios Treg/CD3, CD4T/CD3, CD8T/CD3, CD4T/CD8T, Treg/CD4T, Treg/CD8T
      if("Treg" %in% markers) {
        out = cbind(out, 
                    "Treg/CD3" = out[,"Treg"]/out[,"CD3"],
                    "Treg/CD4T" = out[,"Treg"]/out[,"CD4T"],
                    "Treg/CD8T" = out[,"Treg"]/out[,"CD8T"],
                    "CD4T/CD3" = out[,"CD4T"]/out[,"CD3"],
                    "CD8T/CD3" = out[,"CD8T"]/out[,"CD3"],
                    "CD4T/CD8T" = out[,"CD4T"]/out[,"CD8T"]
        )
      }
      
      list(norm.counts=out, sizes=NA, counts=NA, densities=NA)
    })
    
  })
  names(out) = glue("{biopsies}_Les_avg")
  
  if (iter==1){
    count.tables = c(count.tables, out)
  } else {
    count.tables.2 = c(count.tables.2, out)
  }
}



# save the objects except alldat to a file
# alldat is large and only used for plotting histograms
save(list = setdiff(ls(), "alldat"), file = "R/forRmd.RData")
