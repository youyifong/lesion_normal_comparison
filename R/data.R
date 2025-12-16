library(kyotil)
library(glue)

{ # slide level data
data.dir = "../Chipcytometry_images_from32bit"
key.markers=c("CD3", "CD11c", "CD123", "CD15", "CD56", "CD68")
cell_types=c("CD8T", "CD4T", "Treg")
get.slide.name=function(ss)  sapply(ss, function (s) sub(".*?/(.*)", "\\1", s)) # "DD/DD1_Les" -> DD1_Les
get.initial=function(ss)  sapply(ss, function (s) sub("^([^0-9_]+).*", "\\1", s) ) # "DD1_Les" -> "DD"

depths = c(1000, 500, Inf)
modes = c("DEJ", "Skin"); names(modes)=modes
depth_vals = c("DEJ1000", "DEJ500", "DEJInf", "Skin1000", "Skin500", "SkinInf"); names(depth_vals)=depth_vals

ratio_pairs = list(
  c("Treg", "CD3"),
  c("Treg", "CD4T"),
  c("Treg", "CD8T"),
  c("CD4T", "CD3"),
  c("CD8T", "CD3"),
  c("CD8T", "CD4T")
)

ratio_pairs_labels=sapply(ratio_pairs, function (x) glue("{x[1]}/{x[2]}"))
ylim.ratio = c("Treg/CD3"=50, "Treg/CD4T"=50, "Treg/CD8T"=50, "CD4T/CD3"=100, "CD8T/CD3"=100, "CD8T/CD4T"=200)

# in tissue_size, the depth column is the starting point of the bin.
tissue_size=read.csv(paste0(data.dir, "/tissue_size_by_depth_and_region.csv"))
tissue_size$biopsy = get.slide.name(tissue_size$biopsy)
all.slides = sort(unique(tissue_size$biopsy)); names(all.slides)=all.slides

# for each biopsy, determine if it has no intact region
tmp = mysapply(all.slides, function(s) {
  tmp1 = tissue_size[tissue_size$biopsy==s & tissue_size$region=="intact" & tissue_size$area>0, ]
  tmp2 = tissue_size[tissue_size$biopsy==s & tissue_size$region=="ulcerative" & tissue_size$area>0, ]
  c(no_intact = nrow(tmp1)==0, no_lesion = nrow(tmp2)==0)
}); tmp
lesion_slides = rownames(tmp)[!tmp[,"no_lesion"]]

markers.ls = lapply(all.slides, function(s) {  # s = 'JT_Les'
  tmp = sort(list.files(paste0(data.dir, "/", get.initial(s), "/", s, "/TSP_stitched"), pattern=".*_d2b_regmem.csv"))
  out = sapply(tmp, function(n) sub("_HIL.*.csv", "", sub("_m.*.csv", "", sub(s%.%"_","",n))))
  # move CD4T, CD8T, Treg to the end. by using setdiff, replicates are also lost
  out = c(setdiff(out, cell_types), setdiff(out, setdiff(out, cell_types)))
  names(out)=out; out
})

# remove slides with no markers
all.slides = all.slides[sapply(markers.ls, function(x) length(x)>0)]
# choose one slide for each biopsy
select.slides = c("AS1_Les", "CF1_Les", "CJ1_Les", "DD1_Les",  "EB_Les", "JK1_Les", "JM1_Les", "JT_Les",  "LD_Les", "LG2_Les",  "LM_Les", "LS2_Les", "SK2_Les"); names(select.slides)=get.slide.name(select.slides)
# CF1_Les has no CD68. We used to use CF5_Les CD68 to replace it, but that is not a very good soln. For one thing, CD5_Les has no ulcerative part.
#stopifnot(all(sapply(select.slides, function(b) length(intersect(markers.ls[[b]], key.markers)))==6))# make sure all select.slides have all key markers

biopsies=sapply(select.slides, function(x) get.initial(x)); names(biopsies)=biopsies
avg_slides = glue("{biopsies}_Les_avg")# average norm.counts of all slides from a single biopsy
slides.by.biopsy = lapply(biopsies, function(b) all.slides[startsWith(all.slides,b)])

severity = read.csv("../HSV406/Lesion_HSV_rates_04082025.csv")
severity$initial = sapply(severity$PID, function(x) sub("[0-9]+$", "", x))# add a column for the initial, which are the letter in front of the arabic number

region.labels=c(intact="intact", lesion="ulcerative")
region.labels.title=c(intact="Intact", lesion="Ulcerative")
dist.labels=c("dist2skin", "dist2DEJ"); names(dist.labels)=dist.labels
tissue.types=c("intact","ulcerative"); names(tissue.types)=tissue.types
xlim.ls=list(dist2DEJ=c(-750,1750), dist2skin=c(0,2500))
ylim=list(dist2DEJ=c(0,100), dist2skin=c(0,100))
red=rgb(red=1, blue=0, green=0, alpha=1)
blue=rgb(red=0, blue=1, green=0, alpha=1)

# create output/csv if not exists
if (!dir.exists("output/csv")) dir.create("output/csv", recursive = TRUE)
}



################################################################################
# cell-level data

get.dat=function(s, marker) {
# s="CF1_Les"; dist="dist2skin"; marker="CD68"
    
  # in the beginning, we use CF5 CD68 instead of CF1 CD68 due to quality
  # later, we decide not to do this
  # if (s == "CF1_Les" && marker == "CD68") s = "CF5_Les" 
  
  # human edited and augmented phenotypic
  fn = glue("{data.dir}/{get.initial(s)}/{s}/TSP_stitched/{s}_{marker}_HIL_aug_d2b_regmem.csv")
  if (!file.exists(fn)) {
    # human edited 
    fn = glue("{data.dir}/{get.initial(s)}/{s}/TSP_stitched/{s}_{marker}_HIL_d2b_regmem.csv")
    if (!file.exists(fn)) {
      fn = glue("{data.dir}/{get.initial(s)}/{s}/TSP_stitched/{s}_{marker}_m_c7_d2b_regmem.csv")
      if (!file.exists(fn)) {
        # T cell subset files are named differently
        fn = glue("{data.dir}/{get.initial(s)}/{s}/TSP_stitched/{s}_{marker}_masks_d2b_regmem.csv")
        if (!file.exists(fn)) stop("error: fn not exists: "%.%fn) #return (NULL)
      }
    }
  } else {
    print(fn)
  }
  dat=read.csv(fn)
  
  # find the column with epi_region and intact in the name
  normal_epi_col = grep("epi_region", grep("normal", colnames(dat), value=TRUE), value=TRUE)
  # find the column with dermis_region and intact in the name
  normal_dermis_col = grep("dermis_region", grep("normal", colnames(dat), value=TRUE), value=TRUE)
  # find the column with epi_region and lesion in the name
  lesion_epi_col = grep("epi_region", grep("lesion", colnames(dat), value=TRUE), value=TRUE)
  # find the column with dermis_region and lesion in the name
  lesion_dermis_col = grep("dermis_region", grep("lesion", colnames(dat), value=TRUE), value=TRUE)
  
  # define column intact.epi as any of the normal_epi_col columns equals "True"
  dat$intact.epi = rowSums(dat[, normal_epi_col, drop=F] == "True", na.rm=TRUE) > 0
  # define column intact.dermis as any of the normal_dermis_col columns equals "True"
  dat$intact.dermis = rowSums(dat[, normal_dermis_col, drop=F] == "True", na.rm=TRUE) > 0
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
  
  
  dat$epi = dat$intact.epi | dat$lesion.epi
  dat$dermis = dat$intact.dermis | dat$lesion.dermis

  dat$region = NA
  dat$region[dat$intact.epi | dat$intact.dermis] = "intact"
  dat$region[dat$lesion.epi | dat$lesion.dermis] = "ulcerative"

  # each pixel is 1/2 micron, negative dist in epi
  dat$dist2DEJ = dat$dist2DEJ * ifelse(dat$epi, -1, 1)/2
  dat$dist2skin = dat$dist2skin/2
  dat
}


cell_dat =  lapply (all.slides, function(s) {
          lapply (markers.ls[[s]], function (marker) {
            myprint(s, marker)
            get.dat(s, marker)
          })
          })



################################################################################
## count cells in epi/dermis x intact/lesion
# count_dat comprises a list of slides
# count_dat$AS1_Les comprises DEJ500, DEJ1000, DEJInf, Skin500, Skin1000, SkinInf
# count_dat$AS1_Les$DEJ500 comprises sizes, counts, densities, norm.counts, depth_density

count_dat = lapply(all.slides, function(s) {
  # s="SK1_Les"; depth_val="DEJ1000"
  markers = markers.ls[[s]]; names(markers) = markers
  
  # region_stats is region (EI, EU etc.)-level data
  region_stats = lapply(depth_vals, function(depth_val) {
    myprint(s, depth_val)
    
    # split depth_val into mode and max_depth
    if (startsWith(depth_val, "DEJ")) {
      mode = "DEJ"
      max_depth = as.numeric(sub("DEJ", "", depth_val))
    } else if (startsWith(depth_val, "Skin")) {
      mode = "Skin"
      max_depth = as.numeric(sub("Skin", "", depth_val))
    } else {
      stop("Invalid depth_val")
    }
    
    dist_col = if (mode == "DEJ") "dist2DEJ" else "dist2skin"

    counts = sapply(markers, function(a) {
      if (mode=="DEJ") c(
        EI = nrow(subset(cell_dat[[s]][[a]], get(dist_col) < max_depth & epi    & region == "intact")),
        DI = nrow(subset(cell_dat[[s]][[a]], get(dist_col) < max_depth & dermis & region == "intact")),
        EU = nrow(subset(cell_dat[[s]][[a]], get(dist_col) < max_depth & epi    & region == "ulcerative")),
        DU = nrow(subset(cell_dat[[s]][[a]], get(dist_col) < max_depth & dermis & region == "ulcerative"))
      ) else c(
        I  = nrow(subset(cell_dat[[s]][[a]], get(dist_col) < max_depth & region == "intact")),
        U  = nrow(subset(cell_dat[[s]][[a]], get(dist_col) < max_depth & region == "ulcerative"))
      )
    })
    
    dat_depth = tissue_size[tissue_size$biopsy == s & tissue_size$dist_type == dist_col & tissue_size$depth < max_depth, ]
    dat_depth$area = dat_depth$area / 100 # convert size from 100um x 100um to 1mm x 1mm
    
    sizes = with(dat_depth, 
      if (mode=="DEJ") c(
        EI = sum(area[depth < 0 & region == "intact"]),
        DI = sum(area[depth >= 0 & region == "intact"]),
        EU = sum(area[depth < 0 & region == "ulcerative"]),
        DU = sum(area[depth >= 0 & region == "ulcerative"])
      ) else c(
        I  = sum(area[region == "intact"]),
        U  = sum(area[region == "ulcerative"])
      )
    )
    
    sizes[sizes==0] = NA; sizes 
    counts[is.na(sizes),]=NA; counts
    densities = counts / sizes; mode(densities) <- "integer"; densities
    
    depth_density = lapply(tissue.types, function (t) {
      out = dat_depth[dat_depth$region == t, ]
      if (sum(out$area)==0) return (NULL) else {
        # remove rows with 0 area
        out = subset(out, area>0)
        # if depth are not contiguous, raise an error
        if (any(diff(out$depth) != 100)) {
          print(out)
          stop("depths are not contiguous")
        }
        for (a in markers) {
          tmp = subset(cell_dat[[s]][[a]], get(dist_col) < max_depth & region == t)
          tmp$depth = floor(tmp[[dist_col]] / 100) * 100
          count_by_depth = mytable(tmp$depth)
          out[[a%.%".count"]] = count_by_depth[out$depth %.% ""]
          out[[a%.%".count"]][is.na(out[[a%.%".count"]])] = 0 # NA is introduced in count_by_depth_normal when there are no cells at a depth. 0 is more interpretable
          out[[a%.%".density"]] = out[[a%.%".count"]] / out$area
        }
      }
      
      # Add density ratio at each depth
      for (i in seq_along(ratio_pairs)) {
        m1 <- ratio_pairs[[i]][1]
        m2 <- ratio_pairs[[i]][2]
        if (m1%.%".density" %in% colnames(out) && m2%.%".density" %in% colnames(out)) {
          out = cbind(out, c(out[, m1%.%".density"] / out[, m2%.%".density"]))
        } else {
          out = cbind(out, NA)
        }
        colnames(out)[ncol(out)] = glue("{m1}/{m2}")
      }
      
      out
    }); depth_density
    
    # if depth_val ends with Inf, save each non-null element of depth_density to a csv file
    if (endsWith(depth_val, "Inf")) {
      for (region in names(depth_density)) 
        if (!is.null(depth_density[[region]])) 
          write.csv(depth_density[[region]], file = paste0("output/csv/depth_density_", s, "_", region, "_", depth_val, ".csv"), row.names = FALSE)
    }

    
    # sum over depth_density to get norm.counts 
    norm.counts = list()
    if (mode=="DEJ") {
      # DEJ mode, count by EI,DI,EU,DU
      for (t in c("intact","ulcerative")) {
        if (is.null(depth_density[[t]])) E<-D<-NA else {
          # get density columns in epidermis
          E = depth_density[[t]][depth_density[[t]]$depth<0, grep("density$", colnames(depth_density[[t]]), value=TRUE)]
          if (nrow(E)==0) E=NA else {
            colnames(E) = sub("\\.density$", "", colnames(E))
            E = colSums(E)
          }
          
          # get density columns in dermis
          D = depth_density[[t]][depth_density[[t]]$depth>=0, grep("density$", colnames(depth_density[[t]]), value=TRUE)]
          if (nrow(D)==0) D=NA else {
            colnames(D) = sub("\\.density$", "", colnames(D))
            D = colSums(D)
          }
        }
        norm.counts = append(norm.counts, list(E, D))
      }
      norm.counts=do.call(rbind, norm.counts); rownames(norm.counts) = c("EI","DI","EU","DU")
      
    } else {
      # skin mode, count I and U
      for (t in c("intact","ulcerative")) {
        if (is.null(depth_density[[t]])) ED<-NA else {
          # get density columns 
          ED = depth_density[[t]][, grep("density$", colnames(depth_density[[t]]), value=TRUE)]
          if (nrow(ED)==0) ED=NA else {
            colnames(ED) = sub("\\.density$", "", colnames(ED))
            ED = colSums(ED)
          }
        }
        norm.counts = append(norm.counts, list(ED))
      }
      norm.counts=do.call(rbind, norm.counts); rownames(norm.counts) = c("I","U")
      
    }
    # divide by 100 to get the number of cells in a 100um wide vertical strip
    norm.counts = norm.counts/100
    norm.counts
    
    
    # Add stats for ratios to norm.counts
    for (i in seq_along(ratio_pairs)) {
      m1 <- ratio_pairs[[i]][1]
      m2 <- ratio_pairs[[i]][2]
      if (m1 %in% colnames(norm.counts) && m2 %in% colnames(norm.counts)) {
        norm.counts = cbind(norm.counts, norm.counts[, m1] / norm.counts[, m2])
      } else {
        norm.counts = cbind(norm.counts, NA)
      }
      colnames(norm.counts)[ncol(norm.counts)] = glue("{m1}/{m2}")
      
      # add a column to counts and densities so that they can be rbind with norm.counts
      counts = cbind(counts, NA)
      densities = cbind(densities, NA)
    }
    
    
    list(counts=counts, densities=densities, sizes=sizes, depth_density=depth_density, norm.counts=norm.counts)
  })
  
}) # slides



# add avg.slides 
out = lapply(biopsies, function(b){
      lapply(depth_vals, function(max_depth) { 
      # b="AS"; max_depth="DEJ1000"; s="AS1_Les"
  
  slides = slides.by.biopsy[[b]]
  
  # find the union of markers for slides
  markers=c(); for (s in slides) markers=c(markers, markers.ls[[s]])
  markers=sort(unique(markers))
  markers=markers[!endsWith(markers, "_raw")]      # remove markers ending in raw
  names(markers)=markers
  
  # add a column of NAs to norm.counts for each marker in markers that is not in colnames(norm.counts)
  tab=sapply(slides,simplify="array", function(s){
    norm.counts=count_dat[[s]][[max_depth]]$norm.counts
    missing_markers <- setdiff(markers, colnames(norm.counts))
    if (length(missing_markers) > 0) {
      na_matrix <- matrix(NA, nrow = nrow(norm.counts), ncol = length(missing_markers))
      colnames(na_matrix) <- missing_markers
      norm.counts <- cbind(norm.counts, na_matrix)
    }
    norm.counts[,markers,drop=F]
  })
  out=apply(tab,c(1,2), mean, na.rm=T)
  
  # add ratios 
  # ratios for combined slides := sum 1/sum 2
  # make sure only use slides that have both markers in the ratio
  for (i in seq_along(ratio_pairs)) {
    m1 <- ratio_pairs[[i]][1]
    m2 <- ratio_pairs[[i]][2]
    good.slides = slides[sapply(slides, function(s) m1 %in% markers.ls[[s]] && m2 %in% markers.ls[[s]])]
    tab1=sapply(good.slides,simplify="array", function(s) count_dat[[s]][[max_depth]]$norm.counts[,c(m1,m2)])
    res=apply(tab1,c(1,2), sum)
    out = cbind(out, res[,1]/res[,2] )
    colnames(out)[ncol(out)] = glue("{m1}/{m2}")
  }
  
  
  list(norm.counts=out, sizes=NA, counts=NA, densities=NA) # no depth density for average slide
})
})
names(out) = glue("{biopsies}_Les_avg")

count_dat = c(count_dat, out)




# save the objects except cell_dat to a file
# cell_dat is large and only used for plotting histograms
save(list = setdiff(ls(), "cell_dat"), file = "R/forRmd.RData")
