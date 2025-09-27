################################################################################
# compare several markers in one plot
# 4 plots for each of dist2skin and dist2DEJ:
# ulcerative/intact x DD/CF

dir.exists("output/figures") || dir.create("output/figures", recursive = TRUE)

get.peak = function(x_start, freq) {
  x_start[which.max(freq)] + (x_start[2]-x_start[1])/2
}

for (marker_set in c(1,2,3)) {
for (biopsy in c("select_slides", biopsies)) {
    
  if (marker_set==1) {
    markers = key.markers
  } else if (marker_set==2) {
    markers = c("CD3", cell_types)
  } else {
    markers = ratio_pairs_labels
  }
  # use all markers
  # markers=markers.ls[[s]]; names(markers)=markers
  # use a subset of mostly shared markers
  # markers = c("CD3", "CD15", "CD11c", "CD123", "CD4", "CD8", if(!s%in%c("LG/LG1_Les","SK/SK1_Les")) "CD68")
  names(markers)=markers
  
  if (biopsy=="select_slides") {
    # 13 select slides on one page"
    slides = select.slides
  } else {
    slides = slides.by.biopsy[[biopsy]]
  }
  
  mypng(mfrow=c(length(slides),4), file=glue("output/figures/compare_markers_{biopsy}_set{marker_set}"), 
      width=8*ifelse(marker_set==3,1.15,1), height=4/2*length(slides), mar=c(4,ifelse(marker_set==3,2,0),2,0), oma=c(1,5,0,0))

  dist.bt.peaks = 
  lapply (slides, function(s) {
  sapply (tissue.types, function(tissue.type) { 
  sapply (dist.labels, simplify="array", function(dist) {
  # tissue.type="ulcerative"; dist="dist2skin"; s="CF5_Les"; marker_set=1
    myprint(tissue.type, dist, s, marker_set)
  
    if (marker_set==1 | marker_set==2) {
      dat.ls = lapply(markers, function (a) if (is.null(cell_dat[[s]][[a]])) NULL else subset(cell_dat[[s]][[a]], region==tissue.type)[[dist]])
      dat.ls = dat.ls[sapply(dat.ls, function(x) !is.null(x))]

      if(all(sapply(dat.ls,length)==0)) {
        empty.plot(); title(main=ifelse(s==slides[1],tissue.type,""))
      } else {
        # not all markers may be present
        markers.1 = names(dat.ls)
        # print(seq_along(markers)[match(markers.1, markers)])
        # print(markers.1)
        # print(subset(tissue_size, biopsy==s & dist_type==dist & region==tissue.type))
        # print(ylim[[dist]])
        all.freq = scaled.hist(dat.ls,
           scale.factors=lapply(markers.1, function(a) subset(tissue_size, biopsy==s & dist_type==dist & region==tissue.type)),
           xlim=xlim.ls[[dist]], ylim=ylim[[dist]],
           xlab=dist %.% " (um)", ylab="Number of cells/100x100 um2",
           main=ifelse(s==slides[1],tissue.type,""),
           col=seq_along(markers)[match(markers.1, markers)],
           cex.legend=0.7,
           yaxt=ifelse(tissue.types[1]==tissue.type & dist.labels[1]==dist,"s","n"), # only have y axis in column 1
           legend=markers.1)
      }
      
    } else {
      # market_set 3 is ratio data
      if (dist=="dist2DEJ") {
        dat.tmp = count_dat[[s]]$DEJInf
      } else {
        dat.tmp = count_dat[[s]]$SkinInf
      }
      dat.tmp = dat.tmp$depth_density[[tissue.type]][,c("depth",markers)]
      # convert to percent
      dat.tmp[,-1]=100*dat.tmp[,-1]
      
      if (all(is.na(dat.tmp[,-1]))) {
        empty.plot(); title(main=ifelse(s==slides[1],tissue.type,""))
      } else {
        mymatplot(dat.tmp$depth, dat.tmp[,-1],
                  type="l", lty=1, 
                  col=seq_along(markers), lwd=2,
                  main=ifelse(s==slides[1],tissue.type,""),
                  make.legend=F, 
                  xlab=dist %.% " (um)", xlim=xlim.ls[[dist]],
                  ylab="", ylim=c(0,max(ylim.ratio)), 
                  bty="n", impute.missing.for.line = F)
        # only show legend for markers that are present
        kp=!apply(dat.tmp[,-1],2,function(x) all(is.na(x)))
        mylegend(x=3, legend=markers[kp], lty=1, col=seq_along(markers)[kp], y.intersp=1, text.width=700, cex=0.7)
        
      }
      
    }
    # dividing line between epi and dermis
    if (dist == "dist2DEJ") abline(v=0, lty=2)
    
    # return distance between CD3 and CD15 peaks
    # get.peak(all.freq[,"depth_start"], all.freq[,"CD3"]) - get.peak(all.freq[,"depth_start"], all.freq[,"CD15"])
  
  }) 
  })
  })
  
  # leftmost labeling
  tmp=get.slide.name(slides)
  for (i in 1:length(tmp)) {
    # ylab
    mtext(ifelse(marker_set==3, "percent","Cells/100x100 um2"), side = 2, line = 2., outer = T, at = ((length(tmp)-i)*2 + 1)/(length(tmp)*2)+0., 
          adj = NA, padj = NA, cex = .7, col = NA, font = NA)
    mtext(tmp[i], side = 2, line = 3.5, outer = T, at = ((length(tmp)-i)*2 + 1)/(length(tmp)*2), 
          adj = NA, padj = NA, cex = NA, col = NA, font = NA)
  }

  dev.off()

}
  
}
  
