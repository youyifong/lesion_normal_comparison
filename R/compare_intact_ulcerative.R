################################################################################
# compare DD and CF in one plot
# 8 plots for each of dist2skin and dist2DEJ:
# 4 markers in lesion 
# 4 markers in intact

dir.exists("output/figures") || dir.create("output/figures", recursive = TRUE)

for (dist in dist.labels) {
for (marker_set in c(1,2,3)) {
for (biopsy in c("lesion_slides", "select_slides", biopsies)) {
    
  if (marker_set==1) {
    toplot = key.markers
  } else if (marker_set==2) {
    toplot = c("CD3", cell_types)
  } else {
    toplot = ratio_pairs_labels
  }
  
  if (biopsy=="select_slides") {
    # 13 select slides on one page"
    slides = select.slides
  } else if (biopsy=="lesion_slides") {
    # 9 slides with lesions on one page"
    slides = lesion_slides
  } else {
    slides = slides.by.biopsy[[biopsy]]
  }
      
  mypng(mfrow=c(length(slides),length(toplot)), 
        width=ifelse(marker_set==3,2.3,2)*length(toplot), height=30/13*length(slides), 
        mar=c(4,ifelse(marker_set==3,2,0),2,0), oma=c(1,5,0,0),
        file=glue("output/figures/compare_normal_lesion_{dist}_{biopsy}_set{marker_set}"))
  
  for (s in slides) { 
  for (a in toplot) {
    # dist="dist2skin"; a="CD4T/CD3"; marker_set=3; biopsy="LG"; s="LG3_Les"
  
    if (marker_set==1 | marker_set==2) {
      if (is.null(cell_dat[[s]][[a]])) {
        empty.plot() 
        title(main=ifelse(s==slides[1],a,""))
      } else {
        dat.ls = lapply(tissue.types, function (t) subset(cell_dat[[s]][[a]], region==t)[[dist]])
        # remove empty element
        dat.ls =  dat.ls[sapply(dat.ls, length)>0]
        scaled.hist(dat.ls,
                    scale.factors=lapply(tissue.types, function(t) 
                      subset(tissue_size, biopsy==s & dist_type==dist & region==t)),
                    main=ifelse(s==slides[1],a,""),
                    col=rev(c("red","blue")),
                    legend=names(dat.ls),
                    xlim=xlim.ls[[dist]], xlab=dist %.% " (um)", 
                    ylim=ylim[[dist]], ylab="",
                    yaxt=ifelse(a==toplot[1],"s","n") # only have y axis in column 1
        )
      }      
    } else {
      # market_set 3 is ratio data
      if (dist=="dist2DEJ") {
        dat.tmp = count_dat[[s]]$DEJInf
      } else {
        dat.tmp = count_dat[[s]]$SkinInf
      }
      if (!is.null(dat.tmp$depth_density$ulcerative)) {
        dat.tmp = merge(dat.tmp$depth_density$intact[,c("depth",a)], dat.tmp$depth_density$ulcerative[,c("depth",a)], by="depth", suffixes=c(".I",".U"))
      } else {
        dat.tmp = cbind(dat.tmp$depth_density$intact[,c("depth",a)], NA)
      }
      # convert to percent
      dat.tmp[,2]=100*dat.tmp[,2]; dat.tmp[,3]=100*dat.tmp[,3]; colnames(dat.tmp)[2:3]=c("perc.I","perc.U")
      
      if (all(is.na(dat.tmp[,c("perc.I","perc.U")]))) {
        empty.plot(); title(main=ifelse(s==slides[1],a,""))
      } else {
        mymatplot(dat.tmp$depth, dat.tmp[,c("perc.I","perc.U")],
                  type="l", lty=c(1,1), col=rev(c("red","blue")), lwd=2,
                  main=ifelse(s==slides[1],a,""),
                  legend=tissue.types[!c(all(is.na(dat.tmp$perc.I)), all(is.na(dat.tmp$perc.U)))], legend.x=3,
                  xlab=dist %.% " (um)", xlim=xlim.ls[[dist]],
                  ylab="", ylim=c(0,ylim.ratio[a]), 
                  # yaxt=ifelse(a==toplot[1],"s","n"),
                  bty="n", impute.missing.for.line = F)
      }
    }  
    
    # to divide between epi and dermis
    if (dist == "dist2DEJ") abline(v=0, lty=2)
  
  }
  }

  # leftmost labeling
  tmp=get.slide.name(slides)
  for (i in 1:length(tmp)) {
    # ylab
    mtext(ifelse(marker_set==3, "percent","Cells/100x100 um2"), side = 2, line = 2., outer = T, at = ((length(tmp)-i)*2 + 1)/(length(tmp)*2)+0.00, 
          adj = NA, padj = NA, cex = .7, col = NA, font = NA)
    # slide name
    mtext(tmp[i], side = 2, line = 3.5, outer = T, at = ((length(tmp)-i)*2 + 1)/(length(tmp)*2), 
          adj = NA, padj = NA, cex = NA, col = NA, font = NA)
  }
  
  dev.off()
  
}
}
}

