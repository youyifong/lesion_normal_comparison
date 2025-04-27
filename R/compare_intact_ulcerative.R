################################################################################
# compare DD and CF in one plot
# 8 plots for each of dist2skin and dist2DEJ:
# 4 markers in lesion 
# 4 markers in intact

dir.exists("output/figures") || dir.create("output/figures", recursive = TRUE)

for (dist in dist.labels) {
for (marker_set in c(1,2,3)) {
for (biopsy in c("select_slides", biopsies)) {
    
  if (marker_set==1) {
    toplot = key.markers
  } else if (marker_set==2) {
    toplot = c("CD3", cell_types)
  } else {
    toplot = c("Treg/CD3", "Treg/CD4T", "Treg/CD8T", "CD4T/CD3", "CD8T/CD3", "CD4T/CD8T")
    ylim.ratio = c("Treg/CD3"=50, "Treg/CD4T"=50, "Treg/CD8T"=50, "CD4T/CD3"=100, "CD8T/CD3"=100, "CD4T/CD8T"=500)
  }
  
  if (biopsy=="select_slides") {
    # 13 select slides on one page"
    slides = select.slides
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
      if (is.null(alldat[[s]][[a]])) {
        empty.plot() 
        title(main=ifelse(s==slides[1],a,""))
      } else {
        dat.ls = lapply(tissue.types, function (t) subset(alldat[[s]][[a]], region==t)[[dist]])
        # remove empty element
        dat.ls =  dat.ls[sapply(dat.ls, length)>0]
        scaled.hist(dat.ls,
                    scale.factors=lapply(tissue.types, function(t) 
                      subset(tissue_size, biopsy==s & dist_type==dist & region==t)),
                    main=ifelse(s==slides[1],a,""),
                    col=c("red","blue"),
                    legend=names(dat.ls),
                    xlim=xlim.ls[[dist]], xlab=dist %.% " (um)", 
                    ylim=ylim[[dist]], ylab="",
                    yaxt=ifelse(a==toplot[1],"s","n") # only have y axis in column 1
        )
      }      
    } else {
      # market_set 3 is ratio data
      if (is.null(count.tables[[s]]$DEJInf$density.ratio.by.depth[[a]])) {
        empty.plot() 
        title(main=ifelse(s==slides[1],a,""))
      } else {
        if (dist=="dist2DEJ") {
          dat.tmp = count.tables[[s]]$DEJInf$density.ratio.by.depth[[a]]
        } else {
          dat.tmp = count.tables.2[[s]]$SkinInf$density.ratio.by.depth[[a]]
        }
        # merge dat.tmp$I and dat.tmp$U on depth
        dat.tmp = merge(dat.tmp$I, dat.tmp$U, by="depth", suffixes=c(".I",".U"))
        mymatplot(dat.tmp$depth, dat.tmp[,c("perc.I","perc.U")],
                  type="l", lty=c(1,1), col=c("red","blue"), lwd=2,
                  main=ifelse(s==slides[1],a,""),
                  legend=tissue.types[!c(all(is.na(dat.tmp$perc.I)), all(is.na(dat.tmp$perc.U)))], legend.x=3,
                  xlab=dist %.% " (um)", xlim=xlim.ls[[dist]],
                  ylab="", ylim=c(0,ylim.ratio[a]), 
                  # yaxt=ifelse(a==toplot[1],"s","n"),
                  bty="n", impute.missing.for.line = F
        )
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

