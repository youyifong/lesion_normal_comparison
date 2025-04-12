################################################################################
# compare DD and CF in one plot
# 8 plots for each of dist2skin and dist2DEJ:
# 4 markers in lesion 
# 4 markers in intact

dir.exists("output/figures") || dir.create("output/figures", recursive = TRUE)

# 13 representative slides on one page
for (dist in dist.labels) {

  mypng(mfrow=c(length(select.slides),length(key.markers)), width=12, height=30, mar=c(4,0,2,0), oma=c(1,5,0,0),
        file=paste0("output/figures/compare_normal_lesion_"%.%dist))
  
  for (s in select.slides) { 
  for (a in key.markers) {
      # dist="dist2skin"; a="CD11c"; s="DD"
      
    dat.ls = lapply(tissue.types, function (t) subset(alldat[[s]][[a]], region==t)[[dist]])
    
    scaled.hist(dat.ls,
         scale.factors=lapply(tissue.types, function(t) subset(tissue_size_by_depth_and_region, biopsy==get.initial(s)%.%"/"%.%s & dist_type==dist & region==t)),
         xlim=xlim.ls[[dist]], ylim=ylim[[dist]],
         xlab=dist %.% " (um)", ylab="",
         main=ifelse(s==select.slides[1],a,""),
         col=c("red","blue"),
         legend=if (a==key.markers[length(key.markers)]) tissue.types else NULL,
         yaxt=ifelse(a==key.markers[1],"s","n") # only have y axis in column 1
    )
    
    # to divide between epi and dermis
    if (dist == "dist2DEJ") abline(v=0, lty=2)
      
  }
  }

  # leftmost labeling
  tmp=get.slide.name(select.slides)
  for (i in 1:length(tmp)) {
    # ylab
    mtext("Cells/100x100 um2", side = 2, line = 2., outer = T, at = ((length(tmp)-i)*2 + 1)/(length(tmp)*2)+0.00, 
          adj = NA, padj = NA, cex = .7, col = NA, font = NA)
    # slide name
    mtext(tmp[i], side = 2, line = 3.5, outer = T, at = ((length(tmp)-i)*2 + 1)/(length(tmp)*2), 
          adj = NA, padj = NA, cex = NA, col = NA, font = NA)
  }
  
  dev.off()
  
}



# all slides for one biopsy on one page
for (dist in dist.labels) {
for (biopsy in biopsies) {
  
  slides = slides.by.biopsy[[biopsy]]
  mypng(mfrow=c(length(slides),length(key.markers)), width=12, height=30/13*length(slides), mar=c(4,0,2,0), oma=c(1,5,0,0),
      file=paste0("output/figures/compare_normal_lesion_",dist,"_",biopsy))
  
    for (s in slides) {
    for (a in key.markers) {
    # dist="dist2skin"; a="CD11c"; s="DD"
      if (is.null(alldat[[s]][[a]])) empty.plot() else {
        dat.ls = lapply(tissue.types, function (t) subset(alldat[[s]][[a]], region==t)[[dist]])
        scaled.hist(dat.ls,
                    scale.factors=lapply(tissue.types, function(t) 
                      subset(tissue_size_by_depth_and_region, biopsy==get.initial(s)%.%"/"%.%s & dist_type==dist & region==t)),
                    xlim=xlim.ls[[dist]], ylim=ylim[[dist]],
                    xlab=dist %.% " (um)", ylab="",
                    main=ifelse(s==slides[1],a,""),
                    col=c("red","blue"),
                    legend=if (a==key.markers[length(key.markers)]) tissue.types else NULL,
                    yaxt=ifelse(a==key.markers[1],"s","n") # only have y axis in column 1
        )
        if (dist == "dist2DEJ") abline(v=0, lty=2) # to divide between epi and dermis
      }
    }
    }
    
    # leftmost labeling
    tmp=get.slide.name(slides)
    for (i in 1:length(tmp)) {
      # ylab
      mtext("Cells/100x100 um2", side = 2, line = 2., outer = T, at = ((length(tmp)-i)*2 + 1)/(length(tmp)*2)+0.00, 
            adj = NA, padj = NA, cex = .7, col = NA, font = NA)
      # slide name
      mtext(tmp[i], side = 2, line = 3.5, outer = T, at = ((length(tmp)-i)*2 + 1)/(length(tmp)*2), 
            adj = NA, padj = NA, cex = NA, col = NA, font = NA)
    }
  
  dev.off()
  
}
}

