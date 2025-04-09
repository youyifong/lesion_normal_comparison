################################################################################
# compare DD and CF in one plot
# 8 plots for each of dist2skin and dist2DEJ:
# 4 markers in lesion 
# 4 markers in intact

dir.exists("output/figures") || dir.create("output/figures", recursive = TRUE)

for (dist in dist.labels) {

  mypng(mfrow=c(length(biopsies),length(key.markers)), width=12, height=30, mar=c(4,0,2,0), oma=c(1,5,0,0),
        file=paste0("output/figures/compare_normal_lesion_"%.%dist))
  
  for (b in biopsies) { 
  for (a in key.markers) {
      # dist="dist2skin"; a="CD11c"; b="DD"
      
    dat.ls = lapply(tissue.types, function (t) subset(alldat[[b]][[a]], region==t)[[dist]])
    
    scaled.hist(dat.ls,
         scale.factors=lapply(tissue.types, function(t) subset(tissue_size_by_depth_and_region, biopsy==b & dist_type==dist & region==t)),
         xlim=xlim.ls[[dist]], ylim=ylim[[dist]],
         xlab=dist %.% " (um)", ylab="",
         main=ifelse(b==biopsies[1],a,""),
         col=c("red","blue"),
         legend=if (a==key.markers[length(key.markers)]) tissue.types else NULL,
         yaxt=ifelse(a==key.markers[1],"s","n") # only have y axis in column 1
    )
    
    # to divide between epi and dermis
    if (dist == "dist2DEJ") abline(v=0, lty=2)
      
  }
  }

  # leftmost labeling
  tmp=get.slide.name(biopsies)
  for (i in 1:length(tmp)) {
    # ylab
    mtext("Cells/100x100 um2", side = 2, line = 2., outer = T, at = ((length(tmp)-i)*2 + 1)/(length(tmp)*2)+0.00, 
          adj = NA, padj = NA, cex = .7, col = NA, font = NA)
    mtext(tmp[i], side = 2, line = 3.5, outer = T, at = ((length(tmp)-i)*2 + 1)/(length(tmp)*2), 
          adj = NA, padj = NA, cex = NA, col = NA, font = NA)
  }
  
  dev.off()
  
}

