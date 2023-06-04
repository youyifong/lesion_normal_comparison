################################################################################
# compare different biopsies in one plot
# CF, DD and JM have both ulcerative and intact regions, others just have intact regions
# separately for dist2skin and dist2DEJ
# separately for ulcerative and intact


# CF, DD and JM only

for (dist in dist.labels) { # dist2skin, dist2DEJ
  xlim=xlim.ls[[dist]]
  myprint(dist)
  
  mypng(mfcol=c(2,length(markers.1)), width=10/4*length(markers.1), height=5, mar=c(4,4,2,1), oma=c(1,3,0,0),
        file=paste0("figures/compare_DD_CF_JM_",dist))
  
  for (a in markers.1) {
    for (tissue.type in tissue.types) { # intact, ulcerative
      # dist="dist2DEJ"; a="CD11c"; tissue.type="intact"
      
      dat.ls = lapply(biopsies.1, function (b) subset(alldat[[b]][[a]], region==tissue.type)[[dist]])
      
      scale.factors=sapply(biopsies.1, function(b) unname(100/tissue.widths[[b]][[tissue.type]]))
      
      smoothed.scaled.hist(dat.ls, bin_width=100, scale.factors=scale.factors,
                           xlim=xlim.ls[[dist]], ylim=ylim,
                           xlab=dist %.% " (um)", ylab="Number of cells/100x100 um2", 
                           main=ifelse(tissue.type==tissue.types[1],a,""),
                           col=1:length(biopsies.1),
                           legend=biopsies.1
      )
      
      # dividing line between epi and dermis
      if (dist == "dist2DEJ") abline(v=0, lty=2)
      
    }
  }
  
  # leftmost labeling
  tmp=tissue.types
  for (i in 1:length(tmp)) {
    mtext(tmp[i], side = 2, line = 1.5, outer = T, at = (i*2-1)/(length(tmp)*2), 
          adj = NA, padj = NA, cex = NA, col = NA, font = NA)
  }
  
  dev.off()
  
}


# All biopsies

for (dist in dist.labels) { # dist2skin, dist2DEJ
  xlim=xlim.ls[[dist]]
  myprint(dist)

  mypng(mfcol=c(1,4), width=10, height=2.6, mar=c(4,4,2,1), oma=c(1,3,0,0),
        file=paste0("figures/compare_biopsies_intact_",dist))
  
  for (a in markers.1) {
    tissue.type = "intact"
    
    # dist="dist2skin"; a="CD15"; tissue.type="intact"
    
    dat.ls = lapply(biopsies, function (b) subset(alldat[[b]][[a]], region=="intact")[[dist]])
    
    scale.factors=sapply(biopsies, function(b) unname(100/tissue.widths[[b]][tissue.type]))

    smoothed.scaled.hist(dat.ls, bin_width=100, scale.factors=scale.factors,
                         xlim=xlim.ls[[dist]], ylim=ylim,
                         xlab=dist %.% " (um)", ylab="Number of cells/100x100 um2", 
                         main=a,
                         col=1:length(biopsies),
                         legend=biopsies
    )
    
    # to divide between epi and dermis
    if (dist == "dist2DEJ") abline(v=0, lty=2)
    
  }
  
  tmp="Intact"
  for (i in 1:length(tmp)) {
    mtext(tmp[i], side = 2, line = 1.5, outer = T, at = ((length(tmp)-i)*2 + 1)/(length(tmp)*2), 
          adj = NA, padj = NA, cex = NA, col = NA, font = NA)
  }
  
  dev.off()
  
}



