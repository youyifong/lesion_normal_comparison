################################################################################
# compare DD and CF in one plot
# 8 plots for each of dist2skin and dist2DEJ:
# 4 markers in lesion 
# 4 markers in intact

for (dist in dist.labels) {

  mypng(mfrow=c(2,4), width=10, height=5, mar=c(4,4,2,1), oma=c(1,3,0,0),
        file=paste0("figures/compare_normal_lesion_"%.%dist))
  
  for (s in sample.labels) { 
  for (a in markers) {
      # dist="dist2skin"; a="CD3"; s="DD"
      
      dat.ls = alldat[[s]][[dist]][[a]] 
      dat.ls = lapply(tissue.types, function (t) dat.ls[dat.ls[[t]],"dist"])
      
      scale.factors=sapply(tissue.types, function(tissue.type) unname(100/tissue.widths[[tissue.type]][s]))
          
      smoothed.scaled.hist(dat.ls, bin_width=100, scale.factors=scale.factors,
                           xlim=xlim.ls[[dist]], ylim=c(0,150),
                           xlab=dist %.% " (um)", ylab="Number of cells/100x100 um2", 
                           main=ifelse(s=="DD",a,""),
                           col=c("red","blue"),
                           legend=c("ulcerative","intact")
                          )
      
      # to divide between epi and dermis
      if (dist == "dist2DEJ") abline(v=0, lty=2)
      
  }
  }

  mtext(sample.labels[1], side = 2, line = 1.5, outer = T, at = .75, adj = NA, padj = NA, cex = NA, col = NA, font = NA)
  mtext(sample.labels[2], side = 2, line = 1.5, outer = T, at = .25, adj = NA, padj = NA, cex = NA, col = NA, font = NA)
  
  dev.off()
  
}

