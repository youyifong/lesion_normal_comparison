################################################################################
# compare different biopsies in one plot
# CF, DD and JM have both ulcerative and intact regions, others just have intact regions
# separately for dist2skin and dist2DEJ
# separately for ulcerative and intact

for (dist in dist.labels) { # dist2skin, dist2DEJ
  xlim=xlim.ls[[dist]]
  myprint(dist)

  myfigure(mfcol=c(2,4), width=10, height=5, mar=c(4,4,2,1), oma=c(1,3,0,0))
  
  for (a in markers) {
  for (tissue.type in tissue.types) { # intact, ulcerative
    
    # dist="dist2skin"; a="CD15"; tissue.type="ulcerative"
    
    if (tissue.type=="intact") {
      dat.ls = lapply(biopsies, function (t) subset(alldat[[t]][[dist]][[a]], intact)$dist)
      scale.factors=sapply(biopsies, function(t) unname(100/tissue.widths[[t]][[tissue.type]]))
    } else {
      dat.ls = lapply(biopsies.1, function (t) subset(alldat[[t]][[dist]][[a]], ulcerative)$dist)
      scale.factors=sapply(biopsies.1, function(t) unname(100/tissue.widths[[t]][[tissue.type]]))
    }
    
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
  mtext("Ulcerative", side = 2, line = 1.5, outer = T, at = .75, adj = NA, padj = NA, cex = NA, col = NA, font = NA)
  mtext("Intact", side = 2, line = 1.5, outer = T, at = .25, adj = NA, padj = NA, cex = NA, col = NA, font = NA)

  mydev.off(ext="png",file=paste0("figures/compare_DD_CF_JM_",dist,"_",ifelse(idx==1,"pdf","hist")))
  
}



