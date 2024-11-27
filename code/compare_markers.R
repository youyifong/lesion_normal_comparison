################################################################################
# compare several markers in one plot
# 4 plots for each of dist2skin and dist2DEJ:
# ulcerative/intact x DD/CF

dir.exists("report/figures") || dir.create("report/figures", recursive = TRUE)

get.peak = function(density_est) {
  max_density_index <- which.max(density_est$y)
  peak_location <- density_est$x[max_density_index]
  return (peak_location)
}


mypng(mfrow=c(length(biopsies),4), file="report/figures/compare_markers", 
      width=10, height=4/2*length(biopsies), mar=c(4,4,2,1), oma=c(1,3,0,0))
      
dist.bt.peaks = 
lapply (biopsies, function(s) {
sapply (dist.labels, simplify="array", function(dist) {
sapply (tissue.types, function(tissue.type) { 
  # tissue.type="ulcerative"; dist="dist2skin"; s="AS"; 
  myprint(tissue.type, dist, s)

  if (!s %in% biopsies.1 & tissue.type=="ulcerative") {
    empty.plot()
    NA
    
  } else {
    markers=markers.ls[[s]]; names(markers)=markers
    
    dat.ls = lapply(markers, function (a) subset(alldat[[s]][[a]], region==tissue.type)[[dist]])
    
    all.densities = scaled.hist(dat.ls,
         scale.factors=lapply(markers, function(a) subset(tissue_size_by_depth_and_region, biopsy==s & dist_type==dist & region==tissue.type)),
         xlim=xlim.ls[[dist]], ylim=ylim[[dist]],
         xlab=dist %.% " (um)", ylab="Number of cells/100x100 um2",
         main=ifelse(s==biopsies[1],tissue.type,""),
         col=seq_along(markers),
         legend=markers,
         cex.legend=0.7
    )
    
    # dividing line between epi and dermis
    if (dist == "dist2DEJ") abline(v=0, lty=2)
    
    # return distance between CD3 and CD15 peaks
    get.peak(all.densities$CD3) - get.peak(all.densities$CD15)
  }
}) 
})
})


# leftmost labeling
tmp=biopsies
for (i in 1:length(tmp)) {
  mtext(tmp[i], side = 2, line = 1.5, outer = T, at = ((length(tmp)-i)*2 + 1)/(length(tmp)*2), 
        adj = NA, padj = NA, cex = NA, col = NA, font = NA)
}

dev.off()


