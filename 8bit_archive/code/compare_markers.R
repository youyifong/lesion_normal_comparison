################################################################################
# compare several markers in one plot
# 4 plots for each of dist2skin and dist2DEJ:
# ulcerative/intact x DD/CF

dir.exists("report/figures") || dir.create("report/figures", recursive = TRUE)

get.peak = function(x_start, freq) {
  x_start[which.max(freq)] + (x_start[2]-x_start[1])/2
}


mypng(mfrow=c(length(biopsies),4), file="report/figures/compare_markers", 
      width=10, height=4/2*length(biopsies), mar=c(4,0,2,0), oma=c(1,5,0,0))
      
dist.bt.peaks = 
lapply (biopsies, function(s) {
sapply (tissue.types, function(tissue.type) { 
sapply (dist.labels, simplify="array", function(dist) {
# tissue.type="ulcerative"; dist="dist2skin"; s="CF_Les"; 
  myprint(tissue.type, dist, s)

  if (!s %in% biopsies.1 & tissue.type=="ulcerative") {
    empty.plot()
    NA
    
  } else {
    # use all markers
    # markers=markers.ls[[s]]; names(markers)=markers
    # use a subset of mostly shared markers
    markers = c("CD3", "CD15", "CD11c", "CD123", "CD4", "CD8a", if(s!="CF_Les") "CD68"); names(markers)=markers
    
    dat.ls = lapply(markers, function (a) subset(alldat[[s]][[a]], region==tissue.type)[[dist]])
    
    all.freq = scaled.hist(dat.ls,
         scale.factors=lapply(markers, function(a) subset(tissue_size_by_depth_and_region, biopsy==s & dist_type==dist & region==tissue.type)),
         xlim=xlim.ls[[dist]], ylim=ylim[[dist]],
         xlab=dist %.% " (um)", ylab="Number of cells/100x100 um2",
         main=ifelse(s==biopsies[1],tissue.type,""),
         col=seq_along(markers),
         legend=markers,
         cex.legend=0.7,
         yaxt=ifelse(tissue.types[1]==tissue.type & dist.labels[1]==dist,"s","n") # only have y axis in column 1
    )
    
    # dividing line between epi and dermis
    if (dist == "dist2DEJ") abline(v=0, lty=2)
    
    # return distance between CD3 and CD15 peaks
    get.peak(all.freq[,"depth_start"], all.freq[,"CD3"]) - get.peak(all.freq[,"depth_start"], all.freq[,"CD15"])
  }
}) 
})
})


# leftmost labeling
tmp=biopsies
for (i in 1:length(tmp)) {
  mtext("Cells/100x100 um2", side = 2, line = 2., outer = T, at = ((length(tmp)-i)*2 + 1)/(length(tmp)*2)+0., 
        adj = NA, padj = NA, cex = .7, col = NA, font = NA)
  mtext(tmp[i], side = 2, line = 3.5, outer = T, at = ((length(tmp)-i)*2 + 1)/(length(tmp)*2), 
        adj = NA, padj = NA, cex = NA, col = NA, font = NA)
}

dev.off()


