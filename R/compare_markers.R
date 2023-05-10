################################################################################
# compare several markers in one plot
# 4 plots for each of dist2skin and dist2DEJ:
# lesion/normal x DD/CF

get.peak = function(density_est) {
  max_density_index <- which.max(density_est$y)
  peak_location <- density_est$x[max_density_index]
  return (peak_location)
}


mypng(mfrow=c(2,4), file="figures/compare_markers", 
      width=10, height=5, mar=c(4,4,2,1), oma=c(1,3,0,0))
      
dist.bt.peaks = lapply (tissue.types, function(tissue.type) { 
  sapply (dist.labels, function(dist) {
  sapply (sample.labels, function(s) {
    # tissue.type="lesion"; dist="dist2skin"; s="DD"; 
    myprint(tissue.type, dist, s)
    
    xlim=xlim.ls[[dist]]
    markers=markers.ls[[s]]; names(markers)=markers
    
    dat.ls=lapply(markers, function (l) alldat[[s]][[dist]][[l]] )
    dat.ls = lapply(dat.ls, function (dat) dat[dat[[tissue.type]], "dist"])

    # compare.markers(dat.ls, 
    #                 xlim=xlim, 
    #                 main=s%.%" Les, "%.%region.labels.title["normal"], 
    #                 xlab=dist %.% " (um)",
    #                 markers=markers
    # )    
    
    # plot smoothed histogram
    
    # compute breaks with 100 um width
    bin_width <- 100
    min_break <- floor(min(unlist(dat.ls)) / bin_width) * bin_width
    max_break <- ceiling(max(unlist(dat.ls)) / bin_width) * bin_width
    breaks <- seq(min_break, max_break, by = bin_width)
    
    # get histogram info without plotting
    hist_info <- lapply(markers, function (a) hist(dat.ls[[a]], breaks = breaks, plot = FALSE) )

    # normalize the histogram counts to 100 um width
    for (a in markers) {
      hist_info[[a]]$counts <- hist_info[[a]]$counts / tissue.widths[[tissue.type]][s] * 100
    }

    # estimate density curve
    densities <- lapply(markers, function (a) density(dat.ls[[a]]) )

    for (a in markers) {
      # Scale the height of the density curve to match the histogram's frequency
      densities[[a]]$y <- densities[[a]]$y * diff(hist_info[[a]]$mids[1:2]) * length(dat.ls[[a]])
      # normalize the density curve height to 100 um width
      densities[[a]]$y <- densities[[a]]$y / tissue.widths[[tissue.type]][s] * 100
    }
    
    ylim=c(0,150)
    
    plot(hist_info[[1]], freq = TRUE, 
         border="white",  col="gray", 
         xlab=dist %.% " (um)", ylab="Number of cells/100x100 um2", 
         xlim=xlim, ylim=ylim,
         main=ifelse(tissue.type=="lesion",s,"") )
    
    # plot density curves, which can be viewed as smoothed histograms
    for (a in markers) {
      lines(densities[[a]], col = match(a, markers))
    }
    # to divide between epi and dermis
    if (dist == "dist2DEJ") abline(v=0, lty=2)
    
    mylegend(x=3, legend=markers, lty=1, col=1:length(markers), y.intersp=1, text.width=700)
    
    
    # return distance between CD3 and CD15 peaks
    get.peak(densities$CD3) - get.peak(densities$CD15)
    
}) 
})
})

mtext("Ulcerative", side = 2, line = 1.5, outer = T, at = .75, adj = NA, padj = NA, cex = NA, col = NA, font = NA)
mtext("Intact", side = 2, line = 1.5, outer = T, at = .25, adj = NA, padj = NA, cex = NA, col = NA, font = NA)
dev.off()


