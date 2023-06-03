################################################################################
# compare DD and CF in one plot
# 8 plots for each of dist2skin and dist2DEJ:
# 4 markers in lesion 
# 4 markers in intact

if (!file.exists("../figures")) dir.create("../figures")

for (dist in dist.labels) {
  xlim=xlim.ls[[dist]]
  myprint(dist)
  
  for (idx in 2) { # 1: density, 2: smoothed histogram
    mypng(mfcol=c(2,4), width=10, height=5, mar=c(4,4,2,1), oma=c(1,3,0,0),
          file=paste0("../figures/compare_DD_CF_JM_",dist,"_",ifelse(idx==1,"pdf","hist")))
    
    for (a in markers) {
      
      for (tissue.type in tissue.types) { 
        # dist="dist2skin"; a="CD15"; tissue.type="lesion"
        
        if (tissue.type=="normal") {
          dat.dd=subset(alldat[["DD"]][[dist]][[a]], normal)$dist
          dat.cf=subset(alldat[["CF"]][[dist]][[a]], normal)$dist
          dat.jm=subset(alldat[["JM"]][[dist]][[a]], normal)$dist
        } else {
          dat.dd=subset(alldat[["DD"]][[dist]][[a]], lesion)$dist
          dat.cf=subset(alldat[["CF"]][[dist]][[a]], lesion)$dist
          dat.jm=subset(alldat[["JM"]][[dist]][[a]], lesion)$dist
        }
        
        if (idx==1) {
          # plot density
          ylim = range(density(dat.dd)$y, 
                       density(dat.cf)$y,
                       density(dat.jm)$y)
          
          plot(density(dat.dd), col=1, xlab=dist %.% " (um)", ylab="Density", ylim=ylim, xlim=xlim,
               main=a)
          lines(density(dat.cf), col=2)
          lines(density(dat.jm), col=2)
          mylegend(x=3, legend=sample.labels.1, lty=1, col=1:3, y.intersp=1, text.width=700)
          
        } else {
          # plot smoothed histogram
          
          # compute breaks with 100 um width
          bin_width <- 100
          min_break <- floor  (min(dat.dd,dat.cf,dat.jm) / bin_width) * bin_width
          max_break <- ceiling(max(dat.dd,dat.cf,dat.jm) / bin_width) * bin_width
          breaks <- seq(min_break, max_break, by = bin_width)
          
          # get histogram info without plotting
          hist_info.dd <- hist(dat.dd, breaks = breaks, plot = FALSE) 
          hist_info.cf <- hist(dat.cf, breaks = breaks, plot = FALSE) 
          hist_info.jm <- hist(dat.jm, breaks = breaks, plot = FALSE) 
          
          # normalize the histogram counts to 100 um width
          hist_info.dd$counts <- hist_info.dd$counts / tissue.widths[["DD"]][tissue.type] * 100
          hist_info.cf$counts <- hist_info.cf$counts / tissue.widths[["CF"]][tissue.type] * 100
          hist_info.jm$counts <- hist_info.jm$counts / tissue.widths[["JM"]][tissue.type] * 100
          
          # estimate density curve
          densities.dd <- density(dat.dd)
          densities.cf <- density(dat.cf)
          densities.jm <- density(dat.jm)
          
          # Scale the height of the density curve to match the histogram's frequency
          densities.dd$y <- densities.dd$y * diff(hist_info.dd$mids[1:2]) * length(dat.dd)
          densities.cf$y <- densities.cf$y * diff(hist_info.cf$mids[1:2]) * length(dat.cf)
          densities.jm$y <- densities.jm$y * diff(hist_info.jm$mids[1:2]) * length(dat.jm)
          
          # normalize the density curve height to 100 um width
          densities.dd$y <- densities.dd$y / tissue.widths[["DD"]][tissue.type] * 100
          densities.cf$y <- densities.cf$y / tissue.widths[["CF"]][tissue.type] * 100
          densities.jm$y <- densities.jm$y / tissue.widths[["JM"]][tissue.type] * 100
          
          # ylim = range(hist_info.dd$counts, hist_info.cf$counts)
          ylim=c(0,150)
          
          plot(hist_info.dd, freq = TRUE, 
               border="white",  col=rgb(red=1, blue=0, green=0, alpha=.1), 
               xlab=dist %.% " (um)", ylab="Number of cells/100x100 um2", 
               xlim=xlim, ylim=ylim,
               main=ifelse(tissue.type=="lesion",a,"") )
          
          # plot density curves, which can be viewed as smoothed histograms
          lines(densities.dd, col = "red")
          lines(densities.cf, col = "blue")
          lines(densities.jm, col = "green")
          # to divide between epi and dermis
          if (dist == "dist2DEJ") abline(v=0, lty=2)
          
          mylegend(x=3, legend=sample.labels.1, lty=1, col=c("red","blue","green"), y.intersp=1.5, text.width=700)
          
        }
      }
    }
    mtext("Ulcerative", side = 2, line = 1.5, outer = T, at = .75, adj = NA, padj = NA, cex = NA, col = NA, font = NA)
    mtext("Intact", side = 2, line = 1.5, outer = T, at = .25, adj = NA, padj = NA, cex = NA, col = NA, font = NA)
    dev.off()
    
  }
}



