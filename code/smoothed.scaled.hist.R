smoothed.scaled.hist = function(dat.ls, bin_width, scale.factors, cols=NULL, legend=NULL, cex.legend=1, ...) {
  
  bin_width <- 100
  min_break <- floor(min(unlist(dat.ls)) / bin_width) * bin_width
  max_break <- ceiling(max(unlist(dat.ls)) / bin_width) * bin_width
  breaks <- seq(min_break, max_break, by = bin_width)
  
  labels=names(dat.ls)
  names(labels)=labels
  
  if (is.null(col)) cols=1:length(dat.ls)
  
  # get histogram info without plotting
  hist_info <- lapply(labels, function (a) hist(dat.ls[[a]], breaks = breaks, plot = FALSE) )
  for (a in labels) {
    # normalize the histogram counts to 100 um width
    hist_info[[a]]$counts <- hist_info[[a]]$counts * scale.factors[a]
  }
  
  # estimate density curve
  densities <- lapply(labels, function (a) density(dat.ls[[a]]) )
  for (a in labels) {
    # Scale the height of the density curve to match the histogram's frequency
    densities[[a]]$y <- densities[[a]]$y * diff(hist_info[[a]]$mids[1:2]) * length(dat.ls[[a]])
    # normalize the density curve height to 100 um width
    densities[[a]]$y <- densities[[a]]$y * scale.factors[a]
  }
  
  hist.col = col2rgb(cols[1])
  hist.col = rgb(hist.col[1], hist.col[2], hist.col[3], alpha=255*.15, maxColorValue=255)

  plot(hist_info[[1]], freq = TRUE, border="white",  col=hist.col, ...)
  for (i in 1:length(labels)) lines(densities[[i]], col = cols[i])
  
  mylegend(x=3, legend=if(is.null(legend)) labels else legend, lty=1, col=cols, y.intersp=1, text.width=700, cex=cex.legend)
  
  densities
}