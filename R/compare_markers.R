################################################################################
# compare several markers in one plot
# 4 plots for each of dist2skin and dist2DEJ:
# ulcerative/intact x DD/CF

dir.exists("output/figures") || dir.create("output/figures", recursive = TRUE)

get.peak = function(x_start, freq) {
  x_start[which.max(freq)] + (x_start[2]-x_start[1])/2
}


mypng(mfrow=c(length(select.slides),4), file="output/figures/compare_markers", 
      width=10, height=4/2*length(select.slides), mar=c(4,0,2,0), oma=c(1,5,0,0))
        
  dist.bt.peaks = 
  lapply (select.slides, function(s) {
  sapply (tissue.types, function(tissue.type) { 
  sapply (dist.labels, simplify="array", function(dist) {
  # tissue.type="intact"; dist="dist2skin"; s="AS/AS1_Les"; 
    myprint(tissue.type, dist, s)
  
    # use all markers
    # markers=markers.ls[[s]]; names(markers)=markers
    # use a subset of mostly shared markers
    # markers = c("CD3", "CD15", "CD11c", "CD123", "CD4", "CD8", if(!s%in%c("LG/LG1_Les","SK/SK1_Les")) "CD68")
    markers = key.markers
    names(markers)=markers
    
    dat.ls = lapply(markers, function (a) subset(alldat[[s]][[a]], region==tissue.type)[[dist]])
    
    all.freq = scaled.hist(dat.ls,
         scale.factors=lapply(markers, function(a) subset(tissue_size_by_depth_and_region, biopsy==get.initial(s)%.%"/"%.%s & dist_type==dist & region==tissue.type)),
         xlim=xlim.ls[[dist]], ylim=ylim[[dist]],
         xlab=dist %.% " (um)", ylab="Number of cells/100x100 um2",
         main=ifelse(s==select.slides[1],tissue.type,""),
         col=seq_along(markers),
         legend=markers,
         cex.legend=0.7,
         yaxt=ifelse(tissue.types[1]==tissue.type & dist.labels[1]==dist,"s","n") # only have y axis in column 1
    )
    
    # dividing line between epi and dermis
    if (dist == "dist2DEJ") abline(v=0, lty=2)
    
    # return distance between CD3 and CD15 peaks
    # get.peak(all.freq[,"depth_start"], all.freq[,"CD3"]) - get.peak(all.freq[,"depth_start"], all.freq[,"CD15"])
  
  }) 
  })
  })
  
  
  # leftmost labeling
  tmp=get.slide.name(select.slides)
  for (i in 1:length(tmp)) {
    # ylab
    mtext("Cells/100x100 um2", side = 2, line = 2., outer = T, at = ((length(tmp)-i)*2 + 1)/(length(tmp)*2)+0., 
          adj = NA, padj = NA, cex = .7, col = NA, font = NA)
    mtext(tmp[i], side = 2, line = 3.5, outer = T, at = ((length(tmp)-i)*2 + 1)/(length(tmp)*2), 
          adj = NA, padj = NA, cex = NA, col = NA, font = NA)
  }

dev.off()



for (biopsy in biopsies) {
  slides = slides.by.biopsy[[biopsy]]
  mypng(mfrow=c(length(slides),4), file="output/figures/compare_markers_"%.%biopsy, width=8, height=4/2*length(slides), mar=c(4,0,2,0), oma=c(1,5,0,0))
  
    lapply (slides, function(s) {
      sapply (tissue.types, function(tissue.type) { 
        sapply (dist.labels, simplify="array", function(dist) {
          # tissue.type="intact"; dist="dist2skin"; s="AS_Les"; 
          myprint(tissue.type, dist, s)
          
          # use all markers
          # markers=markers.ls[[s]]; names(markers)=markers
          # use a subset of mostly shared markers
          # markers = c("CD3", "CD15", "CD11c", "CD123", "CD4", "CD8", if(!s%in%c("LG/LG1_Les","SK/SK1_Les")) "CD68")
          markers = key.markers
          names(markers)=markers
          
          dat.ls = lapply(markers, function (a) {
            if (is.null(alldat[[s]][[a]])) NULL else subset(alldat[[s]][[a]], region==tissue.type)[[dist]]
          })
          dat.ls=dat.ls[sapply(dat.ls, function(x) !is.null(x))]
          
          all.freq = scaled.hist(dat.ls,
                                 scale.factors=lapply(markers, function(a) subset(tissue_size_by_depth_and_region, biopsy==get.initial(s)%.%"/"%.%s & dist_type==dist & region==tissue.type)),
                                 xlim=xlim.ls[[dist]], ylim=ylim[[dist]],
                                 xlab=dist %.% " (um)", ylab="Number of cells/100x100 um2",
                                 main=ifelse(s==slides[1],tissue.type,""),
                                 col=seq_along(markers),
                                 legend=markers,
                                 cex.legend=0.7,
                                 yaxt=ifelse(tissue.types[1]==tissue.type & dist.labels[1]==dist,"s","n") # only have y axis in column 1
          )
          
          # dividing line between epi and dermis
          if (dist == "dist2DEJ") abline(v=0, lty=2)
          
          # return distance between CD3 and CD15 peaks
          # get.peak(all.freq[,"depth_start"], all.freq[,"CD3"]) - get.peak(all.freq[,"depth_start"], all.freq[,"CD15"])
          
        }) 
      })
    })
  
  
  # leftmost labeling
  tmp=get.slide.name(slides)
  for (i in 1:length(tmp)) {
    # ylab
    mtext("Cells/100x100 um2", side = 2, line = 2., outer = T, at = ((length(tmp)-i)*2 + 1)/(length(tmp)*2)+0., 
          adj = NA, padj = NA, cex = .7, col = NA, font = NA)
    mtext(tmp[i], side = 2, line = 3.5, outer = T, at = ((length(tmp)-i)*2 + 1)/(length(tmp)*2), 
          adj = NA, padj = NA, cex = NA, col = NA, font = NA)
  }
  
  dev.off()

}


################################################################################
# table combining select.slides
# 
# for (tissue.type in tissue.types) {
#   tab=rbind(
#     round(sapply(select.slides.1, function(b) dist.bt.peaks[[b]]["dist2skin",tissue.type])),
#     round(sapply(select.slides.1, function(b) dist.bt.peaks[[b]]["dist2DEJ",tissue.type])),
#     cbind(c(t(count.tables$DD[[tissue.type]])), 
#           c(t(count.tables$CF[[tissue.type]])),
#           c(t(count.tables$JM[[tissue.type]])))
#   )
#   
#   rownames(tab)=c(
#     "CD3 to CD15 (dist2skin)",
#     "CD3 to CD15 (dist2DEJ)",
#     c(t(outer(key.markers, c("epi","dermis","e/d"), paste)))
#   )
#   tab
#   
#   # add ratio cols
#   tab=cbind(tab, 
#             "CF/DD" = formatDouble(as.numeric(tab[,2])/as.numeric(tab[,1]),2,remove.leading0 = F),
#             "JM/DD" = formatDouble(as.numeric(tab[,3])/as.numeric(tab[,1]),2,remove.leading0 = F))
#   tab
#   
#   mytex(tab, file=glue("output/tables/{tissue.type}_comparison"), 
#         add.to.row=list(as.list(seq(3,by=3,length=4)),
#                         rep("\\hline\n",4)
#         )
#   )  
# }




