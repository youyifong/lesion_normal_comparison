library(kyotil)

region.labels=c(normal="intact", lesion="ulcerative")
region.labels.title=c(normal="Intact", lesion="Ulcerative")

red=rgb(red=1, blue=0, green=0, alpha=1)
blue=rgb(red=0, blue=1, green=0, alpha=1)

sample.labels=c("DD","CF"); names(sample.labels)=sample.labels
dist.labels=c("dist2skin", "dist2DEJ"); names(dist.labels)=dist.labels

markers=c("CD3","CD15","CD123","CD11c"); names(markers)=markers
# CD68 did not pass QC for CF
markers.ls=list(DD=c(markers,CD68="CD68"), CF=markers)

xlim.ls=list(dist2DEJ=c(-750,1750), dist2skin=c(0,2500))


get.dat=function(filename, sample.label, dist) {
    dat=read.csv("data/"%.%filename)
    
    if (sample.label=="DD") {
      dat$normal.epi =    dat$in_DD_Les_epi_region_1_normal=="True" 
      dat$normal.dermis = dat$in_DD_Les_dermis_region_1_normal=="True" 
      dat$lesion.epi =    dat$in_DD_Les_epi_region_2_lesion_alternate=="True" 
      dat$lesion.dermis = dat$in_DD_Les_dermis_region_2_lesion_alternate=="True" 

      # each pixel is 1/2 micron, negative dist in epi
      dat$dist2skin = dat$dist2DD_Les_skin_boundary/2
      dat$dist2DEJ  = dat$dist2DD_Les_DEJviaDAPI.WGA_alternate/2

    } else if (sample.label=="CF") {
      dat$normal.epi =    dat$in_CF_Les_epi_region_1_normal=="True" | dat$in_CF_Les_epi_region_3_normal=="True" 
      dat$normal.dermis = dat$in_CF_Les_dermis_region_1_normal=="True" | dat$in_CF_Les_dermis_region_3_normal=="True" 
      dat$lesion.epi =    dat$in_CF_Les_epi_region_2_lesion=="True" 
      dat$lesion.dermis = dat$in_CF_Les_dermis_region_2_lesion=="True" 
      
      # each pixel is 1/2 micron, negative dist in epi
      dat$dist2skin = dat$dist2CF_Les_skin_boundary/2
      dat$dist2DEJ  = dat$dist2CF_Les_DEJ_v2/2
    
    }
    
    dat$epi = dat$normal.epi | dat$lesion.epi
    dat$dermis = dat$normal.dermis | dat$lesion.dermis
    dat$normal = dat$normal.epi | dat$normal.dermis
    dat$lesion = dat$lesion.epi | dat$lesion.dermis
    
    if (dist=="dist2DEJ") {
      dat$dist = dat$dist2DEJ  * ifelse(dat$epi, -1, 1)
    } else {
      dat$dist = dat$dist2skin 
    }
    dat
}


alldat = lapply (sample.labels, function(s) 
  lapply (dist.labels, function(dist) 
    lapply (markers.ls[[s]], function (l) 
      get.dat(s%.%"_Les_"%.%l%.%"_masks_d2b_regmem.csv", s, dist)
)))




compare.regions=function(dat, legend, ...) {
    hist(dat[[dist]][dat$lesion], border="red", col=NULL, xlab=dist %.% " (um)", ylab="Number of cells", breaks=20
        , ylim=range(c(hist(dat[[dist]][dat$lesion], breaks=20, plot=F)$counts, hist(dat[[dist]][dat$normal], breaks=20, plot=F)$counts))
        , ...)
    hist(dat[[dist]][dat$normal], border="blue", col=NULL, breaks=20, add=T)
    mylegend(x=3, legend=legend, lty=1, col=c("blue","red"), y.intersp=1.5, text.width=700)
    
    plot(density(dat[[dist]][dat$lesion]), col="red", xlab=dist %.% " (um)", ylab="Density"
        , ylim=range(c(density(dat[[dist]][dat$lesion])$y, density(dat[[dist]][dat$normal])$y))
        , ...)
    lines(density(dat[[dist]][dat$in_Les_normal_side]), col="blue")
    mylegend(x=3, legend=legend, lty=1, col=c("blue","red"), y.intersp=1, text.width=700)
}


compare.markers=function(dist.ls, legend, ...) {
    cnts.all = sapply(dist.ls, function (dist) hist(dist, breaks=20, plot=F)$counts)
    ylim = range(cnts.all)    
    
    hist(dist.ls[[1]], border=1, ylab="Number of cells", breaks=20, col=NULL, ylim=ylim, ...)
    for (i in 2:length(dist.ls)) hist(dist.ls[[i]], border=i, breaks=20, col=NULL, add=T)
    mylegend(x=3, legend=legend, lty=1, col=1:length(dist.ls), y.intersp=1, text.width=700)
    
    y.all = sapply(dist.ls, function (dist) density(dist)$y)
    ylim = range(y.all)    
    
    plot(density(dist.ls[[1]]), col=1, ylab="Density", ylim=ylim, ...)
    for (i in 2:length(dist.ls)) lines(density(dist.ls[[i]]), col=i)
    mylegend(x=3, legend=legend, lty=1, col=1:i, y.intersp=1, text.width=700)
    
}



#################################################################################
# compare distribution of distances to a line (skin or DEJ) between lesion and intact

dist.bt.peaks = sapply (sample.labels, function(s) {
 sapply (dist.labels, function(dist) {
  # s="DD"; dist="dist2skin"
  myprint(s, dist)
  
  xlim=xlim.ls[[dist]]
  legend=markers.ls[[s]]; names(legend)=legend
  
  dat.ls=lapply(legend, function (l) alldat[[s]][[dist]][[l]] )
  dist.ls.normal = lapply(dat.ls, function (dat) dat[dat$normal, "dist"])
  dist.ls.lesion = lapply(dat.ls, function (dat) dat[dat$lesion, "dist"])

  mypng(mfcol=c(2,2), file="figures/"%.%s%.%"_"%.%dist%.%"_compare_markers")
    compare.markers(dist.ls.normal, legend=legend, xlim=xlim, main=s%.%" Les, "%.%region.labels.title["normal"], xlab=dist %.% " (um)")    
    compare.markers(dist.ls.lesion, legend=legend, xlim=xlim, main=s%.%" Les, "%.%region.labels.title["lesion"], xlab=dist %.% " (um)")    
  dev.off()
  
  
  ## tabulate cell counts in epi/dermis x normal/lesion
  if (dist=="dist2DEJ"){
    
    counts = sapply(legend, function (l) 
      c(nrow(subset(dat.ls[[l]], epi & normal)),
        nrow(subset(dat.ls[[l]], dermis & normal)),
        nrow(subset(dat.ls[[l]], epi & lesion)),
        nrow(subset(dat.ls[[l]], dermis & lesion)))
    )

    # counts/area
    tmp = read.csv("data/"%.%s%.%"_region_size.csv")
    tmp$epi = contain(tmp$name, "epi")
    tmp$dermis = contain(tmp$name, "dermis")
    tmp$normal = contain(tmp$name, "normal")
    tmp$lesion = contain(tmp$name, "lesion")
    sizes = with(tmp, c(sum(size[epi&normal]), 
                        sum(size[dermis & normal]), 
                        sum(size[epi & lesion]), 
                        sum(size[dermis & lesion])
                       ))
    
    counts2=counts/sizes
    counts2=counts2*4*1e6 #convert to per mm2
    mode(counts2) <- "integer"
    tab=t(rbind(counts, counts2)) 
    mode(tab) <- "character" # convert to string
    tab=data.frame(tab)
    # add area
    tab=rbind(tab, area=formatDouble(c(sizes/4/1e6, rep(NA,4)), digits = 2, remove.leading0=F))
      
    mytex(data.frame(tab), file.name="tables/"%.%s%.%"_counts", digit=0,  include.colnames = F, align=c("l","r","r","r","r|","r","r","r","r"),
          col.headers="\\hline\n 
             \\multicolumn{1}{c}{} & \\multicolumn{4}{c|}{Counts} & \\multicolumn{4}{c}{Counts/mm2} \\\\ 
             \\multicolumn{1}{c}{} & \\multicolumn{2}{c}{normal} & \\multicolumn{2}{c|}{lesion} & \\multicolumn{2}{c}{normal} & \\multicolumn{2}{c}{lesion}  \\\\ 
             \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{epi} & dermis & \\multicolumn{1}{c}{epi} & \\multicolumn{1}{c|}{dermis} & \\multicolumn{1}{c}{epi} & dermis & \\multicolumn{1}{c}{epi} & dermis  \\\\ 
             \\hline
        ")
    
  }# end if dist
  

  
  # return distance between CD3 and CD15 peaks
  get.peak = function(density_est) {
    max_density_index <- which.max(density_est$y)
    peak_location <- density_est$x[max_density_index]
    return (peak_location)
  }

  get.peak(density(dist.ls.lesion$CD3)) - get.peak(density(dist.ls.lesion$CD15))

  
 }) # end dist loop
})# end s loop
  

lesion.width=c(DD=5723/2, CF=2404/2) # computed from analysis.py
tab=rbind(dist.bt.peaks, "lesion width"=lesion.width)
rownames(tab)[1:2]=c("space between peaks (skin)", "space between peaks (DEJ)")
mytex(tab, file="tables/peaks_distance_lesion_size", digits=0)



################################################################################
# compare DD and CF in 8 figures
# one for each of the 4 markers in lesion and one for each of the 4 markers in intact

for (dist in dist.labels) {
  xlim=xlim.ls[[dist]]
  
  mypng(mfcol=c(2,4), file=paste0("figures/compare_DD_CF_",dist))
  for (a in markers) {
    for (i in 1:2) {
      if (i==1) {
        dat.dd=subset(alldat[["DD"]][[dist]][[a]], normal)
        dat.cf=subset(alldat[["CF"]][[dist]][[a]], normal)
      } else {
        dat.dd=subset(alldat[["DD"]][[dist]][[a]], lesion)
        dat.cf=subset(alldat[["CF"]][[dist]][[a]], lesion)
      }
      
      ylim = range(density(dat.dd$dist)$y, 
                   density(dat.cf$dist)$y)
      
      plot(density(dat.dd$dist), col=1, xlab=dist %.% " (um)", ylab="Density", ylim=ylim, xlim=xlim,
           main=paste0("Les, ", region.labels.title[i], ", ", a))
      lines(density(dat.cf$dist), col=2)
      mylegend(x=3, legend=sample.labels, lty=1, col=1:2, y.intersp=1, text.width=700)
      
      # hist(dat.dd$dist, border="red",  col=NULL, breaks=20, 
      #      xlab=dist %.% " (um)", ylab="Number of cells", ylim=ylim
      #      # , xlim=xlim, main=s%.%" Les, "%.%region.labels.title["normal"], legend=legend
      #      )
      # hist(dat.cf$dist, border="blue", col=NULL, breaks=20, add=T)
      # mylegend(x=3, legend=sample.labels, lty=1, col=c("blue","red"), y.intersp=1.5, text.width=700)
      
    }
  }
  dev.off()
  
}



# ### compare normal vs lesion tissues ####
# # generate one png per marker
# 
# for (s in c("DD","CF")) {
#   if (s=="DD") markers=c(        "CD3", "CD4", "CD8a",          "CD15", "CD45", "CD68", "CD69", "CD123", "CD279")
#   if (s=="CF") markers=c("CD1a", "CD3", "CD4", "CD8a", "CD11c", "CD15", "CD45",                 "CD123", "CD279")
#   for(m in markers) {
#     #s="CF"; m="CD1a"
#     myfigure(mfcol=c(1,2))
#     compare.regions(get.dat(s%.%"_Les_"%.%m%.%"_masks_d2b_regmem.csv", s), xlim=xlim, main=s%.%" Les, "%.%m, legend=region.labels)    
#     mydev.off(file="figures/"%.%s%.%"_dist2skin_compare_regions_"%.%m, ext="png")
#   }
# }


