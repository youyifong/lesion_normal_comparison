library(kyotil)

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
      dat$dist2skin=dat$dist2skin * ifelse(dat$epi, -1, 1)
      dat$dist2DEJ =dat$dist2DEJ  * ifelse(dat$epi, -1, 1)
    }
    dat
}

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
    
    hist(dist.ls[[1]], border=1, xlab=dist %.% " (um)", ylab="Number of cells", breaks=20, col=NULL, ylim=ylim, ...)
    for (i in 2:length(dist.ls)) hist(dist.ls[[i]], border=i, breaks=20, col=NULL, add=T)
    mylegend(x=3, legend=legend, lty=1, col=1:length(dist.ls), y.intersp=1, text.width=700)
    
    y.all = sapply(dist.ls, function (dist) density(dist)$y)
    ylim = range(y.all)    
    
    plot(density(dist.ls[[1]]), col=1, xlab=dist %.% " (um)", ylab="Density", ylim=ylim, ...)
    for (i in 2:length(dist.ls)) lines(density(dist.ls[[i]]), col=i)
    mylegend(x=3, legend=legend, lty=1, col=1:i, y.intersp=1, text.width=700)
    
}


region.labels=c(normal="intact", lesion="ulcerative")
region.labels.title=c(normal="Intact", lesion="Ulcerative")

red=rgb(red=1, blue=0, green=0, alpha=1)
blue=rgb(red=0, blue=1, green=0, alpha=1)

sample.labels=c("DD","CF"); names(sample.labels)=sample.labels
dist.labels=c("dist2skin", "dist2DEJ"); names(dist.labels)=dist.labels

dist.bt.peaks = sapply (sample.labels, function(s) {
 sapply (dist.labels, function(dist) {
  # s="DD"; dist="dist2skin"
  myprint(s, dist)
  
  if(dist=="dist2DEJ") {
    xlim=c(-750,1750)
  } else if (dist=="dist2skin") {
    xlim=c(0,2500)
  }
  
  if (s=="DD") {
    legend=c("CD3","CD15","CD123","CD11c","CD68")
  } else if (s=="CF") {
    legend=c("CD3","CD15","CD123","CD11c") # CD68 did not pass QC
    # lesion.width=/2 # computed from analysis.py
  }
  names(legend)=legend
  
  alldat=lapply(legend, function (l) get.dat(s%.%"_Les_"%.%l%.%"_masks_d2b_regmem.csv", s, dist) )
  
  ## compare distribution of markers
  
  dist.ls.normal = lapply(alldat, function (dat) subset(dat, normal)[[dist]])
  dist.ls.lesion = lapply(alldat, function (dat) subset(dat, lesion)[[dist]])
  
  mypng(mfcol=c(2,2), file="figures/"%.%s%.%"_"%.%dist%.%"_compare_markers")
    compare.markers(dist.ls.normal, xlim=xlim, main=s%.%" Les, "%.%region.labels.title["normal"], legend=legend)    
    compare.markers(dist.ls.lesion, xlim=xlim, main=s%.%" Les, "%.%region.labels.title["lesion"], legend=legend)    
  dev.off()
  
  if (dist=="dist2DEJ"){
    ## tabulate cell counts in epi/dermis x normal/lesion
    
    counts = sapply(legend, function (l) 
      c(nrow(subset(alldat[[l]], epi & normal)),
        nrow(subset(alldat[[l]], dermis & normal)),
        nrow(subset(alldat[[l]], epi & lesion)),
        nrow(subset(alldat[[l]], dermis & lesion)))
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
rownames(tab)[1:2]=c("interval between peaks (skin)", "interval between peaks (DEJ)")
mytex(tab, file="tables/peaks_distance_lesion_size", digits=0)



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


