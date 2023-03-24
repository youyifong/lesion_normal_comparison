library(kyotil)

region.labels=c(normal="intact epi", lesion="ulcerative")

red=rgb(red=1, blue=0, green=0, alpha=1)
blue=rgb(red=0, blue=1, green=0, alpha=1)


get.dat=function(filename) {
    dat=read.csv("data/"%.%filename)
    dat$in_DD_Les_lesion_side = dat$in_DD_Les_lesion_side=="True"
    dat$in_DD_Les_normal_side = dat$in_DD_Les_normal_side=="True"
    
    # each pixel is 1/2 micron
    dat$dist2boundary=dat$dist2boundary/2
    
    mean(dat$in_DD_Les_normal_side)
    mean(dat$in_DD_Les_lesion_side)
    
    with(dat, table(in_DD_Les_normal_side, in_DD_Les_lesion_side))
    
    dat
}

compare.regions=function(dat, legend, ...) {
    hist(dat$dist2boundary[dat$in_DD_Les_lesion_side], border="red", col=NULL, xlab="Distance to skin (um)", ylab="Number of cells", breaks=20, 
        ...)
    hist(dat$dist2boundary[dat$in_DD_Les_normal_side], border="blue", col=NULL, breaks=20, add=T)
    mylegend(x=3, legend=legend, lty=1, col=c("blue","red"), y.intersp=1.5, text.width=700)
    
    plot(density(dat$dist2boundary[dat$in_DD_Les_lesion_side]), col="red", xlab="Distance to skin (um)", ylab="Density",  
        ylim=range(c(density(dat$dist2boundary[dat$in_DD_Les_lesion_side])$y, density(dat$dist2boundary[dat$in_DD_Les_normal_side])$y)),
        ...)
    lines(density(dat$dist2boundary[dat$in_DD_Les_normal_side]), col="blue")
    mylegend(x=3, legend=legend, lty=1, col=c("blue","red"), y.intersp=1, text.width=700)
}


compare.markers=function(dat.m1, dat.m2, legend, ...) {
    hist(dat.m1$dist2boundary, border=blue, xlab="Distance to skin (um)", ylab="Number of cells", breaks=20, col=NULL,
        , ylim=range(c(hist(dat.m1$dist2boundary, plot=F)$counts, hist(dat.m2$dist2boundary, plot=F)$counts)),
        , ...)
    hist(dat.m2$dist2boundary, border=red, breaks=20, col=NULL, add=T)
    mylegend(x=3, legend=legend, lty=1, col=c("blue","red"), y.intersp=1.5, text.width=700)
    
    plot(density(dat.m1$dist2boundary), col=blue, xlab="Distance to skin (um)", ylab="Density" 
        , ylim=range(c(density(dat.m1$dist2boundary)$y, density(dat.m2$dist2boundary)$y))
        , ...)
    lines(density(dat.m2$dist2boundary), col=red)
    mylegend(x=3, legend=legend, lty=1, col=c("blue","red"), y.intersp=1, text.width=700)
}


dat.cd4=get.dat("DD_Les_CD4-PE_stitched_masks_d2b_regmem.csv")
dat.cd15=get.dat("DD_Les_CD15_stitched_masks_d2b_regmem.csv")


myfigure(mfcol=c(2,2))
    xlim=c(0,2500)
    compare.regions(dat.cd4, xlim=xlim, main="DD Les, CD4", legend=region.labels)    
    compare.regions(dat.cd15, xlim=xlim, main="DD Les, CD15", legend=region.labels)
mydev.off(file="figures/dist2skin_compare_region", ext="png")


myfigure(mfcol=c(2,2))
    xlim=c(0,2500)
    compare.markers(subset(dat.cd4, in_DD_Les_normal_side), subset(dat.cd15, in_DD_Les_normal_side), xlim=xlim, main="DD Les, "%.%region.labels["normal"], legend=c("CD4","CD15"))    
    compare.markers(subset(dat.cd4, in_DD_Les_lesion_side), subset(dat.cd15, in_DD_Les_lesion_side), xlim=xlim, main="DD Les, "%.%region.labels["lesion"], legend=c("CD4","CD15"))    
mydev.off(file="figures/dist2skin_compare_cd4_cd15", ext="png")
