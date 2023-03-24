library(kyotil)

compare.dist=function(filename, ...) {
    dat=read.csv("data/"%.%filename)
    dat$in_DD_Les_lesion_side = dat$in_DD_Les_lesion_side=="True"
    dat$in_DD_Les_normal_side = dat$in_DD_Les_normal_side=="True"
    
    # each pixel is 1/2 micron
    dat$dist2boundary=dat$dist2boundary/2
    
    mean(dat$in_DD_Les_normal_side)
    mean(dat$in_DD_Les_lesion_side)
    
    with(dat, table(in_DD_Les_normal_side, in_DD_Les_lesion_side))
    
    hist(dat$dist2boundary[dat$in_DD_Les_lesion_side], col="red", xlab="Distance to skin (um)", ylab="Number of cells", main=concatList(strsplit(filename, "_")[[1]][1:3], "_"), breaks=20, ...)
    hist(dat$dist2boundary[dat$in_DD_Les_normal_side], col="blue", breaks=20, add=T)
    mylegend(x=3, legend=c("normal","lesion"), lty=1, col=c("blue","red"), y.intersp=1.5, text.width=700)
}

myfigure(mfcol=c(1,2))
    xlim=c(0,2500)
    compare.dist("DD_Les_CD4-PE_stitched_masks_d2b_regmem.csv", xlim=xlim)    
    compare.dist("DD_Les_CD15_stitched_masks_d2b_regmem.csv", xlim=xlim)
mydev.off(file="figures/dist2skin_comparison", ext="png")
