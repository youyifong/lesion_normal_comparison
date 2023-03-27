library(kyotil)

get.dat=function(filename, sample.label) {
    dat=read.csv("data/"%.%filename)
    dat$in_Les_lesion_side = dat[["in_"%.%sample.label%.%"_Les_lesion_side"]]=="True"
    dat$in_Les_normal_side = dat[["in_"%.%sample.label%.%"_Les_normal_side"]]=="True"
    
    # each pixel is 1/2 micron
    dat$dist2boundary=dat$dist2boundary/2
    
    mean(dat$in_Les_normal_side)
    mean(dat$in_Les_lesion_side)
    
    with(dat, table(in_Les_normal_side, in_Les_lesion_side))
    
    dat
}

compare.regions=function(dat, legend, ...) {
    hist(dat$dist2boundary[dat$in_Les_lesion_side], border="red", col=NULL, xlab="Distance to skin (um)", ylab="Number of cells", breaks=20
        , ylim=range(c(hist(dat$dist2boundary[dat$in_Les_lesion_side], breaks=20, plot=F)$counts, hist(dat$dist2boundary[dat$in_Les_normal_side], breaks=20, plot=F)$counts))
        , ...)
    hist(dat$dist2boundary[dat$in_Les_normal_side], border="blue", col=NULL, breaks=20, add=T)
    mylegend(x=3, legend=legend, lty=1, col=c("blue","red"), y.intersp=1.5, text.width=700)
    
    plot(density(dat$dist2boundary[dat$in_Les_lesion_side]), col="red", xlab="Distance to skin (um)", ylab="Density"
        , ylim=range(c(density(dat$dist2boundary[dat$in_Les_lesion_side])$y, density(dat$dist2boundary[dat$in_Les_normal_side])$y))
        , ...)
    lines(density(dat$dist2boundary[dat$in_Les_normal_side]), col="blue")
    mylegend(x=3, legend=legend, lty=1, col=c("blue","red"), y.intersp=1, text.width=700)
}


compare.markers=function(dist.ls, legend, ...) {
    cnts.all = sapply(dist.ls, function (dist) hist(dist, breaks=20, plot=F)$counts)
    ylim = range(cnts.all)    
    
    hist(dist.ls[[1]], border=1, xlab="Distance to skin (um)", ylab="Number of cells", breaks=20, col=NULL, ylim=ylim, ...)
    for (i in 2:length(dist.ls)) hist(dist.ls[[i]], border=i, breaks=20, col=NULL, add=T)
    mylegend(x=3, legend=legend, lty=1, col=1:length(dist.ls), y.intersp=1, text.width=700)
    
    y.all = sapply(dist.ls, function (dist) density(dist)$y)
    ylim = range(y.all)    
    
    plot(density(dist.ls[[1]]), col=1, xlab="Distance to skin (um)", ylab="Density", ylim=ylim, ...)
    for (i in 2:length(dist.ls)) lines(density(dist.ls[[i]]), col=i)
    mylegend(x=3, legend=legend, lty=1, col=1:i, y.intersp=1, text.width=700)
}


region.labels=c(normal="intact epi", lesion="ulcerative")
region.labels.title=c(normal="Intact Epidermis", lesion="Ulcerative")

red=rgb(red=1, blue=0, green=0, alpha=1)
blue=rgb(red=0, blue=1, green=0, alpha=1)

xlim=c(0,2500)


### compare normal vs lesion tissues ####

# generate one png per marker

for (s in c("DD","CF")) {
    if (s=="DD") markers=c(        "CD3", "CD4", "CD8a",          "CD15", "CD45", "CD68", "CD69", "CD123", "CD279")
    if (s=="CF") markers=c("CD1a", "CD3", "CD4", "CD8a", "CD11c", "CD15", "CD45",                 "CD123", "CD279")
    for(m in markers) {
        #s="CF"; m="CD1a"
        myfigure(mfcol=c(1,2))
            compare.regions(get.dat(s%.%"_Les_"%.%m%.%"_masks_d2b_regmem.csv", s), xlim=xlim, main=s%.%" Les, "%.%m, legend=region.labels)    
        mydev.off(file="figures/"%.%s%.%"_dist2skin_compare_regions_"%.%m, ext="png")
    }
}




#### compare markers ####

## DD
legend=c("CD3","CD15","CD123","CD11c","CD68")

dat.cd3=get.dat("DD_Les_CD3_masks_d2b_regmem.csv", "DD")
dat.cd15=get.dat("DD_Les_CD15_masks_d2b_regmem.csv", "DD")
dat.cd123=get.dat("DD_Les_CD123_masks_d2b_regmem.csv", "DD")
dat.cd11c=get.dat("DD_Les_CD11c_masks_d2b_regmem.csv", "DD")
dat.cd68=get.dat("DD_Les_CD68_masks_d2b_regmem.csv", "DD")

dist.ls.normal = list(subset(dat.cd3, in_Les_normal_side)$dist2boundary, 
                      subset(dat.cd15, in_Les_normal_side)$dist2boundary, 
                      subset(dat.cd123, in_Les_normal_side)$dist2boundary, 
                      subset(dat.cd11c, in_Les_normal_side)$dist2boundary,
                      subset(dat.cd68, in_Les_normal_side)$dist2boundary)
dist.ls.lesion = list(subset(dat.cd3, in_Les_lesion_side)$dist2boundary, 
                      subset(dat.cd15, in_Les_lesion_side)$dist2boundary, 
                      subset(dat.cd123, in_Les_lesion_side)$dist2boundary, 
                      subset(dat.cd11c, in_Les_lesion_side)$dist2boundary,
                      subset(dat.cd68, in_Les_lesion_side)$dist2boundary)

myfigure(mfcol=c(2,2))
    compare.markers(dist.ls.normal, xlim=xlim, main="DD Les, "%.%region.labels.title["normal"], legend=legend)    
    compare.markers(dist.ls.lesion, xlim=xlim, main="DD Les, "%.%region.labels.title["lesion"], legend=legend)    
mydev.off(file="figures/DD_dist2skin_compare_markers", ext="png")


## CF
legend=c("CD3","CD15","CD123","CD11c") # CD68 did not pass QC

dat.cd3=get.dat("CF_Les_CD3_masks_d2b_regmem.csv", "CF")
dat.cd15=get.dat("CF_Les_CD15_masks_d2b_regmem.csv", "CF")
dat.cd123=get.dat("CF_Les_CD123_masks_d2b_regmem.csv", "CF")
dat.cd11c=get.dat("CF_Les_CD11c_masks_d2b_regmem.csv", "CF")

dist.ls.normal = list(subset(dat.cd3, in_Les_normal_side)$dist2boundary, 
                      subset(dat.cd15, in_Les_normal_side)$dist2boundary, 
                      subset(dat.cd123, in_Les_normal_side)$dist2boundary, 
                      subset(dat.cd11c, in_Les_normal_side)$dist2boundary)
dist.ls.lesion = list(subset(dat.cd3, in_Les_lesion_side)$dist2boundary, 
                      subset(dat.cd15, in_Les_lesion_side)$dist2boundary, 
                      subset(dat.cd123, in_Les_lesion_side)$dist2boundary, 
                      subset(dat.cd11c, in_Les_lesion_side)$dist2boundary)
    
myfigure(mfcol=c(2,2))
    compare.markers(dist.ls.normal, xlim=xlim, main="CF Les, "%.%region.labels.title["normal"], legend=legend)    
    compare.markers(dist.ls.lesion, xlim=xlim, main="CF Les, "%.%region.labels.title["lesion"], legend=legend)    
mydev.off(file="figures/CF_dist2skin_compare_markers", ext="png")
