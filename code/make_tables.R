################################################################################
## tabulate cell counts in epi/dermis x normal/lesion

count.tables = lapply (sample.labels, function(s) {
  # s="DD"
  dist="dist2DEJ"
  
  xlim=xlim.ls[[dist]]
  markers=markers.ls[[s]]; names(markers)=markers
  
  dat.ls=lapply(markers, function (l) alldat[[s]][[dist]][[l]] )
  dist.ls.normal = lapply(dat.ls, function (dat) dat[dat$normal, "dist"])
  dist.ls.lesion = lapply(dat.ls, function (dat) dat[dat$lesion, "dist"])
  
  
  counts = sapply(markers, function (l) 
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
  
  ret = counts2[3:4,1:4]
  mode(ret) = "character"
  cbind(t(ret), formatDouble(counts2[3,1:4]/counts2[4,1:4],2, remove.leading0 = F))
}) 


################################################################################
# table to compare DD and CF

tab=rbind(
  round(tissue.widths$lesion), 
  round(tissue.widths$normal), 
  round(t(dist.bt.peaks$lesion)),
  cbind(c(t(count.tables$DD)), c(t(count.tables$CF)))
)
rownames(tab)=c(
  "lesion width",
  "normal width",
  "separation between CD3/CD15 peaks (ref: skin)", 
  "separation between CD3/CD15 peaks (ref: DEJ)", 
  c(t(outer(markers, c("epi","dermis","e/d"), paste)))
)

# add ratio col
tab=cbind(tab, "DD/CF" = formatDouble(as.numeric(tab[,1])/as.numeric(tab[,2]),2,remove.leading0 = F))
tab

mytex(tab, file="tables/DD_CF_comparison", 
      add.to.row=list(as.list(seq(4,by=3,length=4)),
                      rep("\\hline\n",4)
      )
)      




