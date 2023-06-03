################################################################################
## tabulate cell counts in epi/dermis x normal/lesion

count.tables = lapply (biopsies, function(s) {
  # s="DD"
  
  dist="dist2DEJ"
  
  markers=markers.ls[[s]]; names(markers)=markers
  
  counts = sapply(markers, function (a) 
    c(nrow(subset(alldat[[s]][[a]], epi &    region=="intact")),
      nrow(subset(alldat[[s]][[a]], dermis & region=="intact")),
      nrow(subset(alldat[[s]][[a]], epi &    region=="ulcerative")),
      nrow(subset(alldat[[s]][[a]], dermis & region=="ulcerative")))
  )
  
  dat=read.csv(glue("{data.dir}/{s}_Les_stitched_alignedtoCD45_cropped/{s}_Les_region_area_length.csv"))
  
  # remove rows that are lengths
  dat = dat[!contain(dat$name, "length"), ]
  
  dat$epi = contain(dat$name, "epi")
  dat$dermis = contain(dat$name, "dermis")
  dat$normal = contain(dat$name, "normal")
  dat$lesion = contain(dat$name, "lesion")
  
  sizes = with(dat, c(sum(size[epi&normal]), 
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
  print(tab)
  
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
  round(sapply(biopsies.1, function(b) tissue.widths[[b]]['ulcerative'])), 
  round(sapply(biopsies.1, function(b) tissue.widths[[b]]['intact'])), 
  
  round(sapply(biopsies.1, function(b) dist.bt.peaks[[b]]["ulcerative","dist2skin"])),
  round(sapply(biopsies.1, function(b) dist.bt.peaks[[b]]["ulcerative","dist2DEJ"])),
  
  cbind(c(t(count.tables$DD)), 
        c(t(count.tables$CF)),
        c(t(count.tables$JM)))
)

rownames(tab)=c(
  "lesion width",
  "normal width",
  "separation between CD3/CD15 peaks (ref: skin)",
  "separation between CD3/CD15 peaks (ref: DEJ)",
  c(t(outer(markers, c("epi","dermis","e/d"), paste)))
)
tab

# add ratio cols
tab=cbind(tab, 
          "CF/DD" = formatDouble(as.numeric(tab[,2])/as.numeric(tab[,1]),2,remove.leading0 = F),
          "JM/DD" = formatDouble(as.numeric(tab[,3])/as.numeric(tab[,1]),2,remove.leading0 = F))
tab

mytex(tab, file="tables/DD_CF_JM_comparison", 
      add.to.row=list(as.list(seq(4,by=3,length=4)),
                      rep("\\hline\n",4)
      )
)      




