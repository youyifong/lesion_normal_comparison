################################################################################
## tabulate cell counts in epi/dermis x normal/lesion

dir.exists("report/tables") || dir.create("report/tables", recursive = TRUE)

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
  
  dat=read.csv(glue("{data.dir}/{s}_stitched_alignedtoCD45_cropped/{s}_region_area_length.csv"))
  
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
  
  ret = counts2[,markers.1]
  mode(ret) = "character"
  list(ulcerative=cbind(t(ret[3:4,]), 
                        formatDouble(counts2[3,1:4]/counts2[4,1:4],2, remove.leading0 = F)
        ),
       intact=cbind(t(ret[1:2,]), 
                    formatDouble(counts2[1,1:4]/counts2[2,1:4],2, remove.leading0 = F) 
       )
  )
}) 


################################################################################
# table combining biopsies

for (tissue.type in tissue.types) {
  tab=rbind(
    round(sapply(biopsies.1, function(b) unname(tissue.widths[[b]][tissue.type]))),
    round(sapply(biopsies.1, function(b) dist.bt.peaks[[b]][tissue.type,"dist2skin"])),
    round(sapply(biopsies.1, function(b) dist.bt.peaks[[b]][tissue.type,"dist2DEJ"])),
    cbind(c(t(count.tables$DD[[tissue.type]])), 
          c(t(count.tables$CF[[tissue.type]])),
          c(t(count.tables$JM[[tissue.type]])))
  )
  
  rownames(tab)=c(
    "tissue width",
    "CD3 to CD15 (dist2skin)",
    "CD3 to CD15 (dist2DEJ)",
    c(t(outer(markers.1, c("epi","dermis","e/d"), paste)))
  )
  tab
  
  # add ratio cols
  tab=cbind(tab, 
            "CF/DD" = formatDouble(as.numeric(tab[,2])/as.numeric(tab[,1]),2,remove.leading0 = F),
            "JM/DD" = formatDouble(as.numeric(tab[,3])/as.numeric(tab[,1]),2,remove.leading0 = F))
  tab
  
  mytex(tab, file=glue("tables/{tissue.type}_comparison"), 
        add.to.row=list(as.list(seq(3,by=3,length=4)),
                        rep("\\hline\n",4)
        )
  )  
}
    




