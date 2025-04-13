################################################################################
## tabulate cell counts in epi/dermis x intact/ulcerative

dir.exists("output/tables") || dir.create("output/tables", recursive = TRUE)

# count cells to different depths
for (s in all.slides) {
# s="DD/DD1_Les"
  
  s1 = sub(".*?/(.*)", "\\1", s)
  
  markers=markers.ls[[s]]; names(markers)=markers
  
  for (max_depth in c(500, 1000)) { # um 
    myprint(s, max_depth)

    tab=with (count.tables[[s]][["DEJ"%.%max_depth]], t(rbind(counts, densities, norm.densities)))
    mode(tab) <- "character" # convert to string
    tab=data.frame(tab)
    # add area
    tab=rbind(tab, area=formatDouble(c(count.tables[[s]][["DEJ"%.%max_depth]]$sizes, rep(NA,8)), digits = 2, remove.leading0=F))
    print(tab)
    
    mytex(data.frame(tab), file.name="output/tables/"%.%s1%.%"_counts_"%.%max_depth, digit=0,  include.colnames = F, align=c("l","r","r","r","r|","r","r","r","r|","r","r","r","r"),
          floating=T, caption=escapeUnderline(s1)%.%". All measurements are up to "%.%max_depth%.%" um below DEJ. Area is measured in mm2. Density in count per mm2. normalized density refers to the average of positive density measurements taken at 100um DEJ depth intervals, i.e. each depth is weighted equally instead of based on tissue shape.",
          col.headers="\\hline\n 
           \\multicolumn{1}{c}{} & \\multicolumn{4}{c|}{Counts} & \\multicolumn{4}{c|}{Density} & \\multicolumn{4}{c}{normalized Density} \\\\ 
           \\multicolumn{1}{c}{} & \\multicolumn{2}{c}{intact} & \\multicolumn{2}{c|}{ulcerative} & \\multicolumn{2}{c}{intact} & \\multicolumn{2}{c|}{ulcerative} & \\multicolumn{2}{c}{intact} & \\multicolumn{2}{c}{ulcerative}  \\\\ 
           \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{epi} & dermis & \\multicolumn{1}{c}{epi} & \\multicolumn{1}{c|}{dermis} & \\multicolumn{1}{c}{epi} & dermis & \\multicolumn{1}{c}{epi} & \\multicolumn{1}{c|}{dermis} & \\multicolumn{1}{c}{epi} & dermis & \\multicolumn{1}{c}{epi} & dermis  \\\\ 
           \\hline
      ")
    
  }

}

