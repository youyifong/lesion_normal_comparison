################################################################################
## tabulate cell counts in epi/dermis x intact/ulcerative

dir.exists("output/tables") || dir.create("output/tables", recursive = TRUE)

format_mixed_rounding <- function(x) {
  sapply(x, function(val) {
    if (is.nan(val)) NA else if (is.na(val)) NA else {
      if (val %% 1 == 0) {
        formatC(val, format = "f", digits = 0)  # integer: 0 decimal digits
      } else {
        formatC(val, format = "f", digits = 2)  # non-integer: 2 decimal digits
      }
    }
  })
}

# count cells to different depths
for (s in all.slides) {
# s="DD/DD1_Les"
  
  s1 = sub(".*?/(.*)", "\\1", s)
  
  markers=markers.ls[[s]]; names(markers)=markers
  
  for (max_depth in c(500, 1000)) { # um 
    myprint(s, max_depth)

    # put norm.counts first so that the row names are complete
    tab=with (count.tables[[s]][["DEJ"%.%max_depth]], t(rbind(norm.counts, densities, counts)))
    tab = apply(tab, 2, format_mixed_rounding)
    tab=data.frame(tab)
    # add area
    tab=rbind(tab, "Size of biopsy"=formatDouble(c(rep(NA,8), count.tables[[s]][["DEJ"%.%max_depth]]$sizes), digits = 2, remove.leading0=F))
    print(tab)
    
    mytex(data.frame(tab), file.name="output/tables/"%.%s1%.%"_counts_"%.%max_depth, digit=0,  include.colnames = F, align=c("l","r","r","r","r|","r","r","r","r|","r","r","r","r"),
          floating=T, caption=escapeUnderline(s1)%.%". All measurements are up to "%.%max_depth%.%" um below DEJ. The unit of area is mm2; the unit of density is count/mm2. Normalized count is cell count in a vertical strip 100um wide, i.e. each depth is weighted equally instead of variable depending on the tissue shape.",
          col.headers="\\hline\n 
           \\multicolumn{1}{c}{} & \\multicolumn{4}{c|}{Normalized Counts} & \\multicolumn{4}{c|}{Density} & \\multicolumn{4}{c}{Counts} \\\\ 
           \\multicolumn{1}{c}{} & \\multicolumn{2}{c}{intact} & \\multicolumn{2}{c|}{ulcerative} & \\multicolumn{2}{c}{intact} & \\multicolumn{2}{c|}{ulcerative} & \\multicolumn{2}{c}{intact} & \\multicolumn{2}{c}{ulcerative}  \\\\ 
           \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{epi} & dermis & \\multicolumn{1}{c}{epi} & \\multicolumn{1}{c|}{dermis} & \\multicolumn{1}{c}{epi} & dermis & \\multicolumn{1}{c}{epi} & \\multicolumn{1}{c|}{dermis} & \\multicolumn{1}{c}{epi} & dermis & \\multicolumn{1}{c}{epi} & dermis  \\\\ 
           \\hline
      ")
    
  }

}

