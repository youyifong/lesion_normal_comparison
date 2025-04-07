################################################################################
## tabulate cell counts in epi/dermis x normal/lesion

dir.exists("output/tables") || dir.create("output/tables", recursive = TRUE)

tissue_size = read.csv(glue("{data.dir}/tissue_size_by_depth_and_region.csv"))

# count cells to different depths
count.tables = lapply (biopsies, function(s) {
# s="DD/DD1_Les"
  
  s1 = sub(".*?/(.*)", "\\1", s)
  
  markers=markers.ls[[s]]; names(markers)=markers
  
  for (max_depth in c(500, 1000)) { # um 
    myprint(s, max_depth)
    
    counts = sapply(markers, function (a) 
      c(nrow(subset(alldat[[s]][[a]], dist2DEJ<max_depth & epi &    region=="intact")),
        nrow(subset(alldat[[s]][[a]], dist2DEJ<max_depth & dermis & region=="intact")),
        nrow(subset(alldat[[s]][[a]], dist2DEJ<max_depth & epi &    region=="ulcerative")),
        nrow(subset(alldat[[s]][[a]], dist2DEJ<max_depth & dermis & region=="ulcerative")))
    )
    
    
    # compute sizes of epi/dermis, normal/lesion
    dat = tissue_size[tissue_size$biopsy==s & tissue_size$dist_type=="dist2DEJ",] # max_depth is relative to DEJ
    dat = dat[dat$depth<max_depth,] # < b/c depth is the beginning of the bin
    sizes = with(dat, c(sum(area[depth<0 & region=="intact"]),
                        sum(area[depth>0 & region=="intact"]),
                        sum(area[depth<0 & region=="ulcerative"]),
                        sum(area[depth>0 & region=="ulcerative"])))
    # convert size from 100um x 100um to 1mm x 1mm
    sizes = sizes/100
    
    # deprecated because it does not handle a depth limit  
    # # get size from XX_region_area_length.csv
    # dat=read.csv(glue("{data.dir}/{s}/TSP_stitched/{s1}_region_area_length.csv"))
    # dat = dat[!contain(dat$name, "length"), ] # remove rows that are lengths
    # dat$epi = contain(dat$name, "epi")
    # dat$dermis = contain(dat$name, "dermis")
    # dat$normal = contain(dat$name, "normal")
    # dat$lesion = contain(dat$name, "lesion")
    # sizes = with(dat, c(sum(size[epi&normal]), 
    #                     sum(size[dermis & normal]), 
    #                     sum(size[epi & lesion]), 
    #                     sum(size[dermis & lesion])
    # ))
    # sizes = sizes/4/1e6 # convert to mm
    
    counts2=counts/sizes
    mode(counts2) <- "integer"
    
    ## compute depth-normalized counts/mm2
    
    # compute sizes of epi/dermis, normal/lesion
    dat = tissue_size[tissue_size$biopsy==s & tissue_size$dist_type=="dist2DEJ",] # max_depth is relative to DEJ
    dat = dat[dat$depth<max_depth,] # use < instead of <= b/c depth is the beginning of the bin
    dat$area = dat$area/100    # convert size from 100um x 100um to 1mm x 1mm
    size_epi_normal = dat[dat$depth<0 & dat$region=="intact",]
    size_epi_lesion = dat[dat$depth<0 & dat$region=="ulcerative",]
    size_dermis_normal = dat[dat$depth>0 & dat$region=="intact",]
    size_dermis_lesion = dat[dat$depth>0 & dat$region=="ulcerative",]
    
    # first compute depth-specific density
    counts3 = sapply(markers, function (a) {
      # epi normal
      epi_normal = subset(alldat[[s]][[a]], dist2DEJ<max_depth & epi &    region=="intact")
      # add a depth column computed as the floor of dist2DEJ/100
      epi_normal$depth = floor(epi_normal$dist2DEJ/100)*100
      count_by_depth = mytable(epi_normal$depth)
      # add count_by_depth column to size_epi_normal by match the depth in size_epi_normal and count_by_depth
      size_epi_normal$count_by_depth = count_by_depth[size_epi_normal$depth%.%""]
      epi_normal_density = mean(size_epi_normal$count_by_depth/size_epi_normal$area, na.rm=TRUE)
      
      # dermis normal
      dermis_normal = subset(alldat[[s]][[a]], dist2DEJ<max_depth & dermis & region=="intact")
      # add a depth column computed as the floor of dist2DEJ/100
      dermis_normal$depth = floor(dermis_normal$dist2DEJ/100)*100
      count_by_depth = mytable(dermis_normal$depth)
      # add count_by_depth column to size_dermis_normal by match the depth in size_dermis_normal and count_by_depth
      size_dermis_normal$count_by_depth = count_by_depth[size_dermis_normal$depth%.%""]
      dermis_normal_density = mean(size_dermis_normal$count_by_depth/size_dermis_normal$area, na.rm=TRUE)
      
      # epi lesion
      epi_lesion = subset(alldat[[s]][[a]], dist2DEJ<max_depth & epi &    region=="ulcerative")
      # add a depth column computed as the floor of dist2DEJ/100
      epi_lesion$depth = floor(epi_lesion$dist2DEJ/100)*100
      count_by_depth = mytable(epi_lesion$depth)
      # add count_by_depth column to size_epi_lesion by match the depth in size_epi_lesion and count_by_depth
      size_epi_lesion$count_by_depth = count_by_depth[size_epi_lesion$depth%.%""]
      epi_lesion_density = mean(size_epi_lesion$count_by_depth/size_epi_lesion$area, na.rm=TRUE)
      
      # dermis lesion
      dermis_lesion = subset(alldat[[s]][[a]], dist2DEJ<max_depth & dermis & region=="ulcerative")
      # add a depth column computed as the floor of dist2DEJ/100
      dermis_lesion$depth = floor(dermis_lesion$dist2DEJ/100)*100
      count_by_depth = mytable(dermis_lesion$depth)
      # add count_by_depth column to size_dermis_lesion by match the depth in size_dermis_lesion and count_by_depth
      size_dermis_lesion$count_by_depth = count_by_depth[size_dermis_lesion$depth%.%""]
      dermis_lesion_density = mean(size_dermis_lesion$count_by_depth/size_dermis_lesion$area, na.rm=TRUE)
      
      # return the counts
      c(epi_normal_density, dermis_normal_density, epi_lesion_density, dermis_lesion_density)
    })
    mode(counts3) <- "integer"
    
    
    tab=t(rbind(counts, counts2, counts3))
    mode(tab) <- "character" # convert to string
    tab=data.frame(tab)
    # add area
    tab=rbind(tab, area=formatDouble(c(sizes, rep(NA,8)), digits = 2, remove.leading0=F))
    print(tab)
    
    mytex(data.frame(tab), file.name="output/tables/"%.%s1%.%"_counts_"%.%max_depth, digit=0,  include.colnames = F, align=c("l","r","r","r","r|","r","r","r","r|","r","r","r","r"),
          floating=T, caption=escapeUnderline(s1)%.%". All measurements are up to "%.%max_depth%.%" um below DEJ. Area is measured in mm2. Density in count per mm2. Regularized density refers to the average of positive density measurements taken at 100um DEJ depth intervals, i.e. each depth is weighted equally instead of based on tissue shape.",
          col.headers="\\hline\n 
           \\multicolumn{1}{c}{} & \\multicolumn{4}{c|}{Counts} & \\multicolumn{4}{c|}{Density} & \\multicolumn{4}{c}{Regularized Density} \\\\ 
           \\multicolumn{1}{c}{} & \\multicolumn{2}{c}{normal} & \\multicolumn{2}{c|}{lesion} & \\multicolumn{2}{c}{normal} & \\multicolumn{2}{c|}{lesion} & \\multicolumn{2}{c}{normal} & \\multicolumn{2}{c}{lesion}  \\\\ 
           \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{epi} & dermis & \\multicolumn{1}{c}{epi} & \\multicolumn{1}{c|}{dermis} & \\multicolumn{1}{c}{epi} & dermis & \\multicolumn{1}{c}{epi} & \\multicolumn{1}{c|}{dermis} & \\multicolumn{1}{c}{epi} & dermis & \\multicolumn{1}{c}{epi} & dermis  \\\\ 
           \\hline
      ")
    
    # ret = counts3[,markers.1]
    # mode(ret) = "character"
    # list(ulcerative=cbind(t(ret[3:4,]), 
    #                       formatDouble(counts2[3,1:4]/counts2[4,1:4],2, remove.leading0 = F)
    #       ),
    #      intact=cbind(t(ret[1:2,]), 
    #                   formatDouble(counts2[1,1:4]/counts2[2,1:4],2, remove.leading0 = F) 
    #      )
    # )
  }
}) 


