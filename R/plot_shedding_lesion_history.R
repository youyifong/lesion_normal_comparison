# this file is run from local drive and the data resides from the X drive

# JoinFile_20190719_JZ_KP.xlsx
# Column K is the PCR data of the entire genital sweep, independent of the ulcer lesion (I &J), trying to get an overall shedding status in the genital area.
# Column I is genital_lesion
# Column L is the PCR data for the lesion site from which we collected the biopsy.
# Columns N & O recorded a new lesion site and PCR data.
# Column P counts the days for positive PCR; it could come from different lesion sites, as long as that day is HSV PCR positive.
# Column Q based on the patient’s diary, recorded in Column J and N.
# Column R combines both data from P and Q to count for HSV+ activity for that day. We used this to calculate HSV activity rate.

{
if (.Platform$OS.type == "windows") {
  base_dir = 'X:/fast/fong_y/images_from_K'
} else {
  base_dir = '..'
}
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggnewscale)  # for multiple color scales
library(readxl)
library(kyotil)
library(ggplot2)
library(dplyr)
library(tidyr)
dir.create("output", showWarnings = FALSE)
dir.create("output/by_ptid", showWarnings = FALSE)
}

{
## cannot link ID from this source to ptid im the JoinFile
# # read "all days & les epi" sheet from HSV406_HSV406M_comprehensive_data_04172024_2.xlsx
# tmp <- read_excel("X:/fast/fong_y/images_from_K/HSV406/HSV406_HSV406M_comprehensive_data_04172024_2.xlsx", sheet = "all days & les epi")
# tmp = tmp[,1]
# tmp = unlist(c(tmp))
# ptids_fullstudy = tmp[(1+last(which(is.na(tmp)))):length(tmp)]
# print(len(ptids_fullstudy))

# take out 13555 for HIV+
IDs.to.rm = c(13484, 14864, 14862, 14259, 14981, 13369, 13373, 15396, 15402, 13555)

# specify the file path (replace with your actual file path)
file_path <- base_dir%.%"/HSV406/JoinFile_20190719_JZ_KP.xlsx"

# read the data from the "Shedding Rate" sheet
shedding_data <- read_excel(file_path, sheet = "Shedding Rate")

# view the first few rows
str(shedding_data)
# print the number of unique IDs
length(unique(shedding_data$ID))
# print the unique values of site_lesion and new_gen_lesion
unique(shedding_data$site_lesion) # column J
unique(shedding_data$new_gen_lesion) # column N
# check BUT to BUTT in shedding_data$new_gen_lesion
shedding_data$new_gen_lesion[shedding_data$new_gen_lesion == "BUT"] <- "BUTT"
# remove data for IDs in IDs.to.rm
shedding_data <- shedding_data[!shedding_data$ID %in% IDs.to.rm, ]
length(unique(shedding_data$ID))

# standardize column J lesion 
# site groups:
# "BUTT"  
# "PERI"  "PRNL" 
# "LABIA"  "MONS" Female genital area
# "GEN" 
# "LEG"

# column J
mytable(shedding_data$site_lesion)
shedding_data$site_lesion <- case_when(
  shedding_data$site_lesion %in% c("BUTT") ~ "BUTT",
  shedding_data$site_lesion %in% c("PERI", "PRNL") ~ "PRNL",
  shedding_data$site_lesion %in% c("LABIA", "MONS") ~ "FEM_GEN",
  shedding_data$site_lesion %in% c("GEN") ~ "GEN",
  shedding_data$site_lesion %in% c("LEG") ~ "LEG"
)
mytable(shedding_data$site_lesion)

colG = shedding_data[[which(LETTERS == "G")]]
colK = shedding_data[[which(LETTERS == "K")]]
colL = shedding_data[[which(LETTERS == "L")]]
colO = shedding_data[[which(LETTERS == "O")]]
MaxGKLO = pmax(colG, colK, colL, colO, na.rm = TRUE)



# column N
mytable(shedding_data$new_gen_lesion)
shedding_data$new_gen_lesion <- case_when(
  shedding_data$new_gen_lesion %in% c("BUTT") ~ "BUTT",
  shedding_data$new_gen_lesion %in% c("PERI", "PRNL") ~ "PRNL",
  shedding_data$new_gen_lesion %in% c("LABIA", "MONS") ~ "FEM_GEN",
  shedding_data$new_gen_lesion %in% c("GEN") ~ "GEN",
  shedding_data$new_gen_lesion %in% c("LEG") ~ "LEG"
)
mytable(shedding_data$new_gen_lesion)

# used to sort ptids 
activity_data <- read_excel(file_path, sheet = "HSV Activity Site")
# keep only two columns: ID, With All Day
activity_data1 <- activity_data[, c("ID", "With All Day")]
# remove rows with NA
activity_data1 <- activity_data1[!is.na(activity_data1$`With All Day`), ]
str(activity_data1)
# order activity_data1 by With All Day in descending order
activity_data1 <- activity_data1[order(-activity_data1$`With All Day`), ]
# remove data for IDs in IDs.to.rm
activity_data1 <- activity_data1[!activity_data1$ID %in% IDs.to.rm, ]
print(activity_data1, n=100)

}



{
# plot the overall shedding status in the genital area (column K)

value=round (shedding_data[[which(LETTERS == "K")]])
shedding_plot_data <- subset(shedding_data, !is.na(value))
value=value[!is.na(value)]
mytable(value)

# make a plot where x axis is day. Each row corresponds to one ID: plot a dot with varying degrees of red that correspond to the value of shedding_data[[which(LETTERS == "K")]], and a white dot if it is NA
ggplot(shedding_plot_data, aes(x = day, y = factor(ID, levels = activity_data1$ID))) +
  geom_point(aes(color = value), size = 3, shape = 16) +
  scale_color_gradient(low = "white", high = "red") +
  labs(x = "Day", y = "ID", color = "Genital Lesion") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 6)) +
  theme(axis.title.y = element_blank()) + xlim(0, NA)
ggsave("output/clinical_outcome_profile_col_K.png", width = 10, height = 8)


# plot lesions

# primary lesion site

value2=factor(shedding_data[[which(LETTERS == "J")]])
shedding_plot_data2 <- subset(shedding_data, !is.na(value2))
value2=value[!is.na(value2)]
mytable(value2)

ggplot(shedding_plot_data2, aes(x = day, y = factor(ID))) +
  geom_point(aes(color = as.factor(value2)), size = 3, shape = 16, na.rm = TRUE) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Day", y = "ID", color = "Genital Lesion") +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 6),
    axis.title.y = element_blank()
  ) +
  xlim(0, NA)

ggsave("output/clinical_outcome_profile_col_J.png", width = 10, height = 8)
}




# combine copy number and lesion in the same plot
# two iterations, first all ptids, second only ptids in correlates analysis
for (iter in 1:2) {
  if (iter==1) {
    IDs=activity_data1$ID
  } else {
    # 22 ptids for correlates analysis
    IDs=c(15052, 14655, 14995, 15094, 13581, 9149, 15206, 14999, 14989, 14866, 15017, 13279, 10110, 15368, 14341, 14939, 15328, 14883, 15146, 14900, 14954, 14888)
    # sort IDs in the order they appear in activity_data1$ID
    IDs <- activity_data1$ID[sort(match(IDs, activity_data1$ID))]
  }
  
  
  # Extract values
  # value  <- round(shedding_data[[which(LETTERS == "K")]])
  value  <- MaxGKLO
  value2 <- factor(shedding_data[[which(LETTERS == "J")]])
  value3 <- factor(shedding_data[[which(LETTERS == "N")]])
  
  # Prepare datasets
  df_value <- shedding_data %>%
    mutate(
      variable = "value",
      val_num = value,
      val_chr = NA_character_,
      y_pos = as.numeric(factor(ID, levels = IDs))
    )
  
  df_value2 <- shedding_data %>%
    mutate(
      variable = "value2",
      val_num = NA_real_,
      val_chr = as.character(value2),
      y_pos = as.numeric(factor(ID, levels = IDs))
    )
  
  df_value3 <- shedding_data %>%
    mutate(
      variable = "value3",
      val_num = NA_real_,
      val_chr = as.character(value3),
      y_pos = as.numeric(factor(ID, levels = IDs))
    )
  
    
  # Combine
  shedding_plot_data_combined <- bind_rows(df_value, df_value2, df_value3)
  
  # CATEGORICAL AS FOREGROUND SYMBOLS (no background tiles)
  df_cat <- subset(shedding_plot_data_combined, variable == "value2")
  df_cat$val_chr <- factor(df_cat$val_chr, levels=c("BUTT", "PRNL", "FEM_GEN", "GEN", "LEG"))  # ensure factor
  df_cat3 <- subset(shedding_plot_data_combined, variable == "value3")
  df_cat3$val_chr <- factor(df_cat3$val_chr, levels=c("BUTT", "PRNL", "FEM_GEN", "GEN", "LEG"))  # ensure factor

  # column N first, then J, so that primary site will not be overridden by new site;
  # open and line glyphs are see-through, so a day with two sites shows both
  df_sites <- rbind(df_cat3, df_cat)
  # ggplot cannot mix a one-character shape with named shapes in one layer, so the dagger gets its own layers
  is_dagger <- df_sites$val_chr %in% "LEG"

  # shapes chosen to stay distinct when small (circle vs square was not); the dagger is a font glyph,
  # so stroke has no effect on it and it needs a larger size
  site_shapes <- c(BUTT="circle open", PRNL="plus", FEM_GEN="triangle open", GEN="triangle down open", LEG="\u2020")
  site_sizes  <- c(BUTT=.6, PRNL=.6, FEM_GEN=.6, GEN=.6, LEG=1.4)
  site_strokes <- c(BUTT=.25, PRNL=.5, FEM_GEN=.25, GEN=.25, LEG=.5)  # thinner outlines for circle and triangles
  
  
  ggplot() +
    
  # ------------------------------------------------------------------
  # Foreground dots for numeric variable (value) 
  geom_point(
    data = subset(shedding_plot_data_combined, variable == "value"),
    aes(x = day, y = y_pos, color = val_num),
    size = 3, shape = 16, na.rm = TRUE
  ) +
  scale_color_gradient(
    low = "white", high = "firebrick4", na.value = "transparent",
    name = "Copy number (log10)",
    guide = guide_colourbar(theme = theme(legend.key.width = unit(8, "lines")))  # long bar in the horizontal legend
  ) +

  # ------------------------------------------------------------------
  # Categorical variable plotted as black glyphs; the copy-number color shows through the open shapes
  geom_point(
    data = df_sites[!is_dagger, ],
    aes(x = day, y = y_pos, shape = val_chr, size = val_chr, stroke = val_chr),
    colour = "black",
    na.rm = TRUE
  ) +
  geom_point(
    data = df_sites[is_dagger, ],
    aes(x = day, y = y_pos, shape = val_chr, size = val_chr, stroke = val_chr),
    colour = "black",
    na.rm = TRUE
  ) +
  scale_shape_manual(
    values = site_shapes,
    name   = "Lesion site",
    drop   = FALSE,
    na.translate = FALSE   # <--- HIDE NA from the legend
  ) +
  scale_size_manual(
    values = site_sizes,
    name   = "Lesion site",  # same name as the shape scale, so the scales merge into one legend
    drop   = FALSE,
    na.translate = FALSE
  ) +
  scale_discrete_manual(
    aesthetics = "stroke",
    values = site_strokes,
    name   = "Lesion site",
    drop   = FALSE,
    na.translate = FALSE
  ) +
  
  ggnewscale::new_scale_color() +  # keep this (now switching cleanly to numeric layer)
  
  # ------------------------------------------------------------------
  # Layout and theme
  labs(x = "Day", y = "ID") +
    # theme_minimal() +
    theme_gray() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",   # copy-number bar and lesion-site key side by side under the panel
      legend.title.position = "top",
      axis.text.y = element_text(size = 6),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    scale_y_continuous(
      breaks = seq_along(IDs),
      labels = IDs,
      expand = expansion(mult = 0.02)
    ) +
    scale_x_continuous(
      limits = c(1, NA),                         # Force axis to start at 1
      breaks = function(x) seq(1, max(x), by=10), # Optional: Force ticks at every day (1, 2, 3...)
      # expand = expansion(mult = c(0, 0.05))      # Removes the empty padding before Day 1
    )  
  ggsave(glue("output/clinical_outcome_profile_combined{ifelse(iter==2, '_22', '')}.png"), width = 10, height = 8*ifelse(iter==1,1,22/37))
  
}



# plot for each ptid

unique_IDs <- unique(shedding_data$ID)
out = mysapply (unique_IDs, function(ptid) {
  
  ptid_data <- shedding_data[shedding_data$ID == ptid, ]
  ptid_data$MaxGKLO = pmax(
    ptid_data[[which(LETTERS == "K")]],
    ptid_data[[which(LETTERS == "L")]],
    ptid_data[[which(LETTERS == "O")]],
        na.rm = TRUE
      )
    
  
  # Extract values
  value4 <- factor(ptid_data[[which(LETTERS == "J")]])# site lesion
  value35<-  round(ptid_data[[which(LETTERS == "L")]]) # pcr bx lesion, real num
  value3 <-  round(ptid_data[[which(LETTERS == "K")]]) # genital sweep, real num
  value2 <- factor(ptid_data[[which(LETTERS == "N")]])# new lesion
  value1 <-  round(ptid_data[['MaxGKLO']]) # real num
  
  pcr_positivity = mean(value1>0, na.rm=TRUE)
  len_study = nrow(subset(ptid_data, day>=1))
  lesion_positivity = sum(!is.na(value4) | !is.na(value2))/len_study
  
  # Prepare datasets
  df_value4 <- ptid_data %>%
    mutate(
      variable = "value4",
      val_num = NA_real_,
      val_chr = as.character(value4),
      y_pos = 5
    )
  
  df_value35 <- ptid_data %>%
    mutate(
      variable = "value35",
      val_num = value35,
      val_chr = as.character(value35),
      y_pos = 4
    )
  
  df_value3 <- ptid_data %>%
    mutate(
      variable = "value3",
      val_num = value3,
      val_chr = as.character(value3),
      y_pos = 3
    )
  
  df_value2 <- ptid_data %>%
    mutate(
      variable = "value2",
      val_num = NA_real_,
      val_chr = as.character(value2),
      y_pos = 2
    )
  
  df_value1 <- ptid_data %>%
    mutate(
      variable = "value1",
      val_num = value1,
      val_chr = as.character(value1),
      y_pos = 1
    )
  
  # Combine
  shedding_plot_data_combined <- bind_rows(df_value1, df_value2, df_value3, df_value35, df_value4)
  
  # CATEGORICAL AS HOLLOW FOREGROUND SYMBOLS (no background tiles)
  df_cat <- subset(shedding_plot_data_combined, variable == "value2")
  df_cat$val_chr <- factor(df_cat$val_chr, levels=c("BUTT", "PRNL", "FEM_GEN", "GEN", "LEG"))  # ensure factor
  df_cat3 <- subset(shedding_plot_data_combined, variable == "value4")
  df_cat3$val_chr <- factor(df_cat3$val_chr, levels=c("BUTT", "PRNL", "FEM_GEN", "GEN", "LEG"))  # ensure factor
  
  # pick hollow shapes (21–25 support outline + fill)
  shape_vals <- c(21, 22, 24, 25, 23, 21)
  
  
  
  #### -------------------------------------------------------------------------------------

  p = ggplot() +
    
  # Foreground dots for numeric variable (value) 
    
  geom_point(
    data = subset(shedding_plot_data_combined, variable == "value1"),
    aes(x = day, y = y_pos, color = val_num),
    size = 3, shape = 16, na.rm = TRUE
  ) +
  # scale_color_gradient(
  #   low = "white", high = "firebrick4", na.value = "transparent",
  #   name = "Copy number (log10)"
  # ) +
    
  geom_point(
    data = subset(shedding_plot_data_combined, variable == "value3"),
    aes(x = day, y = y_pos, color = val_num),
    size = 3, shape = 16, na.rm = TRUE
  ) +
  # scale_color_gradient(
  #   low = "white", high = "firebrick4", na.value = "transparent",
  #   name = "Copy number (log10)"
  # ) +
  
  geom_point(
    data = subset(shedding_plot_data_combined, variable == "value35"),
    aes(x = day, y = y_pos, color = val_num),
    size = 3, shape = 16, na.rm = TRUE
  ) +
  scale_color_gradient(
    low = "white", high = "firebrick4", na.value = "transparent",
    name = "Copy number (log10)"
  ) +
    
    
  # Categorical variable now plotted as SHAPE-only foreground symbols
    
  # column N first, so that primary site will not be overridden by new site
  geom_point(
    data = df_cat3,
    aes(x = day, y = y_pos, shape = val_chr),
    size = 2,
    stroke =1,
    fill = NA,          # <-- hollow centers (transparent)
    colour = "black",   # <-- outline color fixed (not mapped)
    na.rm = TRUE
  ) +
  geom_point(
    data = df_cat,
    aes(x = day, y = y_pos, shape = val_chr),
    size = 2,
    stroke =1,
    fill = NA,          # <-- hollow centers (transparent)
    colour = "black",   # <-- outline color fixed (not mapped)
    na.rm = TRUE
  ) +
  scale_shape_manual(
    values = shape_vals,
    name   = "",
    drop   = FALSE,
    na.translate = FALSE   # <--- HIDE NA from the legend
  ) +
  
  # ggnewscale::new_scale_color() +  # keep this (now switching cleanly to numeric layer)
  
  
  # Layout and theme
  labs(x = "Day", 
       y = NULL,
       title = ptid
      ) +
    # theme_minimal() +
    theme_gray() +
    theme(
      plot.title = element_text(hjust = 0, size = 14),
      legend.position = "right",
      legend.text  = element_text(size = 7),  
      legend.title = element_text(size = 8, face = "bold"),
      axis.text.y  = element_text(size = 8, color = "black"),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(b = 10, t=10, unit = "pt"),
      panel.grid.major.y = element_blank()
    ) +
    scale_y_continuous(
      breaks = 5:1, 
      labels = c("site lesion", "PCR bx lesion", "PCR genital sweep", "new gen lesion", "max KLO")
    ) +    
    scale_x_continuous(
      limits = c(1, 76),                         # Force axis to start at 1
      breaks = function(x) seq(1, max(x), by=10), # Optional: Force ticks at every day (1, 2, 3...)
      # expand = expansion(mult = c(0, 0.05))      # Removes the empty padding before Day 1
    )  
  
  ggsave(paste0("output/by_ptid/", ptid, ".png"), plot = p, width = 8, height = 3.1)
  
  # return 
  c(ptid, pcr_positivity, lesion_positivity) 
})


colnames(out) = c("ptid", "PCR.positivity.rate", "lesion.positivity.rate")
out


severity = read.csv(base_dir %.% "/HSV406/Lesion_HSV_rates_04082025.csv")
severity$ptid = substr(severity$PID, 3, 100)

final_out = merge (out, severity, by.x="ptid", by.y="ptid", all.x=TRUE)

write.csv(final_out, base_dir%.%"/HSV406/Lesion_HSV_rates_04082025_plus.csv", row.names = FALSE)


mypdf("output/correlations_activity_rates", width=6, height=6)
  mypairs(final_out[,c("HSV.activity.rate", "PCR.positivity.rate", "lesion.positivity.rate")])
dev.off()
