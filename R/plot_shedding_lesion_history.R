{
library(readxl)
library(kyotil)
library(ggplot2)
library(dplyr)
library(tidyr)

## cannot link ID from this source to ptid im the JoinFile
# # read "all days & les epi" sheet from HSV406_HSV406M_comprehensive_data_04172024_2.xlsx
# tmp <- read_excel("X:/fast/fong_y/images_from_K/HSV406/HSV406_HSV406M_comprehensive_data_04172024_2.xlsx", sheet = "all days & les epi")
# tmp = tmp[,1]
# tmp = unlist(c(tmp))
# ptids_fullstudy = tmp[(1+last(which(is.na(tmp)))):length(tmp)]
# print(len(ptids_fullstudy))

IDs.to.rm = c(13484, 14864, 14862, 14259, 14981, 13369, 13373, 15396, 15402)

# specify the file path (replace with your actual file path)
file_path <- "X:/fast/fong_y/images_from_K/HSV406/JoinFile_20190719_JZ_KP.xlsx"

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

# Column K is the PCR data of the entire genital sweep, independent of the ulcer lesion (I &J), trying to get an overall shedding status in the genital area.
# Column I is genital_lesion
# Column L is the PCR data for the lesion site from which we collected the biopsy.
# Columns N & O recorded a new lesion site and PCR data.
# Column P counts the days for positive PCR; it could come from different lesion sites, as long as that day is HSV PCR positive.
# Column Q based on the patient’s diary, recorded in Column J and N.
# Column R combines both data from P and Q to count for HSV+ activity for that day. We used this to calculate HSV activity rate.


################################################################################
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


################################################################################
# plot lesions

# primary lesion site

value2=factor(shedding_data[[which(LETTERS == "J")]])
shedding_plot_data2 <- subset(shedding_data, !is.na(value2))
value2=value[!is.na(value2)]
mytable(value2)

# make a plot where x axis is day. Each row corresponds to one ID: plot a dot with varying degrees of red that correspond to the value of shedding_data[[which(LETTERS == "K")]], and a white dot if it is NA
ggplot(shedding_plot_data2, aes(x = day, y = factor(ID))) +
  geom_point(aes(color = value2), size = 3, shape = 16, na.rm = TRUE) +
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


################################################################################
# combine copy number and lesion
{
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(ggnewscale)  # for multiple color scales
  
  # Extract values
  value  <- round(shedding_data[[which(LETTERS == "K")]])
  value2 <- factor(shedding_data[[which(LETTERS == "J")]])
  value3 <- factor(shedding_data[[which(LETTERS == "N")]])
  
  # Prepare datasets
  df_value <- shedding_data %>%
    mutate(
      variable = "value",
      val_num = value,
      val_chr = NA_character_,
      y_pos = as.numeric(factor(ID, levels = activity_data1$ID))
    )
  
  df_value2 <- shedding_data %>%
    mutate(
      variable = "value2",
      val_num = NA_real_,
      val_chr = as.character(value2),
      y_pos = as.numeric(factor(ID, levels = activity_data1$ID))
    )
  
  df_value3 <- shedding_data %>%
    mutate(
      variable = "value3",
      val_num = NA_real_,
      val_chr = as.character(value3),
      y_pos = as.numeric(factor(ID, levels = activity_data1$ID))
    )
  
  # Combine
  shedding_plot_data_combined <- bind_rows(df_value, df_value2, df_value3)
  
  # CATEGORICAL AS HOLLOW FOREGROUND SYMBOLS (no background tiles)
  df_cat <- subset(shedding_plot_data_combined, variable == "value2")
  df_cat$val_chr <- factor(df_cat$val_chr, levels=c("BUTT", "PRNL", "FEM_GEN", "GEN", "LEG"))  # ensure factor
  df_cat3 <- subset(shedding_plot_data_combined, variable == "value3")
  df_cat3$val_chr <- factor(df_cat3$val_chr, levels=c("BUTT", "PRNL", "FEM_GEN", "GEN", "LEG"))  # ensure factor
  
  # pick hollow shapes (21–25 support outline + fill)
  shape_vals <- c(21, 22, 24, 25, 23, 21)
  
  
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
    name = "Copy number (log10)"
  ) +
    
  #   # ------------------------------------------------------------------
  # # Background rectangles for categorical variable (value2)
  # geom_tile(
  #   data = subset(shedding_plot_data_combined, variable == "value2"),
  #   aes(x = day, y = y_pos, fill = val_chr),
  #   width = 1, height = 0.8, na.rm = TRUE
  # ) +
  #   scale_fill_manual(
  #     values = RColorBrewer::brewer.pal(6, "Set1"),
  #     na.value = "transparent",
  #     name = "Lesion site"
  #   ) +
  # 
  #   ggnewscale::new_scale_color() +  # reset for next color scale
  
  
  # ------------------------------------------------------------------
  # Categorical variable now plotted as SHAPE-only foreground symbols
  # column N first, so that primary site will not be overridden by new site
  geom_point(
    data = df_cat3,
    aes(x = day, y = y_pos, shape = val_chr),
    size = .6,
    stroke =.5,
    fill = NA,          # <-- hollow centers (transparent)
    colour = "black",   # <-- outline color fixed (not mapped)
    na.rm = TRUE
  ) +
  geom_point(
    data = df_cat,
    aes(x = day, y = y_pos, shape = val_chr),
    size = .6,
    stroke =.5,
    fill = NA,          # <-- hollow centers (transparent)
    colour = "black",   # <-- outline color fixed (not mapped)
    na.rm = TRUE
  ) +
  scale_shape_manual(
    values = shape_vals,
    name   = "Lesion site",
    drop   = FALSE,
    na.translate = FALSE   # <--- HIDE NA from the legend
  ) +
  
  ggnewscale::new_scale_color() +  # keep this (now switching cleanly to numeric layer)
  
  # ------------------------------------------------------------------
  # Layout and theme
  labs(x = "Day", y = "ID") +
    theme_minimal() +
    theme(
      legend.position = "right",
      axis.text.y = element_text(size = 6),
      axis.title.y = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    scale_y_continuous(
      breaks = seq_along(activity_data1$ID),
      labels = activity_data1$ID,
      expand = expansion(mult = 0.02)
    ) +
    xlim(0, NA)
  
  ggsave("output/clinical_outcome_profile_combined.png", width = 10, height = 8)
  
}


# add column N lesions