library(readxl)
library(kyotil)
library(ggplot2)

# specify the file path (replace with your actual file path)
file_path <- "X:/fast/fong_y/images_from_K/HSV406/JoinFile_20190719_JZ_KP.xlsx"

# read the data from the "Shedding Rate" sheet
shedding_data <- read_excel(file_path, sheet = "Shedding Rate")

# view the first few rows
str(shedding_data)

# print the number of unique IDs
length(unique(shedding_data$ID))

# print the unique values of site_lesion and new_gen_lesion
unique(shedding_data$site_lesion)
unique(shedding_data$new_gen_lesion)

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
mytable(value)

# make a plot where x axis is day. Each row corresponds to one ID: plot a dot with varying degrees of red that correspond to the value of shedding_data[[which(LETTERS == "K")]], and a white dot if it is NA
ggplot(shedding_data, aes(x = day, y = factor(ID))) +
  geom_point(aes(color = value), size = 3, shape = 16) +
  scale_color_gradient(low = "white", high = "red", na.value = "gray") +
  labs(x = "Day", y = "ID", color = "Genital Lesion") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 6)) +
  theme(axis.title.y = element_blank()) + xlim(0, NA)
ggsave("output/clinical_outcome_profile_col_K.png", width = 10, height = 8)


################################################################################
# show all four episodes. use shades of color to reflect HSV copy number, and different colors for each episode 

# first, show the primary lesion shedding status

value=round (shedding_data[[which(LETTERS == "L")]])
mytable(value)

# make a plot where x axis is day. Each row corresponds to one ID: plot a dot with varying degrees of red that correspond to the value of shedding_data[[which(LETTERS == "K")]], and a white dot if it is NA
ggplot(shedding_data, aes(x = day, y = factor(ID))) +
  geom_point(aes(color = value), size = 3, shape = 16) +
  scale_color_gradient(low = "white", high = "red", na.value = "gray") +
  labs(x = "Day", y = "ID", color = "Genital Lesion") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 6)) +
  theme(axis.title.y = element_blank()) + xlim(0, NA)
ggsave("output/clinical_outcome_profile_col_L.png", width = 10, height = 8)



value=round (shedding_data[[which(LETTERS == "L")]])
mytable(value)


# make a plot where x axis is day. Each row corresponds to one ID: plot a dot with varying degrees of red that correspond to the value of shedding_data[[which(LETTERS == "K")]], and a white dot if it is NA
ggplot(shedding_data, aes(x = day, y = factor(ID))) +
  geom_point(aes(color = value), size = 3, shape = 16) +
  scale_color_gradient(low = "white", high = "red", na.value = "gray") +
  labs(x = "Day", y = "ID", color = "Genital Lesion") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 6)) +
  theme(axis.title.y = element_blank()) + xlim(0, NA)
ggsave("output/clinical_outcome_profile_col_K.png", width = 10, height = 8)


# show column N and O and only for those rows where either is not NA
value=cbind (shedding_data$ID, shedding_data[[which(LETTERS == "N")]], shedding_data[[which(LETTERS == "O")]])
value[is.na(value[,2]) & !is.na(value[,3]), ]


