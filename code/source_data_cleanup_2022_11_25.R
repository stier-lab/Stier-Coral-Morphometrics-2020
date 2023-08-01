rm(list=ls())

# Load libraries
library(here)
library(tidyverse)
library(vegan)
library(reshape)
library(conflicted)
library(cowplot)
library(ez)
library(nls2)

conflict_prefer("rename","dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("select","dplyr")
conflict_prefer("plot_grid","cowplot")
conflict_prefer("get_legend", "cowplot")


# Load data 
branch_width <- read.csv(here("morphometric_data","branchwidth_data.csv")) 
coral_pg <- read.csv(here("morphometric_data","photogrammetry_data_v3_2020_9_10.csv")) 
coral_field <- read.csv(here("morphometric_data","field_experiment_colony_measurements_moorea_summer2019.csv"))
updated_cafi <- read.csv(here("cafi_data","cafi_data_w_taxonomy_summer2019_2020_5_21.csv"))

# Data Preparation and Cleaning

# CAFI data
updated_cafi2 <- updated_cafi %>% filter(str_detect(coral_id, "^FE")) %>%
  select(master_sort, coral_id, code, type, search_term, lowest_level, phylum, genus, species, general_notes)

# Calculate CAFI richness, abundance, and diversity (shannon weiner) for each coral 

cafi_summarized2 <- group_by(updated_cafi2, coral_id) %>%
  summarise(num_cafi = n(), cafi_richness = length(unique(code)), cafi_present = paste(sort(unique(code)), collapse = ";"))
# cafi_summarized2$sw <- updated_cafi2 %>% 
#   count(code, coral_id = coral_id) %>% 
#   spread(code,n) %>% 
#   mutate_all(list(~tidyr::replace_na(.,0))) %>% 
#   select(-coral_id) %>% 
#   diversity(index = "shannon")

# Clean field data
coral_field2 <- coral_field %>%
  rename(branch = branch_width) 


# Clean photogrammetry morphometric data
coral_pg$volume_pg <- as.numeric(as.character(coral_pg$volume_pg)) # Convert factor to numeric
coral_pg <- coral_pg %>% filter(use == "y") %>%   mutate(volume_pg=volume_pg*10^6, #convert m^3 to cm^3
                                                         max_hull_volume=max_hull_volume*10^6, #convert m^3 to cm^3
                                                         max_hull_surface_area=max_hull_surface_area*10^4,#convert m^2 to cm^2
                                                         surface_area=surface_area*10^4, #convert m^2 to cm^2
                                                         height_pg=height_pg*100, #convert m to cm
                                                         length_pg=length_pg*100, #convert m to cm
                                                         width_pg=width_pg*100) #convert m to cm

# Clean photogrammetry branch width data and take averages of branch distance for each coral
branch_width$branch_distance_mm <- as.numeric(as.character(branch_width$branch_distance_mm)) # Convert factor to numeric
branch_w_summarized <- group_by(branch_width, coral_id) %>%
  summarise(avg_w_cm = mean(branch_distance_mm)/10, #take average branch distance and convert mm to cm
            measurements = n(), 
            locations = paste(sort(unique(location)), collapse = ";"))

# Add photogrammetry ellipsoid volume

coral_pg <- coral_pg %>% mutate(ellipsoid_pg = (4/3)*pi*(height_pg/2)*(length_pg/2)*(width_pg/2))

#Create primary data frame for analysis

coral_dim <- merge(coral_pg, branch_w_summarized, by = "coral_id") %>%
  merge(coral_field2, by = "coral_id") %>%
  rename(max_hull_SA = max_hull_surface_area,
         SA = surface_area)

# Swap length to be the bigger linear dimension, width to be the smaller one

coral_dim <- coral_dim %>% mutate(length_pg_temp = if_else(length_pg > width_pg, length_pg, width_pg),
                                  width_pg_temp = if_else(length_pg < width_pg, length_pg, width_pg),
                                  length_field_temp = if_else(length_field > width_field, length_field, width_field),
                                  width_field_temp = if_else(length_field < width_field, length_field, width_field)) %>% 
  select(-length_pg, -width_pg, -length_field, -width_field) %>%
  rename(length_pg = length_pg_temp, width_pg = width_pg_temp,
         length_field = length_field_temp, width_field = width_field_temp)


coral_dim$interstitial_space <- coral_dim$max_hull_volume - coral_dim$volume_pg #calculate available space by subtracting software estimated volume from convex hull volume

coral_dim$SAV <- coral_dim$SA / coral_dim$volume_pg #calculate surface area to volume relationship

coral_dim$convexity <- coral_dim$volume_pg / coral_dim$max_hull_volume #calculate proportion occupied, high ratios indicate less free space between branches

coral_dim$packing <- coral_dim$SA / coral_dim$max_hull_SA #how much of an objects surface area is situated internally

coral_dim$sphericity <- ((pi^(1/3))*((6*coral_dim$volume_pg)^(2/3)))/coral_dim$SA #calculate sphericity or how close the object is to a sphere

coral_dim$diff_est <- coral_dim$volume_field-coral_dim$volume_pg #actual overestimation

coral_dim$prop_est <- (coral_dim$volume_field - coral_dim$volume_pg)/coral_dim$volume_pg #proportion overestimated 

coral_dim <- coral_dim %>% select(-"use", -position, -photos_not_aligned, -measurements, -locations) %>% 
  rename(notes_pg = notes.x, notes_field = notes.y)

coral_dim$branch_conv <- ifelse(coral_dim$convexity>=0.5, "tight","wide") #classifies wide and tight branching coral based on convexity

# Merge coral and CAFI data
cafi_coral <- left_join(coral_dim, cafi_summarized2, by ="coral_id") %>% filter(cafi == "empty")

# Growth Data
december_manual <- read.csv(here("morphometric_data", "field_experiment_colony_measurements_moorea_december2019_v2.csv"))
december_photo <- read.csv(here("morphometric_data","maatea_experiment_photo_measurements_december_JC_2020_7_8.csv"))

december_manual$coral_id <- as.character(december_manual$coral_id)
december_photo$coral_id <- as.character(december_photo$coral_id)

december_photo <- december_photo %>% filter(coral_id %in% coral_dim$coral_id, volume_pg != "NA")

#Swap width and length

december_photo <- december_photo %>% mutate(length_pg_temp = if_else(length_pg > width_pg, length_pg, width_pg),
                                  width_pg_temp = if_else(length_pg < width_pg, length_pg, width_pg)) %>% 
  select(-length_pg, -width_pg) %>%
  rename(length_pg = length_pg_temp, width_pg = width_pg_temp)

# Add August (initial) data to December DF
december_photo$vol_pg_aug <- coral_dim$volume_pg[match(december_photo$coral_id,coral_dim$coral_id)] # Photogrammetry volume
december_photo$sa_aug <-  coral_dim$SA[match(december_photo$coral_id,coral_dim$coral_id)] # PG surface area
december_photo$vol_hull_aug <- coral_dim$max_hull_volume[match(december_photo$coral_id,coral_dim$coral_id)] # Convex Hull
december_photo$pg_ellipse_aug <- coral_dim$ellipsoid_pg[match(december_photo$coral_id,coral_dim$coral_id)] # Photogrammetry ellipsoid volume
december_photo$length_aug <- coral_dim$length_pg[match(december_photo$coral_id,coral_dim$coral_id)] # PG length
december_photo$width_aug <- coral_dim$width_pg[match(december_photo$coral_id,coral_dim$coral_id)] # PG width
december_photo$height_aug <- coral_dim$height_pg[match(december_photo$coral_id,coral_dim$coral_id)] # PG height

december_photo <- december_photo %>% mutate(volume_pg = volume_pg*10^6, max_hull_volume = max_hull_volume*10^6,
                                            ellipsoid_pg = 4/3*pi*(height_pg/2)*(length_pg/2)*(width_pg/2)) # Unit conversion, december PG ellipsoid calculation
december_photo <- december_photo %>% mutate(vol_growth_pg = volume_pg - vol_pg_aug,
                                            prop_growth_pg = vol_growth_pg/vol_pg_aug*100,
                                            hull_growth = max_hull_volume - vol_hull_aug,
                                            prop_hull_growth = hull_growth/vol_hull_aug*100,
                                            ellipsoid_growth_pg = ellipsoid_pg - pg_ellipse_aug,
                                            prop_ellipsoid_growth_pg = ellipsoid_growth_pg/pg_ellipse_aug*100) # Proportional hull, volume, and ellipsoid growth calculations


december_manual <- december_manual %>% filter(coral_id %in% december_photo$coral_id) %>%
  mutate(vol_growth_man = volume_est_dec - volume_est) %>% mutate(prop_growth_man = vol_growth_man/volume_est*100) # Growth calculations

aug_and_dec <- left_join(december_photo, december_manual, by = "coral_id")

# Create stacked data frame for visualization of growth estimates by method
manual_dec_trim <- december_manual %>% select(coral_id, volume = volume_est_dec, volume_aug = volume_est, growth = vol_growth_man, prop_growth = prop_growth_man) %>% mutate(method = "Ellipsoid")

photo_dec_trim <- december_photo %>% select(coral_id, volume = volume_pg, volume_aug = vol_pg_aug, growth = vol_growth_pg, prop_growth = prop_growth_pg) %>% mutate(method = "Skeleton")

photo_hull_trim <- december_photo %>% select(coral_id, volume = max_hull_volume, volume_aug = vol_hull_aug, growth = hull_growth, prop_growth = prop_hull_growth) %>% mutate(method = "Convex Hull")

ellipsoid_pg_trim <- december_photo %>% select(coral_id, volume = ellipsoid_pg, volume_aug = pg_ellipse_aug, growth = ellipsoid_growth_pg, prop_growth = prop_ellipsoid_growth_pg) %>% mutate(method = "PG Ellipsoid")

december_stacked <- rbind(manual_dec_trim, ellipsoid_pg_trim, photo_hull_trim, photo_dec_trim) %>% drop_na()

# Create summary data frame for visualization of growth estimates by method
december_summary_stats <- december_stacked %>% group_by(method) %>% summarize(
  mean_growth = mean(growth, na.rm = TRUE), mean_prop_growth = mean(prop_growth, na.rm = TRUE),
  se_vol = sd(growth, na.rm = TRUE)/sqrt(n()), se_prop = sd(prop_growth, na.rm = TRUE)/sqrt(n()))

december_summary_stats$method <- factor(december_summary_stats$method, levels = c("Ellipsoid", "PG Ellipsoid", "Convex Hull", "Skeleton"))
december_stacked$method <- factor(december_stacked$method, levels = c("Ellipsoid", "PG Ellipsoid", "Convex Hull", "Skeleton"))


#Cleanup

remove(branch_w_summarized, cafi_summarized2, branch_width, coral_field, coral_field2, updated_cafi, coral_pg, updated_cafi2)
