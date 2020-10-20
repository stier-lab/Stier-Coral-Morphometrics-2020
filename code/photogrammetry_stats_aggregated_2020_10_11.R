rm(list=ls())

# Load libraries
library(here)
library(tidyverse)
library(vegan)
library(reshape)
library(conflicted)
library(cowplot)

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

# CAFI data - Inverts only
updated_cafi2 <- updated_cafi %>% filter(str_detect(coral_id, "^FE")) %>%
  select(master_sort, coral_id, code, type, search_term, lowest_level, phylum, genus, species, general_notes)

# Calculate CAFI richness, abundance, and diversity (shannon weiner) for each coral (inverts)

cafi_summarized2 <- group_by(updated_cafi2, coral_id) %>%
  summarise(num_cafi = n(), cafi_richness = length(unique(code)), cafi_present = paste(sort(unique(code)), collapse = ";"))
cafi_summarized2$sw <- updated_cafi2 %>% 
  count(code, coral_id = coral_id) %>% 
  spread(code,n) %>% 
  mutate_all(list(~tidyr::replace_na(.,0))) %>% 
  select(-coral_id) %>% 
  diversity(index = "shannon")

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

coral_dim <- coral_dim %>% select(-"use", -position, -total_photos, -photos_not_aligned, -measurements, -locations) %>% 
  rename(notes_pg = notes.x, notes_field = notes.y)

# Undid filtering of removal corals only.. may need to reintroduce step at some later point

coral_dim$branch_conv <- ifelse(coral_dim$convexity>=0.5, "tight","wide") #classifies wide and tight branching coral based on convexity

# Merge coral and CAFI data
cafi_coral <- left_join(coral_dim, cafi_summarized2, by ="coral_id") %>% filter(cafi == "empty")

#Create trimmed dataset (not sure why this df is necessary.. maybe clear later on?)
dim_bind <- coral_dim %>% select(coral_id, ends_with("pg"), ends_with("field"))

#Cleanup

remove(branch_w_summarized, cafi_summarized2, branch_width, coral_field, coral_field2, updated_cafi, coral_pg, updated_cafi2)


### Regressions

# Ellipsoid - Skeleton vs skeletal volume

lm_diff_vs_vol <- lm(diff_est ~ volume_pg, data = coral_dim)
summary(lm_diff_vs_vol)

# Ellipsoid - Skeleton Prop vs skeletal volume

lm_prop_vs_vol <- lm(prop_est ~ volume_pg, data = coral_dim)
summary(lm_prop_vs_vol)

# Convex Hull vs Ellipse

cor_conv_vs_ellipse <- cor.test(coral_dim$max_hull_volume, coral_dim$volume_field,
                                method = "spearman")
cor_conv_vs_ellipse

# Ellipse vs PG volume

cor_vol_vs_ellipse <- cor.test(coral_dim$volume_pg, coral_dim$volume_field,
                                method = "spearman")
cor_vol_vs_ellipse

# Growth ellipse vs Growth PG

cor_growth_man_vs_pg <- cor.test(december_manual$vol_growth_man, december_photo$vol_growth_pg,
                                 method = "spearman")
cor_growth_man_vs_pg

# Growth Prop ellipse vs PG

cor_prop_growth_man_vs_pg <- cor.test(december_manual$prop_growth_man, december_photo$prop_growth_pg,
                                 method = "spearman")
cor_prop_growth_man_vs_pg
