rm(list=ls())

# Load libraries
library(here)
library(tidyverse)
library(vegan)
library(reshape)
library(conflicted)
library(cowplot)
library(ez)

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

coral_dim <- coral_dim %>% select(-"use", -position, -photos_not_aligned, -measurements, -locations) %>% 
  rename(notes_pg = notes.x, notes_field = notes.y)

# Undid filtering of removal corals only.. may need to reintroduce step at some later point

coral_dim$branch_conv <- ifelse(coral_dim$convexity>=0.5, "tight","wide") #classifies wide and tight branching coral based on convexity

# Merge coral and CAFI data
cafi_coral <- left_join(coral_dim, cafi_summarized2, by ="coral_id") %>% filter(cafi == "empty")

# Growth Data
december_manual <- read.csv(here("morphometric_data", "field_experiment_colony_measurements_moorea_december2019_v2.csv"))
december_photo <- read.csv(here("morphometric_data","maatea_experiment_photo_measurements_december_JC_2020_7_8.csv"))

december_manual$coral_id <- as.character(december_manual$coral_id)
december_photo$coral_id <- as.character(december_photo$coral_id)

december_photo <- december_photo %>% filter(coral_id %in% coral_dim$coral_id, volume_pg != "NA")

december_photo$vol_pg_aug <- coral_dim$volume_pg[match(december_photo$coral_id,coral_dim$coral_id)]
december_photo$sa_aug <-  coral_dim$SA[match(december_photo$coral_id,coral_dim$coral_id)]
december_photo$vol_hull_aug <- coral_dim$max_hull_volume[match(december_photo$coral_id,coral_dim$coral_id)]

december_photo <- december_photo %>% mutate(volume_pg = volume_pg*10^6, max_hull_volume = max_hull_volume*10^6)
december_photo <- december_photo %>% mutate(vol_growth_pg = volume_pg - vol_pg_aug,
                                            prop_growth_pg = vol_growth_pg/vol_pg_aug*100,
                                            hull_growth = max_hull_volume - vol_hull_aug,
                                            prop_hull_growth = hull_growth/vol_hull_aug*100)


december_manual <- december_manual %>% filter(coral_id %in% december_photo$coral_id) %>%
  mutate(vol_growth_man = volume_est_dec - volume_est) %>% mutate(prop_growth_man = vol_growth_man/volume_est*100)

aug_and_dec <- left_join(december_photo, december_manual, by = "coral_id")

manual_dec_trim <- december_manual %>% select(coral_id, volume = volume_est_dec, volume_aug = volume_est, growth = vol_growth_man, prop_growth = prop_growth_man) %>% mutate(method = "Ellipsoid")

photo_dec_trim <- december_photo %>% select(coral_id, volume = volume_pg, volume_aug = vol_pg_aug, growth = vol_growth_pg, prop_growth = prop_growth_pg) %>% mutate(method = "Skeleton")

photo_hull_trim <- december_photo %>% select(coral_id, volume = max_hull_volume, volume_aug = vol_hull_aug, growth = hull_growth, prop_growth = prop_hull_growth) %>% mutate(method = "Convex Hull")

december_stacked <- rbind(manual_dec_trim,photo_dec_trim, photo_hull_trim) %>% drop_na()
december_summary_stats <- december_stacked %>% group_by(method) %>% summarize(
  mean_growth = mean(growth, na.rm = TRUE), mean_prop_growth = mean(prop_growth, na.rm = TRUE),
  se_vol = sd(growth, na.rm = TRUE)/sqrt(n()), se_prop = sd(prop_growth, na.rm = TRUE)/sqrt(n()))

december_summary_stats$method <- factor(december_summary_stats$method, levels = c("Ellipsoid", "Convex Hull", "Skeleton"))
december_stacked$method <- factor(december_stacked$method, levels = c("Ellipsoid", "Convex Hull", "Skeleton"))


#Cleanup

remove(branch_w_summarized, cafi_summarized2, branch_width, coral_field, coral_field2, updated_cafi, coral_pg, updated_cafi2)


# Summary Stats

num_photos <- summary(c(coral_dim$total_photos, december_photo$total_photos))
num_photos
### Regressions

# Ellipsoid - Skeleton vs skeletal volume

lm_diff_vs_vol <- lm(diff_est ~ volume_pg, data = coral_dim)
summary(lm_diff_vs_vol)

# Ellipsoid - Skeleton Prop vs skeletal volume

lm_prop_vs_vol <- lm(prop_est ~ volume_pg, data = coral_dim)
summary(lm_prop_vs_vol)

# Convex Hull vs Ellipse

cor_conv_vs_ellipse <- cor.test(coral_dim$max_hull_volume, coral_dim$volume_field)
cor_conv_vs_ellipse

# Ellipse vs PG volume

cor_vol_vs_ellipse <- cor.test(coral_dim$volume_pg, coral_dim$volume_field)
cor_vol_vs_ellipse

# Growth ellipse vs Growth PG

cor_growth_man_vs_pg <- cor.test(december_manual$vol_growth_man, december_photo$vol_growth_pg)
cor_growth_man_vs_pg

# Growth Prop ellipse vs PG

cor_prop_growth_man_vs_pg <- cor.test(december_manual$prop_growth_man, december_photo$prop_growth_pg)
cor_prop_growth_man_vs_pg
 
# T-test: Volume and proportional growth manual vs. photo

t_vol_growth <- t.test(december_manual$vol_growth_man, december_photo$vol_growth_pg)
t_vol_growth

sd(december_manual$vol_growth_man)
sd(december_photo$vol_growth_pg)
sd(december_photo$hull_growth, na.rm = TRUE)

t_prop_growth <- t.test(december_manual$prop_growth_man, december_photo$prop_growth_pg, paired = TRUE)
t_prop_growth

t_prop_growth <- t.test(december_manual$prop_growth_man, december_photo$prop_growth_pg, paired = TRUE)
t_prop_growth

t_prop_growth_hull <- t.test(december_manual$prop_growth_man, december_photo$prop_hull_growth, paired = TRUE)
t_prop_growth_hull

sd(december_manual$prop_growth_man)
sd(december_photo$prop_growth_pg)
sd(december_photo$prop_hull_growth, na.rm = TRUE)

# ANOVA for growth

december_stacked <- december_stacked %>% filter(coral_id != "FE-POC16", coral_id != "FE-POC56")

growth_anova <- anova_test(data = december_stacked, dv = prop_growth, wid = coral_id, within = method)

get_anova_table(growth_anova)

# Regressions of volume vs growth

man_vol_v_growth <- lm(vol_growth_man ~ volume_est, data = december_manual)
summary(man_vol_v_growth)


pg_vol_v_growth <- lm(vol_growth_pg ~ vol_pg_aug, data = december_photo)
summary(pg_vol_v_growth)

hull_vol_v_growth <- lm(hull_growth ~ vol_hull_aug, data = december_photo)
summary(hull_vol_v_growth)

man_vol_v_growth_prop <- lm(prop_growth_man ~ volume_est, data = december_manual)
summary(man_vol_v_growth_prop)


pg_vol_v_growth_prop <- lm(prop_growth_pg ~ vol_pg_aug, data = december_photo)
summary(pg_vol_v_growth_prop)

hull_vol_v_growth_prop <- lm(prop_hull_growth ~ vol_hull_aug, data = december_photo)
summary(hull_vol_v_growth_prop)

# Residual analysis stats
resid_conv_stats <- lm(resid_diff~convexity, data = cafi_coral)
summary(resid_conv_stats)
  
# Length, width, height comparisons

coral_dim <- coral_dim %>% mutate(height_diff = height_pg-height_field, length_diff = length_pg - length_field, width_diff = width_pg - width_field)

summary(coral_dim$width_diff)
summary(coral_dim$length_diff)
summary(coral_dim$height_diff)

sd(coral_dim$width_diff)
sd(coral_dim$length_diff)
sd(coral_dim$height_diff)

