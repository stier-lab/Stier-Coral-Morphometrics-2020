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
december_photo$vol_hull_aug <- coral_dim$max_hull_volume[match(december_photo$coral_id,coral_dim$coral_id)]

december_photo <- december_photo %>% mutate(volume_pg = volume_pg*10^6, max_hull_volume = max_hull_volume*10^6)
december_photo <- december_photo %>% mutate(vol_growth_pg = volume_pg - vol_pg_aug,
                                            prop_growth_pg = vol_growth_pg/vol_pg_aug*100,
                                            hull_growth = max_hull_volume - vol_hull_aug,
                                            prop_hull_growth = hull_growth/vol_hull_aug)

december_manual <- december_manual %>% filter(coral_id %in% december_photo$coral_id) %>%
  mutate(vol_growth_man = volume_est_dec - volume_est) %>% mutate(prop_growth_man = vol_growth_man/volume_est*100)
aug_and_dec <- left_join(december_photo, december_manual, by = "coral_id")


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

t_vol_growth <- t.test(december_manual$vol_growth_man, december_photo$vol_growth_pg, paired = TRUE)
t_vol_growth

sd(december_manual$vol_growth_man)
sd(december_photo$vol_growth_pg)
sd(december_photo$hull_growth, na.rm = TRUE)

t_prop_growth <- t.test(december_manual$prop_growth_man, december_photo$prop_growth_pg, paired = TRUE)
t_prop_growth

# Regressions of volume vs growth

man_vol_v_growth <- lm(vol_growth_man ~ volume_est, data = december_manual)
summary(man_vol_v_growth)


pg_vol_v_growth <- lm(vol_growth_pg ~ vol_pg_aug, data = december_photo)
summary(pg_vol_v_growth)

hull_vol_v_growth <- lm(hull_growth ~ vol_hull_aug, data = december_photo)
summary(hull_vol_v_growth)
# Average Model Error

# coral_pg <- coral_pg %>% mutate(scalebar_num = if_else(is.na(coral_pg$scalebar_error_z), 2, 3))
# coral_dim <- coral_dim %>% mutate(average_error = scalebar_error_x + scalebar_error_y +
#                                     if_else(is.na(scalebar_error_z), 0, scalebar_error_z)/
#                                     scalebar_num)
# 

all_error <- c(coral_pg$scalebar_error_x, coral_pg$scalebar_error_y, coral_pg$scalebar_error_z,
               december_photo$scalebar_error_x, december_photo$scalebar_error_y, december_photo$scalebar_error_z)

mean(abs(all_error), na.rm = TRUE)
# Mean = 0.0002215286m (.22mm)

# Relationships between manual and photogrammetry from December models

# convex <- ggplot(aug_and_dec, aes(x=max_hull_volume, y = volume_est_dec))+ 
#   geom_point()+
#   geom_smooth(method=lm, color = "coral")+
#   scale_y_continuous(breaks=seq(0,max(coral_dim$volume_field), by = 5000))+
#   scale_x_continuous(breaks=seq(0,max(coral_dim$max_hull_volume, na.rm = TRUE)+1500, by = 5000))+
#   theme_classic()+
#   labs(y= expression(paste("Ellipsoid Volume (cm"^{3},")")), 
#        x= expression(paste("Convex Hull Volume (cm"^{3},")"))) +
#   geom_abline(intercept = 0, slope = 1, size = 1, lty = "dashed", color='black') +
#   theme(axis.text.x = element_text(angle = 305, vjust=0, hjust = .3))#x intercept changed to be forced through 0
# 
# height <- ggplot(aug_and_dec, aes(x=height_pg, y=height_dec))+ 
#   geom_point()+
#   geom_smooth(method=lm, color = "coral")+
#   theme_classic()+
#   labs(y= "Manual Height (cm)", x= "Photogrammetry Height (cm)")+
#   geom_abline(intercept = 0, slope = 1, size = 1, lty = "dashed", color='black') #x intercept changed to be forced through 0
# 
# #Length
# length <- ggplot(aug_and_dec, aes(x=length_pg, y=length_dec))+ 
#   geom_point()+
#   geom_smooth(method=lm, color = "coral")+
#   theme_classic()+
#   labs(y= "Manual Length (cm)", x= "Photogrammetry Length (cm)")+
#   geom_abline(intercept = 0, slope = 1, size = 1, lty = "dashed", color='black') #x intercept changed to be forced through 0
# 
# #Width
# width <- ggplot(aug_and_dec, aes(x=width_pg, y=width_dec))+ 
#   geom_point()+
#   geom_smooth(method=lm, color = "coral")+
#   theme_classic()+
#   labs(y= "Manual With (cm)", x= "Photogrammetry Width (cm)")+
#   geom_abline(intercept = 0, slope = 1, size = 1, lty = "dashed", color='black') #x intercept changed to be forced through 0
# 
# fig1_dec <- cowplot::plot_grid(convex, length, width, height, ncol = 2, nrow = 2, 
#                            align = "vh", labels = "AUTO")

# Residual analysis

# cafi_coral <- cafi_coral %>% filter(coral_id != "FE-POC01")

lm_invert_field <- lm(num_cafi~volume_field, data=cafi_coral)
lm_invert_pg <- lm(num_cafi~volume_pg, data=cafi_coral)
lm_richness_field <- lm(cafi_richness~volume_field, data=cafi_coral)
lm_richness_pg <- lm(cafi_richness~volume_pg, data=cafi_coral)
lm_invert_hull <- lm(num_cafi~max_hull_volume, data=cafi_coral)
lm_richness_hull <- lm(cafi_richness~max_hull_volume, data=cafi_coral)

cafi_coral$ellipse_abun_resid <- residuals(lm_invert_field)
cafi_coral$hull_abun_resid <- residuals(lm_invert_hull)
cafi_coral$abs_ellipse_resid <- abs(residuals(lm_invert_field))
cafi_coral$abs_hull_resid <- abs(residuals(lm_richness_field))
cafi_coral$skel_abun_resid <- residuals(lm_invert_pg)

cafi_coral$resid_diff <- cafi_coral$ellipse_abun_resid - cafi_coral$skel_abun_resid

resid_conv_stats <- lm(resid_diff~convexity, data = cafi_coral)
summary(resid_conv_stats)
  
formula1 <- y~x

resid_diff_plot <- ggplot(cafi_coral, aes(x=convexity, y = resid_diff))+ 
  geom_point()+
  geom_smooth(method=lm, formula = formula1)+
  theme_classic()+
  labs(y="Difference in Residuals (Ellipse - Skeleton ~ Abundance)", x= "Convexity") +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "bottom",
               formula = formula1, parse = TRUE, size = 5)

ellipse_conv <- ggplot(cafi_coral, aes(x=convexity, y = ellipse_abun_resid))+ 
  geom_point()+
  geom_smooth(method=lm, formula = formula1)+
  theme_classic()+
  labs(y="Residual (ellipse vs abundance)", x= "Convexity") +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "bottom",
               formula = formula1, parse = TRUE, size = 5)

# hull_conv <- ggplot(cafi_coral, aes(x=convexity, y=hull_abun_resid))+
#   geom_point()+
#   geom_smooth(method=lm, formula = formula1)+
#   theme_classic()+
#   labs(y="Residual (hull vs abundance)", x= "Convexity") +
#   stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")),
#                label.x.npc = "right", label.y.npc = "bottom",
#                formula = formula1, parse = TRUE, size = 5)

ellipse_conv_abs <- ggplot(cafi_coral, aes(x=convexity, y=abs_ellipse_resid))+
  geom_point()+
  geom_smooth(method=lm, formula = formula1)+
  theme_classic()+
  labs(y="Absolute Residual (ellipse vs abundance)", x= "Convexity") +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")),
               label.x.npc = "right", label.y.npc = "bottom",
               formula = formula1, parse = TRUE, size = 5)

# hull_conv_abs <- ggplot(cafi_coral, aes(x=convexity, y=abs_hull_resid))+
#   geom_point()+
#   geom_smooth(method=lm, formula = formula1)+
#   theme_classic()+
#   labs(y="Absolute Residual (hull vs abundance)", x= "Convexity") +
#   stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")),
#                label.x.npc = "right", label.y.npc = "bottom",
#                formula = formula1, parse = TRUE, size = 5)

# plot_grid(ellipse_conv, hull_conv, ellipse_conv_abs, hull_conv_abs,
#           ncol = 2, nrow = 2)

ellipse_prop <- ggplot(cafi_coral, aes(x=prop_est, y=ellipse_abun_resid))+
  geom_point()+
  geom_smooth(method=lm, formula = formula1)+
  theme_classic()+
  labs(y="Residual (ellipse vs abundance)", x= "Proportion of Ellipse Overestimation") +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")),
               label.x.npc = "right", label.y.npc = "bottom",
               formula = formula1, parse = TRUE, size = 5)

# hull_prop <- ggplot(cafi_coral, aes(x=prop_est, y=hull_abun_resid))+
#   geom_point()+
#   geom_smooth(method=lm, formula = formula1)+
#   theme_classic()+
#   labs(y="Residual (hull vs abundance)", x= "Proportion of Ellipse Overestimation") +
#   stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")),
#                label.x.npc = "right", label.y.npc = "bottom",
#                formula = formula1, parse = TRUE, size = 5)

ellipse_prop_abs <- ggplot(cafi_coral, aes(x=prop_est, y=abs_ellipse_resid))+
  geom_point()+
  geom_smooth(method=lm, formula = formula1)+
  theme_classic()+
  labs(y="Absolute Residual (ellipse vs abundance)", x= "Proportion of Ellipse Overestimation") +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")),
               label.x.npc = "right", label.y.npc = "bottom",
               formula = formula1, parse = TRUE, size = 5)

# hull_prop_abs <- ggplot(cafi_coral, aes(x=prop_est, y=abs_hull_resid))+
#   geom_point()+
#   geom_smooth(method=lm, formula = formula1)+
#   theme_classic()+
#   labs(y="Absolute Residual (hull vs abundance)", x= "Proportion of Ellipse Overestimation") +
#   stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")),
#                label.x.npc = "right", label.y.npc = "bottom",
#                formula = formula1, parse = TRUE, size = 5)

# plot_grid(ellipse_prop, hull_prop, ellipse_prop_abs, hull_prop_abs,
#           ncol = 2, nrow = 2)

ellipse_diff <- ggplot(cafi_coral, aes(x=diff_est, y=ellipse_abun_resid))+ 
  geom_point()+
  geom_smooth(method=lm, formula = formula1)+
  theme_classic()+
  labs(y="Residual (ellipse vs abundance)", x= "Absolute Ellipse Overestimation") +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "bottom",
               formula = formula1, parse = TRUE, size = 5)

# hull_diff <- ggplot(cafi_coral, aes(x=diff_est, y=hull_abun_resid))+
#   geom_point()+
#   geom_smooth(method=lm, formula = formula1)+
#   theme_classic()+
#   labs(y="Residual (hull vs abundance)", x= "Absolute Ellipse Overestimation") +
#   stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")),
#                label.x.npc = "right", label.y.npc = "bottom",
#                formula = formula1, parse = TRUE, size = 5)

ellipse_diff_abs <- ggplot(cafi_coral, aes(x=diff_est, y=abs_ellipse_resid))+
  geom_point()+
  geom_smooth(method=lm, formula = formula1)+
  theme_classic()+
  labs(y="Absolute Residual (ellipse vs abundance)", x= "Absolute Ellipse Overestimation") +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")),
               label.x.npc = "right", label.y.npc = "bottom",
               formula = formula1, parse = TRUE, size = 5)

# hull_diff_abs <- ggplot(cafi_coral, aes(x=diff_est, y=abs_hull_resid))+
#   geom_point()+
#   geom_smooth(method=lm, formula = formula1)+
#   theme_classic()+
#   labs(y="Absolute Residual (hull vs abundance)", x= "Absolute Ellipse Overestimation") +
#   stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")),
#                label.x.npc = "right", label.y.npc = "bottom",
#                formula = formula1, parse = TRUE, size = 5)

# plot_grid(ellipse_diff, hull_diff, ellipse_diff_abs, hull_diff_abs,
#           ncol = 2, nrow = 2)

residual_plot <- plot_grid(ellipse_conv, ellipse_prop, ellipse_diff,
                           ellipse_conv_abs, ellipse_prop_abs, ellipse_diff_abs,
                           nrow = 2, ncol = 3)

# Error Propagation


