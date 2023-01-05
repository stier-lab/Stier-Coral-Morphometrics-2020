rm(list = ls())

library(here)

source(here("code/source_data_cleanup_2022_11_25.R"))

library(rstatix)

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
 
# Growth Prop ellipse vs hull

cor_test_df <- december_photo %>% inner_join(december_manual, by = "coral_id") %>% select(contains("prop"), coral_id) %>% 
  filter(prop_hull_growth != "NA")

cor_prop_growth_hull_vs_pg <- cor.test(cor_test_df$prop_hull_growth, cor_test_df$prop_growth_pg)
cor_prop_growth_hull_vs_pg


# Growth Prop hull vs skeleton

cor_prop_growth_man_vs_hull <- cor.test(cor_test_df$prop_growth_man, cor_test_df$prop_hull_growth)
cor_prop_growth_man_vs_hull

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

growth_anova <- anova_test(data = december_stacked, dv = prop_growth, wid = coral_id, within = method, detailed = TRUE)

get_anova_table(growth_anova, correction = "GG")

aov(prop_growth ~ method, data = december_stacked) %>% tukey_hsd()
games_howell_test(december_stacked, prop_growth ~ method, detailed = TRUE)


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


## AIC Analysis

# Abun v Ellipsoid

lm_abun_man <- lm(num_cafi~volume_field, data = cafi_coral)
nls_abun_man <- nls(num_cafi ~ a*volume_field^beta, cafi_coral, start = list(a = 1.35, beta = 0.32))
summary(nls_abun_man)

aic_abun_man <- AIC(lm_abun_man,nls_abun_man)

# Abun v Skeleton

lm_abun_skel <- lm(num_cafi~volume_pg, data = cafi_coral)
nls_abun_skel <- nls(num_cafi ~ a*volume_pg^beta, cafi_coral, start = list(a = 1, beta = 1))
summary(nls_abun_skel)

aic_abun_skel <- AIC(lm_abun_skel,nls_abun_skel)

# Rich v Ellipsoid

lm_rich_man <- lm(cafi_richness~volume_field, data = cafi_coral)
nls_rich_man <- nls(cafi_richness ~ a*volume_field^beta, cafi_coral, start = list(a = 1, beta = 1))

aic_rich_man <- AIC(lm_rich_man,nls_rich_man)

# Rich v Skeleton

lm_rich_skel <- lm(cafi_richness~volume_pg, data = cafi_coral)
nls_rich_skel <- nls(cafi_richness ~ a*volume_pg^beta, cafi_coral, start = list(a = 1, beta = 1))

aic_rich_skel <- AIC(lm_rich_skel,nls_rich_skel)

# Rich v Hull

lm_rich_hull <- lm(cafi_richness~max_hull_volume, data = cafi_coral)
nls_rich_hull <- nls(cafi_richness ~ a*max_hull_volume^beta, cafi_coral, start = list(a = 1, beta = 1))

aic_rich_hull <- AIC(lm_rich_hull,nls_rich_hull)

# Abun v Hull

lm_abun_hull <- lm(num_cafi~max_hull_volume, data = cafi_coral)
nls_abun_hull <- nls(num_cafi ~ a*max_hull_volume^beta, cafi_coral, start = list(a = 1, beta = 1))

aic_abun_hull <- AIC(lm_abun_hull,nls_abun_hull)

# RMSE values

qpcR::RMSE(nls_abun_man)
qpcR::RMSE(nls_abun_skel)
qpcR::RMSE(nls_abun_hull)
qpcR::RMSE(nls_rich_man)
qpcR::RMSE(nls_rich_skel)
qpcR::RMSE(nls_rich_hull)

confint(nls_abun_man, level = 0.95)

nonnest2::vuongtest(nls_abun_hull, nls_abun_skel)

# Table test

reg_sum <- function(input2, nls_x) {
  temp <- data.frame(matrix(ncol = 11))
  temp_cols <- c("Variable", "a", "SE(a)", "t(a)", "p(a)", "b", "SE(b)", "t(b)", "p(b)", "AIC","RMSE")
  colnames(temp) <- temp_cols
  temp$Variable <- paste(input2)
  temp$a <- signif(as.numeric(coefficients(nls_x)[1]),3)
  temp$"SE(a)" <- signif(summary(nls_x)$coefficients[1,2],3)
  temp$"t(a)" <- signif(summary(nls_x)$coefficients[1,3],3)
  temp$"p(a)" <- signif(summary(nls_x)$coefficients[1,4],3)
  temp$b <- signif(as.numeric(coefficients(nls_x)[2]),3)
  temp$"SE(b)" <- signif(summary(nls_x)$coefficients[2,2],3)
  temp$"t(b)" <- signif(summary(nls_x)$coefficients[2,3],3)
  temp$"p(b)" <- signif(summary(nls_x)$coefficients[2,4],3)
  temp$AIC <- AIC(nls_x)
  temp$RMSE <- qpcR::RMSE(nls_x)
  temp
}

test <- reg_sum("Hull", nls_abun_hull)

nls_abun_man$coefficients
summary(nls_abun_man)
coefficients(nls_abun_hull)

signif(coefficients(nls_abun_man[1]),3)

signif(as.numeric(coefficients(nls_abun_man)[1]),3)

# November 2022 revisions

ellipsoid_comp_cor <- cor.test(coral_dim$ellipsoid_pg,coral_dim$volume_field)
ellipsoid_comp_cor

ellipsoid_hull_pg_cor <- cor.test(coral_dim$ellipsoid_pg,coral_dim$max_hull_volume)
ellipsoid_hull_pg_cor

test <- december_photo %>% select(coral_id, contains("prop")) %>% filter(coral_id != "FE-POC44")
sd(test$prop_ellipsoid_growth_pg)
# coral_dim %>%
#   ggplot(aes(ellipsoid_pg, volume_field)) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# coral_dim %>%
#   ggplot(aes(ellipsoid_pg, max_hull_volume)) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# coral_dim %>%
#   ggplot(aes(ellipsoid_pg, volume_field)) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# coral_dim %>%
#   ggplot(aes(max_hull_volume, volume_field)) +
#   geom_point() +
#   geom_smooth(method = "lm")


# Ellipse vs PG Ellipse
cor_ellipsoid_pg_vs_man <- cor.test(coral_dim$ellipsoid_pg, coral_dim$volume_field)
cor_ellipsoid_pg_vs_man

cor_prop_growth_ell_man_pg <- cor.test(cor_test_df$prop_ellipsoid_growth_pg, cor_test_df$prop_growth_man)
cor_prop_growth_ell_man_pg

cor_prop_growth_ellpg_hull <- cor.test(cor_test_df$prop_ellipsoid_growth_pg, cor_test_df$prop_hull_growth)
cor_prop_growth_ellpg_hull 

cor_prop_growth_ellpg_skel <- cor.test(cor_test_df$prop_ellipsoid_growth_pg, cor_test_df$prop_growth_pg)
cor_prop_growth_ellpg_skel

# Visual examination of PG ellipse correlations

aug_and_dec %>% 
  ggplot(aes(prop_ellipsoid_growth_pg, prop_growth_man)) +
  geom_point() +
  geom_smooth(method = "lm")

aug_and_dec %>% 
  ggplot(aes(prop_growth_pg, prop_growth_man)) +
  geom_point() +
  geom_smooth(method = "lm")

aug_and_dec %>% 
  ggplot(aes(prop_hull_growth, prop_growth_man)) +
  geom_point() +
  geom_smooth(method = "lm")

aug_and_dec %>% 
  ggplot(aes(prop_hull_growth, prop_growth_pg)) +
  geom_point() +
  geom_smooth(method = "lm")

aug_and_dec %>% 
  ggplot(aes(prop_hull_growth, prop_ellipsoid_growth_pg)) +
  geom_point() +
  geom_smooth(method = "lm")

aug_and_dec %>% 
  ggplot(aes(prop_growth_pg, prop_ellipsoid_growth_pg)) +
  geom_point() +
  geom_smooth(method = "lm")


# PG Ellipse growth v volume

ellpg_vol_v_growth_prop <- lm(prop_ellipsoid_growth_pg ~ pg_ellipse_aug, data = december_photo)
summary(ellpg_vol_v_growth_prop)

december_photo %>% 
  ggplot(aes(ellipsoid_pg, pg_ellipse_aug)) +
  geom_point() +
  geom_smooth(method = "lm")

december_photo %>% 
  ggplot(aes(length_pg, length_aug)) +
  geom_point() +
  geom_smooth(method = "lm")

december_photo %>% 
  ggplot(aes(width_pg, width_aug)) +
  geom_point() +
  geom_smooth(method = "lm")

december_photo %>% 
  ggplot(aes(height_pg, height_aug)) +
  geom_point() +
  geom_smooth(method = "lm")


december_photo <- december_photo %>% mutate(length_diff = length_pg - length_aug,
                                            width_diff = width_pg - width_aug,
                                            height_diff = height_pg - height_aug)

december_photo %>% 
  ggplot(aes(height_diff, coral_id)) +
  geom_point() +
  geom_smooth(method = "lm")

# CAFI regression - PG ellipsoid

nls_abun_pgell <- nls(num_cafi ~ a*ellipsoid_pg^beta, cafi_coral, start = list(a = 1, beta = 1))
nls_rich_pgell <- nls(cafi_richness ~ a*ellipsoid_pg^beta, cafi_coral, start = list(a = 1, beta = 1))

qpcR::RMSE(nls_abun_pgell)
qpcR::RMSE(nls_rich_pgell)

AIC(nls_abun_pgell, nls_rich_pgell)
