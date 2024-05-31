# Code for Picardi et al., "Fitness consequences of anthropogenic subsidies
# for a partially migratory wading bird"

# Load packages ####

library(tidyverse)
library(survival)

# Load data ####

dat <- readRDS("input/rsf-data.rds")

# Create groups ####

# South Florida, Jacksonville, Southeast
dat <- dat %>%
  mutate(se = !fl & !jax,
         group = case_when(
           jax ~ "Jacksonville",
           soflo ~ "South Florida",
           TRUE ~ "Southeast"
         ))

dat %>%
  select(att_id, group, choice) %>%
  distinct() %>%
  group_by(group, choice) %>%
  tally()

# RSF for binary developed/undeveloped ####

## South Florida ####

mod_bin_soflo <- clogit(formula = used ~
                    urb_agg_2016 +
                    urb_agg_2016 : choice +
                    strata(att_id),
                  data = dat,
                  subset = soflo,
                  weights = weight,
                  method = "approximate")

## Jacksonville ####

mod_bin_jax <- clogit(formula = used ~
                    urb_agg_2016 +
                    strata(att_id),
                  data = dat,
                  subset = jax,
                  weights = weight,
                  method = "approximate")

## Southeast ####

mod_bin_se <- clogit(formula = used ~
                        urb_agg_2016 +
                        strata(att_id),
                      data = dat,
                      subset = se,
                      weights = weight,
                      method = "approximate")

# RSF for distance to development ####

## 2001 ####

### South Florida ####

mod_dist_2001_soflo <- clogit(formula = used ~
                                dist_to_urb_2001_sc +
                                I(dist_to_urb_2001_sc^2) +
                                dist_to_urb_2001_sc : choice +
                                I(dist_to_urb_2001_sc^2) : choice +
                                strata(att_id),
                              data = dat,
                              subset = soflo,
                              weights = weight,
                              method = "approximate")

### Jacksonville ####

mod_dist_2001_jax <- clogit(formula = used ~
                              dist_to_urb_2001_sc +
                              I(dist_to_urb_2001_sc^2) +
                              strata(att_id),
                            data = dat,
                            subset = jax,
                            weights = weight,
                            method = "approximate")

### Southeast ####

mod_dist_2001_se <- clogit(formula = used ~
                             dist_to_urb_2001_sc +
                             I(dist_to_urb_2001_sc^2) +
                             strata(att_id),
                           data = dat,
                           subset = se,
                           weights = weight,
                           method = "approximate")

## 2016 ####

### South Florida ####

mod_dist_2016_soflo <- clogit(formula = used ~
                          dist_to_urb_2016_sc +
                          I(dist_to_urb_2016_sc^2) +
                          dist_to_urb_2016_sc : choice +
                          I(dist_to_urb_2016_sc^2) : choice +
                          strata(att_id),
                          data = dat,
                          subset = soflo,
                          weights = weight,
                        method = "approximate")

### Jacksonville ####

mod_dist_2016_jax <- clogit(formula = used ~
                                dist_to_urb_2016_sc +
                                I(dist_to_urb_2016_sc^2) +
                                strata(att_id),
                              data = dat,
                              subset = jax,
                              weights = weight,
                              method = "approximate")

### Southeast ####

mod_dist_2016_se <- clogit(formula = used ~
                                dist_to_urb_2016_sc +
                                I(dist_to_urb_2016_sc^2) +
                                strata(att_id),
                              data = dat,
                              subset = se,
                              weights = weight,
                              method = "approximate")

## 2010 ####

### South Florida ####

mod_dist_2010_soflo <- clogit(formula = used ~
                                dist_to_urb_2010_sc +
                                I(dist_to_urb_2010_sc^2) +
                                dist_to_urb_2010_sc : choice +
                                I(dist_to_urb_2010_sc^2) : choice +
                                strata(att_id),
                              data = dat,
                              subset = soflo,
                              weights = weight,
                              method = "approximate")

### Jacksonville ####

mod_dist_2010_jax <- clogit(formula = used ~
                              dist_to_urb_2010_sc +
                              I(dist_to_urb_2010_sc^2) +
                              strata(att_id),
                            data = dat,
                            subset = jax,
                            weights = weight,
                            method = "approximate")

### Southeast ####

mod_dist_2010_se <- clogit(formula = used ~
                             dist_to_urb_2010_sc +
                             I(dist_to_urb_2010_sc^2) +
                             strata(att_id),
                           data = dat,
                           subset = se,
                           weights = weight,
                           method = "approximate")

# Model predictions for binary model ####

## South Florida ####

predd_bin <- expand.grid(urb_agg_2016 = c(0, 1),
                         choice = c("migrant", "resident")) %>%
  mutate(att_id = case_when(
    choice == "migrant" ~ dat[dat$choice == "migrant" & dat$soflo, ]$att_id[1],
    choice == "resident" ~ dat[dat$choice == "resident" & dat$soflo, ]$att_id[1]
  )) %>%
  arrange(urb_agg_2016, choice)

refd_bin <- expand.grid(urb_agg_2016 = 0,
                        choice = c("migrant", "resident")) %>%
  mutate(att_id = case_when(
    choice == "migrant" ~ dat[dat$choice == "migrant" & dat$soflo, ]$att_id[1],
    choice == "resident" ~ dat[dat$choice == "resident" & dat$soflo, ]$att_id[1]
  )) %>%
  arrange(urb_agg_2016, choice)

refd_bin <- rbind(refd_bin, refd_bin)

pred_X <- model.matrix(mod_bin_soflo, data = predd_bin)
ref_X <- model.matrix(mod_bin_soflo, data = refd_bin)
B <- coef(mod_bin_soflo)

# Difference
diff_X <- pred_X - ref_X

# Linear predictor for difference
diff_bin <- diff_X %*% B

# Get variance of difference
var_pred <- diag(diff_X %*% vcov(mod_bin_soflo) %*% t(diff_X))

# Standard error for difference
SE <- sqrt(var_pred)

predd_bin_soflo <- predd_bin %>%
  mutate(log_rss = diff_bin,
         se = SE,
         lwr = log_rss - 1.96*se,
         upr = log_rss + 1.96*se)

## Jacksonville ####

predd_bin <- expand.grid(urb_agg_2016 = c(0, 1),
                         choice = "resident") %>%
  mutate(att_id = dat[dat$choice == "resident" & dat$jax, ]$att_id[1]) %>%
  arrange(urb_agg_2016, choice)

refd_bin <- expand.grid(urb_agg_2016 = 0,
                        choice = "resident") %>%
  mutate(att_id = dat[dat$choice == "resident" & dat$jax, ]$att_id[1]) %>%
  arrange(urb_agg_2016, choice)

refd_bin <- rbind(refd_bin, refd_bin)

pred_X <- model.matrix(mod_bin_jax, data = predd_bin)
ref_X <- model.matrix(mod_bin_jax, data = refd_bin)
B <- coef(mod_bin_jax)

# Difference
diff_X <- pred_X - ref_X

# Linear predictor for difference
diff_bin <- diff_X %*% B

# Get variance of difference
var_pred <- diag(diff_X %*% vcov(mod_bin_jax) %*% t(diff_X))

# Standard error for difference
SE <- sqrt(var_pred)

predd_bin_jax <- predd_bin %>%
  mutate(log_rss = diff_bin,
         se = SE,
         lwr = log_rss - 1.96*se,
         upr = log_rss + 1.96*se)

## Southeast ####

predd_bin <- expand.grid(urb_agg_2016 = c(0, 1),
                         choice = "migrant") %>%
  mutate(att_id = dat[dat$choice == "migrant" & dat$se, ]$att_id[1]) %>%
  arrange(urb_agg_2016, choice)

refd_bin <- expand.grid(urb_agg_2016 = 0,
                        choice = "migrant") %>%
  mutate(att_id = dat[dat$choice == "migrant" & dat$se, ]$att_id[1]) %>%
  arrange(urb_agg_2016, choice)

refd_bin <- rbind(refd_bin, refd_bin)

pred_X <- model.matrix(mod_bin_se, data = predd_bin)
ref_X <- model.matrix(mod_bin_se, data = refd_bin)
B <- coef(mod_bin_se)

# Difference
diff_X <- pred_X - ref_X

# Linear predictor for difference
diff_bin <- diff_X %*% B

# Get variance of difference
var_pred <- diag(diff_X %*% vcov(mod_bin_se) %*% t(diff_X))

# Standard error for difference
SE <- sqrt(var_pred)

predd_bin_se <- predd_bin %>%
  mutate(log_rss = diff_bin,
         se = SE,
         lwr = log_rss - 1.96*se,
         upr = log_rss + 1.96*se)

# Model predictions for distance models ####

## 2001 ####

dat %>%
  filter(used == 1 & soflo) %>%
  summarize(max(dist_to_urb_2001, na.rm = TRUE))

mean_dist_2001 <- unique(dat$mean_dist_to_urb_2001)
sd_dist_2001 <- unique(dat$sd_dist_to_urb_2001)

pred_df_dist_2001 <- expand.grid(choice = c("migrant", "resident"),
                                 dist_to_urban_orig = seq(0, 38000, length.out = 100),
                                 att_id = unique(dat[dat$soflo, ]$att_id)[1]) %>%
  mutate(dist_to_urb_2001_sc = (dist_to_urban_orig - mean_dist_2001)/sd_dist_2001,
         dist_to_urban_km = dist_to_urban_orig/1000)

dat %>%
  filter(used == 1) %>%
  summarize(mean(dist_to_urb_2001, na.rm = TRUE))

refd_df_dist_2001 <- expand.grid(
  choice = c("migrant", "resident"),
  dist_to_urban_orig = 38000/2,
  att_id = unique(dat[dat$soflo, ]$att_id)[1]) %>%
  mutate(dist_to_urb_2001_sc = (dist_to_urban_orig - mean_dist_2001)/sd_dist_2001,
         dist_to_urban_km = dist_to_urban_orig/1000)

refd_df_dist_2001_dup <- refd_df_dist_2001

for (i in 1:(nrow(pred_df_dist_2001)/2 - 1)) {
  refd_df_dist_2001_dup <- rbind(refd_df_dist_2001_dup, refd_df_dist_2001)
}

nrow(refd_df_dist_2001_dup) == nrow(pred_df_dist_2001)

pred_X <- model.matrix(mod_dist_2001, data = pred_df_dist_2001)
ref_X <- model.matrix(mod_dist_2001, data = refd_df_dist_2001_dup)
B <- coef(mod_dist_2001)

# Difference
diff_X <- pred_X - ref_X

# Linear predictor for difference
diff_dist_2001 <- diff_X %*% B

# Get variance of difference
var_pred <- diag(diff_X %*% vcov(mod_dist_2001) %*% t(diff_X))

# Standard error for difference
SE <- sqrt(var_pred)

predd_dist_2001 <- pred_df_dist_2001 %>%
  mutate(log_rss = diff_dist_2001,
         se = SE,
         lwr = log_rss - 1.96*se,
         upr = log_rss + 1.96*se)

## 2016 ####

max_distances <- dat %>%
  filter(used == 1) %>%
  group_by(group) %>%
  summarize(xlim = max(dist_to_urb_2016, na.rm = TRUE))

### South Florida ####

predd_dist_2016 <- expand.grid(dist_to_urb_2016 = seq(0,
                                                max_distances[max_distances$group == "South Florida", ]$xlim,
                                                length.out = 100),
                         choice = c("migrant", "resident")) %>%
  mutate(att_id = case_when(
    choice == "migrant" ~ dat[dat$choice == "migrant" & dat$soflo, ]$att_id[1],
    choice == "resident" ~ dat[dat$choice == "resident" & dat$soflo, ]$att_id[1]
  )) %>%
  mutate(dist_to_urb_2016_sc = (dist_to_urb_2016 - unique(dat$mean_dist_to_urb_2016))/unique(dat$sd_dist_to_urb_2016)) %>%
  arrange(choice, dist_to_urb_2016)

refd_dist_2016 <- expand.grid(dist_to_urb_2016 = 0,
                         choice = c("migrant", "resident")) %>%
  mutate(att_id = case_when(
    choice == "migrant" ~ dat[dat$choice == "migrant" & dat$soflo, ]$att_id[1],
    choice == "resident" ~ dat[dat$choice == "resident" & dat$soflo, ]$att_id[1]
  )) %>%
  mutate(dist_to_urb_2016_sc = (dist_to_urb_2016 - unique(dat$mean_dist_to_urb_2016))/unique(dat$sd_dist_to_urb_2016))

refd_dist_2016_dup <- refd_dist_2016

for (i in 1:(nrow(predd_dist_2016)/2 - 1)) {
  refd_dist_2016_dup <- rbind(refd_dist_2016_dup, refd_dist_2016)
}

nrow(refd_dist_2016_dup) == nrow(predd_dist_2016)

refd_dist_2016_dup <- refd_dist_2016_dup %>%
  arrange(choice, dist_to_urb_2016)

pred_X <- model.matrix(mod_dist_2016_soflo, data = predd_dist_2016)
ref_X <- model.matrix(mod_dist_2016_soflo, data = refd_dist_2016_dup)
B <- coef(mod_dist_2016_soflo)

# Difference
diff_X <- pred_X - ref_X

# Linear predictor for difference
diff_dist_2016 <- unname((diff_X %*% B)[,1])

# Get variance of difference
var_pred <- diag(diff_X %*% vcov(mod_dist_2016_soflo) %*% t(diff_X))

# Standard error for difference
SE <- unname(sqrt(var_pred))

predd_dist_2016_soflo <- predd_dist_2016 %>%
  mutate(log_rss = diff_dist_2016,
         se = SE,
         lwr = log_rss - 1.96*se,
         upr = log_rss + 1.96*se)

### Jacksonville ####

predd_dist_2016 <- expand.grid(dist_to_urb_2016 = seq(0,
                                                      max_distances[max_distances$group == "Jacksonville", ]$xlim,
                                                      length.out = 100),
                               choice = "resident") %>%
  mutate(att_id = dat[dat$choice == "resident" & dat$jax, ]$att_id[1]) %>%
  mutate(dist_to_urb_2016_sc = (dist_to_urb_2016 - unique(dat$mean_dist_to_urb_2016))/unique(dat$sd_dist_to_urb_2016)) %>%
  arrange(choice, dist_to_urb_2016)

refd_dist_2016 <- expand.grid(dist_to_urb_2016 = 0,
                              choice = "resident") %>%
  mutate(att_id = dat[dat$choice == "resident" & dat$jax, ]$att_id[1]) %>%
  mutate(dist_to_urb_2016_sc = (dist_to_urb_2016 - unique(dat$mean_dist_to_urb_2016))/unique(dat$sd_dist_to_urb_2016))

refd_dist_2016_dup <- refd_dist_2016

for (i in 1:(nrow(predd_dist_2016) - 1)) {
  refd_dist_2016_dup <- rbind(refd_dist_2016_dup, refd_dist_2016)
}

nrow(refd_dist_2016_dup) == nrow(predd_dist_2016)

refd_dist_2016_dup <- refd_dist_2016_dup %>%
  arrange(choice, dist_to_urb_2016)

pred_X <- model.matrix(mod_dist_2016_jax, data = predd_dist_2016)
ref_X <- model.matrix(mod_dist_2016_jax, data = refd_dist_2016_dup)
B <- coef(mod_dist_2016_jax)

# Difference
diff_X <- pred_X - ref_X

# Linear predictor for difference
diff_dist_2016 <- diff_X %*% B

# Get variance of difference
var_pred <- diag(diff_X %*% vcov(mod_dist_2016_jax) %*% t(diff_X))

# Standard error for difference
SE <- sqrt(var_pred)

predd_dist_2016_jax <- predd_dist_2016 %>%
  mutate(log_rss = diff_dist_2016,
         se = SE,
         lwr = log_rss - 1.96*se,
         upr = log_rss + 1.96*se)

### Southeast ####

predd_dist_2016 <- expand.grid(dist_to_urb_2016 = seq(0,
                                                      max_distances[max_distances$group == "Southeast", ]$xlim,
                                                      length.out = 100),
                               choice = "migrant") %>%
  mutate(att_id = dat[dat$choice == "migrant" & dat$se, ]$att_id[1]) %>%
  mutate(dist_to_urb_2016_sc = (dist_to_urb_2016 - unique(dat$mean_dist_to_urb_2016))/unique(dat$sd_dist_to_urb_2016)) %>%
  arrange(choice, dist_to_urb_2016)

refd_dist_2016 <- expand.grid(dist_to_urb_2016 = 0,
                              choice = "migrant") %>%
  mutate(att_id = dat[dat$choice == "migrant" & dat$se, ]$att_id[1]) %>%
  mutate(dist_to_urb_2016_sc = (dist_to_urb_2016 - unique(dat$mean_dist_to_urb_2016))/unique(dat$sd_dist_to_urb_2016))

refd_dist_2016_dup <- refd_dist_2016

for (i in 1:(nrow(predd_dist_2016) - 1)) {
  refd_dist_2016_dup <- rbind(refd_dist_2016_dup, refd_dist_2016)
}

nrow(refd_dist_2016_dup) == nrow(predd_dist_2016)

refd_dist_2016_dup <- refd_dist_2016_dup %>%
  arrange(dist_to_urb_2016)

pred_X <- model.matrix(mod_dist_2016_se, data = predd_dist_2016)
ref_X <- model.matrix(mod_dist_2016_se, data = refd_dist_2016_dup)
B <- coef(mod_dist_2016_se)

# Difference
diff_X <- pred_X - ref_X

# Linear predictor for difference
diff_dist_2016 <- diff_X %*% B

# Get variance of difference
var_pred <- diag(diff_X %*% vcov(mod_dist_2016_se) %*% t(diff_X))

# Standard error for difference
SE <- sqrt(var_pred)

predd_dist_2016_se <- predd_dist_2016 %>%
  mutate(log_rss = diff_dist_2016,
         se = SE,
         lwr = log_rss - 1.96*se,
         upr = log_rss + 1.96*se)

## 2010 ####

max_distances <- dat %>%
  filter(used == 1) %>%
  group_by(group) %>%
  summarize(xlim = max(dist_to_urb_2010, na.rm = TRUE))

### South Florida ####

predd_dist_2010 <- expand.grid(dist_to_urb_2010 = seq(0,
                                                      max_distances[max_distances$group == "South Florida", ]$xlim,
                                                      length.out = 100),
                               choice = c("migrant", "resident")) %>%
  mutate(att_id = case_when(
    choice == "migrant" ~ dat[dat$choice == "migrant" & dat$soflo, ]$att_id[1],
    choice == "resident" ~ dat[dat$choice == "resident" & dat$soflo, ]$att_id[1]
  )) %>%
  mutate(dist_to_urb_2010_sc = (dist_to_urb_2010 - unique(dat$mean_dist_to_urb_2010))/unique(dat$sd_dist_to_urb_2010)) %>%
  arrange(dist_to_urb_2010, choice)

refd_dist_2010 <- expand.grid(dist_to_urb_2010 = 0,
                              choice = c("migrant", "resident")) %>%
  mutate(att_id = case_when(
    choice == "migrant" ~ dat[dat$choice == "migrant" & dat$soflo, ]$att_id[1],
    choice == "resident" ~ dat[dat$choice == "resident" & dat$soflo, ]$att_id[1]
  )) %>%
  mutate(dist_to_urb_2010_sc = (dist_to_urb_2010 - unique(dat$mean_dist_to_urb_2010))/unique(dat$sd_dist_to_urb_2010)) %>%
  arrange(dist_to_urb_2010, choice)

refd_dist_2010_dup <- refd_dist_2010

for (i in 1:(nrow(predd_dist_2010)/2 - 1)) {
  refd_dist_2010_dup <- rbind(refd_dist_2010_dup, refd_dist_2010)
}

nrow(refd_dist_2010_dup) == nrow(predd_dist_2010)

pred_X <- model.matrix(mod_dist_2010_soflo, data = predd_dist_2010)
ref_X <- model.matrix(mod_dist_2010_soflo, data = refd_dist_2010_dup)
B <- coef(mod_dist_2010_soflo)

# Difference
diff_X <- pred_X - ref_X

# Linear predictor for difference
diff_dist_2010 <- diff_X %*% B

# Get variance of difference
var_pred <- diag(diff_X %*% vcov(mod_dist_2010_soflo) %*% t(diff_X))

# Standard error for difference
SE <- sqrt(var_pred)

predd_dist_2010_soflo <- predd_dist_2010 %>%
  mutate(log_rss = diff_dist_2010,
         se = SE,
         lwr = log_rss - 1.96*se,
         upr = log_rss + 1.96*se)

### Jacksonville ####

predd_dist_2010 <- expand.grid(dist_to_urb_2010 = seq(0,
                                                      max_distances[max_distances$group == "Jacksonville", ]$xlim,
                                                      length.out = 100),
                               choice = "resident") %>%
  mutate(att_id = dat[dat$choice == "resident" & dat$jax, ]$att_id[1]) %>%
  mutate(dist_to_urb_2010_sc = (dist_to_urb_2010 - unique(dat$mean_dist_to_urb_2010))/unique(dat$sd_dist_to_urb_2010)) %>%
  arrange(dist_to_urb_2010, choice)

refd_dist_2010 <- expand.grid(dist_to_urb_2010 = 0,
                              choice = "resident") %>%
  mutate(att_id = dat[dat$choice == "resident" & dat$jax, ]$att_id[1]) %>%
  mutate(dist_to_urb_2010_sc = (dist_to_urb_2010 - unique(dat$mean_dist_to_urb_2010))/unique(dat$sd_dist_to_urb_2010)) %>%
  arrange(dist_to_urb_2010, choice)

refd_dist_2010_dup <- refd_dist_2010

for (i in 1:(nrow(predd_dist_2010) - 1)) {
  refd_dist_2010_dup <- rbind(refd_dist_2010_dup, refd_dist_2010)
}

nrow(refd_dist_2010_dup) == nrow(predd_dist_2010)

pred_X <- model.matrix(mod_dist_2010_jax, data = predd_dist_2010)
ref_X <- model.matrix(mod_dist_2010_jax, data = refd_dist_2010_dup)
B <- coef(mod_dist_2010_jax)

# Difference
diff_X <- pred_X - ref_X

# Linear predictor for difference
diff_dist_2010 <- diff_X %*% B

# Get variance of difference
var_pred <- diag(diff_X %*% vcov(mod_dist_2010_jax) %*% t(diff_X))

# Standard error for difference
SE <- sqrt(var_pred)

predd_dist_2010_jax <- predd_dist_2010 %>%
  mutate(log_rss = diff_dist_2010,
         se = SE,
         lwr = log_rss - 1.96*se,
         upr = log_rss + 1.96*se)

### Southeast ####

predd_dist_2010 <- expand.grid(dist_to_urb_2010 = seq(0,
                                                      max_distances[max_distances$group == "Southeast", ]$xlim,
                                                      length.out = 100),
                               choice = "migrant") %>%
  mutate(att_id = dat[dat$choice == "migrant" & dat$se, ]$att_id[1]) %>%
  mutate(dist_to_urb_2010_sc = (dist_to_urb_2010 - unique(dat$mean_dist_to_urb_2010))/unique(dat$sd_dist_to_urb_2010)) %>%
  arrange(dist_to_urb_2010, choice)

refd_dist_2010 <- expand.grid(dist_to_urb_2010 = 0,
                              choice = "migrant") %>%
  mutate(att_id = dat[dat$choice == "migrant" & dat$se, ]$att_id[1]) %>%
  mutate(dist_to_urb_2010_sc = (dist_to_urb_2010 - unique(dat$mean_dist_to_urb_2010))/unique(dat$sd_dist_to_urb_2010)) %>%
  arrange(dist_to_urb_2010, choice)

refd_dist_2010_dup <- refd_dist_2010

for (i in 1:(nrow(predd_dist_2010) - 1)) {
  refd_dist_2010_dup <- rbind(refd_dist_2010_dup, refd_dist_2010)
}

nrow(refd_dist_2010_dup) == nrow(predd_dist_2010)

pred_X <- model.matrix(mod_dist_2010_se, data = predd_dist_2010)
ref_X <- model.matrix(mod_dist_2010_se, data = refd_dist_2010_dup)
B <- coef(mod_dist_2010_se)

# Difference
diff_X <- pred_X - ref_X

# Linear predictor for difference
diff_dist_2010 <- diff_X %*% B

# Get variance of difference
var_pred <- diag(diff_X %*% vcov(mod_dist_2010_se) %*% t(diff_X))

# Standard error for difference
SE <- sqrt(var_pred)

predd_dist_2010_se <- predd_dist_2010 %>%
  mutate(log_rss = diff_dist_2010,
         se = SE,
         lwr = log_rss - 1.96*se,
         upr = log_rss + 1.96*se)

# Plot model predictions for binary model ####

predd_bin_soflo$group <- "South Florida"
predd_bin_jax$group <- "Jacksonville"
predd_bin_se$group <- "Southeast"

predd_bin_soflo %>%
  bind_rows(predd_bin_jax) %>%
  bind_rows(predd_bin_se) %>%
  filter(urb_agg_2016 == 1) %>%
  ggplot(aes(x = choice,
             y = log_rss,
             col = choice,
             fill = choice)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ group) +
  scale_x_discrete(name = " ") +
  scale_y_continuous(name = "log-Relative selection strength") +
  scale_color_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

# Plot model predictions for distance models ####

## 2001 ####

n_used <- dat %>%
  filter(used == 1) %>%
  nrow()

n_avail <- dat %>%
  filter(used == 0) %>%
  nrow()

mig_used_2001 <- dat %>%
  filter(choice == "migrant" & used == 1) %>%
  mutate(dist_to_urban_km = dist_to_urb_2001/1000)

res_used_2001 <- dat %>%
  filter(choice == "resident" & used == 1) %>%
  mutate(dist_to_urban_km = dist_to_urb_2001/1000)

ggplot(predd_dist_2001, aes(x = dist_to_urban_km, y = log_rss,
                              group = choice, color = choice)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = choice), alpha = 0.2) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0, color = "gray20", linetype = 2) +
  geom_rug(aes(x = dist_to_urban_km, y = NULL), data = mig_used_2001, sides = "b", alpha = 0.2) +
  geom_rug(aes(x = dist_to_urban_km, y = NULL), data = res_used_2001, sides = "t", alpha = 0.2) +
  scale_x_continuous(name = "Distance to urban (km)") +
  scale_y_continuous(name = "Relative selection strength") +
  scale_color_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  theme_bw()

## 2016 ####

mig_used_2016 <- dat %>%
  filter(choice == "migrant" & used == 1) %>%
  mutate(dist_to_urban_km = dist_to_urb_2016/1000)

res_used_2016 <- dat %>%
  filter(choice == "resident" & used == 1) %>%
  mutate(dist_to_urban_km = dist_to_urb_2016/1000)

predd_dist_2016_soflo$group <- "South Florida"
predd_dist_2016_jax$group <- "Jacksonville"
predd_dist_2016_se$group <- "Southeast"

predd_dist_2016_soflo %>%
  bind_rows(predd_dist_2016_jax) %>%
  bind_rows(predd_dist_2016_se) %>%
  mutate(dist_to_urban_km = dist_to_urb_2016/1000) %>%
  ggplot(aes(x = dist_to_urban_km, y = log_rss,
             group = choice, color = choice)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = choice), alpha = 0.2) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0, color = "gray20", linetype = 2) +
  #geom_rug(aes(x = dist_to_urban_km, y = NULL), data = mig_used_2016, sides = "b") +
  #geom_rug(aes(x = dist_to_urban_km, y = NULL), data = res_used_2016, sides = "t") +
  facet_wrap(~ group) +
  scale_x_continuous(name = "Distance to urban (km)") +
  scale_y_continuous(name = "Relative selection strength") +
  scale_color_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  theme_bw()

## 2010 ####

mig_used_2010 <- dat %>%
  filter(choice == "migrant" & used == 1) %>%
  mutate(dist_to_urban_km = dist_to_urb_2010/1000)

res_used_2010 <- dat %>%
  filter(choice == "resident" & used == 1) %>%
  mutate(dist_to_urban_km = dist_to_urb_2010/1000)

predd_dist_2010_soflo$group <- "South Florida"
predd_dist_2010_jax$group <- "Jacksonville"
predd_dist_2010_se$group <- "Southeast"

predd_dist_2010_soflo %>%
  bind_rows(predd_dist_2010_jax) %>%
  bind_rows(predd_dist_2010_se) %>%
  mutate(dist_to_urban_km = dist_to_urb_2010/1000) %>%
  ggplot(aes(x = dist_to_urban_km, y = log_rss,
             group = choice, color = choice)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = choice), alpha = 0.2) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0, color = "gray20", linetype = 2) +
  geom_rug(aes(x = dist_to_urban_km, y = NULL), data = mig_used_2010, sides = "t", alpha = 0.2) +
  geom_rug(aes(x = dist_to_urban_km, y = NULL), data = res_used_2010, sides = "b", alpha = 0.2) +
  facet_wrap(~ group) +
  scale_x_continuous(name = "Distance to urban (km)") +
  scale_y_continuous(name = "Relative selection strength") +
  scale_color_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  theme_bw()
