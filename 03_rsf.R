# Code for Picardi et al., "Fitness consequences of anthropogenic subsidies
# for a partially migratory wading bird"

# Load packages ####

library(tidyverse)
library(survival)

# Load log-RSS prediction function ####

source("predict_log-rss.R")

# Load data ####

dat <- readRDS("input/rsf-data.rds")

# Choose which datasets to use for predictors ####

# Choose which Landfire version
year <- 2016 # 2001 or 2016

# Distance
# Use the distance calculated from the aggregated Landfire raster?
# If FALSE, defaults to the non-aggregated Landfire raster
use_agg <- FALSE

# Binary
# Use the binary urban variable from TIGER?
# If FALSE, defaults to the binary aggregated Landfire raster
use_tiger <- FALSE

# Get the corresponding data:
# Distance
if (use_tiger) {
  dat <- dat %>%
    mutate(
      dist_to_urb = dist_to_urb_2010,
      dist_to_urb_sc = dist_to_urb_2010_sc,
      mean_dist_to_urb = mean_dist_to_urb_2010,
      sd_dist_to_urb = sd_dist_to_urb_2010
    )

} else if (use_agg & year == 2016) {
  dat <- dat %>%
    mutate(
      dist_to_urb = dist_to_urb_agg_2016,
      dist_to_urb_sc = dist_to_urb_agg_2016_sc,
      mean_dist_to_urb = mean_dist_to_urb_agg_2016,
      sd_dist_to_urb = sd_dist_to_urb_agg_2016
    )

} else if (!use_agg & year == 2016) {
  dat <- dat %>%
    mutate(
      dist_to_urb = dist_to_urb_2016,
      dist_to_urb_sc = dist_to_urb_2016_sc,
      mean_dist_to_urb = mean_dist_to_urb_2016,
      sd_dist_to_urb = sd_dist_to_urb_2016
    )

} else if (use_agg & year == 2001) {
  dat <- dat %>%
    mutate(
      dist_to_urb = dist_to_urb_agg_2001,
      dist_to_urb_sc = dist_to_urb_agg_2001_sc,
      mean_dist_to_urb = mean_dist_to_urb_agg_2001,
      sd_dist_to_urb = sd_dist_to_urb_agg_2001
    )

} else if (!use_agg & year == 2001) {
  dat <- dat %>%
    mutate(
      dist_to_urb = dist_to_urb_2001,
      dist_to_urb_sc = dist_to_urb_2001_sc,
      mean_dist_to_urb = mean_dist_to_urb_2001,
      sd_dist_to_urb = sd_dist_to_urb_2001
    )

}

# Binary
if (use_tiger) {
  dat <- dat %>%
    mutate(urb_bin = urb_tiger_2010)

} else if (!use_tiger & year == 2001) {
  dat <- dat %>%
    mutate(urb_bin = urb_agg_2001)

} else if (!use_tiger & year == 2016) {
  dat <- dat %>%
    mutate(urb_bin = urb_agg_2016)

}

# Create regional groups ####

# South Florida, Jacksonville, Southeast
dat <- dat %>%
  mutate(
    se = !soflo & !jax,
    group = case_when(jax ~ "Jacksonville",
                      soflo ~ "South Florida",
                      TRUE ~ "Southeast")
  )

# Sample sizes
dat %>%
  select(att_id, group, choice) %>%
  distinct() %>%
  group_by(group, choice) %>%
  tally()

# RSF for binary developed/undeveloped ####

## South Florida ####

mod_bin_soflo <- clogit(
  formula = used ~
    urb_bin +
    urb_bin:choice +
    strata(att_id),
  data = dat,
  subset = soflo,
  weights = weight,
  method = "approximate"
)

## Jacksonville ####

mod_bin_jax <- clogit(
  formula = used ~
    urb_bin +
    strata(att_id),
  data = dat,
  subset = jax,
  weights = weight,
  method = "approximate"
)

## Southeast ####

mod_bin_se <- clogit(
  formula = used ~
    urb_bin +
    strata(att_id),
  data = dat,
  subset = se,
  weights = weight,
  method = "approximate"
)

# RSF for distance to development ####

### South Florida ####

mod_dist_soflo <- clogit(
  formula = used ~
    dist_to_urb_sc +
    I(dist_to_urb_sc ^ 2) +
    dist_to_urb_sc:choice +
    I(dist_to_urb_sc ^ 2):choice +
    strata(att_id),
  data = dat,
  subset = soflo,
  weights = weight,
  method = "approximate"
)

### Jacksonville ####

mod_dist_jax <- clogit(
  formula = used ~
    dist_to_urb_sc +
    I(dist_to_urb_sc ^ 2) +
    strata(att_id),
  data = dat,
  subset = jax,
  weights = weight,
  method = "approximate"
)

### Southeast ####

mod_dist_se <- clogit(
  formula = used ~
    dist_to_urb_sc +
    I(dist_to_urb_sc ^ 2) +
    strata(att_id),
  data = dat,
  subset = se,
  weights = weight,
  method = "approximate"
)

# Model predictions for binary model ####

## South Florida ####

predd_bin <- expand.grid(urb_bin = c(0, 1),
                         choice = c("migrant", "resident")) %>%
  mutate(att_id = case_when(
    choice == "migrant" ~ dat[dat$choice == "migrant" &
                                dat$soflo,]$att_id[1],
    choice == "resident" ~ dat[dat$choice == "resident" &
                                 dat$soflo,]$att_id[1]
  )) %>%
  arrange(urb_bin, choice)

refd_bin <- expand.grid(urb_bin = 0,
                        choice = c("migrant", "resident")) %>%
  mutate(att_id = case_when(
    choice == "migrant" ~ dat[dat$choice == "migrant" &
                                dat$soflo,]$att_id[1],
    choice == "resident" ~ dat[dat$choice == "resident" &
                                 dat$soflo,]$att_id[1]
  )) %>%
  arrange(urb_bin, choice)

refd_bin <- rbind(refd_bin, refd_bin)

predd_bin_soflo <- pred_log_rss(mod = mod_bin_soflo,
                                x1 = predd_bin,
                                x2 = refd_bin)

## Jacksonville ####

predd_bin <- expand.grid(urb_bin = c(0, 1),
                         choice = "resident") %>%
  mutate(att_id = dat[dat$choice == "resident" &
                        dat$jax,]$att_id[1]) %>%
  arrange(urb_bin, choice)

refd_bin <- expand.grid(urb_bin = 0,
                        choice = "resident") %>%
  mutate(att_id = dat[dat$choice == "resident" &
                        dat$jax,]$att_id[1]) %>%
  arrange(urb_bin, choice)

refd_bin <- rbind(refd_bin, refd_bin)

predd_bin_jax <- pred_log_rss(mod = mod_bin_jax,
                              x1 = predd_bin,
                              x2 = refd_bin)

## Southeast ####

predd_bin <- expand.grid(urb_bin = c(0, 1),
                         choice = "migrant") %>%
  mutate(att_id = dat[dat$choice == "migrant" &
                        dat$se,]$att_id[1]) %>%
  arrange(urb_bin, choice)

refd_bin <- expand.grid(urb_bin = 0,
                        choice = "migrant") %>%
  mutate(att_id = dat[dat$choice == "migrant" &
                        dat$se,]$att_id[1]) %>%
  arrange(urb_bin, choice)

refd_bin <- rbind(refd_bin, refd_bin)

predd_bin_se <- pred_log_rss(mod = mod_bin_soflo,
                             x1 = predd_bin,
                             x2 = refd_bin)

# Model predictions for distance models ####

max_distances <- dat %>%
  filter(used == 1) %>%
  group_by(group) %>%
  summarize(xlim = max(dist_to_urb, na.rm = TRUE))

### South Florida ####

predd_dist <- expand.grid(
  dist_to_urb = seq(0,
                    max_distances[max_distances$group == "South Florida",]$xlim,
                    length.out = 100),
  choice = c("migrant", "resident")
) %>%
  mutate(att_id = case_when(
    choice == "migrant" ~ dat[dat$choice == "migrant" &
                                dat$soflo,]$att_id[1],
    choice == "resident" ~ dat[dat$choice == "resident" &
                                 dat$soflo,]$att_id[1]
  )) %>%
  mutate(dist_to_urb_sc = (dist_to_urb - unique(dat$mean_dist_to_urb)) /
           unique(dat$sd_dist_to_urb)) %>%
  arrange(choice, dist_to_urb)

refd_dist <- expand.grid(dist_to_urb = 0,
                         choice = c("migrant", "resident")) %>%
  mutate(att_id = case_when(
    choice == "migrant" ~ dat[dat$choice == "migrant" &
                                dat$soflo,]$att_id[1],
    choice == "resident" ~ dat[dat$choice == "resident" &
                                 dat$soflo,]$att_id[1]
  )) %>%
  mutate(dist_to_urb_sc = (dist_to_urb - unique(dat$mean_dist_to_urb)) /
           unique(dat$sd_dist_to_urb))

refd_dist_dup <- refd_dist

for (i in 1:(nrow(predd_dist) / 2 - 1)) {
  refd_dist_dup <- rbind(refd_dist_dup, refd_dist)
}

nrow(refd_dist_dup) == nrow(predd_dist)

refd_dist_dup <- refd_dist_dup %>%
  arrange(choice, dist_to_urb)

predd_dist_soflo <- pred_log_rss(mod = mod_dist_soflo,
                                 x1 = predd_dist,
                                 x2 = refd_dist_dup)

### Jacksonville ####

predd_dist <- expand.grid(dist_to_urb = seq(0,
                                            max_distances[max_distances$group == "Jacksonville",]$xlim,
                                            length.out = 100),
                          choice = "resident") %>%
  mutate(att_id = dat[dat$choice == "resident" &
                        dat$jax,]$att_id[1]) %>%
  mutate(dist_to_urb_sc = (dist_to_urb - unique(dat$mean_dist_to_urb)) /
           unique(dat$sd_dist_to_urb)) %>%
  arrange(choice, dist_to_urb)

refd_dist <- expand.grid(dist_to_urb = 0,
                         choice = "resident") %>%
  mutate(att_id = dat[dat$choice == "resident" &
                        dat$jax,]$att_id[1]) %>%
  mutate(dist_to_urb_sc = (dist_to_urb - unique(dat$mean_dist_to_urb)) /
           unique(dat$sd_dist_to_urb))

refd_dist_dup <- refd_dist

for (i in 1:(nrow(predd_dist) - 1)) {
  refd_dist_dup <- rbind(refd_dist_dup, refd_dist)
}

nrow(refd_dist_dup) == nrow(predd_dist)

refd_dist_dup <- refd_dist_dup %>%
  arrange(choice, dist_to_urb)

predd_dist_jax <- pred_log_rss(mod = mod_dist_jax,
                               x1 = predd_dist,
                               x2 = refd_dist_dup)

### Southeast ####

predd_dist <- expand.grid(dist_to_urb = seq(0,
                                            max_distances[max_distances$group == "Southeast",]$xlim,
                                            length.out = 100),
                          choice = "migrant") %>%
  mutate(att_id = dat[dat$choice == "migrant" &
                        dat$se,]$att_id[1]) %>%
  mutate(dist_to_urb_sc = (dist_to_urb - unique(dat$mean_dist_to_urb)) /
           unique(dat$sd_dist_to_urb)) %>%
  arrange(choice, dist_to_urb)

refd_dist <- expand.grid(dist_to_urb = 0,
                         choice = "migrant") %>%
  mutate(att_id = dat[dat$choice == "migrant" &
                        dat$se,]$att_id[1]) %>%
  mutate(dist_to_urb_sc = (dist_to_urb - unique(dat$mean_dist_to_urb)) /
           unique(dat$sd_dist_to_urb))

refd_dist_dup <- refd_dist

for (i in 1:(nrow(predd_dist) - 1)) {
  refd_dist_dup <- rbind(refd_dist_dup, refd_dist)
}

nrow(refd_dist_dup) == nrow(predd_dist)

refd_dist_dup <- refd_dist_dup %>%
  arrange(dist_to_urb)

predd_dist_se <- pred_log_rss(mod = mod_dist_soflo,
                              x1 = predd_dist,
                              x2 = refd_dist_dup)

# Plot model predictions for binary model ####

predd_bin_soflo$group <- "South Florida"
predd_bin_jax$group <- "Jacksonville"
predd_bin_se$group <- "Southeast"

predd_bin_soflo %>%
  bind_rows(predd_bin_jax) %>%
  bind_rows(predd_bin_se) %>%
  filter(urb_bin == 1) %>%
  ggplot(aes(x = choice, y = log_rss, col = choice, fill = choice)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap( ~ group) +
  scale_x_discrete(name = " ") +
  scale_y_continuous(name = "log-Relative selection strength") +
  scale_color_viridis_d("Tactic", labels = c("Migrant", "Resident"),
                        begin = 0.25, end = 0.75) +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"),
                       begin = 0.25, end = 0.75) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

if (use_tiger) {
  ggsave("output/RSF_bin_tiger_2010.tiff",
       width = 180, height = 90, units = "mm",
       compression = "lzw", dpi = 600)
} else {
  ggsave(paste0("output/RSF_bin_agg_", year, ".tiff"),
         width = 180, height = 90, units = "mm",
         compression = "lzw", dpi = 600)
  }

# Plot model predictions for distance models ####

mig_used <- dat %>%
  filter(choice == "migrant" & used == 1) %>%
  mutate(dist_to_urban_km = dist_to_urb / 1000)

res_used <- dat %>%
  filter(choice == "resident" & used == 1) %>%
  mutate(dist_to_urban_km = dist_to_urb / 1000)

predd_dist_soflo$group <- "South Florida"
predd_dist_jax$group <- "Jacksonville"
predd_dist_se$group <- "Southeast"

predd_dist_soflo %>%
  bind_rows(predd_dist_jax) %>%
  bind_rows(predd_dist_se) %>%
  filter(!(
    choice == "resident" &
      dist_to_urb > max(res_used$dist_to_urb, na.rm = T)
  )) %>%
  mutate(dist_to_urban_km = dist_to_urb / 1000) %>%
  #mutate(group = factor(group, levels = c("South Florida", "Jacksonville", "Southeast")))
  ggplot(aes(x = dist_to_urban_km, y = log_rss, group = choice, color = choice)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = choice), alpha = 0.2) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0, color = "gray20", linetype = 2) +
  geom_rug(aes(x = dist_to_urban_km, y = NULL), data = mig_used,
           sides = "t", alpha = 0.2) +
  geom_rug(aes(x = dist_to_urban_km, y = NULL), data = res_used,
           sides = "b", alpha = 0.2) +
  facet_wrap( ~ factor(group)) +
  scale_x_continuous(name = "Distance to urban (km)") +
  scale_y_continuous(name = "log-Relative selection strength") +
  scale_color_viridis_d("Tactic", labels = c("Migrant", "Resident"),
                        begin = 0.25, end = 0.75) +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"),
                       begin = 0.25, end = 0.75) +
  theme_bw()

if (use_tiger) {
  ggsave("output/RSF_dist_tiger_2010.tiff",
         width = 180, height = 90, units = "mm",
         compression = "lzw", dpi = 600)
} else {
  ggsave(paste0("output/RSF_dist_agg", use_agg, "_", year , ".tiff"),
         width = 180, height = 90, units = "mm",
         compression = "lzw", dpi = 600)
}

# Summaries ####

predd_dist_soflo %>%
  bind_rows(predd_dist_jax) %>%
  bind_rows(predd_dist_se) %>%
  group_by(group, choice) %>%
  filter(log_rss > 0) %>%
  summarize(min(dist_to_urb)/1000, max(dist_to_urb)/1000)

max_lRSS <- predd_dist_soflo %>%
  bind_rows(predd_dist_jax) %>%
  bind_rows(predd_dist_se) %>%
  group_by(group, choice) %>%
  summarize(max_log_rss = max(log_rss))

predd_dist_soflo %>%
  bind_rows(predd_dist_jax) %>%
  bind_rows(predd_dist_se) %>%
  left_join(max_lRSS) %>%
  filter(log_rss == max_log_rss) %>%
  mutate(dist_to_urb_km = dist_to_urb/1000)
