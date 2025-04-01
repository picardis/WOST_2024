# Make plots focusing only on SoFlo to use in presentations.

# Load packages ####

library(tidyverse)
library(survival)
library(nestR)

# Load log-RSS prediction function ####

source("predict_log-rss.R")

# Load data ####

dat <- readRDS("input/rsf-data.rds")

# Choose which dataset to use for predictors ####

dat <- dat %>%
    mutate(
      dist_to_urb = dist_to_urb_2016,
      dist_to_urb_sc = dist_to_urb_2016_sc,
      mean_dist_to_urb = mean_dist_to_urb_2016,
      sd_dist_to_urb = sd_dist_to_urb_2016
    )

# Create regional groups ####

# South Florida, Jacksonville, Southeast
dat <- dat %>%
  mutate(
    se = !soflo & !jax,
    group = case_when(jax ~ "Jacksonville",
                      soflo ~ "South Florida",
                      TRUE ~ "Southeast")
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

# Plot model predictions for distance models ####

mig_used <- dat %>%
  filter(choice == "migrant" & used == 1) %>%
  mutate(dist_to_urban_km = dist_to_urb / 1000)

res_used <- dat %>%
  filter(choice == "resident" & used == 1) %>%
  mutate(dist_to_urban_km = dist_to_urb / 1000)

predd_dist_soflo %>%
  filter(!(
    choice == "resident" &
      dist_to_urb > max(res_used$dist_to_urb, na.rm = T)
  )) %>%
  mutate(dist_to_urban_km = dist_to_urb / 1000) %>%
  ggplot(aes(x = dist_to_urban_km, y = log_rss, group = choice, color = choice)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = choice), alpha = 0.2) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0, color = "gray20", linetype = 2) +
  geom_rug(aes(x = dist_to_urban_km, y = NULL), data = mig_used,
           sides = "t", alpha = 0.2) +
  geom_rug(aes(x = dist_to_urban_km, y = NULL), data = res_used,
           sides = "b", alpha = 0.2) +
  scale_x_continuous(name = "Distance to urban (km)") +
  scale_y_continuous(name = "log-Relative selection strength") +
  scale_color_viridis_d("Tactic", labels = c("Migrant", "Resident"),
                        begin = 0.25, end = 0.75) +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"),
                       begin = 0.25, end = 0.75) +
  theme_bw()

ggsave("output/RSF_SoFlo.tiff",
       width = 180, height = 90, units = "mm",
       compression = "lzw", dpi = 600)

# Load survival model ####

post <- readRDS("input/nest_survival_by_region_2024-06-05.rds")

# Load predictors ####

urb_matrix <- readRDS("input/surv_model_predictors.rds") %>%
  mutate(dist_urb_km = avg_dist_urb/1000) %>%
  mutate(region = case_when(jax == 1 ~ "Jacksonville",
                            soflo == 1 ~ "South Florida",
                            se == 1 ~ "Southeast"))

# Load original data ####

rsf_data <- readRDS("input/rsf-data.rds") %>%
  mutate(
    se = !soflo & !jax,
    region = case_when(jax ~ "Jacksonville",
                       soflo ~ "South Florida",
                       TRUE ~ "Southeast")
  )

ws_nests <- readRDS("input/ws_nests_all.rds")

# Create 'visits' and 'fixes'
ws_attempts <- format_attempts(ws_nests, 110)

rownames(ws_attempts$fixes) <- gsub(pattern = "_",
                                    replacement = "-",
                                    x = rownames(ws_attempts$fixes),
                                    fixed = TRUE)
rownames(ws_attempts$visits) <- gsub(pattern = "_",
                                     replacement = "-",
                                     x = rownames(ws_attempts$visits),
                                     fixed = TRUE)

fixes <- ws_attempts$fixes[unique(rsf_data$att_id), ]
visits <- ws_attempts$visits[unique(rsf_data$att_id), ]

# Survival probability phi ####

phis <- data.frame(phi.b0 = as.vector(post$phi.b0),
                   phi.b1 = as.vector(post$phi.b1),
                   phi.b2 = as.vector(post$phi.b2),
                   phi.b3 = as.vector(post$phi.b3),
                   phi.b4 = as.vector(post$phi.b4),
                   phi.b5 = as.vector(post$phi.b5))

max_distances <- urb_matrix %>%
  group_by(region) %>%
  summarize(lim = max(avg_dist_urb))

preds <- data.frame(dist_urb = c(
  seq(0,
      max_distances[max_distances$region == "South Florida", ]$lim,
      length.out = 100),
  seq(0,
      max_distances[max_distances$region == "Southeast", ]$lim,
      length.out = 100),
  seq(0,
      max_distances[max_distances$region == "Jacksonville", ]$lim,
      length.out = 100)
),
soflo = c(rep(1, 100), rep(0, 200)),
se = c(rep(0, 100), rep(1, 100), rep(0, 100)),
jax = c(rep(0, 200), rep(1, 100)))

est <- lapply(1:nrow(phis), function(iter) {

  est <- preds %>%
    mutate(phi = (plogis((phis[iter, ]$phi.b0 +
                            (phis[iter, ]$phi.b1 * dist_urb) +
                            (phis[iter, ]$phi.b2 * se) +
                            (phis[iter, ]$phi.b3 * jax) +
                            (phis[iter, ]$phi.b4 * dist_urb * se) +
                            (phis[iter, ]$phi.b5 * dist_urb * jax)))),
           iter = iter)

  return(est)

})

est <- est %>%
  bind_rows() %>%
  mutate(region = case_when(jax == 1 ~ "Jacksonville",
                            soflo == 1 ~ "South Florida",
                            se == 1 ~ "Southeast"),
         dist_urb_km = dist_urb/1000) %>%
  group_by(dist_urb, dist_urb_km, region) %>%
  summarize(phi_mean = mean(phi),
            phi_lwr = quantile(phi, 0.025),
            phi_upr = quantile(phi, 0.975)) %>%
  filter(region == "South Florida")

ggplot(est, aes(x = dist_urb_km, y = phi_mean)) +
  geom_ribbon(aes(ymin = phi_lwr, ymax = phi_upr), fill = "gray80") +
  geom_line(linewidth = 1.2) +
  geom_rug(aes(x = dist_urb_km, y = NULL),
           data = urb_matrix, sides = "t", alpha = 0.2) +
  scale_x_continuous(name = "Distance to urban (km)") +
  scale_y_continuous(name = "Daily nest survival") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("output/surv_phi_SoFlo.tiff",
       width = 180, height = 90, units = "mm", compression = "lzw", dpi = 600)

# Cumulative survival z ####

# Getting positional indices of each group
mig_soflo <- which(rownames(visits) %in% rsf_data[rsf_data$choice == "migrant" &
                                                    rsf_data$region == "South Florida", ]$att_id)
res_soflo <- which(rownames(visits) %in% rsf_data[rsf_data$choice == "resident" &
                                                    rsf_data$region == "South Florida", ]$att_id)
mig_se <- which(rownames(visits) %in% rsf_data[rsf_data$choice == "migrant" &
                                                 rsf_data$region == "Southeast", ]$att_id)
res_jax <- which(rownames(visits) %in% rsf_data[rsf_data$choice == "resident" &
                                                  rsf_data$region == "Jacksonville", ]$att_id)

combos <- rsf_data %>%
  select(region, choice) %>%
  distinct()

z_preds <- data.frame()

for (i in 1:110) {
  combos$day <- i
  combos <- combos %>%
    mutate(z_mean = case_when(
      region == "South Florida" &
        choice == "migrant" ~ mean(as.vector(post$z[mig_soflo, i, , ])),
      region == "South Florida" &
        choice == "resident" ~ mean(as.vector(post$z[res_soflo, i, , ])),
      region == "Southeast" &
        choice == "migrant" ~ mean(as.vector(post$z[mig_se, i, , ])),
      region == "Jacksonville" &
        choice == "resident" ~ mean(as.vector(post$z[res_jax, i, , ]))
    ),
    z_lwr = case_when(
      region == "South Florida" &
        choice == "migrant" ~ quantile(as.vector(post$z[mig_soflo, i, , ]), 0.025),
      region == "South Florida" &
        choice == "resident" ~ quantile(as.vector(post$z[res_soflo, i, , ]), 0.025),
      region == "Southeast" &
        choice == "migrant" ~ quantile(as.vector(post$z[mig_se, i, , ]), 0.025),
      region == "Jacksonville" &
        choice == "resident" ~ quantile(as.vector(post$z[res_jax, i, , ]), 0.025)
    ),
    z_upr = case_when(
      region == "South Florida" &
        choice == "migrant" ~ quantile(as.vector(post$z[mig_soflo, i, , ]), 0.975),
      region == "South Florida" &
        choice == "resident" ~ quantile(as.vector(post$z[res_soflo, i, , ]), 0.975),
      region == "Southeast" &
        choice == "migrant" ~ quantile(as.vector(post$z[mig_se, i, , ]), 0.975),
      region == "Jacksonville" &
        choice == "resident" ~ quantile(as.vector(post$z[res_jax, i, , ]), 0.975)
    ))
  z_preds <- rbind(z_preds, combos)
}

z_preds %>%
  filter(region == "South Florida") %>%
  ggplot(aes(x = day, y = z_mean, color = choice)) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  labs(x = "Days", y = "Cumulative nest survival") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

ggsave("output/surv_z_SoFlo.tiff",
       width = 180, height = 90, units = "mm", compression = "lzw", dpi = 600)
