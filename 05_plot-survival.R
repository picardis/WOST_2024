# Load packages ####

library(tidyverse)
library(rjags)
library(coda)
library(nestR)

# Load model ####

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

# Diagnostics ####

# This informs us of how much thinning should we use (done for each beta parameter)
autocorr.plot(as.mcmc.list(post$phi.b0), lag.max = 20)
autocorr.plot(as.mcmc.list(post$phi.b1), lag.max = 20)
autocorr.plot(as.mcmc.list(post$phi.b2), lag.max = 20)
autocorr.plot(as.mcmc.list(post$phi.b3), lag.max = 20)
autocorr.plot(as.mcmc.list(post$phi.b4), lag.max = 20)
autocorr.plot(as.mcmc.list(post$phi.b5), lag.max = 20)

traceplot(as.mcmc.list(post$phi.b0))
traceplot(as.mcmc.list(post$phi.b1))
traceplot(as.mcmc.list(post$phi.b2))
traceplot(as.mcmc.list(post$phi.b3))
traceplot(as.mcmc.list(post$phi.b4))
traceplot(as.mcmc.list(post$phi.b5))
traceplot(as.mcmc.list(post$p.b0))
traceplot(as.mcmc.list(post$p.b1))

gelman.diag(post$phi.b0)
gelman.diag(post$phi.b1)
gelman.diag(post$phi.b2)
gelman.diag(post$phi.b3)
gelman.diag(post$phi.b4)
gelman.diag(post$phi.b5)

densplot(as.mcmc.list(post$phi.b0))
abline(v = 0, col = "red")
densplot(as.mcmc.list(post$phi.b1))
abline(v = 0, col = "red")
densplot(as.mcmc.list(post$phi.b2))
abline(v = 0, col = "red")
densplot(as.mcmc.list(post$phi.b3))
abline(v = 0, col = "red")
densplot(as.mcmc.list(post$phi.b4))
abline(v = 0, col = "red")
densplot(as.mcmc.list(post$phi.b5))
abline(v = 0, col = "red")
densplot(as.mcmc.list(post$p.b0))
abline(v = 0, col = "red")
densplot(as.mcmc.list(post$p.b1))
abline(v = 0, col = "red")

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
            phi_upr = quantile(phi, 0.975))

ggplot(est, aes(x = dist_urb_km, y = phi_mean)) +
  geom_ribbon(aes(ymin = phi_lwr, ymax = phi_upr), fill = "gray80") +
  geom_line(linewidth = 1.2) +
  geom_rug(aes(x = dist_urb_km, y = NULL),
           data = urb_matrix, sides = "t", alpha = 0.2) +
  facet_wrap(~ region) +
  scale_x_continuous(name = "Distance to urban (km)") +
  scale_y_continuous(name = "Daily nest survival") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("output/surv_phi.tiff",
       width = 180, height = 90, units = "mm", compression = "lzw", dpi = 600)

# Detection probability p ####

p_df <- data.frame(t = 1:110,
                   p_lwr = apply(post$p, 1, quantile, 0.025),
                   p_mean = apply(post$p, 1, mean),
                   p_upr = apply(post$p, 1, quantile, 0.975))

ggplot(p_df, aes(x = t, y = p_mean)) +
  geom_ribbon(aes(ymin = p_lwr, ymax = p_upr), fill = "tomato", alpha = 0.3) +
  geom_line(linewidth = 0.1, color = "tomato") +
  labs(x = "Day of nesting attempt", y = "Detection probability") +
  theme_bw() +
  coord_cartesian(ylim = c(0, 1))

ggsave("output/surv_p.tiff",
       width = 80, height = 60, units = "mm", compression = "lzw", dpi = 600)

# Cumulative survival z ####

# Testing the object's dimensions and how to subset it
dim(post$z)
dim(post$z[1, , , ]) #individual
dim(post$z[, 1, , ]) # day
dim(post$z[, , 1, ]) # iteration
dim(post$z[, , , 1]) # chain

mean(post$z[, 110, , ])
mean(as.vector(post$z[1:145, 110, , ]))

# Getting positional indices of each group
mig_soflo <- which(rownames(visits) %in% rsf_data[rsf_data$choice == "migrant" &
                                                    rsf_data$region == "South Florida", ]$att_id)
res_soflo <- which(rownames(visits) %in% rsf_data[rsf_data$choice == "resident" &
                                                    rsf_data$region == "South Florida", ]$att_id)
mig_se <- which(rownames(visits) %in% rsf_data[rsf_data$choice == "migrant" &
                                                 rsf_data$region == "Southeast", ]$att_id)
res_jax <- which(rownames(visits) %in% rsf_data[rsf_data$choice == "resident" &
                                                  rsf_data$region == "Jacksonville", ]$att_id)

# z

mean(as.vector(post$z[mig_soflo, , , ]))
mean(as.vector(post$z[res_soflo, , , ]))
mean(as.vector(post$z[mig_se, , , ]))
mean(as.vector(post$z[res_jax, , , ]))

quantile(as.vector(post$z[mig_soflo, , , ]), 0.025)
quantile(as.vector(post$z[res_soflo, , , ]), 0.025)
quantile(as.vector(post$z[mig_se, , , ]), 0.025)
quantile(as.vector(post$z[res_jax, , , ]), 0.025)

quantile(as.vector(post$z[mig_soflo, , , ]), 0.975)
quantile(as.vector(post$z[res_soflo, , , ]), 0.975)
quantile(as.vector(post$z[mig_se, , , ]), 0.975)
quantile(as.vector(post$z[res_jax, , , ]), 0.975)

quantile(as.vector(post$z[mig_soflo, , , ]))
quantile(as.vector(post$z[res_soflo, , , ]))
quantile(as.vector(post$z[mig_se, , , ]))
quantile(as.vector(post$z[res_jax, , , ]))

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

ggplot(z_preds, aes(x = day, y = z_mean, color = choice)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ region) +
  scale_color_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  labs(x = "Days", y = "Cumulative nest survival") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

ggsave("output/surv_z.tiff",
       width = 180, height = 90, units = "mm", compression = "lzw", dpi = 600)

# z 110
mean(as.vector(post$z[mig_soflo, 110, , ]))
mean(as.vector(post$z[res_soflo, 110, , ]))
mean(as.vector(post$z[mig_se, 110, , ]))
mean(as.vector(post$z[res_jax, 110, , ]))

quantile(as.vector(post$z[mig_soflo, 110, , ]), 0.025)
quantile(as.vector(post$z[res_soflo, 110, , ]), 0.025)
quantile(as.vector(post$z[mig_se, 110, , ]), 0.025)
quantile(as.vector(post$z[res_jax, 110, , ]), 0.025)

quantile(as.vector(post$z[mig_soflo, 110, , ]), 0.975)
quantile(as.vector(post$z[res_soflo, 110, , ]), 0.975)
quantile(as.vector(post$z[mig_se, 110, , ]), 0.975)
quantile(as.vector(post$z[res_jax, 110, , ]), 0.975)

quantile(as.vector(post$z[mig_soflo, 110, , ]))
quantile(as.vector(post$z[res_soflo, 110, , ]))
quantile(as.vector(post$z[mig_se, 110, , ]))
quantile(as.vector(post$z[res_jax, 110, , ]))

z110_df <- z_preds %>%
  filter(day == 110)

ggplot(z110_df, aes(x = choice, y = z_mean, color = choice)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = z_lwr, ymax = z_upr), width = 0.3) +
  facet_wrap(~ region) +
  scale_color_viridis_d("Choice", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  labs(x = " ", y = "Nest survival to day 110") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave("output/surv_z110.tiff",
       width = 180, height = 90, units = "mm", compression = "lzw", dpi = 600)

# Trends in survival ####

surv_trends <- data.frame(att_id = rownames(visits))

mig_yrs <- readRDS("input/mig_yrs.rds")
urban_rsf_df <- readRDS("input/rsf-data.rds")

att_info <- urban_rsf_df %>%
  group_by(att_id) %>%
  arrange(date_ymd) %>%
  slice(1) %>%
  mutate(year = year(date_ymd)) %>%
  filter(!is.na(year)) %>%
  select(att_id, year, group, choice) %>%
  distinct() %>%
  filter(att_id %in% surv_trends$att_id)

surv_trends <- left_join(surv_trends, att_info)

class(post$z[, 110, , ])

mean_z <- c()
lwr_z <- c()
upr_z <- c()

for (i in 1:145) {
  mean_z[i] <- quantile(post$z[i, 110, , ], 0.5)
  lwr_z[i] <- quantile(post$z[i, 110, , ], 0.025)
  upr_z[i] <- quantile(post$z[i, 110, , ], 0.975)
}

surv_trends$mean_z <- mean_z
surv_trends$lwr_z <- lwr_z
surv_trends$upr_z <- upr_z

surv_trends %>%
  group_by(year, group, choice) %>%
  summarize(mean_z = mean(mean_z),
            lwr_z = mean(lwr_z),
            upr_z = mean(upr_z)) %>%
  ggplot(aes(x = factor(year), y = mean_z, color = choice)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr_z, ymax = upr_z)) +
  facet_grid(choice ~ group) +
  scale_color_viridis_d("Tactic", labels = c("Migrant", "Resident"),
                        begin = 0.25, end = 0.75) +
  theme_bw() +
  labs(x = "Year", y = "Nest survival to day 110")

ggsave("output/trends_z110.tiff",
       width = 14, height = 7, compression = "lzw", dpi = 600)
