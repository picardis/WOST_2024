# Load packages ####

library(tidyverse)
library(rjags)
library(coda)
library(nestR)

# Load model ####

post <- readRDS("output/nest_survival_by_region.rds")

# Load predictors ####

urb_matrix <- readRDS("output/surv_model_predictors.rds") %>%
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

# Detection probability p ####

p_df <- data.frame(t = 1:110,
                   p_lwr = apply(post$p, 1, quantile, 0.025),
                   p_mean = apply(post$p, 1, mean),
                   p_upr = apply(post$p, 1, quantile, 0.975))

ggplot(p_df, aes(x = t, y = p_mean)) +
  geom_ribbon(aes(ymin = p_lwr, ymax = p_upr), fill = "gray80") +
  geom_line() +
  labs(x = "Day of nesting attempt", y = "Detection probability") +
  theme_bw()

# Cumulative survival z ####

z <- post$z
rownames(z) <- rownames(visits)

mig_soflo <- rsf_data %>%
  dplyr::select(att_id, choice, region) %>%
  distinct() %>%
  filter(choice == "migrant" & region == "South Florida") %>%
  pull(att_id)

res_soflo <- rsf_data %>%
  dplyr::select(att_id, choice, region) %>%
  distinct() %>%
  filter(choice == "resident" & region == "South Florida") %>%
  pull(att_id)

res_jax <- rsf_data %>%
  dplyr::select(att_id, choice, region) %>%
  distinct() %>%
  filter(region == "Jacksonville") %>%
  pull(att_id)

mig_se <- rsf_data %>%
  dplyr::select(att_id, choice, region) %>%
  distinct() %>%
  filter(region == "Southeast") %>%
  pull(att_id)

z_mig_soflo <- z[rownames(z) %in% mig_soflo, , , ]
z_res_soflo <- z[rownames(z) %in% res_soflo, , , ]
z_mig_se <- z[rownames(z) %in% mig_se, , , ]
z_res_jax <- z[rownames(z) %in% res_jax, , , ]
z_mig_soflo_mean <- apply(z_mig_soflo, 2, mean)
z_res_soflo_mean <- apply(z_res_soflo, 2, mean)
z_mig_se_mean <- apply(z_mig_se, 2, mean)
z_res_jax_mean <- apply(z_res_jax, 2, mean)

z_preds <- expand.grid(choice = c("Migrant", "Resident"),
                       region = c("Jacksonville", "South Florida", "Southeast"),
                       day = 1:110) %>%
  filter(!(choice == "Migrant" & region == "Jacksonville") &
           !(choice == "Resident" & region == "Southeast")) %>%
  arrange(region, choice, day) %>%
  mutate(z_mean = c(z_res_jax_mean,
                    z_mig_soflo_mean,
                    z_res_soflo_mean,
                    z_mig_se_mean))

ggplot(z_preds, aes(x = day, y = z_mean, color = choice)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ region) +
  scale_color_viridis_d("Choice", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  labs(x = "Days", y = "Cumulative nest survival") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

# Survival to fledging z110 ####

z110_mig_soflo_mean <- mean(z_mig_soflo[, 110, , ])
z110_res_soflo_mean <- mean(z_res_soflo[, 110, , ])
z110_mig_se_mean <- mean(z_mig_se[, 110, , ])
z110_res_jax_mean <- mean(z_res_jax[, 110, , ])

z110_mig_soflo_lwr <- quantile(z_mig_soflo[, 110, , ], 0.025)
z110_res_soflo_lwr <- quantile(z_res_soflo[, 110, , ], 0.025)
z110_mig_se_lwr <- quantile(z_mig_se[, 110, , ], 0.025)
z110_res_jax_lwr <- quantile(z_res_jax[, 110, , ], 0.025)

z110_mig_soflo_upr <- quantile(z_mig_soflo[, 110, , ], 0.975)
z110_res_soflo_upr <- quantile(z_res_soflo[, 110, , ], 0.975)
z110_mig_se_upr <- quantile(z_mig_se[, 110, , ], 0.975)
z110_res_jax_upr <- quantile(z_res_jax[, 110, , ], 0.975)

z110_mig_soflo <- data.frame(choice = "Migrant",
                             region = "South Florida") %>%
  mutate(z110_mean = z110_mig_soflo_mean,
         z110_lwr = z110_mig_soflo_lwr,
         z110_upr = z110_mig_soflo_upr)

z110_res_soflo <- data.frame(choice = "Resident",
                             region = "South Florida") %>%
  mutate(z110_mean = z110_res_soflo_mean,
         z110_lwr = z110_res_soflo_lwr,
         z110_upr = z110_res_soflo_upr)

z110_mig_se <- data.frame(choice = "Migrant",
                             region = "Southeast") %>%
  mutate(z110_mean = z110_mig_se_mean,
         z110_lwr = z110_mig_se_lwr,
         z110_upr = z110_mig_se_upr)

z110_res_jax <- data.frame(choice = "Resident",
                             region = "Jacksonville") %>%
  mutate(z110_mean = z110_res_jax_mean,
         z110_lwr = z110_res_jax_lwr,
         z110_upr = z110_res_jax_upr)

z110_df <- z110_mig_soflo %>%
  bind_rows(z110_res_soflo) %>%
  bind_rows(z110_mig_se) %>%
  bind_rows(z110_res_jax)

ggplot(z110_df, aes(x = choice, y = z110_mean, color = choice)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = z110_lwr, ymax = z110_upr), width = 0.3) +
  facet_wrap(~ region) +
  scale_color_viridis_d("Choice", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  labs(x = " ", y = "Nest survival to day 110") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
