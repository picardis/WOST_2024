# Load packages ####

library(tidyverse)
library(nestR)
library(rjags)
library(coda)

# Load data ####

rsf_data <- readRDS("input/rsf-data.rds") %>%
  mutate(
    se = !soflo & !jax,
    region = case_when(jax ~ "Jacksonville",
                      soflo ~ "South Florida",
                      TRUE ~ "Southeast")
  )
ws_nests <- readRDS("input/ws_nests_all.rds")

# Nest survival ####

# Model data
model_data <- rsf_data %>%
  filter(used == 1) %>%
  as.data.frame()

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

# Average distance to urban
urb_matrix <- model_data %>%
  filter(loc_id == 0) %>%
  group_by(att_id, se, soflo, jax) %>%
  summarize(avg_dist_urb = mean(dist_to_urb_2016, na.rm = TRUE)) %>%
  mutate(se = as.numeric(se),
         soflo = as.numeric(soflo),
         jax = as.numeric(jax)) %>%
  as.data.frame()

saveRDS(urb_matrix, "input/surv_model_predictors.rds")

# MCMC Parameters
mcmc_params <- list(burn_in = 10000, # 10000 (was 1000)
                    n_chain = 3,
                    thin = 10, # 10 or 20 (was 5)
                    n_adapt = 10000, # 10000 (was 1000)
                    n_iter = 10 * 2000) # thin * 2000

# Path to the JAGS file
jags_file <- "jags_urban_regions"

# Starting values for survival status
s1 <- nestR:::initialize_z(ch = visits)

# Define JAGS model
jags <- jags.model(file = jags_file,
                   data = list("nests" = nrow(visits),
                               "days" = ncol(visits),
                               "dist_urb" = urb_matrix$avg_dist_urb,
                               "soflo" = urb_matrix$soflo,
                               "jax" = urb_matrix$jax,
                               "se" = urb_matrix$se,
                               "gps_fixes" = fixes,
                               "y" = visits),
                   inits = list("z" = s1),
                   n.chain = mcmc_params$n_chain,
                   n.adapt = mcmc_params$n_adapt)

#Run the burn-in
update(object = jags, n.iter = mcmc_params$burn_in)

#Generate posterior samples
post <- jags.samples(model = jags,
                     variable.names = c("phi.b0", "phi.b1", "phi.b2",
                                        "phi.b3", "phi.b4", "phi.b5",
                                        "p.b0", "p.b1", "p",
                                        "z"),
                     n.iter = mcmc_params$n_iter,
                     thin = mcmc_params$thin)

saveRDS(post, "input/nest_survival_by_region_2024-06-05.rds")
