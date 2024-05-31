# Code for Picardi et al., "Fitness consequences of alternative behavioral 
# tactics in response to anthropogenic tradeoffs"

# Load packages ----

library(sf)
library(terra)
library(tidyverse)
library(viridis)
library(MASS)
library(nestR)
library(rjags)
library(coda)

# Load data ----

# Load raster of distance from urban
dist_urb_r <- rast("input/distance_to_urban_raster.tiff")
# Load WOST GPS data for breeding attempts
attempts <- readRDS("../partmigenv/output/nesting_attempts_110.rds")
attempts_raw <- readRDS("../partmigenv/output/nesting_attempts.rds")
# Load data on individual migratory behavior
mig_yrs <- readRDS("../partmigenv/output/mig_yrs.rds")
# Load data on nests
ws_nests <- readRDS("../partmigenv/output/ws_nests_all.rds")

# Process data ----

# Fix burst column to match between data frames
mig_yrs$burst <- paste0(mig_yrs$stork, "-", mig_yrs$year)

# Fix att_id column to match between data frames
ws_nests$nests$att_id <- paste0(ws_nests$nests$burst, "-", 
                                ws_nests$nests$loc_id)
ws_nests$visits$att_id <- paste0(ws_nests$visits$burst, "-", 
                                 ws_nests$visits$loc_id)

# Join raw attempts data with migration data and with environmental data
nesting_data <- left_join(attempts_raw, mig_yrs, by = "burst") %>% 
  dplyr::select(att_id, burst, stork, choice, strategy, flexible, 
                day_of_att, date, datetime, long, lat, loc_id)

# Keep only birds for which I know the migratory choice and mutate date column
nesting_data <- nesting_data %>% 
  filter(!is.na(choice)) %>% 
  mutate(date_ymd = date)

# Join location of the nest
nest_coords <- ws_nests$nests %>% dplyr::select(att_id, 
                                                nest_long = long, 
                                                nest_lat = lat)
nesting_data <- nesting_data %>% left_join(nest_coords, by = "att_id")

# All nests corresponding to the attempts in nesting_data
nest_loc_data <- nesting_data %>% 
  distinct(att_id, choice, nest_long, nest_lat) %>% 
  arrange(desc(nest_lat))

# Load coastline shapefile
coastline <- read_sf(dsn = "input/fl_bnd/fl_bnd.shp", layer = "fl_bnd") %>% 
  filter(SHAPE_FID < 26) %>% 
  dplyr::select(1)

coastline <- st_transform(coastline, 32617)

# Retain used foraging points within FL only
nesting_data <- nesting_data %>% 
  filter(lat < 30.828 & loc_id == 0) 

# Fit distribution to empirical ----

# Fit an exponential distribution to the observed steps
fitexp <- nesting_data %>%
  mutate(dist_to_nest = geosphere::distGeo(p1 = cbind(long, lat), p2 = cbind(nest_long, nest_lat))/1000) %>% 
  pull(dist_to_nest) %>% 
  fitdistr("exponential")

# Plot probability density function
fit_pdf <- list()
fit_pdf$x <- seq(0, 350, by = 0.001)
fit_pdf$y <- dexp(fit_pdf$x, rate = fitexp$estimate)
plot(fit_pdf$x, fit_pdf$y, type = "l")

# Where does the value of the PDF get to 1e-4?
max_dist <- fit_pdf$x[min(which(fit_pdf$y <= 1e-3))]

# For each individual attempt, generate 10x as many available
# points as the used with decreasing probability as you get
# farther from the nest (from a fitted distribution fitted to the empirical)

ids <- unique(nesting_data$att_id)

system.time({
  urban_rsf_data <- data.frame()
  set.seed(1)
  tracker <- 0
  for (i in ids) {
    tracker = tracker + 1
    print(paste(tracker, "of", length(ids)))
    # Get the data for this attempt
    dat <- nesting_data %>% filter(att_id == i) %>% 
      mutate(used = 1) %>% 
      dplyr::select(att_id, burst, stork, choice, date_ymd, datetime, 
                    day_of_att, long, lat, loc_id,
                    nest_long, nest_lat, used) %>% 
      mutate(dist_to_nest = geosphere::distGeo(p1 = cbind(long, lat), p2 = cbind(nest_long, nest_lat))) 
    # Make nest spatial
    nest <- dat %>% 
      dplyr::select(nest_long, nest_lat) %>% 
      unique() %>% 
      st_as_sf(coords = c("nest_long", "nest_lat"), crs = 4326)
    # Convert to UTM
    nest_utm <- st_transform(nest, 32617)
    # Create buffer around nest
    buff <- st_buffer(x = nest_utm, dist = units::set_units(max_dist, km))
    # Clip it to the coastline
    buff_clip <- st_intersection(buff, coastline)
    
    if (length(is.na(buff_clip$geometry)) != 0) {
      #Get the coordinates of the nest
      coords_nest <- data.frame("x" = st_coordinates(nest_utm)[, 1], 
                                "y" = st_coordinates(nest_utm)[, 2])
      n_pts <- nrow(dat)*10
      #Vector of distances from nest (in meters)
      dists <- rexp(n = n_pts * 3, rate = fitexp$estimate) * 1000
      #Vector of angles in radians
      angs <- runif(n = n_pts * 3, min = 0, max = 2*pi)
      #Calculate delta x and delta y
      delta_x <- dists * cos(pi/2 - angs)
      delta_y <- dists * cos(angs)
      #Calculate new x and y coords
      rand <- data.frame(X = coords_nest$x + delta_x, Y = coords_nest$y + delta_y)
      #Make rand spatial
      rand <- st_as_sf(rand, coords = c("X", "Y"))
      st_crs(rand) <- 32617
      #Clip just those within the coastline
      rand <- rand %>% 
        st_as_sf() %>% 
        st_intersection(coastline) %>% 
        #And keep only the desired amount
        sample_n(size = n_pts)
      #Add the other attributes
      rand <- rand %>% 
        st_transform(4326) %>% 
        st_coordinates() %>% 
        as.data.frame() %>% 
        mutate(att_id = i,
               burst = unique(dat$burst),
               stork = unique(dat$stork),
               choice = unique(dat$choice),
               date_ymd = NA,
               datetime = NA,
               day_of_att = NA, 
               long = X,
               lat = Y,
               loc_id = 0,
               nest_long = unique(dat$nest_long),
               nest_lat = unique(dat$nest_lat),
               used = 0) %>% 
        dplyr::select(att_id, burst, stork, choice, date_ymd, datetime,
                      day_of_att, long, lat, loc_id, nest_long, nest_lat, used) %>% 
        mutate(dist_to_nest = geosphere::distGeo(p1 = cbind(long, lat), p2 = cbind(nest_long, nest_lat))) 
      
      dat <- rbind(dat, rand)
      urban_rsf_data <- rbind(urban_rsf_data, dat)
      
    }
  }
  beepr::beep(3)
})

# Check resulting distribution
urban_rsf_data %>% 
  filter(used == 0) %>% 
  mutate(dist_to_nest_km = dist_to_nest/1000) %>% 
  pull(dist_to_nest_km) %>% 
  hist(freq = FALSE, breaks = 40)
lines(x = fit_pdf$x, y = fit_pdf$y, col = "red")

# Map results (takes a while)
# urban_rsf_data %>% 
#   filter(used == 0) %>% 
#   ggplot(aes(x = long, y = lat, group = att_id)) +
#   geom_point() +
#   facet_wrap(~att_id) +
#   geom_point(aes(x = nest_long, y = nest_lat), col = "red")

# Check that everyone was paired (same number of att_id for used and available):
# some birds might have nested outside of FL and we are clipping the buffers to the FL coastline
urban_rsf_data %>% 
  as.data.frame() %>% 
  filter(used==1) %>% 
  group_by(att_id) %>% 
  tally() %>% 
  arrange(n)

urban_rsf_data %>% 
  as.data.frame() %>% 
  filter(used==0) %>% 
  group_by(att_id) %>% 
  tally() %>% 
  arrange(n)

table(urban_rsf_data$used)

# Prep for model fitting ----

# Reproject into UTMs
urban_rsf_data <- urban_rsf_data %>% 
  rename(x = long, y = lat) %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326) %>% 
  st_transform(crs = 32617)

# Calculate distance of each point to urban
urban_rsf_data$dist_to_urban <- terra::extract(dist_urb_r, 
                                        st_coordinates(urban_rsf_data))[, 1]

# Transform to data frame
urban_rsf_df <- as.data.frame(urban_rsf_data) %>% 
  cbind(st_coordinates(urban_rsf_data)) %>% 
  dplyr::select(-geometry) %>% 
  rename(utm_x = X, utm_y = Y)

# Scale and center variables
mean_dist <- mean(urban_rsf_df$dist_to_urban)
sd_dist <- sd(urban_rsf_df$dist_to_urban)

urban_rsf_scaled <- urban_rsf_df %>% 
  mutate(dist_to_urban = (dist_to_urban - mean_dist)/sd_dist,
         weight = case_when(
           used == 1 ~ 1,
           used == 0 ~ 10000
         ))

# Run model ----

mod <- glm(formula = used ~ dist_to_urban*choice + I(dist_to_urban^2)*choice,
              data = urban_rsf_scaled, family = binomial, weights = weight)

# Model predictions ----

urban_rsf_df %>% 
  filter(used == 1) %>% 
  summarize(max(dist_to_urban))

pred_df <- data.frame(
  choice = rep(c("migrant", "resident"), each = 100),
  dist_to_urban_orig = seq(0, 45000, length.out = 100)) %>% 
  mutate(dist_to_urban = (dist_to_urban_orig - mean_dist)/sd_dist, # unscale
         dist_to_urban_km = dist_to_urban_orig/1000)

urban_rsf_df %>% 
  filter(used == 1) %>% 
  summarize(mean(dist_to_urban))

refd <- data.frame(
  choice = c("migrant", "resident"),
  dist_to_urban_orig = 45000/2) %>% 
  mutate(dist_to_urban = (dist_to_urban_orig - mean_dist)/sd_dist, # unscale
         dist_to_urban_km = dist_to_urban_orig/1000)
  
pred <- predict(mod, newdata = pred_df, type = "response", se.fit = TRUE)
ref_pred <- predict(mod, newdata = refd, type = "response", se.fit = TRUE)

rss <- c(pred$fit[1:100]/ref_pred$fit[1], pred$fit[101:200]/ref_pred$fit[2])
pred_df$rss <- rss

lwr_pred <- pred$fit - 1.96 * pred$se.fit
lwr_ref <- ref_pred$fit - 1.96 * ref_pred$se.fit
lwr <- c(lwr_pred[1:100]/lwr_ref[1], lwr_pred[101:200]/lwr_ref[2])
pred_df$lwr <- lwr

upr_pred <- pred$fit + 1.96 * pred$se.fit
upr_ref <- ref_pred$fit + 1.96 * ref_pred$se.fit
upr <- c(upr_pred[1:100]/upr_ref[1], upr_pred[101:200]/upr_ref[2])
pred_df$upr <- upr

n_used <- urban_rsf_df %>% 
  filter(used == 1) %>% 
  nrow()

n_avail <- urban_rsf_df %>% 
  filter(used == 0) %>% 
  nrow()

mig_used <- urban_rsf_df %>% 
  filter(choice == "migrant" & used == 1) %>% 
  mutate(dist_to_urban_km = dist_to_urban/1000)

res_used <- urban_rsf_df %>% 
  filter(choice == "resident" & used == 1) %>% 
  mutate(dist_to_urban_km = dist_to_urban/1000)

ggplot(pred_df, aes(x = dist_to_urban_km, y = rss, col = choice)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = choice), alpha = 0.2) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 1, color = "gray20", linetype = 2) +
  geom_rug(aes(x = dist_to_urban_km, y = NULL), data = mig_used, sides = "b") +
  geom_rug(aes(x = dist_to_urban_km, y = NULL), data = res_used, sides = "t") +
  scale_x_continuous(name = "Distance to urban (km)") +
  scale_y_continuous(name = "Relative selection strength") +
  scale_color_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  theme_bw()

ggsave("output/2023-11-08_RSF.tiff", 
       width = 180, height = 90, units = "mm", compression = "lzw", dpi = 600)

# Nest survival ----

# Calculate distance from urban
dist_urb_r_ll <- terra::project(dist_urb_r, "epsg:4326")

nesting_data <- st_as_sf(nesting_data, 
                         coords = c("long", "lat"),
                         crs = 4326)

nesting_data$dist_to_urban <- terra::extract(dist_urb_r_ll, 
                                      st_coordinates(nesting_data))[, 1]

nesting_data2 <- as.data.frame(nesting_data) %>% 
  cbind(st_coordinates(nesting_data)) %>% 
  dplyr::select(-geometry)

# Select all birds that nest in FL (for which I have distance from urban)
keepers <- nesting_data2 %>% 
  group_by(att_id) %>% 
  summarize(avg_dist_urb = mean(dist_to_urban, na.rm = TRUE)) %>% 
  filter(!is.nan(avg_dist_urb)) %>% 
  # filter only those that are in the RSF data
  filter(att_id %in% unique(urban_rsf_scaled$att_id)) %>% 
  pull(att_id) 

# Model data
model_data <- nesting_data2 %>% 
  as.data.frame() %>% 
  filter(att_id %in% keepers)

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

fixes <- ws_attempts$fixes[keepers,] 
visits <- ws_attempts$visits[keepers,] 

# Distance to urban
urb_matrix <- model_data %>% 
  filter(loc_id == 0) %>% 
  group_by(att_id) %>% 
  summarize(avg_dist_urb = mean(dist_to_urban, na.rm = TRUE)) %>% 
  filter(!is.nan(avg_dist_urb)) %>% 
  as.data.frame()

# MCMC Parameters
mcmc_params <- list(burn_in = 1000,
                    n_chain = 3,
                    thin = 5,
                    n_adapt = 1000,
                    n_iter = 5 * 2000) 

# Path to the JAGS file
jags_file <- "jags_urban"

# Starting values for survival status
s1 <- nestR:::initialize_z(ch = visits)

# Define JAGS model
jags <- jags.model(file = jags_file,
                   data = list("nests" = nrow(visits),
                               "days" = ncol(visits),
                               "dist_urb" = urb_matrix$avg_dist_urb,
                               "gps_fixes" = fixes,
                               "y" = visits),
                   inits = list("z" = s1),
                   n.chain = mcmc_params$n_chain,
                   n.adapt = mcmc_params$n_adapt)

#Run the burn-in
update(object = jags, n.iter = mcmc_params$burn_in)

#Generate posterior samples
post <- jags.samples(model = jags,
                     variable.names = c("phi.b0", "phi.b1",
                                        "p.b0", "p.b1", "p",
                                        "z"),
                     n.iter = mcmc_params$n_iter,
                     thin = mcmc_params$thin)

saveRDS(post, "output/urban_surv_2023-09-01.rds")

# Get parameter values from posterior----

post <- readRDS("output/urban_surv_2023-08-31.rds")

# Initialize parameter values
param <- list()
param$mean <- list()
param$lwr <- list()
param$upr <- list()
# Get parameter values from posterior distribution
param$mean <- lapply(post, mean)
param$lwr <- lapply(post, quantile, 0.025)
param$upr <- lapply(post, quantile, 0.975)

# Diagnostics----

# This informs us of how much thinning should we use (done for each beta parameter)
autocorr.plot(as.mcmc.list(post$phi.b0), lag.max = 20)
autocorr.plot(as.mcmc.list(post$phi.b1), lag.max = 20)

traceplot(as.mcmc.list(post$phi.b0))
traceplot(as.mcmc.list(post$phi.b1))

densplot(as.mcmc.list(post$phi.b0))
abline(v = 0, col = "red")
densplot(as.mcmc.list(post$phi.b1))
abline(v = 0, col = "red")

# Prep data for plotting----

# Based directly on phi
phi.df.param <- data.frame(phi.b0 = as.vector(post$phi.b0),
                           phi.b1 = as.vector(post$phi.b1)) %>% 
  list() %>% 
  rep(100) %>% 
  bind_rows()

phi.df.cov <- data.frame(dist_urb = seq(0, 65, length.out = 100)) %>% 
  list() %>% 
  rep(6000) %>% 
  bind_rows() %>% 
  arrange(dist_urb)

phi.df.all <- bind_cols(phi.df.param, phi.df.cov) %>% 
  mutate(phi = plogis((phi.b0 + (phi.b1 * dist_urb)))) %>% 
  mutate(phi_lwr = quantile(phi, 0.025),
         phi_upr = quantile(phi, 0.975))

# Based on the parameters
phi.df <- data.frame(dist_urb = seq(0, 45000, length.out = 500)) %>% 
  mutate(phi.lwr = plogis((param$lwr$phi.b0 + (param$lwr$phi.b1 * dist_urb))),
         phi.mean = plogis((param$mean$phi.b0 + (param$mean$phi.b1 * dist_urb))),
         phi.upr = plogis((param$upr$phi.b0 + (param$upr$phi.b1 * dist_urb))),
         dist_urb_km = dist_urb/1000
  )

# Based on acutally calculated p
p.df <- data.frame(t = 1:110,
                   p.lwr = apply(post$p, 1, quantile, 0.025),
                   p.mean = apply(post$p, 1, mean),
                   p.upr = apply(post$p, 1, quantile, 0.975))

#Plots----
plot(p.mean ~ t, data = p.df, type = "l", lwd = 1, ylim = c(0.15, 0.75))
lines(p.lwr~t, data = p.df, lty = 2)
lines(p.upr~t, data = p.df, lty = 2)

raw <- model_data %>% 
  mutate(dist_to_urban_km = dist_to_urban/1000)

ggplot(phi.df, aes(x = dist_urb_km, y = phi.mean)) +
  geom_ribbon(aes(ymin = phi.lwr, ymax = phi.upr, alpha = 0), fill = "gray70") +
  geom_line(linewidth = 1.2) +
  geom_rug(aes(x = dist_to_urban_km, y = NULL), data = raw, sides = "t") +
  scale_x_continuous(name = "Distance to urban (km)") +
  scale_y_continuous(name = "Daily nest survival") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("output/2023-09-01_Surv.tiff", 
       width = 180, height = 90, units = "mm", compression = "lzw", dpi = 600)

# Plot survival of mig vs res
z <- post$z
rownames(z) <- rownames(visits)
mig <- nesting_data2 %>% 
  dplyr::select(att_id, choice) %>% 
  distinct() %>% 
  filter(choice == "migrant") %>% 
  pull(att_id)
res <- nesting_data2 %>% 
  dplyr::select(att_id, choice) %>% 
  distinct() %>% 
  filter(choice == "resident") %>% 
  pull(att_id)
z_mig <- z[rownames(z) %in% mig, , , ]
z_res <- z[rownames(z) %in% res, , , ]
z_mig_mean <- apply(z_mig, 2, mean)
z_res_mean <- apply(z_res, 2, mean)

z_migres <- data.frame(choice = rep(c("Migrant", "Residents"), each = 110),
                       day = c(1:110, 1:110),  
                       z = c(z_mig_mean, z_res_mean))

ggplot(z_migres, aes(x = day, y = z, color = choice)) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis_d("Choice", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  labs(x = "Days", y = "Cumulative nest survival") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

ggsave("output/2023-09-01_Surv-curves.tiff", 
       width = 80, height = 60, units = "mm", compression = "lzw", dpi = 600)

# Plots ----

# Color palette
viridis(begin = 0.25, end = 0.75, n = 2)

# Plot 1: map of urban
urban_r <- readRDS("input/urban_rasterized.rds") %>% 
  rast() %>% 
  project("epsg:32617")

ext(urban_r) <- ext(coastline)
st_crs(urban_r)
st_crs(coastline)

urban_r_df <- as.data.frame(urban_r, xy = TRUE)

p1 <- ggplot() +
  geom_raster(data = urban_r_df, aes(x = x, y = y, fill = layer)) +
  geom_sf(data = coastline, fill = NA) +
  coord_sf() +
  scale_fill_gradient(low = "#000000", high = "#000000", na.value = NA) +
  theme_void() +
  theme(legend.position = "none", panel.grid.major = element_line(colour = "transparent"),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Urban areas")

# Plot 2: map of nests (by mig)

nest_locs <- urban_rsf_data %>%
  as.data.frame() %>% 
  dplyr::select(nest_long, nest_lat, choice) %>% 
  distinct() %>%  
  st_as_sf(coords = c("nest_long", "nest_lat"), crs = 4326) %>% 
  st_transform(crs = 32617)

p2 <- ggplot() +
  geom_sf(data = coastline, fill = NA) +
  geom_sf(data = nest_locs, mapping = aes(col = choice), 
          size = 3, alpha = 0.4, shape = 17) +
  coord_sf() +
  scale_color_viridis_d(begin = 0.25, end = 0.75) +
  theme_void() +
  theme(legend.position = "none", panel.grid.major = element_line(colour = "transparent"),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Nest locations")

# Plot 3: map of foraging sites

forag <- urban_rsf_data %>%
  as.data.frame() %>% 
  cbind(st_coordinates(urban_rsf_data)) %>% 
  filter(used == 1) 

p3 <- ggplot() +
  geom_sf(data = coastline, fill = NA) +
  geom_point(data = forag, aes(x = X, y = Y, color = choice), alpha = 0.1) +
  scale_color_viridis_d("Choice", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  coord_sf() +
  theme_void() +
  theme(legend.position = "none", panel.grid.major = element_line(colour = "transparent"),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Foraging locations")

ggpubr::ggarrange(p3, p1, p2, nrow = 1, ncol = 3, 
                  labels = c("(A)", "(B)", "(C)"))

ggsave("output/2023-09-01_Maps.tiff", 
       width = 180, height = 60, units = "mm", dpi = 600, compression = "lzw")

# Additional exploration ----

# Fit curve to empirical distribution of distances between nest and foraging sites
fitexp <- urban_rsf_df %>%
  filter(used == 1) %>% 
  mutate(dist_to_nest = geosphere::distGeo(p1 = cbind(long, lat), p2 = cbind(nest_long, nest_lat))/1000) %>% 
  pull(dist_to_nest) %>% 
  MASS::fitdistr("exponential")

fitgam <- urban_rsf_df %>%
  filter(used == 1) %>% 
  mutate(dist_to_nest = geosphere::distGeo(p1 = cbind(long, lat), p2 = cbind(nest_long, nest_lat))/1000) %>% 
  pull(dist_to_nest) %>% 
  MASS::fitdistr("gamma")

dists <- urban_rsf_df %>%
  filter(used == 1) %>% 
  mutate(dist_to_nest = geosphere::distGeo(p1 = cbind(long, lat), p2 = cbind(nest_long, nest_lat))/1000) %>% 
  pull(dist_to_nest)

ks.test(dists, "pexp", fitexp$estimate)
ks.test(dists, "pgamma", fitgam$estimate)

EnvStats::qqPlot(dists, distribution = "exp", param.list = list(rate = 0.0478787932))
EnvStats::qqPlot(dists, distribution = "gamma", param.list = list(shape = 0.5452323980, scale = 1/0.0261050488))

empir <- urban_rsf_df %>%
  filter(used == 1) %>% 
  mutate(dist_to_nest = geosphere::distGeo(p1 = cbind(long, lat), p2 = cbind(nest_long, nest_lat))/1000) %>% 
  ggplot(aes(x = dist_to_nest)) +
  geom_density() +
  #  stat_function(fun = dgamma, n = 100, args = list(shape = 0.55, rate = 0.03), col = "red") +
  stat_function(fun = dexp, n = 100, args = list(rate = 0.05), col = "tomato", lwd = 1.5, alpha = 0.7) +
  theme_minimal() +
  scale_y_continuous("Density") +
  scale_x_continuous("Foraging trip distance (km)") +
  ggtitle("Observed")

simul <- urban_rsf_df %>%
  filter(used == 0) %>% 
  mutate(dist_to_nest = geosphere::distGeo(p1 = cbind(long, lat), p2 = cbind(nest_long, nest_lat))/1000) %>% 
  ggplot(aes(x = dist_to_nest)) +
  geom_density() +
  stat_function(fun = dexp, n = 100, args = list(rate = 0.05), col = "tomato", lwd = 1.5, alpha = 0.7) +
  theme_minimal() +
  scale_y_continuous("Density") +
  scale_x_continuous("Foraging trip distance (km)") +
  ggtitle("Simulated")

ggpubr::ggarrange(empir, simul, nrow = 1, ncol = 2, 
                  labels = c("(A)", "(B)"))

#ggsave("output/trip_distance_empir_vs_simul.tiff", width = 8, height = 4, dpi = 500, compression = "lzw")

# Mean distance of foraging locations from nest for migrants and residents
urban_rsf_df %>%
  filter(used == 1) %>% 
  mutate(dist_to_nest = geosphere::distGeo(p1 = cbind(long, lat), p2 = cbind(nest_long, nest_lat))) %>% 
  group_by(choice) %>% 
  summarize(mean(dist_to_nest)/1000)

# Mean distance of foraging locations from urban for migrants and residents
urban_rsf_df %>%
  filter(used == 1) %>% 
  group_by(choice) %>% 
  summarize(mean(dist_to_urban/1000))

# Plot of distance of foraging locations to urban for migrants and residents (used vs. available)

mig_labs <- c(migrant = "Migrant",
              resident = "Resident")

avail <- urban_rsf_df %>% 
  filter(used == 0) %>% 
  mutate(dist_to_urban = dist_to_urban/1000)

urban_rsf_df %>% 
  filter(used == 1) %>% # used points
  mutate(dist_to_urban = dist_to_urban/1000) %>% 
  ggplot(aes(x = dist_to_urban, fill = choice, group = choice)) +
  geom_density(data = avail, fill = "gray80") +
  geom_density(alpha = 0.5) +
  facet_wrap(~choice, labeller = labeller(choice = mig_labs)) +
  scale_fill_viridis_d(begin = 0.25, end = 0.75) +
  theme_minimal() +
  scale_y_continuous("Density") +
  scale_x_continuous("Distance to urban development (km)") +
  theme(legend.position = "none")

#ggsave("output/dist_to_urb_used_vs_avail.tiff", width = 8, height = 4, dpi = 500, compression = "lzw")

used_plot <- urban_rsf_df %>% 
  filter(used == 1) %>% # used points
  mutate(dist_to_urban = dist_to_urban/1000) %>% 
  ggplot(aes(x = dist_to_urban, fill = choice, group = choice)) +
  geom_density() +
  facet_wrap(~choice, labeller = labeller(choice = mig_labs)) +
  scale_fill_viridis_d(begin = 0.25, end = 0.75) +
  theme_minimal() +
  scale_y_continuous("Density", limits = c(0, 0.7)) +
  scale_x_continuous("Distance to urban development (km)") +
  theme(legend.position = "none") +
  ggtitle("Used")

avail_plot <- urban_rsf_df %>% 
  filter(used == 0) %>% # available points
  mutate(dist_to_urban = dist_to_urban/1000) %>% 
  ggplot(aes(x = dist_to_urban, fill = choice, group = choice)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~choice, labeller = labeller(choice = mig_labs)) +
  scale_fill_viridis_d(begin = 0.25, end = 0.75) +
  theme_minimal() +
  scale_y_continuous("Density", limits = c(0, 0.7)) +
  scale_x_continuous("Distance to urban development (km)") +
  theme(legend.position = "none") +
  ggtitle("Available")

ggpubr::ggarrange(used_plot, avail_plot, nrow = 2, ncol = 1, 
                  labels = c("(A)", "(B)"))

#ggsave("output/dist_to_urb_used_vs_avail.tiff", width = 8, height = 4, dpi = 500, compression = "lzw")

# Plot availability by choice
urban_avail <- urban_rsf_df %>% 
  filter(used == 0)

ggplot(urban_avail, aes(x = dist_to_urban/1000, fill = choice, group = choice)) +
  geom_density() +
  facet_wrap(~ choice) +
  scale_x_continuous("Distance to urban development (km)") +
  scale_y_continuous("Density") +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  ggtitle("Distance of available points to urban development") +
  theme_minimal()

# Compare use and availability in terms of distance to urban at the attempt level,
# rank based on discrepancy (larger discrepancies first), add line for distance of nest to urban

us <- urban_rsf_df %>% 
  group_by(att_id) %>% 
  filter(used == 1) %>% 
  summarize(mean_dist_used = mean(dist_to_urban)/1000)

av <- urban_rsf_df %>% 
  group_by(att_id) %>% 
  filter(used == 0) %>% 
  summarize(mean_dist_avai = mean(dist_to_urban)/1000)

discrep <- left_join(us, av, by = "att_id") %>% 
  mutate(discrep = abs(mean_dist_avai - mean_dist_used)) %>% 
  arrange(desc(discrep)) %>% 
  pull(att_id)

urban_rsf_df$att_id_reord <- factor(urban_rsf_df$att_id, levels = discrep)

nests <- urban_rsf_df %>% 
  dplyr::select(att_id_reord, nest_long, nest_lat) %>% 
  distinct() 

coordinates(nests) <- nests[,c("nest_long", "nest_lat")]
proj4string(nests) <- CRS("+init=epsg:4326")

nests$nest_to_urban <- raster::extract(dist_urb_r, nests)

nests <- as.data.frame(nests)

ind_plot <- urban_rsf_df %>% 
  filter(used == 1) %>% 
  ggplot(aes(x = dist_to_urban/1000, fill = choice, group = choice)) +
  geom_density(data = urban_rsf_df[urban_rsf_df$used==0,], alpha = 0.5, fill = "gray50") +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = nest_to_urban/1000), data = nests) + 
  facet_wrap(~att_id_reord) +
  scale_fill_viridis_d(begin = 0.25, end = 0.75) +
  theme_void() +
  scale_y_continuous("Density", limits = c(0, 0.5)) +
  scale_x_continuous("Distance to urban development (km)") +
  theme(legend.position = "none")

#ggsave("output/individual_dist_to_urban_use_vs_avail.tiff", plot = ind_plot, width = 15, height = 15, dpi = 400, compression = "lzw")

# Plot distance of nest to urban for migrants and residents

nests %>% 
  mutate(att_id = att_id_reord) %>% 
  left_join(urban_rsf_df[,c("att_id", "choice")], by = "att_id") %>% 
  mutate(nest_to_urban = nest_to_urban/1000) %>% 
  dplyr::select(att_id, nest_to_urban, choice) %>% 
  distinct() %>% 
  ggplot(aes(x = nest_to_urban, fill = choice)) +
  geom_density() +
  facet_wrap(~ choice, labeller = labeller(choice = mig_labs)) +
  scale_x_continuous(name = "Distance of nests to urban development (km)") +
  scale_y_continuous(name = "Density") +
  scale_fill_viridis_d("Choice", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  theme_minimal() +
  theme(legend.position = "none") 

#ggsave("output/dist_nest_to_urban.tiff", width = 8, height = 4, dpi = 500, compression = "lzw")
