# Code for Picardi et al., "Fitness consequences of anthropogenic subsidies
# for a partially migratory wading bird"

# Load packages ####

library(sf)
library(terra)
library(tidyverse)
library(MASS)
library(nestR)
library(survival)

# Load data ####

# Load WOST GPS data for breeding attempts
attempts <- readRDS("input/nesting_attempts_110.rds")
attempts_raw <- readRDS("input/nesting_attempts.rds")
# Load data on individual migratory behavior
mig_yrs <- readRDS("input/mig_yrs.rds")
# Load data on nests
ws_nests <- readRDS("input/ws_nests_all.rds")

# Load GIS layers
dev_2001_se <- rast("processed_layers/development_1km_2001_Landfire_UTM17N.tiff")
dev_2016_se <- rast("processed_layers/development_1km_2016_Landfire_UTM17N.tiff")
dist_to_urb_2001 <- rast("processed_layers/dist-to-development_1km_2001_Landfire_UTM17N.tiff")
dist_to_urb_2016 <- rast("processed_layers/dist-to-development_1km_2016_Landfire_UTM17N.tiff")
dist_to_urb_2010 <- rast("processed_layers/dist-to-development_1km_2010_TIGER_UTM17N.tiff")
dev_2001_agg <- rast("processed_layers/development_aggregated_300m_2001_Landfire_UTM17N.tiff")
dev_2016_agg <- rast("processed_layers/development_aggregated_300m_2016_Landfire_UTM17N.tiff")
dist_to_urb_agg_2001 <- rast("processed_layers/dist-to-development_aggregated_300m_2001_Landfire_UTM17N.tiff")
dist_to_urb_agg_2016 <- rast("processed_layers/dist-to-development_aggregated_300m_2016_Landfire_UTM17N.tiff")
tiger <- st_read("processed_layers/urban_polygons_TIGER_UTM17N.shp")
wmd <- st_read("processed_layers/fl-water-mgmt-districts_UTM17N.shp")
soflo <- st_read("processed_layers/south-florida-buffer_UTM17N.shp")
se <- st_read("processed_layers/southeast-boundary_UTM17N.shp")

# Prepare data ####

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

# Intersect nests with water mgmt districts
nest_sf <- st_as_sf(nest_coords,
                    coords = c("nest_long", "nest_lat"),
                    crs = 4326) %>%
  st_transform(32617)

nest_wmd <- st_intersects(nest_sf, wmd, sparse = FALSE)

nest_wmd <- apply(nest_wmd, 1, function(x) {
  pol <- which(x)
  if(length(pol) == 0) {return(NA)} else {return(pol)}
})

nest_coords$wmd <- nest_wmd

# Intersect with South Florida buffer
nest_buff <- st_intersects(nest_sf, soflo, sparse = FALSE)

nest_coords$soflo <- as.logical(nest_buff)

nesting_data <- nesting_data %>%
  left_join(nest_coords, by = "att_id")

jax_nests <- nesting_data %>%
  filter(choice == "resident" & nest_lat > 30) %>%
  pull(att_id) %>%
  unique()

nesting_data <- nesting_data %>%
  mutate(fl = case_when(
      # This individual's nest falls on the edge of the WMD so even though
      # it's in FL it gets an NA for wmd
      !is.na(wmd) | burst == "429190-2009" ~ TRUE,
      is.na(wmd) ~ FALSE),
    jax = att_id %in% jax_nests)

# All nests corresponding to the attempts in nesting_data
nest_loc_data <- nesting_data %>%
  distinct(att_id, choice, nest_long, nest_lat) %>%
  arrange(desc(nest_lat))

# Retain only points off the nest (foraging points)
nesting_data <- nesting_data %>%
  filter(loc_id == 0)

# Sample available points ####

# Fit an exponential distribution to the observed steps
fitexp <- nesting_data %>%
  mutate(dist_to_nest = geosphere::distGeo(p1 = cbind(long, lat),
                                           p2 = cbind(nest_long, nest_lat))/1000) %>%
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
                    nest_long, nest_lat, used, wmd, soflo, fl, jax) %>%
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
    buff_clip <- st_intersection(buff, se)

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
        st_intersection(se) %>%
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
               used = 0,
               soflo = unique(dat$soflo),
               fl = unique(dat$fl),
               wmd = unique(dat$wmd),
               jax = unique(dat$jax)) %>%
        dplyr::select(att_id, burst, stork, choice, date_ymd, datetime,
                      day_of_att, long, lat, loc_id, nest_long, nest_lat, used, wmd, soflo, fl, jax) %>%
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

# Intersect covariates ####

# Project WOST data to UTM
urban_rsf_data <- urban_rsf_data %>%
  rename(x = long, y = lat) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  st_transform(crs = 32617)

# Intersect with urban rasters
urban_rsf_data$urb_2016 <- terra::extract(dev_2016_se,
                                          st_coordinates(urban_rsf_data))[, 1]
urban_rsf_data$urb_2001 <- terra::extract(dev_2001_se,
                                          st_coordinates(urban_rsf_data))[, 1]

urban_rsf_data$dist_to_urb_2016 <- terra::extract(dist_to_urb_2016,
                                                  st_coordinates(urban_rsf_data))[, 1]
urban_rsf_data$dist_to_urb_2001 <- terra::extract(dist_to_urb_2001,
                                                  st_coordinates(urban_rsf_data))[, 1]
urban_rsf_data$dist_to_urb_2010 <- terra::extract(dist_to_urb_2010,
                                                  st_coordinates(urban_rsf_data))[, 1]

urban_rsf_data$urb_agg_2016 <- terra::extract(dev_2016_agg,
                                              st_coordinates(urban_rsf_data))[, 1]
urban_rsf_data$urb_agg_2001 <- terra::extract(dev_2001_agg,
                                              st_coordinates(urban_rsf_data))[, 1]

urban_rsf_data$urb_tiger_2010 <- as.numeric(st_intersects(urban_rsf_data, tiger))

urban_rsf_data$dist_to_urb_agg_2001 <- terra::extract(dist_to_urb_agg_2001,
                                                      st_coordinates(urban_rsf_data))[, 1]
urban_rsf_data$dist_to_urb_agg_2016 <- terra::extract(dist_to_urb_agg_2016,
                                                      st_coordinates(urban_rsf_data))[, 1]

# Transform to data frame
urban_rsf_df <- as.data.frame(urban_rsf_data) %>%
  cbind(st_coordinates(urban_rsf_data)) %>%
  dplyr::select(-geometry)

urban_rsf_df <- urban_rsf_df %>%
  mutate(weight = case_when(
    used == 1 ~ 1,
    used == 0 ~ 10000
  ),
  urb_tiger_2010 = ifelse(is.na(urb_tiger_2010), 0, urb_tiger_2010))

# Scale and center variables
mean_2016 <- mean(urban_rsf_df$dist_to_urb_2016, na.rm = TRUE)
sd_2016 <- sd(urban_rsf_df$dist_to_urb_2016, na.rm = TRUE)
mean_2001 <- mean(urban_rsf_df$dist_to_urb_2001, na.rm = TRUE)
sd_2001 <- sd(urban_rsf_df$dist_to_urb_2001, na.rm = TRUE)
mean_2010 <- mean(urban_rsf_df$dist_to_urb_2010, na.rm = TRUE)
sd_2010 <- sd(urban_rsf_df$dist_to_urb_2010, na.rm = TRUE)
mean_2016_agg <- mean(urban_rsf_df$dist_to_urb_agg_2016, na.rm = TRUE)
sd_2016_agg <- sd(urban_rsf_df$dist_to_urb_agg_2016, na.rm = TRUE)
mean_2001_agg <- mean(urban_rsf_df$dist_to_urb_agg_2001, na.rm = TRUE)
sd_2001_agg <- sd(urban_rsf_df$dist_to_urb_agg_2001, na.rm = TRUE)

urban_rsf_df <- urban_rsf_df %>%
  mutate(dist_to_urb_2016_sc = (dist_to_urb_2016 - mean_2016)/sd_2016,
         dist_to_urb_2001_sc = (dist_to_urb_2001 - mean_2001)/sd_2001,
         dist_to_urb_2010_sc = (dist_to_urb_2010 - mean_2010)/sd_2010,
         dist_to_urb_agg_2016_sc = (dist_to_urb_agg_2016 - mean_2016_agg)/sd_2016_agg,
         dist_to_urb_agg_2001_sc = (dist_to_urb_agg_2001 - mean_2001_agg)/sd_2001_agg,
         mean_dist_to_urb_2001 = mean_2001,
         mean_dist_to_urb_2016 = mean_2016,
         mean_dist_to_urb_2010 = mean_2010,
         mean_dist_to_urb_agg_2001 = mean_2001_agg,
         mean_dist_to_urb_agg_2016 = mean_2016_agg,
         sd_dist_to_urb_2001 = sd_2001,
         sd_dist_to_urb_2016 = sd_2016,
         sd_dist_to_urb_2010 = sd_2010,
         sd_dist_to_urb_agg_2001 = sd_2001_agg,
         sd_dist_to_urb_agg_2016 = sd_2016_agg)

# Save processed data ####

saveRDS(urban_rsf_df, "input/rsf-data.rds")

# Plots ####

urban_df <- as.data.frame(dev_2016_agg, xy = TRUE)

p1 <- ggplot() +
  geom_raster(data = urban_df, aes(x = x, y = y, fill = EVT_NAME)) +
  geom_sf(data = se, fill = NA) +
  coord_sf() +
  scale_fill_gradient(low = "#FFFFFF", high = "#000000", na.value = NA) +
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
  geom_sf(data = se, fill = NA) +
  geom_sf(data = nest_locs, mapping = aes(col = choice,
                                          shape = choice),
          size = 3, alpha = 0.4) +
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
  geom_sf(data = se, fill = NA) +
  geom_point(data = forag, aes(x = X, y = Y, color = choice), alpha = 0.1) +
  scale_color_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  coord_sf() +
  theme_void() +
  theme(legend.position = "none", panel.grid.major = element_line(colour = "transparent"),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Foraging locations")

ggpubr::ggarrange(p3, p1, p2, nrow = 1, ncol = 3,
                  labels = c("(A)", "(B)", "(C)"))

ggsave("output/mapABC.tiff",
       width = 180, height = 60, units = "mm", dpi = 600, compression = "lzw")
