# Supplementary materials

# Load packages ####

library(tidyverse)
library(sf)
library(terra)
library(patchwork)

# Load data ####

urban_rsf_df <- readRDS("input/rsf-data.rds")
ind_info <- readRDS("input/wost_info.rds") %>%
  janitor::clean_names() %>%
  mutate(satellite_id = as.character(satellite_id))

# Load WOST GPS data for breeding attempts
attempts <- readRDS("input/nesting_attempts_110.rds")
attempts_raw <- readRDS("input/nesting_attempts.rds")
# Load data on individual migratory behavior
mig_yrs <- readRDS("input/mig_yrs.rds")
# Load data on nests
ws_nests <- readRDS("input/ws_nests_all.rds")

# Load layers ####

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

# Sample composition ####

urban_rsf_df %>%
  select(stork) %>%
  distinct() %>%
  left_join(ind_info, by = c("stork" = "satellite_id")) %>%
  select(stork, sex, general_area_of_origin, state_of_origin) %>%
  group_by(sex, general_area_of_origin, state_of_origin) %>%
  tally() %>%
  arrange(general_area_of_origin, desc(sex))

# Definition of regions ####

nests <- urban_rsf_df %>%
  mutate(
    region = case_when(jax ~ "Jacksonville",
                      soflo ~ "South Florida",
                      TRUE ~ "Southeast")
  ) %>%
  select(att_id, nest_long, nest_lat, region) %>%
  distinct()

nests <- st_as_sf(nests, coords = c("nest_long", "nest_lat"), crs = 4326)

ggplot() +
  geom_sf(data = se, fill = "white") +
  geom_sf(data = wmd[3, ]) +
  geom_sf(data = nests, aes(color = region)) +
  theme_bw() +
  labs(color = "Region")

ggsave("output/regional-groups.tiff", compression = "lzw", dpi = 300,
       width = 6, height = 4)

# Expansion of urban areas ####

# See what the proportion of developed/undeveloped was in 2001
urb_yn_2001_se <- table(values(dev_2001_se))

urb_yn_2001_se[2]/sum(urb_yn_2001_se)

# Versus in 2016
urb_yn_2016_se <- table(values(dev_2016_se))

urb_yn_2016_se[2]/sum(urb_yn_2016_se)

urb2001 <- ggplot() +
  geom_sf(data = se, fill = "white") +
  geom_raster(data = as.data.frame(dev_2001_se, xy = TRUE),
              aes(x = x, y = y, fill = factor(VALUE_1))) +
  scale_fill_manual(values = c("transparent", "tomato")) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("2001")

urb2016 <- ggplot() +
  geom_sf(data = se, fill = "white") +
  geom_raster(data = as.data.frame(dev_2016_se, xy = TRUE),
              aes(x = x, y = y, fill = factor(EVT_NAME))) +
  scale_fill_manual(values = c("transparent", "tomato")) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("2016")

urb2001 + urb2016

ggsave("output/urban_2001-2016.tiff", compression = "lzw", dpi = 300,
       width = 8, height = 4)

urb2001 <- ggplot() +
  geom_sf(data = se, fill = "white") +
  geom_raster(data = as.data.frame(dev_2001_agg, xy = TRUE),
              aes(x = x, y = y, fill = factor(VALUE_1))) +
  scale_fill_manual(values = c("transparent", "tomato")) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("2001")

urb2016 <- ggplot() +
  geom_sf(data = se, fill = "white") +
  geom_raster(data = as.data.frame(dev_2016_agg, xy = TRUE),
              aes(x = x, y = y, fill = factor(EVT_NAME))) +
  scale_fill_manual(values = c("transparent", "tomato")) +
  theme_bw() +
  theme(legend.position = "none")  +
  ggtitle("2016")

urb2001 + urb2016

ggsave("output/urban_2001-2016_aggregated.tiff", compression = "lzw", dpi = 300,
       width = 8, height = 4)

# Metrics of urban ####

u1 <- ggplot() +
  geom_sf(data = se, fill = "white") +
  geom_raster(data = as.data.frame(dev_2001_agg, xy = TRUE),
              aes(x = x, y = y, fill = factor(VALUE_1))) +
  scale_fill_manual(values = c("transparent", "tomato")) +
  theme_bw() +
  theme(legend.position = "none")  +
  ggtitle("Landfire 2001 (aggregated 1-km)")

u2 <- ggplot() +
  geom_sf(data = se, fill = "white") +
  geom_raster(data = as.data.frame(dev_2016_agg, xy = TRUE),
              aes(x = x, y = y, fill = factor(EVT_NAME))) +
  scale_fill_manual(values = c("transparent", "tomato")) +
  theme_bw() +
  theme(legend.position = "none")  +
  ggtitle("Landfire 2016 (aggregated 1-km)")

tiger2 <- st_intersection(tiger, se)

u3 <- ggplot() +
  geom_sf(data = se, fill = "white") +
  geom_sf(data = tiger2, fill = "tomato") +
  theme_bw() +
  theme(legend.position = "none")  +
  ggtitle("Tiger 2010")

u4 <- ggplot() +
  geom_sf(data = se, fill = "white") +
  geom_raster(data = as.data.frame(dist_to_urb_2001, xy = TRUE),
              aes(x = x, y = y, fill = VALUE_1)) +
  theme_bw() +
  theme(legend.position = "none")  +
  scale_fill_viridis_c(direction = -1, option = "magma") +
  ggtitle("Distance using Landfire 2001")

u5 <- ggplot() +
  geom_sf(data = se, fill = "white") +
  geom_raster(data = as.data.frame(dist_to_urb_2016, xy = TRUE),
              aes(x = x, y = y, fill = EVT_NAME)) +
  theme_bw() +
  theme(legend.position = "none")  +
  scale_fill_viridis_c(direction = -1, option = "magma") +
  ggtitle("Distance using Landfire 2016")

u6 <- ggplot() +
  geom_sf(data = se, fill = "white") +
  geom_raster(data = as.data.frame(dist_to_urb_2010, xy = TRUE),
              aes(x = x, y = y, fill = layer)) +
  theme_bw() +
  theme(legend.position = "none")  +
  scale_fill_viridis_c(direction = -1, option = "magma") +
  ggtitle("Distance using TIGER 2010")

(u1 + u2 + u3)/ (u4 + u5 + u6)

ggsave("output/urban_metrics.tiff", compression = "lzw", dpi = 300,
       width = 10, height = 10)

# Distance of nests to urban ####

nests_utm <- st_transform(nests, crs = 32617)

nests_utm$dist <- extract(dist_to_urb_2016, nests_utm)[, 2]

urban_rsf_df %>%
  select(att_id, choice, group) %>%
  distinct() %>%
  left_join(nests_utm) %>%
  rename(dist_to_urb = dist) %>%
  ggplot() +
  geom_density(aes(x = dist_to_urb/1000, fill = choice)) +
  facet_grid(group ~ choice) +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"),
                       begin = 0.25, end = 0.75) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Distance of nests to urban development (km)")

ggsave("output/dist_nest_to_urban.tiff", compression = "lzw", dpi = 300,
       width = 6, height = 4)

# Empirical distribution of foraging trips ####

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
  ggtitle("Observed") +
  coord_cartesian(ylim = c(0, 0.06))

simul <- urban_rsf_df %>%
  filter(used == 0) %>%
  mutate(dist_to_nest = geosphere::distGeo(p1 = cbind(long, lat), p2 = cbind(nest_long, nest_lat))/1000) %>%
  ggplot(aes(x = dist_to_nest)) +
  geom_density() +
  stat_function(fun = dexp, n = 100, args = list(rate = 0.05), col = "tomato", lwd = 1.5, alpha = 0.7) +
  theme_minimal() +
  scale_y_continuous("Density") +
  scale_x_continuous("Foraging trip distance (km)") +
  ggtitle("Simulated") +
  coord_cartesian(ylim = c(0, 0.06))

ggpubr::ggarrange(empir, simul, nrow = 1, ncol = 2,
                  labels = c("(A)", "(B)"))

ggsave("output/trip_distance_empir_vs_simul.tiff", width = 8, height = 4, dpi = 500, compression = "lzw")

# Used versus available ####

urb_rsf_ll <- urban_rsf_df
urb_rsf_ll <- st_as_sf(urb_rsf_ll, coords = c("X", "Y"), crs = 32617)
urb_rsf_ll <- st_transform(urb_rsf_ll, 4326)
urb_rsf_ll <- st_coordinates(urb_rsf_ll)
colnames(urb_rsf_ll) <- c("long", "lat")

urban_rsf_df <- cbind(urban_rsf_df, urb_rsf_ll)

# Mean distance of foraging locations from nest for migrants and residents
urban_rsf_df %>%
  filter(used == 1) %>%
  mutate(dist_to_nest = geosphere::distGeo(p1 = cbind(long, lat), p2 = cbind(nest_long, nest_lat))) %>%
  group_by(choice, group) %>%
  summarize(mean(dist_to_nest)/1000)

# Mean distance of foraging locations from urban for migrants and residents
urban_rsf_df %>%
  filter(used == 1) %>%
  group_by(choice, group) %>%
  summarize(mean(dist_to_urb_2016/1000, na.rm = T))

# Plot of distance of foraging locations to urban for migrants and residents (used vs. available)

mig_labs <- c(migrant = "Migrant",
              resident = "Resident")

avail <- urban_rsf_df %>%
  filter(used == 0) %>%
  mutate(dist_to_urban = dist_to_urb_2016/1000)

urban_rsf_df %>%
  filter(used == 1) %>% # used points
  mutate(dist_to_urban = dist_to_urb_2016/1000) %>%
  ggplot(aes(x = dist_to_urban, fill = choice, group = choice)) +
  geom_density(data = avail, fill = "gray80") +
  geom_density(alpha = 0.5) +
  facet_grid(choice ~ group, labeller = labeller(choice = mig_labs)) +
  scale_fill_viridis_d(begin = 0.25, end = 0.75) +
  theme_minimal() +
  scale_y_continuous("Density") +
  scale_x_continuous("Distance to urban development (km)") +
  theme(legend.position = "none")

used_plot <- urban_rsf_df %>%
  filter(used == 1) %>% # used points
  mutate(dist_to_urban = dist_to_urb_2016/1000) %>%
  ggplot(aes(x = dist_to_urban, fill = choice, group = choice)) +
  geom_density() +
  facet_grid(choice~group, labeller = labeller(choice = mig_labs)) +
  scale_fill_viridis_d(begin = 0.25, end = 0.75) +
  theme_minimal() +
  scale_y_continuous("Density", limits = c(0, 0.8)) +
  scale_x_continuous("Distance to urban development (km)") +
  theme(legend.position = "none") +
  ggtitle("Used")

avail_plot <- urban_rsf_df %>%
  filter(used == 0) %>% # available points
  mutate(dist_to_urban = dist_to_urb_2016/1000) %>%
  ggplot(aes(x = dist_to_urban, fill = choice, group = choice)) +
  geom_density(alpha = 0.5) +
  facet_grid(choice~group, labeller = labeller(choice = mig_labs)) +
  scale_fill_viridis_d(begin = 0.25, end = 0.75) +
  theme_minimal() +
  scale_y_continuous("Density", limits = c(0, 0.8)) +
  scale_x_continuous("Distance to urban development (km)") +
  theme(legend.position = "none") +
  ggtitle("Available")

ggpubr::ggarrange(used_plot, avail_plot, nrow = 2, ncol = 1,
                  labels = c("(A)", "(B)"))

ggsave("output/dist_to_urb_used_vs_avail.tiff", width = 10, height = 6,
       dpi = 300, compression = "lzw")

# Plot availability by choice
urban_avail <- urban_rsf_df %>%
  filter(used == 0)

ggplot(urban_avail, aes(x = dist_to_urb_2016/1000, fill = choice, group = choice)) +
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
  summarize(mean_dist_used = mean(dist_to_urb_2016)/1000)

av <- urban_rsf_df %>%
  group_by(att_id) %>%
  filter(used == 0) %>%
  summarize(mean_dist_avai = mean(dist_to_urb_2016)/1000)

discrep <- left_join(us, av, by = "att_id") %>%
  mutate(discrep = abs(mean_dist_avai - mean_dist_used)) %>%
  arrange(desc(discrep)) %>%
  pull(att_id)

urban_rsf_df$att_id_reord <- factor(urban_rsf_df$att_id, levels = discrep)

nests <- urban_rsf_df %>%
  dplyr::select(att_id_reord, nest_long, nest_lat) %>%
  distinct()

# Facultative migrants ####

# Find facultative migrants

fac <- urban_rsf_df %>%
  dplyr::select(choice, stork) %>%
  distinct() %>%
  group_by(stork) %>%
  tally() %>%
  arrange(desc(n)) %>%
  filter(n > 1) %>%
  pull(stork)

# Plot distance to urban in resident vs migrant years

# 2016
fac1 <- urban_rsf_df %>%
  filter(stork == fac[1] & used == 1) %>%
  ggplot(aes(x = dist_to_urb_2016/1000, fill = choice)) +
  geom_density() +
  facet_wrap(~ choice) +
  # facet_grid(stork ~ att_id) +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  theme_bw() +
  labs(x = "Distance of foraging sites from urban areas (km)", title = fac[1])

fac2 <- urban_rsf_df %>%
  filter(stork == fac[2] & used == 1) %>%
  ggplot(aes(x = dist_to_urb_2016/1000, fill = choice)) +
  geom_density() +
  facet_wrap(~ choice) +
  # facet_grid(stork ~ att_id) +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  theme_bw() +
  labs(x = "Distance of foraging sites from urban areas (km)", title = fac[2],
       y = "N")

fac1 + fac2 + plot_layout(guides = "collect")

ggsave("output/facultative.tiff", width = 10, height = 5,
       dpi = 300, compression = "lzw")

# 2001
urban_rsf_df %>%
  filter(stork == fac[1] & used == 1) %>%
  ggplot(aes(x = dist_to_urb_2001, group = att_id, fill = choice)) +
  geom_histogram() +
  facet_wrap(~ choice) +
  # facet_grid(stork ~ att_id) +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75)

urban_rsf_df %>%
  filter(stork == fac[2] & used == 1) %>%
  ggplot(aes(x = dist_to_urb_2001, group = att_id, fill = choice)) +
  geom_histogram() +
  facet_wrap(~ choice) +
  # facet_grid(stork ~ att_id) +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75)

# Hypothesis testing on distance to urban

res_years_stork1 <- urban_rsf_df%>%
  filter(stork == fac[1] & choice == "resident" & used == 1)
res_years_stork2 <- urban_rsf_df%>%
  filter(stork == fac[2] & choice == "resident" & used == 1)
mig_years_stork1 <- urban_rsf_df%>%
  filter(stork == fac[1] & choice == "migrant" & used == 1)
mig_years_stork2 <- urban_rsf_df%>%
  filter(stork == fac[2] & choice == "migrant" & used == 1)

wilcox.test(res_years_stork1$dist_to_urb_2016,
            mig_years_stork1$dist_to_urb_2016,
            alternative = "less") # Reject null
wilcox.test(res_years_stork1$dist_to_urb_2001,
            mig_years_stork1$dist_to_urb_2001,
            alternative = "less") # Reject null
wilcox.test(res_years_stork2$dist_to_urb_2016,
            mig_years_stork2$dist_to_urb_2016,
            alternative = "less") # Fail to reject null
wilcox.test(res_years_stork2$dist_to_urb_2001,
            mig_years_stork2$dist_to_urb_2001,
            alternative = "less") # Fail to reject null
# Both of these are Jacksonville storks

# Nest locations of facultative migrants

nest_fac1 <- urban_rsf_df%>%
  filter(stork == fac[1] & used == 1) %>%
  mutate(year = year(date_ymd)) %>%
  dplyr::select(choice, year, nest_long, nest_lat) %>%
  distinct() %>%
  arrange(year)

nest_fac1 <- st_as_sf(nest_fac1, coords = c("nest_long", "nest_lat"), crs = 4326)

nest_fac2 <- urban_rsf_df%>%
  filter(stork == fac[2] & used == 1) %>%
  mutate(year = year(date_ymd)) %>%
  dplyr::select(choice, year, nest_long, nest_lat) %>%
  distinct() %>%
  arrange(year)

nest_fac2 <- st_as_sf(nest_fac2, coords = c("nest_long", "nest_lat"), crs = 4326)

map_fac1 <- ggplot() +
  geom_sf(data = se) +
  geom_sf(data = nest_fac1, aes(color = choice)) +
  facet_wrap(~ year) +
  scale_color_viridis_d("Tactic", labels = c("Migrant", "Resident"),
                        begin = 0.25, end = 0.75) +
  theme_bw() +
  labs(title = fac[1])

ggsave("output/facultative1_map.tiff", map_fac1, width = 10, height = 10,
       dpi = 300, compression = "lzw")

map_fac2 <- ggplot() +
  geom_sf(data = se) +
  geom_sf(data = nest_fac2, aes(color = choice)) +
  facet_wrap(~ year) +
  scale_color_viridis_d("Tactic", labels = c("Migrant", "Resident"),
                        begin = 0.25, end = 0.75) +
  theme_bw() +
  labs(title = fac[2])

ggsave("output/facultative2_map.tiff", map_fac2, width = 10, height = 10,
       dpi = 300, compression = "lzw")

# Temporal trends ####

# Trends in use of urban areas?

urban_rsf_df %>%
  filter(used == 1) %>%
  mutate(year = year(date_ymd)) %>%
  group_by(year) %>%
  ggplot(aes(x = factor(year), y = dist_to_urb_2016/1000)) +
  geom_boxplot() +
  labs(x = " ", y = "Mean distance to urban areas (km)") +
  theme_bw()

ggsave("output/trend_urban.tiff", width = 5, height = 5,
       dpi = 300, compression = "lzw")

# Trends in migratory status

urban_rsf_df %>%
  filter(used == 1) %>%
  mutate(year = year(date_ymd)) %>%
  select(stork, year, choice) %>%
  distinct() %>%
  group_by(year, choice) %>%
  tally() %>%
  ggplot(aes(x = factor(year), y = n, fill = choice)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  labs(x = "", y = "Number of individuals") +
  theme_bw()

ggsave("output/trend_mig.tiff", width = 5, height = 5,
       dpi = 300, compression = "lzw")

urban_rsf_df %>%
  filter(used == 1) %>%
  mutate(year = year(date_ymd)) %>%
  select(stork, year, choice) %>%
  distinct() %>%
  group_by(year, choice) %>%
  tally() %>%
  ggplot(aes(x = factor(year), y = n, fill = choice)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ choice) +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  labs(x = "", y = "Number of individuals") +
  theme_bw()

# Relationship between migration status and breeding status ####

# Fix burst column to match between data frames
mig_yrs$burst <- paste0(mig_yrs$stork, "-", mig_yrs$year)

# Join raw attempts data with migration data and with environmental data
nesting_data <- left_join(attempts_raw, mig_yrs, by = "burst") %>%
  dplyr::select(att_id, burst, stork, choice, strategy, flexible,
                day_of_att, date, datetime, long, lat, loc_id)

breeders <- nesting_data %>%
  select(burst, choice, strategy) %>%
  distinct() %>%
  filter(!is.na(choice)) %>%
  mutate(breeding_status = "Breeder")

status <- mig_yrs %>%
  left_join(breeders) %>%
  mutate(breeding_status = case_when(
    is.na(breeding_status) ~ "Non-breeder",
    TRUE ~ breeding_status
  ),
  strategy_new = case_when(
    strategy == "con_mig" ~ "Consistent migrants",
    strategy == "con_res" ~ "Consistent residents",
    strategy == "fac_mig" ~ "Facultative migrants"
  ))

status %>%
  group_by(strategy_new, choice, breeding_status) %>%
  tally() %>%
  ggplot(aes(x = breeding_status, y = n, fill = choice)) +
  facet_wrap(~ strategy_new) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  labs(x = "", y = "Number of individual-years") +
  theme_bw()

ggsave("output/breeding_status.tiff", width = 8, height = 5,
       dpi = 300, compression = "lzw")

status %>%
  group_by(strategy_new, choice, breeding_status) %>%
  tally() %>%
  ggplot(aes(x = breeding_status, y = n, fill = choice)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  labs(x = "", y = "Number of individual-years") +
  theme_bw()
