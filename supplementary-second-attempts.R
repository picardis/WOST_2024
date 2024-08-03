# Load packages ####

library(sf)
library(tidyverse)

# Load data ####

# Load WOST GPS data for breeding attempts
attempts <- readRDS("input/nesting_attempts_110.rds")
attempts_raw <- readRDS("input/nesting_attempts.rds")
# Load data on individual migratory behavior
mig_yrs <- readRDS("input/mig_yrs.rds")
# Load data on nests
ws_nests <- readRDS("input/ws_nests_all.rds")

# Load GIS layers
wmd <- st_read("processed_layers/fl-water-mgmt-districts_UTM17N.shp")
se <- st_read("processed_layers/southeast-boundary_UTM17N.shp")
soflo <- st_read("processed_layers/south-florida-buffer_UTM17N.shp")

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

nesting_data <- nesting_data %>%
  left_join(nest_coords, by = "att_id")

nest_soflo <- st_as_sf(nesting_data, coords = c("nest_long", "nest_lat"), crs = 4326) %>%
  st_transform(32617) %>%
  st_intersects(soflo, sparse = FALSE)

nesting_data$soflo <- as.logical(nest_soflo)

jax_nests <- nesting_data %>%
  filter(choice == "resident" & nest_lat > 30) %>%
  pull(att_id) %>%
  unique()

nesting_data <- nesting_data %>%
  mutate(jax = att_id %in% jax_nests,
         se = !jax & !soflo)

# All nests corresponding to the attempts in nesting_data
nest_loc_data <- nesting_data %>%
  distinct(att_id, choice, nest_long, nest_lat) %>%
  arrange(desc(nest_lat))

se_ll <- st_transform(se, 4326)

ggplot() +
  geom_sf(data = se_ll) +
  geom_point(data = nest_loc_data,
             aes(x = nest_long, y = nest_lat, color = choice))

soflo <- nesting_data %>%
  filter(soflo & day_of_att == 1) %>%
  dplyr::select(burst, att_id, choice, date, nest_long, nest_lat) %>%
  distinct() %>%
  mutate(doy = lubridate::yday(date))

se <- nesting_data %>%
  filter(se & day_of_att == 1) %>%
  dplyr::select(burst, att_id, choice, date, nest_long, nest_lat) %>%
  distinct() %>%
  mutate(doy = lubridate::yday(date))

# How many of the birds that nested in the Southeast also had an attempt in
# South FL in the same year?
sum(se$burst %in% soflo$burst)
sum(se$burst %in% soflo$burst)/nrow(se)

comparison <- se %>%
  left_join(soflo, by = "burst", relationship = "many-to-many",
            unmatched = "drop")
# 24% of the birds that nest out of FL have a previous attempt in FL in the
# same year. The others could be birds that go down to assess and don't even
# try or birds that try but fail too early for the attempt to be detected
