# Code for Picardi et al., "Fitness consequences of anthropogenic subsidies
# for a partially migratory wading bird"

# Load packages ####

library(sf)
library(terra)
library(tidyverse)

# State boundaries ####

states <- read_sf("input/State boundaries/cb_2018_us_state_5m.shp")

se <- states[states$NAME %in% c("Florida",
                                "Georgia",
                                "South Carolina",
                                "North Carolina"), ]

se <- st_union(se)

# Development layers (Landfire) ####

## 2001 ####

landfire_2001 <- rast("input/US_105_EVT/Tif/us_105evt.tif")

# Get extent of Southeast in same projection as urban layer
ext_se <- se %>%
  st_transform(crs = crs(landfire_2001)) %>%
  st_bbox() %>%
  ext()

# Crop urban layer to Southeast extent
crop_southeast_2001 <- crop(landfire_2001, ext_se)

# Create template in UTM
template_se <- rast(extent = ext(crop_southeast_2001),
                    resolution = 30,
                    crs = crs(crop_southeast_2001))
template_se <- project(template_se, "epsg:32617")

crop_southeast_2001_30m <- project(crop_southeast_2001, template_se)

# Coarsen resolution from 30m to 1km and project into UTM
template_se_1km <- rast(extent = ext(crop_southeast_2001),
                        resolution = 1000,
                        crs = crs(crop_southeast_2001))
template_se_1km <- project(template_se_1km, "epsg:32617")

crop_southeast_2001_1km <- project(crop_southeast_2001, template_se_1km)

# Load metadata from Landfire
meta_2001 <- read.csv("input/US_105_EVT/CSV_Data/us_105evt.csv")

unique(meta_2001$EVT_Name)[unique(grep("Developed", unique(meta_2001$EVT_Name)))]

# Select development classes
dev_values_2001 <- meta_2001 %>%
  filter(EVT_Name %in% c("Developed-Open Space",
                         "Developed-Low Intensity",
                         "Developed-Medium Intensity",
                         "Developed-High Intensity",
                         "Developed-Roads"
  )) %>%
  pull(Value) %>%
  unique()

# Create development raster
dev_2001 <- crop_southeast_2001_1km
values(dev_2001) <- ifelse(values(crop_southeast_2001_1km) %in% dev_values_2001, 1, 0)

dev_2001_agg <- aggregate(crop_southeast_2001_30m, fact = 10,
                          fun = function(x) {any(x %in% dev_values_2001)})

# Reproject boundary of Southeast in UTM
se_reproj <- se %>%
  st_transform(crs = "epsg:32617")

# Transform to SpatVector
se_reproj <- vect(se_reproj)

# Mask to coastline
dev_2001_se <- mask(dev_2001, se_reproj)
dev_2001_agg_se <- mask(dev_2001_agg, se_reproj)

## 2016 ####

urban_2016 <- rast("input/LF2016_EVT_200_CONUS/Tif/LC16_EVT_200.tif")

crop_southeast_2016 <- crop(urban_2016, ext_se)

crop_southeast_2016_30m <- project(crop_southeast_2016, template_se)
crop_southeast_2016_1km <- project(crop_southeast_2016, template_se_1km)

meta_2016 <- read.csv("input/LF2016_EVT_200_CONUS/CSV_Data/LF16_EVT_200.csv")

unique(meta_2016$EVT_NAME)[unique(grep("Developed", unique(meta_2016$EVT_NAME)))]

dev_values_2016 <- meta_2016 %>%
  filter(EVT_NAME %in% c("Developed-Open Space",
                         "Developed-Low Intensity",
                         "Developed-Medium Intensity",
                         "Developed-High Intensity",
                         "Developed-Roads")) %>%
  pull(VALUE) %>%
  unique()

dev_2016 <- crop_southeast_2016_1km
values(dev_2016) <- ifelse(values(crop_southeast_2016_1km) %in% dev_values_2016, 1, 0)

dev_2016_agg <- aggregate(crop_southeast_2016_30m, fact = 10,
                          fun = function(x) {any(x %in% dev_values_2016)})

dev_2016_se <- mask(dev_2016, se_reproj)
dev_2016_agg_se <- mask(dev_2016_agg, se_reproj)

# TIGER urban footprint ####

urb_poly <- read_sf("input/TIGER Urban Footprint/tl_2019_us_uac10.shp")

urb_poly <- st_crop(urb_poly, se)

urb_poly <- vect(urb_poly)

urb_poly <- project(urb_poly, "epsg:32617")

urb_poly_rast <- rasterize(urb_poly, template_se_1km, background = 0)

urb_poly_rast_2010 <- mask(urb_poly_rast, se_reproj)

urb_poly <- aggregate(urb_poly)

# Calculate distance to urban ####

dist_to_urb_2001 <- distance(dev_2001_se, target = 0, exclude = NA)
dist_to_urb_2016 <- distance(dev_2016_se, target = 0, exclude = NA)

dist_to_urb_2010_tiger <- distance(urb_poly_rast_2010, target = 0, exclude = NA)

dist_to_urb_agg_2001 <- distance(dev_2001_agg_se, target = 0, exclude = NA)
dist_to_urb_agg_2016 <- distance(dev_2016_agg_se, target = 0, exclude = NA)

# Water management districts ####

wmd <- read_sf("input/Water Mgmt Districts/Water_Management_District_Boundaries.shp") %>%
  st_transform(32617)

soflo <- st_buffer(wmd[3, ], dist = units::as_units(130, "km"))

# Southeast boundary ####

se <- st_transform(se, 32617)

# Save processed layers ####

writeRaster(dev_2001_se, "processed_layers/development_1km_2001_Landfire_UTM17N.tiff", overwrite = TRUE)
writeRaster(dev_2016_se, "processed_layers/development_1km_2016_Landfire_UTM17N.tiff", overwrite = TRUE)
writeRaster(dist_to_urb_2001, "processed_layers/dist-to-development_1km_2001_Landfire_UTM17N.tiff", overwrite = TRUE)
writeRaster(dist_to_urb_2016, "processed_layers/dist-to-development_1km_2016_Landfire_UTM17N.tiff", overwrite = TRUE)
writeRaster(dist_to_urb_2010_tiger, "processed_layers/dist-to-development_1km_2010_TIGER_UTM17N.tiff", overwrite = TRUE)
writeRaster(dev_2001_agg_se, "processed_layers/development_aggregated_300m_2001_Landfire_UTM17N.tiff", overwrite = TRUE)
writeRaster(dev_2016_agg_se, "processed_layers/development_aggregated_300m_2016_Landfire_UTM17N.tiff", overwrite = TRUE)
writeRaster(dist_to_urb_agg_2001, "processed_layers/dist-to-development_aggregated_300m_2001_Landfire_UTM17N.tiff", overwrite = TRUE)
writeRaster(dist_to_urb_agg_2016, "processed_layers/dist-to-development_aggregated_300m_2016_Landfire_UTM17N.tiff", overwrite = TRUE)
st_write(wmd, "processed_layers/fl-water-mgmt-districts_UTM17N.shp", append = FALSE)
st_write(soflo, "processed_layers/south-florida-buffer_UTM17N.shp", append = FALSE)
st_write(se, "processed_layers/southeast-boundary_UTM17N.shp", append = FALSE)
writeVector(urb_poly, "processed_layers/urban_polygons_TIGER_UTM17N.shp", overwrite = TRUE)

# Change in development from 2001 to 2016 ####

# See what the proportion of developed/undeveloped was in 2001
urb_yn_2001_se <- table(values(dev_2001_agg_se))

urb_yn_2001_se[2]/sum(urb_yn_2001_se)

# Versus in 2016
urb_yn_2016_se <- table(values(dev_2016_agg_se))

urb_yn_2016_se[2]/sum(urb_yn_2016_se)

# Do the same for Florida only
fl <- states[states$NAME == "Florida", ]
fl <- st_transform(fl, crs(dev_2001_agg_se))

dev_2001_fl <- mask(dev_2001_agg_se, fl)

urb_yn_2001_fl <- table(values(dev_2001_fl))

urb_yn_2001_fl[2]/sum(urb_yn_2001_fl)

# Versus in 2016

dev_2016_fl <- mask(dev_2016_agg_se, fl)

urb_yn_2016_fl <- table(values(dev_2016_fl))

urb_yn_2016_fl[2]/sum(urb_yn_2016_fl)
