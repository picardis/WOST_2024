# Code for Picardi et al., "Fitness consequences of anthropogenic subsidies
# for a partially migratory wading bird"

# Load packages ####

library(tidyverse)

# Load data ####

dat <- readRDS("input/rsf-data.rds")
dat

# Find facultative migrants ####

fac <- dat %>%
  dplyr::select(choice, stork) %>%
  distinct() %>%
  group_by(stork) %>%
  tally() %>%
  arrange(desc(n)) %>%
  filter(n > 1) %>%
  pull(stork)

# Plot distance to urban in resident vs migrant years ####

# 2016
dat %>%
  filter(stork == fac[1] & used == 1) %>%
  ggplot(aes(x = dist_to_urb_2016, group = att_id, fill = choice)) +
  geom_histogram() +
  facet_wrap(~ att_id) +
  # facet_grid(stork ~ att_id) +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75)

dat %>%
  filter(stork == fac[2] & used == 1) %>%
  ggplot(aes(x = dist_to_urb_2016, group = att_id, fill = choice)) +
  geom_histogram() +
  facet_wrap(~ att_id) +
  # facet_grid(stork ~ att_id) +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75)

# 2001
dat %>%
  filter(stork == fac[1] & used == 1) %>%
  ggplot(aes(x = dist_to_urb_2001, group = att_id, fill = choice)) +
  geom_histogram() +
  facet_wrap(~ att_id) +
  # facet_grid(stork ~ att_id) +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75)

dat %>%
  filter(stork == fac[2] & used == 1) %>%
  ggplot(aes(x = dist_to_urb_2001, group = att_id, fill = choice)) +
  geom_histogram() +
  facet_wrap(~ att_id) +
  # facet_grid(stork ~ att_id) +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75)

# Hypothesis testing on distance to urban ####

res_years_stork1 <- dat %>%
  filter(stork == fac[1] & choice == "resident" & used == 1)
res_years_stork2 <- dat %>%
  filter(stork == fac[2] & choice == "resident" & used == 1)
mig_years_stork1 <- dat %>%
  filter(stork == fac[1] & choice == "migrant" & used == 1)
mig_years_stork2 <- dat %>%
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

# Nest locations ####

dat %>%
  filter(stork == fac[1] & used == 1) %>%
  mutate(year = year(date_ymd)) %>%
  dplyr::select(choice, year, nest_long, nest_lat) %>%
  distinct() %>%
  arrange(year)

dat %>%
  filter(stork == fac[2] & used == 1) %>%
  mutate(year = year(date_ymd)) %>%
  dplyr::select(choice, year, nest_long, nest_lat) %>%
  distinct() %>%
  arrange(year)
