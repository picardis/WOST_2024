# Load packages ####

library(tidyverse)

# Load data ####

dat <- readRDS("input/rsf-data.rds")

# Trends in use of urban areas? ####

dat %>%
  filter(used == 1) %>%
  mutate(year = year(date_ymd)) %>%
  group_by(year) %>%
  summarize(mean_duy = mean(dist_to_urb_2016, na.rm = TRUE),
            sd_duy = sd(dist_to_urb_2016, na.rm = TRUE),
            n = n(),
            se = sd_duy/sqrt(n)) %>%
  ggplot(aes(x = factor(year), y = mean_duy)) +
  geom_point() +
  # geom_errorbar(aes(ymin = mean_duy - sd_duy, ymax = mean_duy + sd_duy),
  #               width = 0.3) +
  geom_errorbar(aes(ymin = mean_duy - se, ymax = mean_duy + se),
                width = 0.3) +
  labs(x = " ", y = "Mean distance to urban areas") +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw()

dat %>%
  filter(used == 1) %>%
  mutate(year = year(date_ymd)) %>%
  group_by(year) %>%
  summarize(mean_duy = mean(dist_to_urb_2001, na.rm = TRUE),
            sd_duy = sd(dist_to_urb_2001, na.rm = TRUE),
            n = n(),
            se = sd_duy/sqrt(n)) %>%
  ggplot(aes(x = factor(year), y = mean_duy)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_duy - sd_duy, ymax = mean_duy + sd_duy),
                width = 0.3) +
  labs(x = " ", y = "Mean distance to urban areas") +
  theme_bw()

# Trends in migratory status? ####

dat %>%
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

dat %>%
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

# Trends in nest survival? ####
