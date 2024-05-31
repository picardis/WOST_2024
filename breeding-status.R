# Load packages ####

library(tidyverse)

# Load data ####

# Load WOST GPS data for breeding attempts
attempts_raw <- readRDS("input/nesting_attempts.rds")
# Load data on individual migratory behavior
mig_yrs <- readRDS("input/mig_yrs.rds")
# Load data on nests
ws_nests <- readRDS("input/ws_nests_all.rds")


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

status %>%
  group_by(strategy_new, choice, breeding_status) %>%
  tally() %>%
  ggplot(aes(x = breeding_status, y = n, fill = choice)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d("Tactic", labels = c("Migrant", "Resident"), begin = 0.25, end = 0.75) +
  labs(x = "", y = "Number of individual-years") +
  theme_bw()
