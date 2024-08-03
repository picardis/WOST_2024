# Load packages ####

library(tidyverse)
library(nestR)

# Load WOST data ####

# Load GPS data for adult WOST
gps <- readRDS("input/GPS_WOST_adults.rds")

# Load data on annual individual migratory behavior. This object was generated
# using methods described in Picardi, S., Frederick, P. C., Borkhataria, R. R.,
# & Basille, M. (2020). Partial migration in a subtropical wading bird in the
# southeastern United States. Ecosphere, 11(2), e03054.
mig_yrs <- readRDS("input/mig_yrs.rds")

# Find nesting attempts ####

# Use parameter values described in Picardi, S., Smith, B. J., Boone, M. E.,
# Frederick, P. C., Cecere, J. G., Rubolini, D., ... & Basille, M. (2020).
# Analysis of movement recursions to detect reproductive events and estimate
# their fate in central place foragers. Movement Ecology, 8(1), 1-14.

ws_nests <- find_nests(gps_data = gps,
                       buffer = 40,
                       sea_start = "11-01",
                       sea_end = "10-31",
                       nest_cycle = 110,
                       min_d_fix = 5,
                       min_consec = 14,
                       min_top_att = 79,
                       min_days_att = 1,
                       min_pts = 2,
                       discard_overlapping = TRUE)

# Isolate GPS data within nesting attempts ####

get_attempts <- function(nest_info,
                         nest_cycle = 110) {

  # Create unique attempt identifier
  attempts <- nest_info$nests %>%
    mutate(attempt_id = paste0(burst, "-", loc_id))

  # Initialize output
  output <- data.frame()

  # Loop over attempts
  for (i in 1:nrow(attempts)) {

    # Select current attempt
    att <- attempts[i,]

    # Data on nest revisits
    visits <- nest_info$visits %>%
      filter(burst == att$burst)

    # Cut between attempt start and end of nesting cycle
    visits <- visits %>%
      filter(between(date,
                     att$attempt_start,
                     att$attempt_end)) %>%
      mutate(att_id = att$attempt_id) %>%
      mutate(day_of_att = floor(as.numeric(difftime(as_date(date), att$attempt_start, units = "days"))) + 1) %>%
      dplyr::select(att_id, burst, date, day_of_att, long, lat, loc_id)

    # Append to output
    output <- rbind(output, visits)

  }

  return(output)

}

nest_att <- get_attempts(ws_nests)

# Make date and datetime two distinct columns
nest_att$datetime <- nest_att$date
nest_att$date <- as_date(nest_att$date)

# Save ####

saveRDS(ws_nests, "input/ws_nests_all.rds")
