# Function to select only 1 sample per SiteID, the closest to the median across all SiteIDs,
# and to flag those cases where the sampling date is inside-outside a certain buffer of days


library(dplyr)
library(lubridate)

select_one_sample_and_flag <- function(Data,
                                       site_col = "SiteID",
                                       date_col = "Date",
                                       buffer_days = 30) {
  
  doy_dist <- function(a, b) {
    pmin(abs(a - b), 365 - abs(a - b))
  }
  
  df <- Data %>%
    mutate(
      .Site = .data[[site_col]],
      .Date = as.Date(.data[[date_col]]),
      .Year = year(.Date),
      .DOY = yday(.Date)
    ) %>%
    filter(!is.na(.Site), !is.na(.Date), !is.na(.Year), !is.na(.DOY))
  
  median_doy <- median(df$.DOY, na.rm = TRUE)
  
  df %>%
    mutate(
      median_doy = median_doy,
      dist_to_median_doy = doy_dist(.DOY, median_doy)
    ) %>%
    group_by(.Site, .Year) %>%
    slice_min(dist_to_median_doy, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      within_median_buffer = dist_to_median_doy <= buffer_days
    ) %>%
    select(-.Site, -.Date, -.Year, -.DOY)
}


# # select relevant columns
# Data_sel <- Data %>%
#   select(SiteID, Date) %>% distinct()
# unique_siteid <- unique(Data_sel$SiteID)
# 
# 
# # run function
# result <- select_one_sample_and_flag(
#   Data = Data_sel,
#   site_col = "SiteID",
#   date_col = "Date",
#   buffer_days = 45
# )