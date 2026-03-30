#' Harmonize sampling effort across sites and years
#'
#' Identifies recurring sampling periods at the group level and matches
#' site-level observations to those reference periods in order to retain
#' a comparable sampling design across sites and years.
#'
#' @param Data A data frame containing site, group, and date information.
#' @param site_col Name of the site column.
#' @param group_col Name of the grouping column used to define shared
#'   reference sampling periods.
#' @param date_col Name of the date column.
#' @param min_years Minimum number of years required for a site to be included.
#' @param buffer_months Allowed time window around each reference sampling
#'   period, expressed in months.
#'
#' @return A data frame containing only rows retained after sampling
#'   harmonization, with visit assignments and reference information.
#'
#' @examples
#' result <- harmonize_sampling(
#'   Data = df,
#'   site_col = "SiteID",
#'   group_col = "SiteID", # Or Country etc..
#'   date_col = "Date",
#'   min_years = 3,
#'   buffer_months = 3
#' )
#'
#' @export
#' 

library(dplyr)
library(lubridate)
library(purrr)
library(tidyr)

# ============================================================
# MAIN FUNCTION
# ============================================================
#
# Harmonize sampling effort across sites by identifying recurring
# sampling periods at the group level and then matching site-level
# observations to those reference periods.
#
# The function:
# 1. Builds group-level reference sampling periods
# 2. Determines a common sampling effort across sites and years
# 3. Matches each site's yearly samples to the reference periods
# 4. Retains only samples that fit the harmonized design
harmonize_sampling <- function(Data,
                               site_col = "SiteID",
                               group_col = "Country",
                               date_col = "Date",
                               min_years = 3,
                               buffer_months = 3) {
  
  # Convert the buffer from months to approximate days.
  buffer_days <- round(buffer_months * 30.44)
  
  # Standardize key columns into internal working variables and
  # derive year and day-of-year from the date column.
  df <- Data %>%
    mutate(
      .Site = .data[[site_col]],
      .Group = .data[[group_col]],
      .Date = as.Date(.data[[date_col]]),
      .Year = lubridate::year(.Date),
      .DOY = lubridate::yday(.Date)
    ) %>%
    filter(!is.na(.Site), !is.na(.Group), !is.na(.Date), !is.na(.Year), !is.na(.DOY))
  
  # --------------------------------------------------
  # Step 1. Find common effort and recurrent periods
  # at the group level
  # --------------------------------------------------
  #
  # For each group:
  # - keep only sites with at least min_years of data
  # - determine the strict common effort k, defined as the minimum
  #   number of samples recorded in any site-year combination
  # - identify k recurring sampling periods (day-of-year centers)
  #   that maximize coverage across site-years
  get_group_reference <- function(gdat) {
    
    site_year_n <- gdat %>%
      count(.Site, .Year, name = "n_year")
    
    # Keep only sites with data from at least min_years years.
    valid_sites <- site_year_n %>%
      count(.Site, name = "n_years") %>%
      filter(n_years >= min_years) %>%
      pull(.Site)
    
    gdat2 <- gdat %>%
      filter(.Site %in% valid_sites)
    
    if (nrow(gdat2) == 0) {
      return(NULL)
    }
    
    site_year_n2 <- gdat2 %>%
      count(.Site, .Year, name = "n_year")
    
    # Define strict common effort across the group as the smallest
    # number of samples observed in any site-year.
    k <- min(site_year_n2$n_year)
    
    if (is.na(k) || k < 1) {
      return(NULL)
    }
    
    remaining <- gdat2
    chosen_centers <- numeric(0)
    
    # Iteratively identify k reference periods.
    for (j in seq_len(k)) {
      
      candidate_doys <- sort(unique(remaining$.DOY))
      
      if (length(candidate_doys) == 0) {
        return(NULL)
      }
      
      # Score each possible center day-of-year by:
      # - how many site-years it can match within the buffer
      # - the average distance of those matches
      candidate_scores <- map_dfr(candidate_doys, function(center) {
        
        matched <- remaining %>%
          mutate(dist = abs(.DOY - center)) %>%
          filter(dist <= buffer_days) %>%
          group_by(.Site, .Year) %>%
          slice_min(dist, n = 1, with_ties = FALSE) %>%
          ungroup()
        
        tibble(
          center_doy = center,
          n_site_years = n_distinct(paste(matched$.Site, matched$.Year, sep = "__")),
          mean_dist = if (nrow(matched) > 0) mean(matched$dist) else Inf
        )
      })
      
      # Choose the best center by:
      # 1. maximizing matched site-years
      # 2. minimizing average distance
      # 3. using the earliest day-of-year as a tie-breaker
      best_center <- candidate_scores %>%
        arrange(desc(n_site_years), mean_dist, center_doy) %>%
        slice(1) %>%
        pull(center_doy)
      
      best_matches <- remaining %>%
        mutate(dist = abs(.DOY - best_center)) %>%
        filter(dist <= buffer_days) %>%
        group_by(.Site, .Year) %>%
        slice_min(dist, n = 1, with_ties = FALSE) %>%
        ungroup()
      
      if (nrow(best_matches) == 0) {
        return(NULL)
      }
      
      # Refine the center using the median matched day-of-year.
      refined_center <- median(best_matches$.DOY, na.rm = TRUE)
      chosen_centers <- c(chosen_centers, refined_center)
      
      # Remove matched rows before identifying the next reference period.
      remaining <- remaining %>%
        anti_join(
          best_matches %>% select(.Site, .Year, .Date),
          by = c(".Site", ".Year", ".Date")
        )
    }
    
    # Return one reference table per group.
    ref_tbl <- tibble(
      .Group = unique(gdat$.Group)[1],
      visit_id = seq_along(chosen_centers),
      center_doy = sort(chosen_centers),
      common_effort = k
    )
    
    ref_tbl
  }
  
  # Build reference sampling periods for each group.
  group_refs <- df %>%
    group_split(.Group) %>%
    map_dfr(get_group_reference)
  
  if (nrow(group_refs) == 0) {
    message("No group-level reference periods identified.")
    return(tibble())
  }
  
  # --------------------------------------------------
  # Step 2. Apply the group reference to each site
  # --------------------------------------------------
  #
  # For each site:
  # - use only the reference periods from its group
  # - keep only years with at least k samples
  # - match yearly samples to the group's reference periods
  # - retain only valid one-to-one matches within the time buffer
  process_one_site <- function(sdat) {
    
    this_group <- unique(sdat$.Group)
    ref_tbl <- group_refs %>%
      filter(.Group == this_group) %>%
      arrange(visit_id)
    
    if (nrow(ref_tbl) == 0) {
      return(NULL)
    }
    
    k <- unique(ref_tbl$common_effort)
    
    yearly_n <- sdat %>%
      count(.Year, name = "n_year")
    
    # Require the site to have data from at least min_years years.
    if (nrow(yearly_n) < min_years) {
      return(NULL)
    }
    
    # Keep only years with enough samples to match all reference visits.
    valid_years <- yearly_n %>%
      filter(n_year >= k) %>%
      pull(.Year)
    
    sdat2 <- sdat %>%
      filter(.Year %in% valid_years)
    
    if (length(unique(sdat2$.Year)) < min_years) {
      return(NULL)
    }
    
    out_years <- map(sort(unique(sdat2$.Year)), function(yy) {
      
      yy_dat <- sdat2 %>%
        filter(.Year == yy) %>%
        arrange(.Date)
      
      if (nrow(yy_dat) < k) {
        return(NULL)
      }
      
      # Compute distances between this year's sample dates and the
      # group reference centers.
      dist_mat <- outer(
        yy_dat$.DOY,
        ref_tbl$center_doy,
        function(a, b) abs(a - b)
      )
      
      # Create all sample-to-visit combinations and rank them by closeness.
      pairs <- expand.grid(
        sample_i = seq_len(nrow(yy_dat)),
        visit_id = seq_len(nrow(ref_tbl))
      ) %>%
        as_tibble() %>%
        mutate(dist_to_center = as.vector(dist_mat)) %>%
        arrange(dist_to_center)
      
      chosen_samples <- integer(0)
      chosen_visits <- integer(0)
      matched_rows <- list()
      
      # Greedy matching:
      # assign the closest available sample to each available visit
      # until all visits are matched or no valid pairs remain.
      for (i in seq_len(nrow(pairs))) {
        s <- pairs$sample_i[i]
        v <- pairs$visit_id[i]
        
        if (!(s %in% chosen_samples) && !(v %in% chosen_visits)) {
          chosen_samples <- c(chosen_samples, s)
          chosen_visits <- c(chosen_visits, v)
          
          matched_rows[[length(matched_rows) + 1]] <- tibble(
            sample_i = s,
            visit_id = v,
            dist_to_center = pairs$dist_to_center[i]
          )
        }
        
        if (length(chosen_visits) == k) break
      }
      
      matched <- bind_rows(matched_rows)
      
      if (nrow(matched) < k) {
        return(NULL)
      }
      
      out <- yy_dat %>%
        slice(matched$sample_i) %>%
        mutate(
          visit_id = matched$visit_id,
          dist_to_center = matched$dist_to_center
        ) %>%
        left_join(ref_tbl %>% select(visit_id, center_doy, common_effort),
                  by = "visit_id") %>%
        filter(dist_to_center <= buffer_days)
      
      # Require all reference visits to be represented after filtering.
      if (n_distinct(out$visit_id) < k) {
        return(NULL)
      }
      
      out %>%
        arrange(visit_id)
    })
    
    out <- bind_rows(out_years)
    
    if (nrow(out) == 0) {
      return(NULL)
    }
    
    out
  }
  
  # Apply site-level processing to all sites.
  result <- df %>%
    group_split(.Site) %>%
    map_dfr(process_one_site)
  
  if (nrow(result) == 0) {
    message("No rows retained after harmonization.")
    return(result)
  }
  
  # Restore original column names and remove internal helper columns.
  if (site_col == group_col) {
    
    result %>%
      select(-all_of(site_col), -all_of(date_col)) %>%
      rename(
        !!site_col := .Site,
        !!date_col := .Date
      ) %>%
      select(-any_of(c(".DOY", ".Group")))
    
  } else {
    
    result %>%
      select(-all_of(site_col), -all_of(group_col), -all_of(date_col)) %>%
      rename(
        !!site_col := .Site,
        !!group_col := .Group,
        !!date_col := .Date
      ) %>%
      select(-any_of(c(".DOY")))
  }
}