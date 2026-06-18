library(dplyr)
library(tibble)
library(tidyr)
library(rlang)

rank_cols <- c(
  "phylum", "subphylum", "class", "subclass", "order",
  "family", "subfamily", "genus", "species"
)

parent_order <- c(
  "genus", "subfamily", "family", "order",
  "subclass", "class", "subphylum", "phylum"
)

clean_taxonomy <- function(df, rank_cols) {
  df %>%
    mutate(
      across(all_of(rank_cols), ~ {
        x <- as.character(.x)
        x <- trimws(x)
        x[x %in% c("", "NA", "Na", "na", "NULL", "null")] <- NA_character_
        x
      })
    )
}

lowest_rank <- function(df, rank_cols) {
  x <- as.data.frame(df[, rank_cols, drop = FALSE])
  apply(x, 1, function(z) {
    idx <- which(!is.na(z) & z != "")
    if (length(idx) == 0) NA_character_ else rank_cols[max(idx)]
  })
}

lowest_name <- function(df, rank_cols) {
  x <- as.data.frame(df[, rank_cols, drop = FALSE])
  apply(x, 1, function(z) {
    idx <- which(!is.na(z) & z != "")
    if (length(idx) == 0) NA_character_ else as.character(z[max(idx)])
  })
}

refresh_rank <- function(df) {
  df %>%
    mutate(
      current_rank = lowest_rank(., rank_cols),
      current_name = lowest_name(., rank_cols)
    )
}

collapse_taxa_flags <- function(df, group_cols, abundance_col = "Value") {
  abundance_sym <- sym(abundance_col)
  
  df %>%
    group_by(across(all_of(c(group_cols, rank_cols)))) %>%
    summarise(
      !!abundance_col := sum(!!abundance_sym, na.rm = TRUE),
      flag_sample_coarse_overridden = any(flag_sample_coarse_overridden, na.rm = TRUE),
      flag_site_coarse_overridden = any(flag_site_coarse_overridden, na.rm = TRUE),
      flag_final_rollup = any(flag_final_rollup, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    refresh_rank()
}

add_year_column <- function(df, date_col, year_col = NULL) {
  if (!is.null(year_col)) {
    return(df %>% mutate(.harm_year = .data[[year_col]]))
  }
  
  if (inherits(df[[date_col]], "Date") || inherits(df[[date_col]], "POSIXt")) {
    df %>% mutate(.harm_year = as.integer(format(.data[[date_col]], "%Y")))
  } else {
    df %>%
      mutate(
        .harm_year = suppressWarnings(
          as.integer(format(as.Date(.data[[date_col]]), "%Y"))
        )
      )
  }
}

harmonize_parent_abundance_flag <- function(df,
                                            decision_cols,
                                            parent,
                                            abundance_col = "Value",
                                            stage1_finer_abundance_threshold = 0.9) {
  
  rank_index <- setNames(seq_along(rank_cols), rank_cols)
  parent_idx <- rank_index[[parent]]
  key_cols <- c(decision_cols, rank_cols[1:parent_idx])
  
  df <- refresh_rank(df)
  
  work <- df %>% filter(!is.na(.data[[parent]]))
  untouched <- df %>% filter(is.na(.data[[parent]]))
  
  if (nrow(work) == 0) return(df)
  
  decisions <- work %>%
    group_by(across(all_of(key_cols))) %>%
    reframe({
      current_rank_idx <- rank_index[current_rank]
      
      is_parent <- current_rank == parent
      is_finer <- current_rank_idx > parent_idx
      
      parent_abund <- sum(.data[[abundance_col]][is_parent], na.rm = TRUE)
      finer_abund <- sum(.data[[abundance_col]][is_finer], na.rm = TRUE)
      
      has_parent <- any(is_parent, na.rm = TRUE)
      has_finer <- any(is_finer, na.rm = TRUE)
      
      total_abund <- parent_abund + finer_abund
      prop_finer <- if (total_abund > 0) finer_abund / total_abund else NA_real_
      
      action <- case_when(
        has_parent & has_finer &
          !is.na(prop_finer) &
          prop_finer >= stage1_finer_abundance_threshold ~
          "keep_finer_flag_parent",
        
        has_parent & has_finer &
          !is.na(prop_finer) &
          prop_finer < stage1_finer_abundance_threshold ~
          "collapse_to_parent",
        
        TRUE ~
          "no_change"
      )
      
      tibble(action = action)
    }) %>%
    ungroup()
  
  # Diptera rule: never roll finer Diptera taxa up to order level.
  if (parent == "order") {
    decisions <- decisions %>%
      mutate(
        action = if_else(
          order == "Diptera" & action == "collapse_to_parent",
          "keep_finer_flag_parent",
          action
        )
      )
  }
  
  work2 <- work %>%
    left_join(decisions, by = key_cols) %>%
    mutate(
      current_rank_idx = rank_index[current_rank],
      flag_sample_coarse_overridden =
        flag_sample_coarse_overridden |
        (action == "keep_finer_flag_parent" & current_rank == parent)
    )
  
  deeper_cols <- if (parent_idx < length(rank_cols)) {
    rank_cols[(parent_idx + 1):length(rank_cols)]
  } else {
    character(0)
  }
  
  if (length(deeper_cols) > 0) {
    work2 <- work2 %>%
      mutate(
        across(
          all_of(deeper_cols),
          ~ ifelse(
            action == "collapse_to_parent" & current_rank_idx > parent_idx,
            NA_character_,
            .x
          )
        )
      )
  }
  
  bind_rows(
    work2 %>% select(-action, -current_rank_idx),
    untouched
  ) %>%
    collapse_taxa_flags(
      group_cols = decision_cols,
      abundance_col = abundance_col
    )
}

harmonize_sample_level <- function(df,
                                   sample_cols,
                                   abundance_col = "Value",
                                   stage1_finer_abundance_threshold = 0.9) {
  
  out <- df %>%
    mutate(
      flag_sample_coarse_overridden = FALSE,
      flag_site_coarse_overridden = FALSE,
      flag_final_rollup = FALSE
    ) %>%
    collapse_taxa_flags(
      group_cols = sample_cols,
      abundance_col = abundance_col
    )
  
  for (parent in parent_order) {
    out <- harmonize_parent_abundance_flag(
      out,
      decision_cols = sample_cols,
      parent = parent,
      abundance_col = abundance_col,
      stage1_finer_abundance_threshold = stage1_finer_abundance_threshold
    )
  }
  
  out
}

harmonize_parent_year_frequency_flag <- function(df,
                                                 site_cols,
                                                 date_col,
                                                 parent,
                                                 abundance_col = "Value",
                                                 stage2_finer_abundance_threshold = 0.9,
                                                 stage2_finer_year_threshold = 0.9) {
  
  rank_index <- setNames(seq_along(rank_cols), rank_cols)
  parent_idx <- rank_index[[parent]]
  key_cols <- c(site_cols, rank_cols[1:parent_idx])
  
  df <- refresh_rank(df)
  
  work <- df %>% filter(!is.na(.data[[parent]]))
  untouched <- df %>% filter(is.na(.data[[parent]]))
  
  if (nrow(work) == 0) return(df)
  
  decisions <- work %>%
    group_by(across(all_of(key_cols))) %>%
    reframe({
      current_rank_idx <- rank_index[current_rank]
      
      is_parent <- current_rank == parent & !flag_sample_coarse_overridden
      is_finer <- current_rank_idx > parent_idx
      
      parent_abund <- sum(.data[[abundance_col]][is_parent], na.rm = TRUE)
      finer_abund <- sum(.data[[abundance_col]][is_finer], na.rm = TRUE)
      total_abund <- parent_abund + finer_abund
      
      prop_finer_abund <- if (total_abund > 0) {
        finer_abund / total_abund
      } else {
        NA_real_
      }
      
      has_parent <- any(is_parent, na.rm = TRUE)
      has_finer <- any(is_finer, na.rm = TRUE)
      
      years_parent <- unique(.harm_year[is_parent])
      years_finer <- unique(.harm_year[is_finer])
      years_total <- union(years_parent, years_finer)
      
      prop_finer_years <- if (length(years_total) > 0) {
        length(years_finer) / length(years_total)
      } else {
        NA_real_
      }
      
      action <- case_when(
        !has_parent | !has_finer ~
          "no_change",
        
        parent %in% c("genus", "family") &
          !is.na(prop_finer_abund) &
          prop_finer_abund >= stage2_finer_abundance_threshold ~
          "keep_finer_flag_parent",
        
        parent %in% c("genus", "family") &
          !is.na(prop_finer_abund) &
          prop_finer_abund < stage2_finer_abundance_threshold ~
          "collapse_to_parent",
        
        !(parent %in% c("genus", "family")) &
          !is.na(prop_finer_years) &
          prop_finer_years >= stage2_finer_year_threshold ~
          "keep_finer_flag_parent",
        
        !(parent %in% c("genus", "family")) &
          !is.na(prop_finer_years) &
          prop_finer_years < stage2_finer_year_threshold ~
          "collapse_to_parent",
        
        TRUE ~
          "no_change"
      )
      
      tibble(action = action)
    }) %>%
    ungroup()
  
  # Diptera rule: never roll finer Diptera taxa up to order level.
  if (parent == "order") {
    decisions <- decisions %>%
      mutate(
        action = if_else(
          order == "Diptera" & action == "collapse_to_parent",
          "keep_finer_flag_parent",
          action
        )
      )
  }
  
  work2 <- work %>%
    left_join(decisions, by = key_cols) %>%
    mutate(
      current_rank_idx = rank_index[current_rank],
      
      flag_site_coarse_overridden =
        flag_site_coarse_overridden |
        (action == "keep_finer_flag_parent" & current_rank == parent),
      
      flag_final_rollup =
        flag_final_rollup |
        (action == "collapse_to_parent")
    )
  
  deeper_cols <- if (parent_idx < length(rank_cols)) {
    rank_cols[(parent_idx + 1):length(rank_cols)]
  } else {
    character(0)
  }
  
  if (length(deeper_cols) > 0) {
    work2 <- work2 %>%
      mutate(
        across(
          all_of(deeper_cols),
          ~ ifelse(
            action == "collapse_to_parent" & current_rank_idx > parent_idx,
            NA_character_,
            .x
          )
        )
      )
  }
  
  bind_rows(
    work2 %>% select(-action, -current_rank_idx),
    untouched
  ) %>%
    collapse_taxa_flags(
      group_cols = c(site_cols, date_col, ".harm_year"),
      abundance_col = abundance_col
    )
}

harmonize_site_year_frequency <- function(df,
                                          site_cols,
                                          date_col,
                                          abundance_col = "Value",
                                          stage2_finer_abundance_threshold = 0.9,
                                          stage2_finer_year_threshold = 0.9) {
  
  out <- df
  
  for (parent in parent_order) {
    out <- harmonize_parent_year_frequency_flag(
      out,
      site_cols = site_cols,
      date_col = date_col,
      parent = parent,
      abundance_col = abundance_col,
      stage2_finer_abundance_threshold = stage2_finer_abundance_threshold,
      stage2_finer_year_threshold = stage2_finer_year_threshold
    )
  }
  
  out
}

add_n_finer_site_level <- function(out, site_cols) {
  
  out2 <- out %>%
    mutate(
      .row_id_internal = row_number(),
      harmonized_rank = lowest_rank(., rank_cols),
      harmonized_taxon = lowest_name(., rank_cols),
      rank_idx = match(harmonized_rank, rank_cols),
      flag_coarse_overridden =
        flag_sample_coarse_overridden | flag_site_coarse_overridden,
      drop_candidate =
        flag_coarse_overridden & !flag_final_rollup
    )
  
  flagged <- out2 %>%
    filter(drop_candidate) %>%
    select(
      .row_id_internal,
      all_of(site_cols),
      coarse_rank_idx = rank_idx,
      all_of(rank_cols)
    )
  
  if (nrow(flagged) == 0) {
    return(
      out2 %>%
        mutate(n_finer = NA_integer_) %>%
        select(-.row_id_internal, -rank_idx)
    )
  }
  
  finer <- out2 %>%
    filter(!is.na(rank_idx)) %>%
    select(
      all_of(site_cols),
      finer_rank_idx = rank_idx,
      finer_taxon = harmonized_taxon,
      all_of(rank_cols)
    )
  
  n_finer <- flagged %>%
    left_join(
      finer,
      by = site_cols,
      suffix = c("_coarse", "_finer"),
      relationship = "many-to-many"
    ) %>%
    filter(finer_rank_idx > coarse_rank_idx)
  
  for (rc in rank_cols) {
    coarse_col <- paste0(rc, "_coarse")
    finer_col <- paste0(rc, "_finer")
    
    n_finer <- n_finer %>%
      filter(
        is.na(.data[[coarse_col]]) |
          .data[[coarse_col]] == .data[[finer_col]]
      )
  }
  
  n_finer <- n_finer %>%
    group_by(.row_id_internal) %>%
    summarise(
      n_finer = n_distinct(finer_taxon, na.rm = TRUE),
      .groups = "drop"
    )
  
  out2 %>%
    left_join(n_finer, by = ".row_id_internal") %>%
    mutate(
      n_finer = if_else(
        drop_candidate,
        replace_na(n_finer, 0L),
        NA_integer_
      )
    ) %>%
    select(-.row_id_internal, -rank_idx)
}

assign_unique_finer_to_drop_candidates <- function(out, site_cols) {
  
  out2 <- out %>%
    mutate(
      .row_id_internal = row_number(),
      rank_idx = match(harmonized_rank, rank_cols),
      flag_reassigned_unique_finer = FALSE
    )
  
  flagged <- out2 %>%
    filter(drop_candidate, n_finer == 1L) %>%
    select(
      .row_id_internal,
      all_of(site_cols),
      coarse_rank_idx = rank_idx,
      all_of(rank_cols)
    )
  
  if (nrow(flagged) == 0) {
    return(
      out2 %>%
        mutate(
          drop_candidate_final = drop_candidate
        ) %>%
        select(-.row_id_internal, -rank_idx)
    )
  }
  
  finer <- out2 %>%
    filter(!is.na(rank_idx)) %>%
    select(
      all_of(site_cols),
      finer_rank_idx = rank_idx,
      all_of(rank_cols)
    )
  
  candidates <- flagged %>%
    left_join(
      finer,
      by = site_cols,
      suffix = c("_coarse", "_finer"),
      relationship = "many-to-many"
    ) %>%
    filter(finer_rank_idx > coarse_rank_idx)
  
  for (rc in rank_cols) {
    coarse_col <- paste0(rc, "_coarse")
    finer_col <- paste0(rc, "_finer")
    
    candidates <- candidates %>%
      filter(
        is.na(.data[[coarse_col]]) |
          .data[[coarse_col]] == .data[[finer_col]]
      )
  }
  
  assignment <- candidates %>%
    group_by(.row_id_internal) %>%
    summarise(
      phylum_assigned = first(phylum_finer),
      subphylum_assigned = first(subphylum_finer),
      class_assigned = first(class_finer),
      subclass_assigned = first(subclass_finer),
      order_assigned = first(order_finer),
      family_assigned = first(family_finer),
      subfamily_assigned = first(subfamily_finer),
      genus_assigned = first(genus_finer),
      species_assigned = first(species_finer),
      .groups = "drop"
    )
  
  out2 %>%
    left_join(assignment, by = ".row_id_internal") %>%
    mutate(
      flag_reassigned_unique_finer = !is.na(phylum_assigned),
      
      phylum = coalesce(phylum_assigned, phylum),
      subphylum = coalesce(subphylum_assigned, subphylum),
      class = coalesce(class_assigned, class),
      subclass = coalesce(subclass_assigned, subclass),
      order = coalesce(order_assigned, order),
      family = coalesce(family_assigned, family),
      subfamily = coalesce(subfamily_assigned, subfamily),
      genus = coalesce(genus_assigned, genus),
      species = coalesce(species_assigned, species),
      
      drop_candidate_final =
        drop_candidate & !flag_reassigned_unique_finer
    ) %>%
    select(
      -ends_with("_assigned"),
      -.row_id_internal,
      -rank_idx
    ) %>%
    mutate(
      harmonized_rank = lowest_rank(., rank_cols),
      harmonized_taxon = lowest_name(., rank_cols)
    )
}

harmonize_taxonomy_two_stage <- function(df,
                                         site_cols = "SiteID",
                                         date_col = "Date",
                                         year_col = NULL,
                                         abundance_col = "Value",
                                         stage1_finer_abundance_threshold = 0.9,
                                         stage2_finer_abundance_threshold = 0.9,
                                         stage2_finer_year_threshold = 0.9) {
  
  sample_cols <- c(site_cols, date_col)
  
  df_clean <- df %>%
    clean_taxonomy(rank_cols) %>%
    refresh_rank()
  
  arthropoda_only <- df_clean %>%
    filter(
      current_rank == "phylum",
      current_name == "Arthropoda"
    )
  
  df_for_harmonization <- df_clean %>%
    filter(
      !(current_rank == "phylum" & current_name == "Arthropoda")
    )
  
  sample_harm <- df_for_harmonization %>%
    select(all_of(c(sample_cols, rank_cols, abundance_col))) %>%
    harmonize_sample_level(
      sample_cols = sample_cols,
      abundance_col = abundance_col,
      stage1_finer_abundance_threshold = stage1_finer_abundance_threshold
    ) %>%
    add_year_column(date_col = date_col, year_col = year_col)
  
  site_harm <- sample_harm %>%
    harmonize_site_year_frequency(
      site_cols = site_cols,
      date_col = date_col,
      abundance_col = abundance_col,
      stage2_finer_abundance_threshold = stage2_finer_abundance_threshold,
      stage2_finer_year_threshold = stage2_finer_year_threshold
    )
  
  final_harm <- site_harm %>%
    mutate(
      harmonized_rank = lowest_rank(., rank_cols),
      harmonized_taxon = lowest_name(., rank_cols),
      flag_coarse_overridden =
        flag_sample_coarse_overridden | flag_site_coarse_overridden,
      drop_candidate =
        flag_coarse_overridden & !flag_final_rollup
    ) %>%
    add_n_finer_site_level(site_cols = site_cols) %>%
    assign_unique_finer_to_drop_candidates(site_cols = site_cols)
  
  if (nrow(arthropoda_only) > 0) {
    arthropoda_final <- arthropoda_only %>%
      select(all_of(c(sample_cols, rank_cols, abundance_col))) %>%
      add_year_column(date_col = date_col, year_col = year_col) %>%
      mutate(
        harmonized_rank = lowest_rank(., rank_cols),
        harmonized_taxon = lowest_name(., rank_cols),
        flag_sample_coarse_overridden = FALSE,
        flag_site_coarse_overridden = FALSE,
        flag_final_rollup = FALSE,
        flag_coarse_overridden = TRUE,
        drop_candidate = TRUE,
        n_finer = NA_integer_,
        flag_reassigned_unique_finer = FALSE,
        drop_candidate_final = TRUE
      )
    
    final_harm <- bind_rows(final_harm, arthropoda_final)
  }
  
  final_harm %>%
    arrange(across(all_of(c(site_cols, date_col)))) %>%
    select(
      all_of(site_cols),
      all_of(date_col),
      .harm_year,
      all_of(rank_cols),
      all_of(abundance_col),
      harmonized_rank,
      harmonized_taxon,
      flag_sample_coarse_overridden,
      flag_site_coarse_overridden,
      flag_final_rollup,
      flag_coarse_overridden,
      drop_candidate,
      n_finer,
      flag_reassigned_unique_finer,
      drop_candidate_final
    )
}