library(dplyr)
library(tibble)
library(tidyr)
library(rlang)

rank_cols <- c(
  "phylum", "class", "subclass", "order",
  "family", "subfamily", "genus", "species"
)

clean_taxonomy <- function(df, rank_cols) {
  df %>%
    mutate(
      across(
        all_of(rank_cols),
        ~ {
          x <- as.character(.x)
          x <- trimws(x)
          x[x %in% c("", "NA", "Na", "na", "NULL", "null")] <- NA_character_
          x
        }
      )
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

collapse_taxa <- function(df, group_cols, abundance_col = "Value") {
  abundance_sym <- rlang::sym(abundance_col)
  
  df %>%
    group_by(across(all_of(c(group_cols, rank_cols)))) %>%
    summarise(
      !!abundance_col := sum(!!abundance_sym, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      current_rank = lowest_rank(., rank_cols),
      current_name = lowest_name(., rank_cols)
    )
}

get_observed_lineages <- function(df, group_cols, abundance_col = "Value") {
  df %>%
    collapse_taxa(group_cols = group_cols, abundance_col = abundance_col) %>%
    transmute(
      across(all_of(c(group_cols, rank_cols))),
      observed_rank = current_rank
    )
}

lineage_exists_in_observed <- function(cr,
                                       observed_lineages,
                                       group_cols,
                                       coarse_idx,
                                       coarse_rank) {
  same <- rep(TRUE, nrow(observed_lineages))
  
  for (g in group_cols) {
    same <- same & (observed_lineages[[g]] == cr[[g]])
  }
  
  same <- same & (observed_lineages$observed_rank == coarse_rank)
  
  for (rc in rank_cols[1:coarse_idx]) {
    val <- cr[[rc]][1]
    if (!is.na(val)) {
      same <- same & (observed_lineages[[rc]] == val)
    }
  }
  
  any(same, na.rm = TRUE)
}

harmonize_parent_vs_all_finer <- function(df,
                                          decision_cols,
                                          parent,
                                          abundance_col = "Value",
                                          threshold = NULL) {
  
  rank_index <- setNames(seq_along(rank_cols), rank_cols)
  parent_idx <- rank_index[[parent]]
  
  key_cols <- c(decision_cols, rank_cols[1:parent_idx])
  
  work <- df %>% filter(!is.na(.data[[parent]]))
  untouched <- df %>% filter(is.na(.data[[parent]]))
  
  if (nrow(work) == 0) return(df)
  
  decisions <- work %>%
    group_by(across(all_of(key_cols))) %>%
    reframe(
      {
        current_rank_idx <- rank_index[current_rank]
        
        parent_abund <- sum(.data[[abundance_col]][current_rank == parent], na.rm = TRUE)
        finer_abund  <- sum(.data[[abundance_col]][current_rank_idx > parent_idx], na.rm = TRUE)
        
        has_parent <- any(current_rank == parent, na.rm = TRUE)
        has_finer  <- any(current_rank_idx > parent_idx, na.rm = TRUE)
        
        total_pf <- parent_abund + finer_abund
        prop_finer <- if (total_pf > 0) finer_abund / total_pf else NA_real_
        
        action <- case_when(
          has_parent & has_finer & !is.null(threshold) & prop_finer >= threshold ~
            "keep_finer_drop_parent",
          has_parent & has_finer ~
            "collapse_to_parent",
          TRUE ~
            "no_change"
        )
        
        tibble(action = action)
      }
    ) %>%
    ungroup()
  
  work2 <- work %>%
    left_join(decisions, by = key_cols)
  
  deeper_cols <- rank_cols[(parent_idx + 1):length(rank_cols)]
  if (parent_idx == length(rank_cols)) deeper_cols <- character(0)
  
  if (length(deeper_cols) > 0) {
    work2 <- work2 %>%
      mutate(
        across(
          all_of(deeper_cols),
          ~ ifelse(action == "collapse_to_parent", NA_character_, .x)
        )
      )
  }
  
  work2 <- work2 %>%
    filter(!(action == "keep_finer_drop_parent" & current_rank == parent)) %>%
    select(-action)
  
  bind_rows(work2, untouched) %>%
    collapse_taxa(group_cols = decision_cols, abundance_col = abundance_col)
}

cleanup_missing_intermediates_decision <- function(df,
                                                   decision_cols,
                                                   observed_lineages,
                                                   abundance_col = "Value",
                                                   threshold = NULL) {
  
  rank_index <- setNames(seq_along(rank_cols), rank_cols)
  
  repeat {
    changed <- FALSE
    
    for (coarse_rank in rank_cols[1:(length(rank_cols) - 1)]) {
      
      coarse_idx <- rank_index[[coarse_rank]]
      deeper_cols <- rank_cols[(coarse_idx + 1):length(rank_cols)]
      
      df <- collapse_taxa(df, group_cols = decision_cols, abundance_col = abundance_col)
      
      coarse_rows <- df %>% filter(current_rank == coarse_rank)
      if (nrow(coarse_rows) == 0) next
      
      keep_row <- rep(TRUE, nrow(df))
      
      for (i in seq_len(nrow(coarse_rows))) {
        cr <- coarse_rows[i, ]
        
        if (!lineage_exists_in_observed(
          cr, observed_lineages,
          decision_cols, coarse_idx, coarse_rank
        )) next
        
        same_group <- Reduce(`&`, lapply(decision_cols, function(g) df[[g]] == cr[[g]]))
        
        match_lineage <- same_group
        for (rc in rank_cols[1:coarse_idx]) {
          val <- cr[[rc]][1]
          if (!is.na(val)) {
            match_lineage <- match_lineage & (df[[rc]] == val)
          }
        }
        
        candidates <- which(match_lineage)
        if (length(candidates) <= 1) next
        
        cand_rank_idx <- rank_index[df$current_rank[candidates]]
        finer <- candidates[cand_rank_idx > coarse_idx]
        coarse_same <- candidates[cand_rank_idx == coarse_idx]
        
        if (length(finer) == 0 || length(coarse_same) == 0) next
        
        coarse_abund <- sum(df[[abundance_col]][coarse_same], na.rm = TRUE)
        finer_abund  <- sum(df[[abundance_col]][finer], na.rm = TRUE)
        total_abund  <- coarse_abund + finer_abund
        
        if (!is.null(threshold) &&
            total_abund > 0 &&
            (finer_abund / total_abund) >= threshold) {
          
          keep_row[coarse_same] <- FALSE
          changed <- TRUE
          
        } else {
          if (length(deeper_cols) > 0) {
            df[finer, deeper_cols] <- NA_character_
          }
          changed <- TRUE
        }
      }
      
      df <- df[keep_row, , drop = FALSE] %>%
        collapse_taxa(group_cols = decision_cols, abundance_col = abundance_col)
    }
    
    if (!changed) break
  }
  
  df
}

apply_decision_harmonization_to_rows <- function(df,
                                                 decision_harm,
                                                 decision_cols,
                                                 abundance_col = "Value") {
  
  rank_index <- setNames(seq_along(rank_cols), rank_cols)
  
  out <- df %>%
    mutate(
      row_id = row_number(),
      original_rank = lowest_rank(., rank_cols),
      flag_coarse_taxon = FALSE
    )
  
  decision_harm2 <- decision_harm %>%
    mutate(final_rank_index = rank_index[current_rank])
  
  mapped_list <- vector("list", nrow(out))
  
  for (i in seq_len(nrow(out))) {
    rw <- out[i, , drop = FALSE]
    rw_rank_idx <- rank_index[[rw$original_rank[1]]]
    
    cand <- decision_harm2
    
    for (g in decision_cols) {
      cand <- cand[cand[[g]] == rw[[g]][1], , drop = FALSE]
    }
    
    for (rc in rank_cols) {
      val <- rw[[rc]][1]
      if (!is.na(val)) {
        cand <- cand[is.na(cand[[rc]]) | cand[[rc]] == val, , drop = FALSE]
      }
    }
    
    cand_same_or_coarser <- cand[cand$final_rank_index <= rw_rank_idx, , drop = FALSE]
    
    if (nrow(cand_same_or_coarser) > 0) {
      best <- cand_same_or_coarser[which.max(cand_same_or_coarser$final_rank_index), ]
      rw[, rank_cols] <- best[, rank_cols]
    } else {
      rw$flag_coarse_taxon <- TRUE
    }
    
    mapped_list[[i]] <- rw
  }
  
  bind_rows(mapped_list) %>%
    select(-row_id, -original_rank)
}

add_flag_no_finer <- function(out, output_site_cols, date_col) {
  
  out_with_id <- out %>%
    mutate(.row_id_internal = row_number())
  
  flag_lookup <- out_with_id %>%
    group_by(across(all_of(c(output_site_cols, date_col)))) %>%
    group_modify(~ {
      
      df <- .x
      df$rank_idx <- match(df$harmonized_rank, rank_cols)
      df$flag_no_finer <- FALSE
      
      flagged_rows <- which(df$flag_coarse_taxon)
      
      for (i in flagged_rows) {
        this_idx <- df$rank_idx[i]
        if (is.na(this_idx)) next
        
        possible_finer <- which(df$rank_idx > this_idx)
        
        if (length(possible_finer) == 0) {
          df$flag_no_finer[i] <- TRUE
          next
        }
        
        lineage_cols <- rank_cols[1:this_idx]
        same_lineage <- rep(TRUE, length(possible_finer))
        
        for (rc in lineage_cols) {
          this_val <- df[[rc]][i]
          if (!is.na(this_val)) {
            same_lineage <- same_lineage & df[[rc]][possible_finer] == this_val
          }
        }
        
        finer_taxa <- df$harmonized_taxon[possible_finer][same_lineage]
        df$flag_no_finer[i] <- dplyr::n_distinct(finer_taxa, na.rm = TRUE) == 0L
      }
      
      df %>% select(.row_id_internal, flag_no_finer)
    }) %>%
    ungroup() %>%
    select(.row_id_internal, flag_no_finer)
  
  out_with_id %>%
    left_join(flag_lookup, by = ".row_id_internal") %>%
    mutate(
      flag_coarse_taxon = as.logical(flag_coarse_taxon),
      flag_no_finer = if_else(flag_coarse_taxon, flag_no_finer, FALSE)
    ) %>%
    select(-.row_id_internal)
}

add_n_finer_fast <- function(out, output_site_cols) {
  
  out2 <- out %>%
    mutate(
      .row_id_internal = row_number(),
      rank_idx = match(harmonized_rank, rank_cols)
    )
  
  flagged <- out2 %>%
    filter(flag_coarse_taxon) %>%
    select(
      .row_id_internal,
      all_of(output_site_cols),
      rank_idx,
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
      all_of(output_site_cols),
      finer_rank_idx = rank_idx,
      finer_taxon = harmonized_taxon,
      all_of(rank_cols)
    )
  
  n_finer <- flagged %>%
    left_join(
      finer,
      by = output_site_cols,
      suffix = c("_coarse", "_finer"),
      relationship = "many-to-many"
    ) %>%
    filter(finer_rank_idx > rank_idx)
  
  for (rc in rank_cols) {
    coarse_col <- paste0(rc, "_coarse")
    finer_col  <- paste0(rc, "_finer")
    
    n_finer <- n_finer %>%
      filter(is.na(.data[[coarse_col]]) | .data[[coarse_col]] == .data[[finer_col]])
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
      flag_coarse_taxon = as.logical(flag_coarse_taxon),
      n_finer = if_else(
        flag_coarse_taxon,
        replace_na(n_finer, 0L),
        NA_integer_
      )
    ) %>%
    select(-.row_id_internal, -rank_idx)
}

harmonize_taxonomy <- function(df,
                               decision_cols = "SiteID",
                               output_site_cols = "SiteID",
                               date_col = "Date",
                               abundance_col = "Value",
                               threshold = NULL) {
  
  df_clean <- df %>%
    clean_taxonomy(rank_cols)
  
  observed_lineages <- df_clean %>%
    select(all_of(c(decision_cols, rank_cols, abundance_col))) %>%
    get_observed_lineages(group_cols = decision_cols)
  
  decision_harm <- df_clean %>%
    select(all_of(c(decision_cols, rank_cols, abundance_col))) %>%
    collapse_taxa(group_cols = decision_cols)
  
  for (parent in c("genus", "subfamily", "family", "order", "subclass", "class", "phylum")) {
    decision_harm <- harmonize_parent_vs_all_finer(
      decision_harm,
      decision_cols,
      parent,
      abundance_col,
      threshold
    )
  }
  
  decision_harm <- cleanup_missing_intermediates_decision(
    decision_harm,
    decision_cols,
    observed_lineages,
    abundance_col,
    threshold
  )
  
  mapped_rows <- apply_decision_harmonization_to_rows(
    df_clean %>%
      select(all_of(c(
        decision_cols, output_site_cols, date_col,
        rank_cols, abundance_col
      ))),
    decision_harm,
    decision_cols,
    abundance_col
  )
  
  out <- mapped_rows %>%
    group_by(across(all_of(c(output_site_cols, date_col, rank_cols)))) %>%
    summarise(
      !!abundance_col := sum(.data[[abundance_col]], na.rm = TRUE),
      flag_coarse_taxon = any(flag_coarse_taxon, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      harmonized_rank = lowest_rank(., rank_cols),
      harmonized_taxon = lowest_name(., rank_cols)
    ) %>%
    add_flag_no_finer(output_site_cols, date_col) %>%
    add_n_finer_fast(output_site_cols) %>%
    arrange(across(all_of(c(output_site_cols, date_col)))) %>%
    select(
      all_of(output_site_cols),
      all_of(date_col),
      all_of(rank_cols),
      all_of(abundance_col),
      harmonized_rank,
      harmonized_taxon,
      flag_coarse_taxon,
      flag_no_finer,
      n_finer
    )
  
  out
}