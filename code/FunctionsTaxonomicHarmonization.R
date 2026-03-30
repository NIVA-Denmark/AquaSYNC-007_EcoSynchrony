#' Harmonize taxonomic data across observations
#'
#' Resolves inconsistencies between coarse and fine taxonomic identifications
#' within user-defined groups by reconciling taxonomic ranks and aggregating
#' abundance accordingly.
#'
#' @param df A data frame containing taxonomy and abundance data.
#' @param decision_cols Column(s) used to define groups for harmonization decisions.
#' @param output_site_cols Column(s) used for grouping in the final output.
#' @param date_col Name of the date column.
#' @param abundance_col Name of the abundance column.
#' @param threshold Numeric value between 0 and 1, or NULL. If provided,
#'   finer taxa are retained when they contribute at least this proportion
#'   of total abundance within a lineage.
#'
#' @return A data frame with harmonized taxonomy and aggregated abundance.
#'
#' @examples
#' result <- harmonize_taxonomy(
#'   df = Data_taxonomy,
#'   decision_cols = "SiteID",   # or "Country", etc.
#'   output_site_cols = "SiteID",
#'   date_col = "Date",
#'   abundance_col = "Value",
#'   threshold = 0.90
#' )
#'
#' @export


library(dplyr)
library(tibble)

# Taxonomic ranks ordered from broadest to most specific.
rank_cols <- c("phylum", "class", "subclass", "order", "family", "genus", "species")

# ============================================================
# HELPERS
# ============================================================

# Clean taxonomy columns by:
# - converting values to character
# - trimming extra whitespace
# - converting empty strings and common missing-value strings to NA
clean_taxonomy <- function(df, rank_cols) {
  df %>%
    mutate(
      across(
        all_of(rank_cols),
        ~{
          x <- as.character(.x)
          x <- trimws(x)
          x[x %in% c("", "NA", "Na", "na", "NULL", "null")] <- NA_character_
          x
        }
      )
    )
}

# For each row, return the lowest (most specific) taxonomic rank
# that has a non-missing value.
lowest_rank <- function(df, rank_cols) {
  x <- as.data.frame(df[, rank_cols, drop = FALSE])
  
  apply(x, 1, function(z) {
    idx <- which(!is.na(z) & z != "")
    if (length(idx) == 0) NA_character_ else rank_cols[max(idx)]
  })
}

# For each row, return the name at the lowest (most specific)
# non-missing taxonomic rank.
lowest_name <- function(df, rank_cols) {
  x <- as.data.frame(df[, rank_cols, drop = FALSE])
  
  apply(x, 1, function(z) {
    idx <- which(!is.na(z) & z != "")
    if (length(idx) == 0) NA_character_ else as.character(z[max(idx)])
  })
}

# Collapse duplicate taxa by summing abundance within:
# - the grouping columns
# - the full taxonomic lineage
#
# After collapsing, add:
# - current_rank: deepest non-missing rank in each row
# - current_name: name at that deepest rank
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

# Build a table of observed lineages from the data.
# This is later used to prevent the cleanup step from creating
# lineage combinations that were never actually observed.
get_observed_lineages <- function(df, group_cols, abundance_col = "Value") {
  df %>%
    collapse_taxa(group_cols = group_cols, abundance_col = abundance_col) %>%
    transmute(
      across(all_of(c(group_cols, rank_cols))),
      observed_rank = current_rank
    )
}

# Check whether a coarse lineage exists among the originally observed lineages.
#
# This is used during cleanup to decide whether a coarse taxon
# should be compared against finer taxa beneath it.
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

# ============================================================
# STEP 1: Harmonize each parent rank against ALL finer ranks
# ============================================================
#
# For a given parent rank (for example, family), compare rows identified
# at that rank against rows identified at any finer rank below it
# (for example, genus or species) within the same lineage.
#
# Rules:
# - If both parent-level and finer-level taxa are present:
#   - collapse everything to the parent rank
#   - OR, if a threshold is provided and finer abundance reaches that
#     threshold proportion, keep finer taxa and drop the parent row
# - Otherwise, keep the data unchanged
harmonize_parent_vs_all_finer <- function(df,
                                          decision_cols,
                                          parent,
                                          abundance_col = "Value",
                                          threshold = NULL) {
  
  rank_index <- setNames(seq_along(rank_cols), rank_cols)
  parent_idx <- rank_index[[parent]]
  
  key_cols <- c(decision_cols, rank_cols[1:parent_idx])
  
  work <- df %>%
    filter(!is.na(.data[[parent]]))
  
  untouched <- df %>%
    filter(is.na(.data[[parent]]))
  
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
          has_parent & has_finer & !is.null(threshold) & prop_finer >= threshold ~ "keep_finer_drop_parent",
          has_parent & has_finer ~ "collapse_to_parent",
          TRUE ~ "no_change"
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

# ============================================================
# STEP 2: Clean up cases with missing intermediate ranks
# Also compare a coarse rank against ALL finer ranks below it
# ============================================================
#
# This step repeatedly checks whether coarse taxa should remain coarse
# or be reconciled with finer taxa beneath them, especially when some
# intermediate ranks are missing.
#
# Rules:
# - If finer taxa dominate above the threshold, drop the coarse row
# - Otherwise, collapse finer rows upward to the coarse rank
# - Repeat until no more changes occur
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
      
      coarse_rows <- df %>%
        filter(current_rank == coarse_rank)
      
      if (nrow(coarse_rows) == 0) next
      
      keep_row <- rep(TRUE, nrow(df))
      
      for (i in seq_len(nrow(coarse_rows))) {
        cr <- coarse_rows[i, ]
        
        if (!lineage_exists_in_observed(
          cr = cr,
          observed_lineages = observed_lineages,
          group_cols = decision_cols,
          coarse_idx = coarse_idx,
          coarse_rank = coarse_rank
        )) {
          next
        }
        
        same_group <- rep(TRUE, nrow(df))
        for (g in decision_cols) {
          same_group <- same_group & (df[[g]] == cr[[g]])
        }
        
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
        
        if (!is.null(threshold) && total_abund > 0 && (finer_abund / total_abund) >= threshold) {
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

# ============================================================
# STEP 3: Map each original row to the deepest compatible
# retained lineage in decision_harm
# ============================================================
#
# After the harmonized taxonomy has been determined at the decision-group
# level, this step maps each original row back to the best retained lineage.
#
# Rules:
# - Keep only harmonized taxa in the same decision group
# - Keep only harmonized taxa compatible with the row's lineage
# - Prefer the deepest retained lineage that is the same rank or coarser
# - If only finer retained taxa exist below a coarse original row, drop that row
apply_decision_harmonization_to_rows <- function(df,
                                                 decision_harm,
                                                 decision_cols,
                                                 abundance_col = "Value") {
  
  rank_index <- setNames(seq_along(rank_cols), rank_cols)
  
  out <- df %>%
    mutate(
      row_id = row_number(),
      original_rank = lowest_rank(., rank_cols)
    )
  
  decision_harm2 <- decision_harm %>%
    mutate(final_rank_index = rank_index[current_rank])
  
  mapped_list <- vector("list", nrow(out))
  
  for (i in seq_len(nrow(out))) {
    rw <- out[i, , drop = FALSE]
    
    rw_rank <- rw$original_rank[1]
    rw_rank_idx <- rank_index[[rw_rank]]
    
    cand <- decision_harm2
    
    # Restrict to the same decision group.
    for (g in decision_cols) {
      cand <- cand[cand[[g]] == rw[[g]][1], , drop = FALSE]
    }
    
    if (nrow(cand) == 0) {
      mapped_list[[i]] <- rw
      next
    }
    
    # Keep only retained taxa compatible with the row lineage.
    for (rc in rank_cols) {
      row_val <- rw[[rc]][1]
      
      if (!is.na(row_val)) {
        cand <- cand[
          is.na(cand[[rc]]) | cand[[rc]] == row_val,
          ,
          drop = FALSE
        ]
      }
    }
    
    if (nrow(cand) == 0) {
      mapped_list[[i]] <- rw
      next
    }
    
    cand_same_or_coarser <- cand[cand$final_rank_index <= rw_rank_idx, , drop = FALSE]
    
    if (nrow(cand_same_or_coarser) > 0) {
      best_idx <- which.max(cand_same_or_coarser$final_rank_index)
      best <- cand_same_or_coarser[best_idx, , drop = FALSE]
      
      rw[, rank_cols] <- best[, rank_cols]
      mapped_list[[i]] <- rw
      next
    }
    
    # If only finer retained taxa exist under this coarse row,
    # the row is not compatible with the retained output and is dropped.
    mapped_list[[i]] <- NULL
  }
  
  bind_rows(mapped_list) %>%
    select(-row_id, -original_rank)
}

# ============================================================
# MAIN FUNCTION
# ============================================================
#
# Harmonize taxonomy across observations while preserving output grouping.
#
# Arguments:
# - df: input data frame
# - decision_cols: columns used to decide harmonization rules
# - output_site_cols: columns used in the final grouped output
# - date_col: date column preserved in output
# - abundance_col: abundance/count column
# - threshold:
#     * NULL: always collapse mixed coarse/fine taxa to the coarser rank
#     * numeric in [0, 1]: if finer taxa contribute at least this proportion
#       of abundance, keep finer taxa and drop the coarse parent row
#
# Returns:
# A harmonized data frame aggregated by output_site_cols, date_col,
# and taxonomic lineage, with:
# - harmonized_rank
# - harmonized_taxon
harmonize_taxonomy <- function(df,
                               decision_cols = "SiteID",
                               output_site_cols = "SiteID",
                               date_col = "Date",
                               abundance_col = "Value",
                               threshold = NULL) {
  
  required_cols <- unique(c(decision_cols, output_site_cols, date_col, rank_cols, abundance_col))
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Standardize taxonomy values before harmonization.
  df_clean <- df %>%
    clean_taxonomy(rank_cols)
  
  # Process from more specific parent levels upward toward broader levels.
  parent_ranks <- c("genus", "family", "order", "subclass", "class", "phylum")
  
  # Save the set of originally observed lineages for later validation.
  observed_lineages <- df_clean %>%
    select(all_of(c(decision_cols, rank_cols, abundance_col))) %>%
    get_observed_lineages(group_cols = decision_cols, abundance_col = abundance_col)
  
  # Start from collapsed taxonomy within each decision group.
  decision_harm <- df_clean %>%
    select(all_of(c(decision_cols, rank_cols, abundance_col))) %>%
    collapse_taxa(group_cols = decision_cols, abundance_col = abundance_col)
  
  # Step 1: harmonize each parent rank against all finer ranks below it.
  for (parent in parent_ranks) {
    decision_harm <- harmonize_parent_vs_all_finer(
      df = decision_harm,
      decision_cols = decision_cols,
      parent = parent,
      abundance_col = abundance_col,
      threshold = threshold
    )
  }
  
  # Step 2: clean up remaining inconsistencies involving missing intermediate ranks.
  decision_harm <- cleanup_missing_intermediates_decision(
    df = decision_harm,
    decision_cols = decision_cols,
    observed_lineages = observed_lineages,
    abundance_col = abundance_col,
    threshold = threshold
  )
  
  # Step 3: map original rows back onto the retained harmonized taxonomy.
  mapped_rows <- apply_decision_harmonization_to_rows(
    df = df_clean %>%
      select(all_of(c(decision_cols, output_site_cols, date_col, rank_cols, abundance_col))),
    decision_harm = decision_harm,
    decision_cols = decision_cols,
    abundance_col = abundance_col
  )
  
  # Final aggregation and labeling.
  out <- mapped_rows %>%
    group_by(across(all_of(c(output_site_cols, date_col, rank_cols)))) %>%
    summarise(
      !!abundance_col := sum(.data[[abundance_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      harmonized_rank = lowest_rank(., rank_cols),
      harmonized_taxon = lowest_name(., rank_cols)
    ) %>%
    arrange(
      across(all_of(c(output_site_cols, date_col))),
      phylum, class, subclass, order, family, genus, species
    )
  
  out
}