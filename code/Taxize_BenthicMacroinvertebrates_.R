library(taxize)
library(arrow)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

# Step 1 -------------------------------------------------------------------

Data1 <- read_parquet("Combined_dataset/BenthicMacroinvertebrates.parquet")
Data2 <- read_parquet("Combined_dataset/BenthicMacroinvertebrates_part2.parquet")
Data  <- rbind(Data1, Data2)

# cleaning
TaxaList <- Data %>%
  select(Taxon) %>%
  distinct() %>%
  mutate(
    Taxon_clean = Taxon,
    
    # fix uppercase taxon names at the beginning of the string
    Taxon_clean = str_replace(
      as.character(Taxon_clean),
      "^[A-Z]+",
      ~ str_to_title(str_to_lower(.x))
    ),
    
    # capitalize names that are entirely lowercase
    Taxon_clean = if_else(
      Taxon_clean == tolower(Taxon_clean),
      str_to_title(Taxon_clean),
      Taxon_clean
    ),
    
    # remove dot between genus and species
    Taxon_clean = gsub("(?<=[A-Za-z])\\.(?=[A-Za-z])", " ", Taxon_clean, perl = TRUE),
    
    # specific corrections
    Taxon_clean = case_when(
      Taxon_clean == "Valvata." ~ "Valvata sp.",
      Taxon_clean == "Boreobdella verrucata" ~ "Glossiphonia verrucata",
      Taxon_clean == "Hydropsyhce angustipennis" ~ "Hydropsyche angustipennis",
      Taxon_clean == "Gyraulus." ~ "Gyraulus sp.",
      Taxon_clean == "Amphnemoura sulcicollis" ~ "Amphinemoura sulcicollis",
      Taxon_clean == "Philopotomidae" ~ "Philopotamidae",
      Taxon_clean == "Paraleutra" ~ "Paraleuctra",
      Taxon_clean == "Oligochaetae" ~ "Oligochaeta",
      Taxon_clean == "Elminthidae" ~ "Elmidae",
      Taxon_clean %in% c("Elmis aena Lv.", "Elmis aena Ad.", "Elmis aena lv.") ~ "Elmis aenea",
      TRUE ~ Taxon_clean
    )
  ) %>%
  distinct()

ToBeMatched <- TaxaList

# Verify names with Global Names
ver <- gna_verifier(ToBeMatched$Taxon_clean)

matched_names <- ToBeMatched %>%
  mutate(
    matched_name = ver$matchedCanonicalFull,
    matchType    = ver$matchType
  ) %>%
  mutate(
    matched_name = if_else(
      matchType %in% c("Exact", "Fuzzy"),
      trimws(matched_name),
      NA_character_
    )
  )

pick_id <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    NA_character_
  } else {
    as.character(x[[1]])
  }
}

pick_rank <- function(df, rank_name) {
  out <- df$name[tolower(df$rank) == tolower(rank_name)]
  if (length(out) == 0) NA_character_ else out[1]
}

# Query GBIF only once per unique matched_name
unique_names <- matched_names %>%
  distinct(matched_name) %>%
  filter(!is.na(matched_name)) %>%
  mutate(
    gbif_id = map(
      matched_name,
      \(nm) {
        tryCatch(
          get_gbifid(
            nm,
            method = "backbone",
            ask = FALSE,
            messages = FALSE
          ),
          error = function(e) NA
        )
      }
    )
  ) %>%
  mutate(
    gbif_id_chr = map_chr(gbif_id, pick_id)
  )

# Get GBIF classification once per unique matched_name
unique_names <- unique_names %>%
  mutate(
    classif = map(
      gbif_id,
      \(id) {
        if (is.null(id) || length(id) == 0 || all(is.na(id))) {
          return(NULL)
        }
        
        out <- tryCatch(
          classification(id, db = "gbif"),
          error = function(e) NULL
        )
        
        if (is.null(out) || length(out) == 0) {
          NULL
        } else if (is.list(out) && !is.data.frame(out)) {
          out[[1]]
        } else {
          out
        }
      }
    )
  )

# Build taxonomy lookup
tax_lookup <- map2_dfr(
  unique_names$classif,
  seq_len(nrow(unique_names)),
  \(x, i) {
    if (is.null(x) || !is.data.frame(x) || NROW(x) == 0) {
      return(tibble(
        matched_name  = unique_names$matched_name[i],
        gbif_id       = unique_names$gbif_id_chr[i],
        accepted_name = NA_character_,
        matched_rank  = NA_character_,
        species       = NA_character_,
        genus         = NA_character_,
        family        = NA_character_,
        order         = NA_character_,
        class         = NA_character_,
        phylum        = NA_character_
      ))
    }
    
    x <- x %>% mutate(rank = tolower(rank))
    
    tibble(
      matched_name  = unique_names$matched_name[i],
      gbif_id       = unique_names$gbif_id_chr[i],
      accepted_name = dplyr::last(x$name),
      matched_rank  = dplyr::last(x$rank),
      species       = pick_rank(x, "species"),
      genus         = pick_rank(x, "genus"),
      family        = pick_rank(x, "family"),
      order         = pick_rank(x, "order"),
      class         = pick_rank(x, "class"),
      phylum        = pick_rank(x, "phylum")
    )
  }
)

tax_table <- matched_names %>%
  left_join(tax_lookup, by = "matched_name")

# write.csv(tax_table, "TaxonomicHarm/BenthicMacroinvertebrate/Harmonization_table_GBIF.csv", row.names = FALSE)

# Unresolved - Step 2 ------------------------------------------------------

tax_table_final <- read.table(
  "TaxonomicHarm/BenthicMacroinvertebrate/Harmonization_table_GBIF.csv",
  h = TRUE, sep = ",", encoding = "UTF-8"
)

unresolved <- tax_table_final %>%
  filter(is.na(phylum)) 

tax_na <- unresolved %>%
  filter(is.na(gbif_id)) %>%
  distinct(Taxon, Taxon_clean, matched_name, matchType)

na_names <- tax_na %>%
  filter(!is.na(matched_name)) %>%
  distinct(matched_name) %>%
  pull(matched_name)

manual_ids <- vector("list", length(na_names))

for (i in seq_along(na_names)) {
  nm <- na_names[i]
  cat("\n", i, "/", length(na_names), ":", nm, "\n")
  
  manual_ids[[i]] <- tryCatch(
    get_gbifid(
      nm,
      method = "backbone",
      ask = TRUE,
      messages = FALSE
    ),
    error = function(e) {
      message("Failed for ", nm, ": ", conditionMessage(e))
      NA
    }
  )
}

manual_choices <- tibble(
  matched_name = na_names,
  gbif_id = vapply(
    manual_ids,
    function(x) {
      if (is.null(x) || length(x) == 0 || all(is.na(x))) {
        NA_character_
      } else {
        as.character(x[[1]])
      }
    },
    character(1)
  )
)

# specific manual decision in step 2
manual_choices <- manual_choices %>%
  mutate(gbif_id = if_else(matched_name == "Stagnicola", NA_character_, gbif_id))

manual_lookup <- dplyr::bind_rows(
  lapply(seq_len(nrow(manual_choices)), function(i) {
    id <- manual_choices$gbif_id[i]
    
    x <- if (is.na(id) || id == "") {
      NULL
    } else {
      tryCatch(
        {
          out <- classification(id, db = "gbif")
          if (is.null(out) || length(out) == 0) {
            NULL
          } else if (is.list(out) && !is.data.frame(out)) {
            out[[1]]
          } else {
            out
          }
        },
        error = function(e) NULL
      )
    }
    
    if (is.null(x) || !is.data.frame(x) || NROW(x) == 0) {
      return(tibble(
        matched_name  = manual_choices$matched_name[i],
        gbif_id       = manual_choices$gbif_id[i],
        accepted_name = NA_character_,
        matched_rank  = NA_character_,
        species       = NA_character_,
        genus         = NA_character_,
        family        = NA_character_,
        order         = NA_character_,
        class         = NA_character_,
        phylum        = NA_character_
      ))
    }
    
    x <- x %>% mutate(rank = tolower(rank))
    
    tibble(
      matched_name  = manual_choices$matched_name[i],
      gbif_id       = manual_choices$gbif_id[i],
      accepted_name = dplyr::last(x$name),
      matched_rank  = dplyr::last(x$rank),
      species       = pick_rank(x, "species"),
      genus         = pick_rank(x, "genus"),
      family        = pick_rank(x, "family"),
      order         = pick_rank(x, "order"),
      class         = pick_rank(x, "class"),
      phylum        = pick_rank(x, "phylum")
    )
  })
)

unresolved_2 <- unresolved %>%
  select(-gbif_id, -accepted_name, -matched_rank,
         -species, -genus, -family, -order, -class, -phylum) %>%
  left_join(manual_lookup, by = "matched_name")

resolved_step2 <- unresolved_2 %>%
  filter(!is.na(matched_rank)) %>%
  filter(matched_name != "Branchiopoda")

# write.csv(resolved_step2, "TaxonomicHarm/BenthicMacroinvertebrate/resolved_step2.csv", row.names = FALSE)

# Unresolved - Step 3 ------------------------------------------------------
# same manual GBIF workflow as Step 2
# only new part kept here:

resolved_step2 <- read.table("TaxonomicHarm/BenthicMacroinvertebrate/resolved_step2.csv", h = TRUE, sep = ",", encoding = "UTF-8")

unresolved2 <- unresolved %>%
  filter(!Taxon_clean %in% resolved_step2$Taxon_clean)

# run same manual GBIF resolution block as in Step 2 on unresolved2
# output: resolved_step3

# write.csv(resolved_step3, "TaxonomicHarm/BenthicMacroinvertebrate/resolved_step3.csv", row.names = FALSE)

# Unresolved - Step 4 ------------------------------------------------------
# same manual GBIF workflow as Step 2
# only new part kept here:

resolved_step3 <- read.table("TaxonomicHarm/BenthicMacroinvertebrate/resolved_step3.csv", h = TRUE, sep = ",", encoding = "UTF-8")

resolved_step2_3 <- rbind(resolved_step2, resolved_step3)

unresolved2 <- unresolved %>%
  filter(!Taxon_clean %in% resolved_step2_3$Taxon_clean)

unresolved2[unresolved2$Taxon == "Sialis scotti", "matched_name"] <- "Sialis scotti"

# run same manual GBIF resolution block as in Step 2 on unresolved2
# output: resolved_step4

# Unresolved - Step 5 ------------------------------------------------------
# same manual GBIF workflow as Step 2
# only new part kept here:

resolved_step4 <- read.table("TaxonomicHarm/BenthicMacroinvertebrate/resolved_step4.csv", h = TRUE, sep = ",", encoding = "UTF-8")

resolved_step2_4 <- rbind(resolved_step2, resolved_step3, resolved_step4)

unresolved2 <- unresolved %>%
  filter(!Taxon_clean %in% resolved_step2_4$Taxon_clean)

unresolved2[unresolved2$Taxon == "Heteroptera", "matched_name"] <- "Hemiptera"

# run same manual GBIF resolution block as in Step 2 on unresolved2
# output: resolved_step5

# Unresolved - Step 6 ------------------------------------------------------
# same manual GBIF workflow as Step 2
# only new part kept here:

resolved_step5 <- read.table("TaxonomicHarm/BenthicMacroinvertebrate/resolved_step5.csv", h = TRUE, sep = ",", encoding = "UTF-8")

resolved_step2_5 <- rbind(resolved_step2, resolved_step3, resolved_step4, resolved_step5)

unresolved2 <- unresolved %>%
  filter(!Taxon_clean %in% resolved_step2_5$Taxon_clean)

unresolved2[unresolved2$Taxon == "Cladocera", "matched_name"]    <- "Diplostraca"
unresolved2[unresolved2$Taxon == "Prodiamesinae", "matched_name"] <- "Chironomidae"

# run same manual GBIF resolution block as in Step 2 on unresolved2
# output: resolved_step6

# Unresolved - Step 7 ------------------------------------------------------
# same manual GBIF workflow as Step 2
# only new part kept here:

resolved_step6 <- read.table("TaxonomicHarm/BenthicMacroinvertebrate/resolved_step6.csv", h = TRUE, sep = ",", encoding = "UTF-8")

resolved_step2_6 <- rbind(resolved_step2, resolved_step3, resolved_step4, resolved_step5, resolved_step6)

unresolved2 <- unresolved %>%
  filter(!Taxon_clean %in% resolved_step2_6$Taxon_clean)

unresolved2[unresolved2$Taxon == "Gresnia", "matched_name"] <- "Gresnia"
unresolved2[unresolved2$Taxon == "Branchiopoda.Cladocera", "matched_name"] <- "Diplostraca"

# run same manual GBIF resolution block as in Step 2 on unresolved2
# output: resolved_step7

# Unresolved - Step 8 ------------------------------------------------------
# here the workflow changes: use NBN instead of GBIF

resolved_step7 <- read.table("TaxonomicHarm/BenthicMacroinvertebrate/resolved_step7.csv", h = TRUE, sep = ",", encoding = "UTF-8")

resolved_step2_7 <- rbind(resolved_step2, resolved_step3, resolved_step4, resolved_step5, resolved_step6, resolved_step7)

unresolved2 <- unresolved %>%
  filter(!Taxon_clean %in% resolved_step2_7$Taxon_clean)

unresolved2 <- unresolved2 %>%
  mutate(
    matched_name = case_when(
      Taxon_clean == "Hydropsyche/Hydroptila" ~ "Trichoptera",
      Taxon_clean == "Basommatophora"         ~ "Gastropoda",
      Taxon_clean == "Chaetopterygini"        ~ "Limnephilidae",
      Taxon_clean == "Pseudochironomini"      ~ "Chironominae",
      TRUE ~ matched_name
    )
  )

# run same manual resolution workflow, but using NBN:
# - get_nbnid(...)
# - classification(..., db = "nbn")
# - tax_rank(..., db = "nbn")
# output: resolved_step8

# Unresolved - Step 9 / Final merge ---------------------------------------

resolved_step8 <- read.table("TaxonomicHarm/BenthicMacroinvertebrate/resolved_step8.csv", h = TRUE, sep = ",", encoding = "UTF-8")

resolved_step2_8 <- bind_rows(
  resolved_step2, resolved_step3, resolved_step4,
  resolved_step5, resolved_step6, resolved_step7, resolved_step8
) %>%
  unique()

tax_table_final_updated <- tax_table_final %>%
  filter(!Taxon %in% resolved_step2_8$Taxon) %>%
  bind_rows(resolved_step2_8)

tax_table_final_updated <- tax_table_final_updated %>%
  relocate(nbn_id, .after = gbif_id) %>%
  arrange(phylum, class, order, family, genus, species)

write.csv(
  tax_table_final_updated,
  "TaxonomicHarm/BenthicMacroinvertebrate/Harmonization_table_BenthicMacroinvertebrates.csv",
  row.names = FALSE
)

# Merge to main dataset ----------------------------------------------------

Data1 <- read_parquet("Combined_dataset/BenthicMacroinvertebrates.parquet")
Data2 <- read_parquet("Combined_dataset/BenthicMacroinvertebrates_part2.parquet")
Data  <- rbind(Data1, Data2)

HarmonizationTable <- read.table(
  "TaxonomicHarm/BenthicMacroinvertebrate/Harmonization_BenthicMacroinvertebrates.csv",
  h = TRUE, sep = ",", encoding = "UTF-8"
)

Data_taxonomy <- Data %>%
  left_join(
    HarmonizationTable %>%
      select(
        Taxon, Taxon_clean, accepted_name,
        taxonomy_rank = matched_rank,
        species, genus, family, order, class, phylum
      ),
    by = "Taxon"
  )

Data_taxonomy <- Data_taxonomy %>%
  select(
    BioticGroup:Taxon,
    Taxon_clean, accepted_name, taxonomy_rank,
    species, genus, family, order, class, phylum,
    Value,
    ValueType, Unit, FishPresence,
    ReferenceCondition, source_file
  )

# write_parquet(Data_taxonomy, "TaxonomicHarm/BenthicMacroinvertebrate/BenthicMacroinvertebrates_withTaxonomy.parquet")

