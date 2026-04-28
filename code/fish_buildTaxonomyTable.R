## Make fish taxonomy table for AquaSync analyses

rm(list=ls())

library(taxize)
library(arrow)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Load data and helper functions, some initial checking -----------------------------------------

dat.raw <- read_parquet("../Fish_all.parquet")

#option to load benthic macroinverts to check formatting consistency
#bm1 <- read_parquet("../BenthicMacroinvertebrates_all.parquet")
#bm2 <- read_parquet("../BenthicMacroinvertebrates_all_withTaxonomy.parquet")

dat.raw$Taxon <- as.character(dat.raw$Taxon)


### Build taxonomy table -------------------------------------------------------------------------

# For species listed with common name, assign a scientific name
dat.raw$Taxon_orig <- dat.raw$Taxon #keep original
dat.raw$Taxon[dat.raw$Taxon=="Arctic grayling"] <- "Thymallus arcticus"
dat.raw$Taxon[dat.raw$Taxon=="Brook_trout"] <- "Salvelinus fontinalis"
dat.raw$Taxon[dat.raw$Taxon=="Char"] <- "Salvelinus alpinus"
dat.raw$Taxon[dat.raw$Taxon=="Chinook"] <- "Oncorhynchus tshawytscha"
dat.raw$Taxon[dat.raw$Taxon=="Coho"] <- "Oncorhynchus kisutch"
dat.raw$Taxon[dat.raw$Taxon=="Dolly Varden"] <- "Salvelinus malma"
dat.raw$Taxon[dat.raw$Taxon=="Perch"] <- "Perca fluviatilis"
dat.raw$Taxon[dat.raw$Taxon=="Phoxinus"] <- "Phoxinus phoxinus"
dat.raw$Taxon[dat.raw$Taxon=="slimy sculpin"] <- "Cottus cognatus"
dat.raw$Taxon[dat.raw$Taxon=="Stickleback"] <- "Gasterosteidae"
dat.raw$Taxon[dat.raw$Taxon=="Trout"] <- "Salmo trutta"

dat.raw <- dat.raw[dat.raw$Taxon != "Unidentified Unidentified",]


# cleaning
TaxaList <- dat.raw %>%
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
    
    # # specific corrections
    # Taxon_clean = case_when(
    #   Taxon_clean == "Valvata." ~ "Valvata sp.",
    #   Taxon_clean == "Boreobdella verrucata" ~ "Glossiphonia verrucata",
    #   Taxon_clean == "Hydropsyhce angustipennis" ~ "Hydropsyche angustipennis",
    #   Taxon_clean == "Gyraulus." ~ "Gyraulus sp.",
    #   Taxon_clean == "Amphnemoura sulcicollis" ~ "Amphinemoura sulcicollis",
    #   Taxon_clean == "Philopotomidae" ~ "Philopotamidae",
    #   Taxon_clean == "Paraleutra" ~ "Paraleuctra",
    #   Taxon_clean == "Oligochaetae" ~ "Oligochaeta",
    #   Taxon_clean == "Elminthidae" ~ "Elmidae",
    #   Taxon_clean %in% c("Elmis aena Lv.", "Elmis aena Ad.", "Elmis aena lv.") ~ "Elmis aenea",
    #   TRUE ~ Taxon_clean
    # )
  ) %>%
  distinct()




## match taxa names with online databases ---------------------------------------------------------

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
            ask = TRUE,
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

tax_table$subclass <- NA_character_

tax_table <- tax_table[,c(1:11,14,12,13)]


write.csv(tax_table, "../fish_taxonomyTable_gbif.csv", row.names=FALSE)
