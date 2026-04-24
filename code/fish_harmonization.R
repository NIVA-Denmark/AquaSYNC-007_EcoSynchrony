## Explore data and harmonize taxonomy for AquaSync Fish Data

rm(list=ls())

library(taxize)
library(arrow)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#TODO check whether ITIS works better for fish (e.g. has taxonomic class)


#source("FunctionsTaxonomicHarmonization.R")

dat.raw <- read_parquet("../Fish_all.parquet")

#bm1 <- read_parquet("../BenthicMacroinvertebrates_all.parquet")
#bm2 <- read_parquet("../BenthicMacroinvertebrates_all_withTaxonomy.parquet")


table(dat.raw$Lake.river)
table(dat.raw$Sampled.habitat)
table(dat.raw$Country)
table(dat.raw$SiteName)
table(dat.raw$Taxon)
table(dat.raw$Taxon[dat.raw$Country=="USA"])
table(dat.raw$Taxon[dat.raw$Country=="Norway"])

dat.raw$Taxon <- as.character(dat.raw$Taxon)

## USA data have arctic grayling, chinook, coho, dolly varden, slimy sculpin, and stickleback
## Question for data providers: 
##    Are stickleback one species given a general name, or are they multiple spp? (not many here)

## Norway data have brook trout, char, cottus, perch, phoxinus, and trout
## Question for data providers: 
##    Are trout and brook trout the same? There's not many "brook trout" so maybe lump? Are trout multiple spp?
##    Are char all artic char? or a group of species?
##    Are cottus (freshwater sculpins) one species or multiple?
##    Are perch one species or multiple? (not many here)
##    Are phoxinus (minnows) one species or multiple?



## Taxonomy Harmonization -------------------------------------------------------------------------

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


write.csv(tax_table, "../fish_taxonomyTable_gbif.csv", row.names=FALSE)






# ## Sampling through time - USA --------------------------------------------------------
# dat.usa <- dat.raw[dat.raw$Country=="USA",]
# spp <- unique(as.character(dat.usa$Taxon))
# 
# #sampling through time
# par(mar=c(4.1,8,2,1))
# plot(NA,NA, ylim=c(0.5, length(spp)+0.5), xlim=range(dat.usa$Date),
#      xaxt="n", yaxt="n", xlab="Time", ylab="")
# axis(2, at=1:length(spp), labels=spp, las=2)
# axis(1, at=pretty(dat.usa$Date), labels=as.character(as.Date(pretty(dat.usa$Date))))
# abline(h=1:length(spp))
# for(ii in 1:length(spp)){
#   dat.sub <- dat.usa[dat.usa$Taxon==spp[ii],]
#   points(dat.sub$Date, rep(ii, nrow(dat.sub)), pch=16, col="blue")
# }
# mtext("USA Fishes")
# 
# #sampling by location
# dat.usa$SiteName <- factor(dat.usa$SiteName)
# 
# 
# par(mar=c(12,8,2,1))
# plot(NA,NA, ylim=c(0.5, length(spp)+0.5), xlim=c(1,length(unique(dat.usa$SiteName))),
#      xaxt="n", yaxt="n", xlab="", ylab="")
# axis(2, at=1:length(spp), labels=spp, las=2)
# axis(1, at=1:length(unique(dat.usa$SiteName)), labels=unique(dat.usa$SiteName), las=2)
# abline(h=1:length(spp))
# abline(v=1:length(unique(dat.usa$SiteName)))
# for(ii in 1:length(spp)){
#   dat.sub <- dat.usa[dat.usa$Taxon==spp[ii],]
#   points(as.numeric(dat.sub$SiteName), rep(ii, nrow(dat.sub)), pch=16, col="blue")
# }
# mtext("USA Fishes")
# 
# 
# 
# ## Sampling through time - Norway --------------------------------------------------------
# dat.nor <- dat.raw[dat.raw$Country=="Norway",]
# spp <- unique(as.character(dat.nor$Taxon))
# 
# par(mar=c(4.1,8,2,1))
# plot(NA,NA, ylim=c(0.5, length(spp)+0.5), xlim=range(dat.nor$Date),
#      xaxt="n", yaxt="n", xlab="Time", ylab="")
# axis(2, at=1:length(spp), labels=spp, las=2)
# axis(1, at=pretty(dat.nor$Date), labels=as.character(as.Date(pretty(dat.nor$Date))))
# abline(h=1:length(spp))
# for(ii in 1:length(spp)){
#   dat.sub <- dat.nor[dat.nor$Taxon==spp[ii],]
#   points(dat.sub$Date, rep(ii, nrow(dat.sub)), pch=16, col="blue")
# }
# mtext("Norway Fishes")
# 
# 
# dat.nor$SiteName <- factor(dat.nor$SiteName)
# 
# par(mar=c(12,8,2,1))
# plot(NA,NA, ylim=c(0.5, length(spp)+0.5), xlim=c(1,length(unique(dat.nor$SiteName))),
#      xaxt="n", yaxt="n", xlab="", ylab="")
# axis(2, at=1:length(spp), labels=spp, las=2)
# axis(1, at=1:length(unique(dat.nor$SiteName)), labels=unique(dat.nor$SiteName), las=2)
# abline(h=1:length(spp))
# abline(v=1:length(unique(dat.nor$SiteName)))
# for(ii in 1:length(spp)){
#   dat.sub <- dat.nor[dat.nor$Taxon==spp[ii],]
#   points(as.numeric(dat.sub$SiteName), rep(ii, nrow(dat.sub)), pch=16, col="blue")
# }
# mtext("Norway Fishes")
# 
