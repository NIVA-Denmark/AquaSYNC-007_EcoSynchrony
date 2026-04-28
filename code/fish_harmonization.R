## Explore data and harmonize taxonomy for AquaSync Fish Data

rm(list=ls())

library(taxize)
library(arrow)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## TODO data with ValueType=="Weight" appear to be weights of individual fish, need to sum by site and sampling event
## TODO is ValueType=="Total weight" the same as ValueType=="Biomass"? 
## TODO For Salmo salar, "Total weight" there are often two rows, maybe for adult and juvenile and the lifestage info has been lost?
## - could also be M/F, get identifying info and talk with Jen
## TODO how to use length data? A lot of sites have both length data and another kind of data, but not all.
##  - It looks like 16 north american sites only have length data.
##  - For these sites, is sampling scheme such that we can use counts of numbers of fish as an abundance metric?
##  - Can we assume effort is consistent enough over time to use data that aren't normalized by effort?

### Load data and helper functions, some initial checking -----------------------------------------

source("FunctionsTaxonomicHarmonization.R")
source("FunctionsTemporalHarmonization.R")

dat.raw <- read_parquet("../Fish_all.parquet")

#option to load benthic macroinverts to check formatting consistency
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

### ID cases where there are 2 rows for Salmo salar -----------------------------------------------

# sub1 <- dat.raw[dat.raw$Taxon=="Salmo salar" & dat.raw$ValueType=="Total weight",]
# 
# site_x_date <- unique(cbind(as.character(sub1$SiteID), as.character(sub1$Date)))
# 
# write.csv(sub1, "../SalmoSalar_2weights_full.csv", row.names=FALSE)
# write.csv(site_x_date, "../SalmoSalar_2weights_SiteIDxDate.csv", row.names=FALSE)


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



### Add taxonomy to data and export parquet file --------------------------------------------------

dat_taxa <- left_join(dat.raw, tax_table)
write_parquet(dat_taxa, "../fish_all_withTaxonomy.parquet")



### Exclude or otherwise deal with ValueTypes that will not analyzed ------------------------------
#origN <- nrow(dat_taxa)
#dat_taxa <- dat_taxa[!grepl("length", dat_taxa$ValueType),]
dat_taxa <- dat_taxa[dat_taxa$ValueType != "PresenceAbsence",] #possibly omit for richness-based analyses, but only 12 occurrences 
#table(as.character(dat_taxa$ValueType))


### Aggregate individual-level data to counts and total biomass -----------------------------------



### Temporal harmonization ------------------------------------------------------------------------

# Keep only the columns needed for temporal harmonization.
Data_temporal <- dat_taxa %>%
  select(SiteID, Date) %>%
  distinct()

# Harmonize sampling dates across site-level time series.
Data_temporal_harmonized <- harmonize_sampling(
  Data = Data_temporal,
  site_col = "SiteID",
  group_col = "SiteID",   # or "Country", "climate", etc.
  date_col = "Date",
  min_years = 6,
  buffer_months = 3
)

length(unique(dat.raw$SiteID)) #212 sites
length(unique(Data_temporal_harmonized$SiteID)) #211 sites--only one gets dropped.

Data_temporal_harmonized$sampleID <- paste(Data_temporal_harmonized$SiteID, Data_temporal_harmonized$Date)
dat_taxa_timesub <- dat_taxa
dat_taxa_timesub$sampleID <- paste(dat_taxa_timesub$SiteID, dat_taxa_timesub$Date)
dat_taxa_timesub <- dat_taxa_timesub[dat_taxa_timesub$sampleID %in% Data_temporal_harmonized$sampleID,]

### Taxonomic harmonization -----------------------------------------------------------------------

# Keep the columns required for taxonomic harmonization.
Data_taxonomy <- dat_taxa_timesub %>%
  select(
    SiteID, Date, Taxon_clean,
    species, genus, family, order, subclass, class, phylum,
    Value
  )

# Keep sample-level metadata to reattach after harmonization.
Data_sample_level <- dat_taxa_timesub %>%
  select(
    BioticGroup, Lake.river, Sampled.habitat,
    Country, SiteID, SiteName,
    WaterBody, Lon, Lat, Day, Month,
    Year, Date,
    ValueType, Unit,
    FishPresence, ReferenceCondition, climate, source_file
  ) %>%
  distinct()

# Harmonize taxonomy 
Data_taxonomy_harmonized <- harmonize_taxonomy(
  df = Data_taxonomy,
  decision_cols = "SiteID",   # or "Country", "climate", etc.
  output_site_cols = "SiteID",
  date_col = "Date",
  abundance_col = "Value",
  threshold = 0.90
)

# Reattach sample-level metadata 


Data_harmonized <- Data_taxonomy_harmonized %>%
  left_join(Data_sample_level, by = c("SiteID", "Date")) %>%
  select(
    any_of(names(dat_taxa)),  # keep original column order where possible
    everything()          # then append any additional columns
  )

table(dat.raw$ValueType)

#trying to resolve warning messages
Data_taxonomy_harmonized$SiteID[1133]; Data_taxonomy_harmonized$Date[1133]


check <- Data_sample_level[Data_sample_level$SiteID=="Könkämäeno 1" & Data_sample_level$Date=="1986-09-05",]

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
