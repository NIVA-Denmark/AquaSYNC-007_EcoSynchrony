
library(arrow)
library(dplyr)
library(stringr)

# Merge taxonomy to main dataset ---------------------------------------------------

Data <- read_parquet( "Combined_dataset/all/BenthicMacroinvertebrates_all.parquet")

HarmonizationTable <- read.table("TaxonomicHarm/BenthicMacroinvertebrate/BenthicMacroinvertebrates_Taxonomy_table_qualityChecked.csv",h=T, sep=",", encoding="UTF-8")

length(unique(Data$Taxon))
length(unique(HarmonizationTable$Taxon))

Data_taxonomy <- Data %>%
  left_join(HarmonizationTable %>%
              select(Taxon,Taxon_clean, accepted_name, taxonomy_rank= matched_rank, species,genus,subfamily, family,order,subclass, class,subphylum,phylum), by="Taxon")
dim(Data)
dim(Data_taxonomy)
Data_taxonomy <- Data_taxonomy %>%
  select(BioticGroup:Taxon,
         Taxon_clean, accepted_name, taxonomy_rank,     
         species,genus,subfamily, family,order,subclass,class,subphylum,phylum,
         Value,
         ValueType, Unit, FishPresence,
         ReferenceCondition, source_file, climate)


# change taxon name and rank based on taxonomy info
rank_cols <- c(
  "phylum", "subphylum", "class", "subclass", "order",
  "family", "subfamily", "genus", "species"
)

Data_taxonomy <- Data_taxonomy %>%
  mutate(
    across(
      all_of(rank_cols),
      ~ {
        x <- as.character(.x)
        x <- trimws(x)
        x[x %in% c("", "NA", "Na", "na", "NULL", "null")] <- NA_character_
        x
      }
    ),
    
    Rank.harm = apply(
      as.data.frame(across(all_of(rank_cols))),
      1,
      function(z) {
        idx <- which(!is.na(z) & z != "")
        if (length(idx) == 0) NA_character_ else rank_cols[max(idx)]
      }
    ),
    
    Taxon.harm = apply(
      as.data.frame(across(all_of(rank_cols))),
      1,
      function(z) {
        idx <- which(!is.na(z) & z != "")
        if (length(idx) == 0) NA_character_ else as.character(z[max(idx)])
      }
    )
  ) %>%
  select(
    -Taxon,
    -Taxon_clean,
    -accepted_name,
    -taxonomy_rank
  ) %>%
  relocate(Taxon.harm, Rank.harm, .before = species) %>%
  rename(Taxon = Taxon.harm)




write_parquet(Data_taxonomy, "TaxonomicHarm/BenthicMacroinvertebrate/BenthicMacroinvertebrates_withTaxonomy.parquet")


# Compute averages for cases when > 1 sample per sampling occasion  ------------------------------------------------


Data <- read_parquet( "TaxonomicHarm/BenthicMacroinvertebrate/BenthicMacroinvertebrates_withTaxonomy.parquet")

N_samples <- Data %>%
  group_by(Country, SiteID, Date, source_file) %>%
  summarise(nsampl = length(unique(SampleID))) %>%
  ungroup()

N_samples %>% 
  filter(nsampl>1) %>% 
  select(Country, SiteID, source_file, Date, nsampl) %>% 
  distinct() %>% print(n=260)


View(Data %>%
  filter(SiteID =="BEL01_c"))



# correct Finnish site (multiple Sample IDs have typos): Finland S045 Aquasync_Macroinvertebrates_Finland_Aroviita.xlsx

Data %>%
  filter(SiteID =="S045" &  Date=="2024-09-03") %>%
  select(SiteID, Date, SampleID)
Data <- Data %>%
  mutate(
    SampleID = if_else(
      SiteID == "S045" & Date == "2024-09-03",
      "S045_20240903",
      SampleID
    )
  )

# compute average abundances across samples for the other cases: 
N_samples <- Data %>%
  group_by(Country, SiteID, Date, source_file) %>%
  summarise(nsampl = length(unique(SampleID))) %>%
  ungroup()
N_samples %>% 
  filter(nsampl>1) %>% 
  select(Country, SiteID, source_file, Date, nsampl) %>% 
  distinct() %>% print(n=140)



library(dplyr)
library(tidyr)


tax_cols <- c(
  "Taxon", "Rank.harm",
  "species", "genus", "subfamily", "family", "order",
  "subclass", "class", "subphylum", "phylum"
)

site_cols <- c("Country", "SiteID", "Date", "source_file")

meta_cols <- c(
  "BioticGroup", "Lake.river", "Sampled.habitat",
  "SiteName", "WaterBody", "Lon", "Lat",
  "Day", "Month", "Year",
  "ValueType", "Unit",
  "FishPresence", "ReferenceCondition",
  "climate"
)

N_samples <- Data %>%
  group_by(across(all_of(site_cols))) %>%
  summarise(
    nsampl = n_distinct(SampleID),
    original_SampleID = first(SampleID),
    across(all_of(meta_cols), first),
    .groups = "drop"
  )

Data_avg <- Data %>%
  group_by(across(all_of(c(site_cols, tax_cols)))) %>%
  summarise(
    Value = sum(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(N_samples, by = site_cols) %>%
  mutate(
    Value = Value / nsampl,
    SampleID = if_else(
      nsampl == 1,
      original_SampleID,
      paste0(SiteID, "_", format(Date, "%Y%m%d"))
    )
  ) %>%
  select(
    all_of(meta_cols),
    all_of(site_cols),
    SampleID,
    all_of(tax_cols),
    Value,
    nsampl
  )


# check
Data %>%
  filter(SiteID == "BEA01_a",
         Date == "2009-10-01") %>%
  distinct(SampleID)
Data_avg %>%
  filter(SiteID == "BEA01_a",
         Date == "2009-10-01") %>%
  distinct(SampleID)
Data %>%
  filter(SiteID == "BEA01_a",
         Date == "2009-10-01",
         family == "Chironomidae") %>%
  select(SampleID, family, Value)
Data_avg %>%
  filter(SiteID == "BEA01_a",
         Date == "2009-10-01",
         family == "Chironomidae") %>%
  select(SampleID, family, Value)

Data_avg %>%
  filter(SiteID == "126",
         Date == "2000-09-12",
         family == "Chironomidae") %>%
  select(SampleID, family, Value)


Data_old <- Data
Data <- Data_avg

Data_old %>%
  select(SiteID) %>%
  distinct()
Data %>%
  select(SiteID) %>%
  distinct()

Data_old %>%
  select(SiteID, Date) %>%
  distinct()
Data%>%
  select(SiteID, Date) %>%
  distinct()


# check if there are repeated Taxon in SiteID-Dte combinations
Data %>%
  count(SiteID, Date, Taxon, name = "n") %>%
  filter(n > 1)
Data %>%
  count(across(all_of(c(site_cols, "Taxon"))), name = "n") %>%
  filter(n > 1)
Data %>%
  group_by(across(all_of(c(site_cols, "Taxon")))) %>%
  filter(n() > 1) %>%
  arrange(SiteID, Date, Taxon)
Data %>%
  count(SiteID, Date, SampleID, Taxon, name = "n") %>%
  filter(n > 1)


# Harmonize sampling effort - Select 1 sample per year ------------------------------------------------
source("FunctionsTemporalHarmonization_v2.R")


length(unique(Data$SiteID)) # 739


# Filter original data
# e.g. discard ValueType == Presence/Absence   ?
Data %>%
  select(ValueType, Unit) %>%
  distinct()    
Data <-  Data %>%
  filter(ValueType!= "Presence/Absence")

length(unique(Data$SiteID)) # 738


# e.g. discard taxa not resolved: without phylum 
Data <- Data %>%
  filter(!is.na(phylum))


# remove profundal and sublittoral samples
Data <- Data %>%
  filter(!Sampled.habitat %in% c("Profundal", "Sublittoral")) 
Data_old %>%
  filter(Sampled.habitat %in% c("Profundal", "Sublittoral")) %>%
  select(SiteID) %>% distinct()
length(unique(Data$SiteID))  # 580
summary(factor(Data$Sampled.habitat))

# select relevant columns
Data_sel <- Data %>%
  select(SiteID, Date) %>% distinct()
unique_siteid <- unique(Data_sel$SiteID)

# run harm function
Sel_samples <- select_one_sample_and_flag(
  Data = Data_sel,
  site_col = "SiteID",
  date_col = "Date",
  buffer_days = 60
)

# check
Sel_samples %>%
  mutate(.Year= year(Date)) %>%
  group_by(SiteID,.Year) %>%
  summarise(n_samples = length(unique(Date))) %>% 
  filter(n_samples>1) %>% print(n=300)

# subset selected samples
Data_selectedSamples <- Data %>%
  semi_join(Sel_samples, by = c("SiteID", "Date"))
  
Data_selectedSamples %>%
  group_by(SiteID,Year) %>%
  summarise(n_samples = length(unique(Date))) %>% 
  filter(n_samples>1) 
  
length(unique(Data_selectedSamples$SiteID))
# check if there are repeated Taxon in SiteID-Dte combinations
Data_selectedSamples %>%
  count(SiteID, Date,  Taxon, name = "n") %>%
  filter(n > 1)
Data_selectedSamples %>%
  count(across(all_of(c(site_cols, "Taxon"))), name = "n") %>%
  filter(n > 1)
Data_selectedSamples %>%
  group_by(across(all_of(c(site_cols, "Taxon")))) %>%
  filter(n() > 1) %>%
  arrange(SiteID, Date, Taxon)
Data_selectedSamples %>%
  count(SiteID, Date, SampleID, Taxon, name = "n") %>%
  filter(n > 1)

# Prepare for taxonomic harmonization -----------------------------------------------------------
source("FunctionsTaxonomicHarmonization_v5.R") # # 2 step procedure

Data <- Data_selectedSamples


# e.g. exclude zooplankton, fish, algae etc..
summary(factor(Data$phylum))
Data <- Data %>%
  filter(!phylum %in% c("Chordata",  "Cnidaria",              
                        "Nematoda", "Nematomorpha","Platyhelminthes", 
                        "Porifera", "Nemertea"))   # , "Tardigrada"
summary(factor(Data$class))
Data <- Data %>%
  filter(!class %in% c("Branchiopoda", "Copepoda", "Collembola", "Ostracoda"))   #, "Ostracoda","Charophyceae","Anostraca", 
summary(factor(Data$order))
Data <- Data %>%
  filter(!order %in% c("Trombidiformes", "Sarcoptiformes", "Arguloida"))  


# Select relevant columns
Data_taxonomy <- Data %>%
  select(SiteID, Date, Taxon, species, genus, subfamily, family, order, subclass, class, subphylum, phylum, Value)

Data_SampleLevel <- Data %>%
   select(BioticGroup, Lake.river,Sampled.habitat,
                    Country,SiteID,SiteName,        
                    WaterBody,Lon,Lat, Day,Month,
                    Year,Date,
                    ValueType,Unit,             
                    FishPresence,ReferenceCondition, source_file, climate ) %>%
  distinct()

sites <- unique(Data_taxonomy$SiteID)

testdata <- Data_taxonomy[Data_taxonomy$SiteID %in% c("S183", "Vid","Sol", "Vol", "BNP-BOW-11"),]

# check if there are repeated Taxon in SiteID-Date combinations
Data_taxonomy %>%
  count(SiteID, Date,  Taxon, name = "n") %>%
  filter(n > 1)

# genus level upwards
#Data_taxonomy$species <- NA




# manual roll ups ---------------------------------------------------------
manual_rollup_rules <- tribble(
  ~Country,  ~SiteID,             ~target_rank, ~target_taxon,
  "Canada",  "ASH01_b",           "family",     "Chironomidae",
  "Canada",  "BNP-Johnston3",     "family",     "Chironomidae",
  "Canada",  "BOI001",            "family",     "Nemouridae",
  "Canada",  "BOI002",            "family",     "Nemouridae",
  "Canada",  "CLD01_a",           "family",     "Naididae",
  "Canada",  "CLD02_a",           "family",     "Naididae",
  "Canada",  "CLD02_a",           "family",     "Baetidae",
  "Canada",  "CYT-L4-0.8",        "family",     "Chironomidae",
  "Canada",  "FIRMID02",          "family",     "Chironomidae",
  "Canada",  "JPRIFFDN",          "family",     "Baetidae",
  "Canada",  "K-K-07-42",         "family",     "Nemouridae",
  "Canada",  "M2A-1",             "family",     "Baetidae",
  "Canada",  "MCKRIFF04",         "family",     "Chironomidae",
  "Canada",  "QUI01",             "family",     "Naididae",
  "Canada",  "S-2",               "family",     "Ephemerellidae",
  "Canada",  "SAL05",             "family",     "Naididae",
  "Sweden",  "SE660713-165636",   "genus",      "Baetis",
  "Canada",  "STBRIFF11",         "family",     "Brachycentridae",
  "Canada",  "STBRIFF20",         "family",     "Naididae",
  "Canada",  "STBRIFF21",         "family",     "Baetidae",
  "Canada",  "UFR07",             "family",     "Ephemerellidae",
  "Canada",  "WIG03",             "family",     "Heptageniidae",
  "Canada",  "YPS-107",           "family",     "Heptageniidae",
  "Canada",  "ALO01",             "family",     "Chironomidae",
  "Canada",  "BIRTRIB10",         "family",     "Chironomidae",
  "Canada",  "BNP-BOW-6B",        "family",     "Chironomidae",
  "Canada",  "COW03",             "family",     "Chironomidae",
  "Canada",  "DOVRIFF05",         "family",     "Chironomidae",
  "Canada",  "ELLSRIFF03A",       "family",     "Chironomidae",
  "Canada",  "ELLSRIFF07",        "family",     "Chironomidae",
  "Canada",  "ELLSRIFF09",        "family",     "Chironomidae",
  "Canada",  "ELLSRIFF12",        "family",     "Chironomidae",
  "Canada",  "FIRMID",            "family",     "Chironomidae",
  "Canada",  "FIRMID03",        "family",     "Chironomidae",
  "Canada",  "HRSRIFF01",        "family",     "Chironomidae",
  "Canada",  "JPRIFFDN",        "family",     "Chironomidae",
  "Canada",  "M6-3",            "family",     "Chironomidae",
  "Canada",  "MCKRIFF02",        "family",     "Chironomidae",
  "Canada",  "MCKRIFF08",        "family",     "Chironomidae",
  "Canada",  "MCKRIFF09",        "family",     "Chironomidae",
  "Canada",  "MCKRIFF12",        "family",     "Chironomidae",
  "Canada",  "MCKRIFF15",        "family",     "Chironomidae",
  "Canada",  "NAL-NIAG-01",        "family",     "Chironomidae",
  "Canada",  "SAN01_b",        "family",     "Chironomidae",
  "Canada",  "STBRIFF03",        "family",     "Chironomidae",
  "Canada",  "STBRIFF07",        "family",     "Chironomidae",
  "Canada",  "STBRIFF09",        "family",     "Chironomidae",
  "Canada",  "STBRIFF16",        "family",     "Chironomidae", 
  "Canada",  "STBRIFF20",        "family",     "Chironomidae",
  "Canada",  "CLD01_a",           "family",     "Baetidae",
  "Canada",  "COW03",            "family",     "Baetidae",
  "Canada",  "ELK01_a",         "family",     "Baetidae",
  "Canada",  "ELLSRIFF02",         "family",     "Baetidae",
  "Canada",  "EQU02",         "family",     "Baetidae",
  "Canada",  "GRA01",         "family",     "Baetidae",
  "Canada",  "GRA02",         "family",     "Baetidae",
  "Canada",  "HRSRIFF02",         "family",     "Baetidae",
  "Canada",  "KOS01",         "family",     "Baetidae",
  "Canada",  "LIZ001",         "family",     "Baetidae",
  "Canada",  "MCKRIFF05",         "family",     "Baetidae",
  "Canada",  "MCKRIFF09",         "family",     "Baetidae",
  "Canada",  "NWL-TSUL-01",         "family",     "Baetidae",
  "Canada",  "STBRIFF09",         "family",     "Baetidae",
  "Canada",  "STBRIFF10",         "family",     "Baetidae",
  "Iceland",  "Vesturdalsá",         "family",     "Limnephilidae",
  "Iceland",  "Elliðaár",         "genus",     "Diamesa",
  "Canada",  "ALX003",         "family",     "Ephemerellidae",
  "Canada",  "B-B-17-03",         "family",     "Ephemerellidae",
  "Canada",  "HRS01",         "family",     "Ephemerellidae",
  "Canada",  "LIZ001",         "family",     "Ephemerellidae"
  
  
)

apply_manual_rollups_raw <- function(df, rules, value_col = "Value") {
  
  rank_index <- setNames(seq_along(rank_cols), rank_cols)
  
  out <- df %>%
    mutate(manual_rollup = FALSE)
  
  use_country <- "Country" %in% names(out) && "Country" %in% names(rules)
  
  for (i in seq_len(nrow(rules))) {
    
    target_rank_i <- rules$target_rank[i]
    target_taxon_i <- rules$target_taxon[i]
    target_idx <- rank_index[[target_rank_i]]
    
    deeper_cols <- if (target_idx < length(rank_cols)) {
      rank_cols[(target_idx + 1):length(rank_cols)]
    } else {
      character(0)
    }
    
    match_rule <-
      out$SiteID == rules$SiteID[i] &
      !is.na(out[[target_rank_i]]) &
      out[[target_rank_i]] == target_taxon_i
    
    if (use_country) {
      match_rule <- match_rule &
        as.character(out$Country) == rules$Country[i]
    }
    
    match_rule[is.na(match_rule)] <- FALSE
    
    if (length(deeper_cols) > 0) {
      out[match_rule, deeper_cols] <- NA_character_
    }
    
    if ("Taxon" %in% names(out)) {
      out$Taxon[match_rule] <- target_taxon_i
    }
    
    out$manual_rollup[match_rule] <- TRUE
  }
  
  group_cols <- setdiff(names(out), c(value_col, "manual_rollup"))
  
  out %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      !!value_col := sum(.data[[value_col]], na.rm = TRUE),
      manual_rollup = any(manual_rollup),
      .groups = "drop"
    )
}

Data_taxonomy2 <- Data_taxonomy %>%
  apply_manual_rollups_raw(manual_rollup_rules)

Data_taxonomy2 %>% 
  filter(SiteID == "ASH01_b", family == "Chironomidae") %>%
  distinct(Taxon, family, subfamily, genus, species, order, class, manual_rollup)
Data_taxonomy %>% 
  filter(SiteID== "ASH01_b", family == "Chironomidae") %>%
  select(Taxon) %>% distinct()

Data_taxonomy2 %>% 
  filter(SiteID == "M2A-1", family == "Baetidae") %>%
  distinct(Taxon, family, subfamily, genus, species, order, class, manual_rollup)
Data_taxonomy %>% 
  filter(SiteID== "M2A-1", family == "Baetidae") %>%
  select(Taxon) %>% distinct()


Data_taxonomy2 %>% 
  filter(SiteID == "SE660713-165636", genus == "Baetis") %>%
  distinct(Taxon, family, subfamily, genus, species, order, class, manual_rollup)
Data_taxonomy %>% 
  filter(SiteID== "SE660713-165636", genus == "Baetis") %>%
  select(Taxon) %>% distinct()

# check
check_manual_rollups <- function(before, after, rules) {
  
  out <- list()
  
  for (i in seq_len(nrow(rules))) {
    
    site_i <- rules$SiteID[i]
    target_rank_i <- rules$target_rank[i]
    target_taxon_i <- rules$target_taxon[i]
    
    before_i <- before %>%
      filter(
        SiteID == site_i,
        !is.na(.data[[target_rank_i]]),
        .data[[target_rank_i]] == target_taxon_i
      ) %>%
      distinct(
        SiteID,
        target_rank = target_rank_i,
        target_taxon = target_taxon_i,
        dataset = "before",
        Taxon,
        family, subfamily, genus, species,
        order, class
      )
    
    after_i <- after %>%
      filter(
        SiteID == site_i,
        !is.na(.data[[target_rank_i]]),
        .data[[target_rank_i]] == target_taxon_i
      ) %>%
      distinct(
        SiteID,
        target_rank = target_rank_i,
        target_taxon = target_taxon_i,
        dataset = "after",
        Taxon,
        family, subfamily, genus, species,
        order, class,
        manual_rollup
      )
    
    out[[i]] <- bind_rows(
      before_i %>% mutate(manual_rollup = NA),
      after_i
    )
  }
  
  bind_rows(out)
}

manual_check <- check_manual_rollups(
  before = Data_taxonomy,
  after = Data_taxonomy2,
  rules = manual_rollup_rules
)
manual_check %>%
  arrange(SiteID, target_taxon, dataset, Taxon)

manual_check %>%
  group_by(SiteID, target_rank, target_taxon, dataset) %>%
  summarise(
    n_taxa = n_distinct(Taxon),
    taxa = paste(sort(unique(Taxon)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(SiteID, target_taxon, dataset) %>%  print(n=200)


sum(Data_taxonomy2$Value)
sum(Data_taxonomy$Value)
sum(Data$Value)

# check if there are repeated Taxon in SiteID-Date combinations
Data_taxonomy2 %>%
  count(SiteID, Date,  Taxon, name = "n") %>%
  filter(n > 1)

# Run harm function -------------------------------------------------------


# run function
# 2 steps :
Data_harm <- harmonize_taxonomy_two_stage(
  Data_taxonomy2,
  site_cols = "SiteID",
  date_col = "Date",
  abundance_col = "Value",
  stage1_finer_abundance_threshold = 0.8,   # stage1_finer_abundance_threshold = X: Operates at sample level: If finer taxa represent at least x% of abundance in the SAMPLE, it keeps the finer taxa and flags the coarse records.Otherwise it rolls up to the coarse level.
  stage2_finer_abundance_threshold = 0.8,   # stage2_finer_abundance_threshold = X: If finer taxa represent at least x% of abundance, it keeps the finer taxa and flags the coarse records.If finer taxa represent less than x%, it rolls up to the coarse genus or family.
  stage2_finer_year_threshold = 0.8         # stage2_finer_year_threshold = X: If finer taxa occur in at least x% of the years where that lineage is observed, it keeps the finer taxa and flags the coarse records. Otherwise it rolls up to the coarse level.
)

# check if there are repeated Taxon in SiteID-Date combinations
Data_harm %>%
  count(SiteID, Date,  harmonized_taxon, name = "n") %>%
  filter(n > 1)

# combine repeated taxa
Data_harm_final <- Data_harm %>%
  group_by(
    SiteID, Date, .harm_year,
    across(all_of(rank_cols)),
    harmonized_rank,
    harmonized_taxon
  ) %>%
  summarise(
    Value = sum(Value, na.rm = TRUE),
    flag_sample_coarse_overridden = any(flag_sample_coarse_overridden, na.rm = TRUE),
    flag_site_coarse_overridden = any(flag_site_coarse_overridden, na.rm = TRUE),
    flag_final_rollup = any(flag_final_rollup, na.rm = TRUE),
    flag_coarse_overridden = any(flag_coarse_overridden, na.rm = TRUE),
    drop_candidate = any(drop_candidate, na.rm = TRUE),
    n_finer = if (all(is.na(n_finer))) NA_integer_ else max(n_finer, na.rm = TRUE),
    flag_reassigned_unique_finer = any(flag_reassigned_unique_finer, na.rm = TRUE),
    drop_candidate_final = any(drop_candidate_final, na.rm = TRUE),
    .groups = "drop"
  )

# check if there are mixed drop_candidate_final per duplicated taxa
Data_harm %>%
  group_by(SiteID, Date, harmonized_taxon) %>%
  summarise(
    n_drop_final = n_distinct(drop_candidate_final),
    .groups = "drop"
  ) %>%
  filter(n_drop_final > 1) # must be 0


# check if there are repeated Taxon in SiteID-Date combinations
Data_harm_final %>%
  count(SiteID, Date,  harmonized_taxon, name = "n") %>%
  filter(n > 1)


Data_harm %>%
  filter(SiteID == "ALX001",
         Date == as.Date("2023-09-27"),
         harmonized_taxon == "Glossosoma") %>%
  select(
    SiteID, Date, Value,
    harmonized_rank, harmonized_taxon,
    all_of(rank_cols),
    flag_sample_coarse_overridden,
    flag_site_coarse_overridden,
    flag_final_rollup,
    drop_candidate,
    n_finer,
    flag_reassigned_unique_finer,
    drop_candidate_final
  )
Data_harm_final %>%
  filter(SiteID == "ALX001",
         Date == as.Date("2023-09-27"),
         harmonized_taxon == "Glossosoma") %>%
  select(
    SiteID, Date, Value,
    harmonized_rank, harmonized_taxon,
    all_of(rank_cols),
    flag_sample_coarse_overridden,
    flag_site_coarse_overridden,
    flag_final_rollup,
    drop_candidate,
    n_finer,
    flag_reassigned_unique_finer,
    drop_candidate_final
  )


# checks
(Data_harm_final %>% filter(family =="Leptoceridae") %>% print(n=50))

Data_harm_final %>% filter(n_finer==1)
(Data_harm_final %>% filter(flag_sample_coarse_overridden==T ))

(Data_harm_final %>% 
       filter(SiteID == "231")%>% 
       filter(family =="Perlodidae") %>% print(n=50))
(Data_harm_final %>% 
       filter(SiteID == "126")%>% 
       filter(genus =="Rhyacophila") %>% print(n=50))
(Data_harm_final %>% 
       filter(SiteID == "126")%>% 
       filter(family =="Polycentropodidae") %>% print(n=50))

(Data_harm_final %>% 
            filter(SiteID == "BNP-Johnston2")%>% 
            filter(order =="Plecoptera") %>% print(n=50))

Data_harm_final %>% 
  filter(SiteID=="S183") %>%
  group_by(SiteID) %>%
  summarise(Abund = sum(Value))

testdata %>% 
  filter(SiteID=="S183") %>%
  group_by(SiteID) %>%
  summarise(Abund = sum(Value))

# check abund & richness loss
Data_harm_final %>%
  summarise(
    total_abundance = sum(Value, na.rm = TRUE),
    abundance_removed = sum(
      Value[drop_candidate_final],
      na.rm = TRUE
    ),
    prop_removed =
      abundance_removed / total_abundance
  )

Data_harm_final %>%
  summarise(
    total_taxa = n_distinct(harmonized_taxon),
    taxa_removed = n_distinct(
      harmonized_taxon[drop_candidate_final]
    )
  )
Data_harm_final%>%
  filter(drop_candidate_final=="FALSE") %>%
  group_by(harmonized_rank) %>%
  summarise(n.taxa = length(unique(harmonized_taxon)),
            Abund = sum(Value))
dim(Data_harm_final)

# re-match site-specific columns
Data_harm2 <- Data_harm_final %>%
  left_join(Data_SampleLevel, by = c("SiteID", "Date")) %>%
  select(
    any_of(names(Data)),  # columns in the same order as Data
    everything()          # then keep remaining columns
  )
dim(Data_harm2)
dim(Data)
length(unique(Data$SiteID))
length(unique(Data_harm2$SiteID))

Data_harm2 <- Data_harm2 %>%
  select(BioticGroup, Country, Lake.river,Sampled.habitat, SiteID, 
         SiteName, WaterBody, Lon, Lat, Day, Month, Year,
         Date,
         Taxon.harm = harmonized_taxon, Value,ValueType,Unit, FishPresence, ReferenceCondition,
         source_file, climate,
         Rank.harm = harmonized_rank, species:phylum, drop_candidate_final)


# check if there are repeated Taxon in SiteID-Date combinations
Data_harm2 %>%
  count(SiteID, Date,  Taxon.harm, name = "n") %>%
  filter(n > 1)


write_parquet(Data_harm2, "TaxonomicHarm/BenthicMacroinvertebrate/BenthicMacroinvertebrates_HarmonizedTaxonomy_2steps_tr80.parquet")

sum(Data$Value)
sum(Data_harm2$Value)


# check proportion  deleted ind ----------------------------------




Data_harm2 <- read_parquet("TaxonomicHarm/BenthicMacroinvertebrate/BenthicMacroinvertebrates_HarmonizedTaxonomy_2steps_v5_tr80.parquet")

# check if there are repeated Taxon in SiteID-Date combinations
Data_harm2 %>%
  count(SiteID, Date,  Taxon.harm, name = "n") %>%
  filter(n > 1)
Data_harm2 %>%
  count(across(all_of(c(site_cols, "Taxon.harm"))), name = "n") %>%
  filter(n > 1)
Data_harm2 %>%
  group_by(across(all_of(c(site_cols, "Taxon.harm")))) %>%
  filter(n() > 1) %>%
  arrange(SiteID, Date, Taxon.harm)

#write.csv(Data_harm2, "BenthicMacroinvertebrates_HarmonizedTaxonomy_2steps_v5.csv")

Data_harm2 %>%
  filter(drop_candidate_final == TRUE) %>%
  select(Taxon.harm, Value) %>%
  arrange(desc(Value)) %>% print(n=200)


# check chironomidae

chir_sites <- Data_harm2 %>%
  filter(
    drop_candidate_final == TRUE,
    Taxon.harm == "Chironomidae",
    Value > 300
  ) %>%
  distinct(SiteID)

Data_harm2 %>%
  semi_join(chir_sites, by = "SiteID") %>%
  filter(family == "Chironomidae") %>%
  select(SiteID, Date, Taxon.harm, Value, drop_candidate_final) %>%
  arrange(SiteID, Date, Taxon.harm) %>% print(n=3000)

bae_sites <- Data_harm2 %>%
  filter(
    drop_candidate_final == TRUE,
    Taxon.harm == "Baetidae",
    Value > 300
  ) %>%
  distinct(SiteID)
Data_harm2 %>%
  semi_join(bae_sites, by = "SiteID") %>%
  filter(family == "Baetidae") %>%
  select(SiteID, Date, Taxon.harm, Value, drop_candidate_final) %>%
  arrange(SiteID, Date, Taxon.harm) %>% print(n=3000)

Limn_sites <- Data_harm2 %>%
  filter(
    drop_candidate_final == TRUE,
    Taxon.harm == "Limnephilidae",
    Value > 300
  ) %>%
  distinct(SiteID)
Data_harm2 %>%
  semi_join(Limn_sites, by = "SiteID") %>%
  filter(family == "Limnephilidae") %>%
  select(Country, SiteID, Date, Taxon.harm, Value, drop_candidate_final) %>%
  arrange(Country, SiteID, Date, Taxon.harm) %>% print(n=3000)


Dia_sites <- Data_harm2 %>%
  filter(
    drop_candidate_final == TRUE,
    Taxon.harm == "Diamesa",
    Value > 300
  ) %>%
  distinct(SiteID)
Data_harm2 %>%
  semi_join(Dia_sites, by = "SiteID") %>%
  filter(genus == "Diamesa") %>%
  select(Country, SiteID, Date, Taxon.harm, Value, drop_candidate_final) %>%
  arrange(Country, SiteID, Date, Taxon.harm) %>% print(n=3000)

Ephem_sites <- Data_harm2 %>%
  filter(
    drop_candidate_final == TRUE,
    Taxon.harm == "Ephemerellidae",
    Value > 300
  ) %>%
  distinct(SiteID)
Data_harm2 %>%
  semi_join(Ephem_sites, by = "SiteID") %>%
  filter(family == "Ephemerellidae") %>%
  select(Country, SiteID, Date, Taxon.harm, Value, drop_candidate_final) %>%
  arrange(Country, SiteID, Date, Taxon.harm) %>% print(n=3000)
# Summaries ---------------------------------------------------------------



#
Data.summarySiteLevel <- Data_harm2 %>%
  mutate(drop_candidate_final = coalesce(drop_candidate_final, FALSE)) %>%
  group_by(Country, Sampled.habitat, SiteID) %>%
  summarise(
    totAbund = sum(Value, na.rm = TRUE),
    totRichness = n_distinct(Taxon.harm, na.rm = TRUE),
    
    flaggedAbund = sum(
      if_else(drop_candidate_final, Value, 0),
      na.rm = TRUE
    ),
    
    flaggedRichness = n_distinct(
      Taxon.harm[drop_candidate_final],
      na.rm = TRUE
    ),
    
    ratio_Abund_flagged = if_else(
      totAbund > 0,
      flaggedAbund / totAbund,
      NA_real_
    ),
    
    ratio_Richness_flagged = if_else(
      totRichness > 0,
      flaggedRichness / totRichness,
      NA_real_
    ),
    
    .groups = "drop"
  )

hist(Data.summarySiteLevel$ratio_Abund_flagged, breaks = 10)
summary(Data.summarySiteLevel$ratio_Abund_flagged)


Data.summarySampleLevel <- Data_harm2 %>%
  mutate(drop_candidate_final = coalesce(drop_candidate_final, FALSE)) %>%
  group_by(Country, Sampled.habitat, SiteID, Date) %>%
  summarise(
    totAbund = sum(Value, na.rm = TRUE),
    totRichness = n_distinct(Taxon.harm, na.rm = TRUE),
    
    flaggedAbund = sum(
      if_else(drop_candidate_final, Value, 0),
      na.rm = TRUE
    ),
    
    flaggedRichness = n_distinct(
      Taxon.harm[drop_candidate_final],
      na.rm = TRUE
    ),
    
    ratio_Abund_flagged = if_else(
      totAbund > 0,
      flaggedAbund / totAbund,
      NA_real_
    ),
    
    ratio_Richness_flagged = if_else(
      totRichness > 0,
      flaggedRichness / totRichness,
      NA_real_
    ),
    
    .groups = "drop"
  )

hist(Data.summarySampleLevel$ratio_Abund_flagged, breaks = 10)
summary(Data.summarySampleLevel$ratio_Abund_flagged)
Data.summarySampleLevel %>%
  filter(ratio_Abund_flagged>0.50) %>% print(n=30)
summary(Data.summarySampleLevel$ratio_Richness_flagged)


Data.summarySampleLevel %>%
  filter(ratio_Richness_flagged>0.50) %>% print(n=100)

library(ggplot2)
ggplot(Data.summarySampleLevel, aes(y=ratio_Abund_flagged, x=SiteID)) +
  geom_point()+
  ggtitle("Sample level")+
  ylab("Ratio abundance that would be lost")
ggplot(Data.summarySiteLevel, aes(y=ratio_Abund_flagged, x=SiteID)) +
  coord_cartesian(ylim = c(0, 0.1))+
  geom_point()+
  ggtitle("Site level")+
  ylab("Ratio abundance that would be lost")

library(patchwork)
# Plot 1
p1 <- ggplot(Data.summarySampleLevel, aes(y = ratio_Abund_flagged, x = SiteID)) +
  geom_point() +
  ggtitle("Sample level") +
  ylab("Ratio abundance that would be lost")

# Plot 2
p2 <- ggplot(Data.summarySiteLevel, aes(y = ratio_Abund_flagged, x = SiteID)) +
  coord_cartesian(ylim = c(0, 0.1)) +
  geom_point() +
  ggtitle("Site level") +
  ylab("Ratio abundance that would be lost")

# Combine
combined_plot <- p1 / p2

# Save as TIFF
ggsave("Combined_plot_genusLevel_2steps_v2.tiff",
       plot = combined_plot,
       width = 20, height = 8,
       dpi = 300,
       compression = "lzw")

AA <- Data_harm2 %>%
  filter(SiteID=="126")
write.csv2(AA, "Check_harm_126.csv")

Data_harm%>%
  filter(drop_candidate_final=="FALSE") %>%
  group_by(harmonized_rank) %>%
  summarise(n.taxa = length(unique(harmonized_taxon)),
            Abund = sum(Value))


Data_harm2 %>%
  filter(candida)
  filter(Rank.harm =="class") %>%
  select(Taxon.harm) %>% distinct()

Data_harm2 %>%
  filter(drop_candidate_final=="FALSE") %>%
  filter(Rank.harm =="phylum") %>%
  select(Taxon.harm) %>% distinct()

Data_harm2 %>%
  filter(drop_candidate_final=="FALSE") %>%
  filter(Rank.harm =="subphylum") %>%
  select(Taxon.harm) %>% distinct()

Data_harm2 %>%
  filter(drop_candidate_final=="FALSE") %>%
  filter(Rank.harm =="class") %>%
  select(Taxon.harm) %>% distinct()

Data_harm2 %>%
  filter(drop_candidate_final=="FALSE") %>%
  filter(Rank.harm =="subclass") %>%
  select(Taxon.harm) %>% distinct()

Data_harm2 %>%
  filter(drop_candidate_final=="FALSE") %>%
  filter(Rank.harm =="order") %>%
  select(Taxon.harm) %>% distinct()



