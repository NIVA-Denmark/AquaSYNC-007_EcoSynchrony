## Explore data and harmonize taxonomy for AquaSync Fish Data

rm(list=ls())

library(taxize)
library(arrow)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



### Load data and helper functions, some initial checking -----------------------------------------

source("FunctionsTaxonomicHarmonization.R")
source("FunctionsTemporalHarmonization.R")

dat_raw <- read_parquet("../Fish_all.parquet")
tax_table <- read.csv("../fish_taxonomyTable_gbif.csv")

#option to load benthic macroinverts to check formatting consistency
#bm1 <- read_parquet("../BenthicMacroinvertebrates_all.parquet")
#bm2 <- read_parquet("../BenthicMacroinvertebrates_all_withTaxonomy.parquet")

dat_raw$Taxon <- as.character(dat_raw$Taxon)
dat_raw$ValueType <- as.character(dat_raw$ValueType)

# For species listed with common name, assign a scientific name
dat_raw$Taxon_orig <- dat_raw$Taxon #keep original
dat_raw$Taxon[dat_raw$Taxon=="Arctic grayling"] <- "Thymallus arcticus"
dat_raw$Taxon[dat_raw$Taxon=="Brook_trout"] <- "Salvelinus fontinalis"
dat_raw$Taxon[dat_raw$Taxon=="Char"] <- "Salvelinus alpinus"
dat_raw$Taxon[dat_raw$Taxon=="Chinook"] <- "Oncorhynchus tshawytscha"
dat_raw$Taxon[dat_raw$Taxon=="Coho"] <- "Oncorhynchus kisutch"
dat_raw$Taxon[dat_raw$Taxon=="Dolly Varden"] <- "Salvelinus malma"
dat_raw$Taxon[dat_raw$Taxon=="Perch"] <- "Perca fluviatilis"
dat_raw$Taxon[dat_raw$Taxon=="Phoxinus"] <- "Phoxinus phoxinus"
dat_raw$Taxon[dat_raw$Taxon=="slimy sculpin"] <- "Cottus cognatus"
dat_raw$Taxon[dat_raw$Taxon=="Stickleback"] <- "Gasterosteidae"
dat_raw$Taxon[dat_raw$Taxon=="Trout"] <- "Salmo trutta"

dat_raw <- dat_raw[dat_raw$Taxon != "Unidentified Unidentified",]

dat_raw$Unit <- trimws(dat_raw$Unit)

#fix some value types

dat_raw$ValueType[grepl("fork length", dat_raw$ValueType, ignore.case=TRUE)] <- "Fork length"
dat_raw$ValueType[grepl("total length", dat_raw$ValueType, ignore.case=TRUE)] <- "Total length"

### Exclude or otherwise deal with ValueTypes that will not analyzed ------------------------------
#origN <- nrow(dat_taxa)
#dat_taxa <- dat_taxa[!grepl("length", dat_taxa$ValueType),]
dat_raw <- dat_raw[dat_raw$ValueType != "PresenceAbsence",] #possibly include for richness-based analyses, but only 12 occurrences 
#table(as.character(dat_taxa$ValueType))


### Aggregating certain value types as appropriate ------------------------------------------------

orig_cols <- colnames(dat_raw)

## For Salmo salar, "Total weight" there are sometimes two rows for wild/hatchery. Combine.

dat_totalWeight_ss <- dat_raw[dat_raw$ValueType=="Total weight" & dat_raw$Taxon == "Salmo salar",]

dat_totalWeight_ss_sum <- dat_totalWeight_ss %>%
  group_by(SiteID, Taxon, Date) %>%
  summarize(
    BioticGroup = first(BioticGroup),
    Lake.river = first(Lake.river),
    Sampled.habitat = first(Sampled.habitat),
    Country = first(Country),
    SiteName = first(SiteName),
    WaterBody = first(WaterBody),
    Lon = first(Lon),
    Lat = first(Lat),
    SampleID = first(SampleID),
    Day = first(Day),
    Month = first(Month),
    Year = first(Year),
    Value = sum(Value),
    ValueType = "Total weight",
    Unit = first(Unit),
    FishPresence = "Yes",
    ReferenceCondition = first(ReferenceCondition),
    source_file = first(source_file),
    climate = first(climate),
    Taxon_orig = first(Taxon_orig)
  )

dat_totalWeight_ss_sum <- dat_totalWeight_ss_sum[,match(orig_cols,colnames(dat_totalWeight_ss_sum))] 

## Data with ValueType=="Weight" are individual fish, sum by site and sampling event

dat_weight <- dat_raw[dat_raw$ValueType=="Weight",]

dat_weight_sum <- dat_weight %>%
  group_by(SiteID, Taxon, Date) %>%
  summarize(
    BioticGroup = first(BioticGroup),
    Lake.river = first(Lake.river),
    Sampled.habitat = first(Sampled.habitat),
    Country = first(Country),
    SiteName = first(SiteName),
    WaterBody = first(WaterBody),
    Lon = first(Lon),
    Lat = first(Lat),
    SampleID = first(SampleID),
    Day = first(Day),
    Month = first(Month),
    Year = first(Year),
    Value = sum(Value),
    ValueType = "Total weight",
    Unit = first(Unit),
    FishPresence = "Yes",
    ReferenceCondition = first(ReferenceCondition),
    source_file = first(source_file),
    climate = first(climate),
    Taxon_orig = first(Taxon_orig)
  )

dat_weight_sum <- dat_weight_sum[,match(orig_cols,colnames(dat_weight_sum))] 


## Data with ValueType=="Length" are individual fish, turn into count by site, sampling date, and taxon

dat_length <- dat_raw[grepl("length", dat_raw$ValueType, ignore.case=TRUE),]

dat_length_count <- dat_length %>%
  group_by(SiteID, Taxon, Date) %>%
  summarize(
    BioticGroup = first(BioticGroup),
    Lake.river = first(Lake.river),
    Sampled.habitat = first(Sampled.habitat),
    Country = first(Country),
    SiteName = first(SiteName),
    WaterBody = first(WaterBody),
    Lon = first(Lon),
    Lat = first(Lat),
    SampleID = first(SampleID),
    Day = first(Day),
    Month = first(Month),
    Year = first(Year),
    Value = n(),
    ValueType = "Count",
    Unit = "ind",
    FishPresence = "Yes",
    ReferenceCondition = first(ReferenceCondition),
    source_file = first(source_file),
    climate = first(climate),
    Taxon_orig = first(Taxon_orig)
    )

dat_length_count <- dat_length_count[,match(orig_cols,colnames(dat_length_count))]



## Remove aggregated observations from raw data table and then add the aggregated data on the bottom

dat_agg <- dat_raw
dat_agg <- dat_agg[!(dat_agg$ValueType=="Total weight" & dat_agg$Taxon == "Salmo salar"),]
dat_agg <- dat_agg[dat_agg$ValueType!="Weight",]
dat_agg <- dat_agg[!grepl("length", dat_agg$ValueType, ignore.case=TRUE),]
dat_agg <- rbind(dat_agg, dat_totalWeight_ss_sum, dat_weight_sum, dat_length_count)

#check number of rows
#nrow(dat_raw) - nrow(dat_length) - nrow(dat_weight) - nrow(dat_totalWeight_ss) + nrow(dat_totalWeight_ss_sum) + nrow(dat_weight_sum) + nrow(dat_length_count)

rm(dat_length, dat_weight, dat_totalWeight_ss, dat_length_count, dat_totalWeight_ss_sum, dat_weight_sum)

### Select data type to analyze, as needed --------------------------------------------------------

## Check for multiple data types per site

multitype_data <- NULL

for(site in unique(dat_agg$SiteID)){
  
  tmp <- dat_agg[dat_agg$SiteID==site,]
  if(length(unique(as.character(tmp$ValueType)))>1){
    
    tmp <- tmp %>%
      group_by(SiteID, ValueType) %>%
      summarize(Count = n())
    
    multitype_data <- rbind(multitype_data, tmp)
    
  }
}


## Select which to use, favoring the data type with the longest record, or when equal the hierarchy:

multitype_toUse <- NULL

for(site in unique(multitype_data$SiteID)){
   tmp <- multitype_data[multitype_data$SiteID==site,]
   if(length(unique(tmp$Count))>1){
     multitype_toUse <- rbind(multitype_toUse, tmp[tmp$Count==max(tmp$Count),1:2])
   }
   else{
     if("CPUE" %in% tmp$ValueType & "Count" %in% tmp$ValueType){
       multitype_toUse <- rbind(multitype_toUse, tmp[tmp$ValueType=="CPUE",1:2])
     }
     if("Total weight" %in% tmp$ValueType & "Count" %in% tmp$ValueType){
       multitype_toUse <- rbind(multitype_toUse, tmp[tmp$ValueType=="Total weight",1:2])
     }
     if("Biomass" %in% tmp$ValueType & "Density" %in% tmp$ValueType){
       multitype_toUse <- rbind(multitype_toUse, tmp[tmp$ValueType=="Biomass",1:2])
     }
   }
}

#length(unique(multitype_data$SiteID)) == nrow(multitype_toUse)
#any(duplicated(multitype_toUse$SiteID))

#Identify and drop rows representing data types we will not use
droprow <- NULL

for(ii in 1:nrow(multitype_toUse)){
  droprow <- c(droprow,
               which(dat_agg$SiteID==multitype_toUse$SiteID[ii] & dat_agg$ValueType!=multitype_toUse$ValueType[ii]))
}

dat_agg <- dat_agg[-droprow,]


### Add taxonomy to data and export parquet file --------------------------------------------------

dat_taxa <- left_join(dat_agg, tax_table)
#write_parquet(dat_taxa, "../fish_all_withTaxonomy_individualsAggregated.parquet")




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

length(unique(dat_raw$SiteID)) #211 sites
length(unique(Data_temporal_harmonized$SiteID)) #210 sites--one gets dropped.

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


# #trying to understand warning messages <--- the causes for these were resolved upsteam and nolonger occur as of 2026-04-29
# Data_taxonomy_harmonized$SiteID[1133]; Data_taxonomy_harmonized$Date[1133]
# check <- Data_sample_level[Data_sample_level$SiteID=="Könkämäeno 1" & Data_sample_level$Date=="1986-09-05",]
# 
# Data_sample_level$SiteID[5714]; Data_sample_level$Date[5714]
# check <- Data_taxonomy_harmonized[Data_taxonomy_harmonized$SiteID=="Anxiety Ridge" & Data_taxonomy_harmonized$Date=="2016-08-05",]
# #Warnings arise because of some sites have multiple units of measure -- need to select the best one and drop others
# 
# check <- Data_harmonized[Data_harmonized$SiteID=="Könkämäeno 1",]
