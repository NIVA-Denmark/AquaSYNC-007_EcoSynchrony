# Example workflow for temporal and taxonomic harmonization
#
# This script:
# 1. Harmonizes sampling dates across sites
# 2. Harmonizes taxonomy within sites
#
# Note:
# Temporal and taxonomic harmonization are shown here as separate examples.
# If needed, the temporally harmonized output can be used to subset the
# dataset before taxonomic harmonization.

library(arrow)
library(dplyr)
library(ggplot2)

# Load input data ===========================================================

Data <- read_parquet("BenthicMacroinvertebrates_all_withTaxonomy.parquet")

# 1. Temporal harmonization ============================================================

source("FunctionsTemporalHarmonization.R")

# Keep only the columns needed for temporal harmonization.
Data_temporal <- Data %>%
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

# 2. Taxonomic harmonization ============================================================

source("FunctionsTaxonomicHarmonization.R")

# Remove records that are not suitable for abundance-based
# taxonomic harmonization.
Data_filtered <- Data %>%
  # Exclude presence/absence records.
  filter(ValueType != "Presence/Absence") %>%
  # Exclude unresolved taxa with missing phylum.
  filter(!is.na(phylum)) %>%
  # Exclude non-target groups such as fish, zooplankton, and algae.
  filter(!phylum %in% c("Chordata", "Branchiopoda", "Copepoda", "Ostracoda", "Charophyceae"))

# Optional check of available measurement types and units.
Data_filtered %>%
  select(ValueType, Unit) %>%
  distinct()

# Keep the columns required for taxonomic harmonization.
Data_taxonomy <- Data_filtered %>%
  select(
    SiteID, Date, Taxon_clean,
    species, genus, family, order, subclass, class, phylum,
    Value
  )

# Keep sample-level metadata to reattach after harmonization.
Data_sample_level <- Data_filtered %>%
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
    any_of(names(Data)),  # keep original column order where possible
    everything()          # then append any additional columns
  )