library(dplyr)

Data <- read_parquet( "Combined_dataset/all/Zooplankton_all.parquet")
HarmonizationTable <- read.table("TaxonomicHarm/Zooplankton/Zooplankton_Taxonomy_table_QualityChecked.csv",h=T, sep=",", encoding="UTF-8")

HarmonizationTable <- HarmonizationTable %>%
  select(submittedName, species,genus,family,order, subclass,class,phylum)

dim(Data)
Data <- Data %>%
  left_join(HarmonizationTable, by=c("Taxon" = "submittedName"))



Data %>%
  filter(is.na(phylum)) %>% 
  select(Taxon) %>% unique()


write_parquet( Data, "TaxonomicHarm/Zooplankton/Zooplankton_all_withTaxonomy.parquet")

library(arrow)
library(dplyr)
source("FunctionsTaxonomicHarmonization_v3.R")

Data <- read_parquet( "TaxonomicHarm/Zooplankton/Zooplankton_all_withTaxonomy.parquet")


summary(Data$Stage)

Data %>%
  select(ValueType, Unit) %>%
  distinct()  


Data %>%
  filter(Stage %in% c("Copepodite","copepodite+adult", "Nauplii"))%>%
  select(Taxon,Stage) %>% distinct() %>% print(n=50)


# Filter original data
# e.g. remove larval stage
Data <- Data %>%
  filter(!Stage %in% c("Copepodite","copepodite+adult", "Nauplii"))
  
# e.g. discard taxa not resolved: without phylum 
Data %>%
  filter(is.na(phylum)) %>%
  select(Taxon) %>% distinct()

Data <- Data %>%
  filter(!is.na(phylum))

# e.g. exclude Rotifers
Data <- Data %>%
  filter(!phylum %in% c("Rotifera"))



#
# what to do with dataset that do not specify larval stage???
#


# Select relevant columns
Data_taxonomy <- Data %>%
  select(SiteID, Date, Taxon, species, genus, family, order, subclass, class, phylum, Value)


# Aggregate 
Data_taxonomy <- Data_taxonomy %>%
  group_by(SiteID, Date, Taxon, species, genus, family, order, subclass, class, phylum) %>%
  summarise(Value= sum(Value))%>% ungroup()

Data_SampleLevel <- Data %>%
   select(BioticGroup, Lake.river,Sampled.habitat,
                    Country,SiteID,SiteName,        
                    WaterBody,Lon,Lat, Day,Month,
                    Year,Date,
                    ValueType,Unit,             
                    FishPresence,ReferenceCondition, climate, source_file ) %>%
  distinct()

sites <- unique(Data_taxonomy$SiteID)

# run function
Data_harm <- harmonize_taxonomy(
  df = Data_taxonomy, # [Data_taxonomy$SiteID %in% sites[1:3],],
  decision_cols = "SiteID",
  output_site_cols = "SiteID",
  date_col = "Date",
  abundance_col = "Value",
  threshold = 0.90
)


# re-match site-specific columns
Data_harm2 <- Data_harm %>%
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
  select(BioticGroup, Country, Lake.river,Sampled.habitat, SiteID : Date,
         Taxon.harm = harmonized_taxon, Value:climate,
         Rank.harm = harmonized_rank, species:phylum,, flag_coarse_taxon, flag_no_finer, n_finer)


write_parquet(Data_harm2, "TaxonomicHarm/Zooplankton/Zooplankton_HarmonizedTaxonomy_SiteLevel.parquet")
write.csv2(Data_harm2, "TaxonomicHarm/Zooplankton/Zooplankton_HarmonizedTaxonomy_SiteLevel.csv")

