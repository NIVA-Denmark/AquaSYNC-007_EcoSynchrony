## Analyze community synchrony in fish

library(arrow)
library(codyn)


rm(list=ls())

dat_long <- read_parquet("../fish_harmonized.parquet")

## make data frame of site information, save for later
site_info <- dat_long %>%
  select(!(species:harmonized_taxon)) %>%
  select(!c(Day:Date)) %>%
  select(!Value) %>%
  distinct()

## Do analyses of community synchrony across all times --------------------------------------------

sites <- unique(as.character(dat_long$SiteID))
loreau <- rep(NA, length(sites))
gross <- rep(NA, length(sites))
ntaxa <- rep(NA, length(sites))

minspan <- 6 #minimum span of years to run analysis
maxpmiss <- 1/minspan #maximum proportion of missing years
mintaxa <- 2 #minimum number of species

for(ii in 1:length(sites)){
  
  if(ii %in% c(184, 191)){next} #Ignore sites with multiple obs/year under old harmonization script
  
  tmp <- dat_long[dat_long$SiteID==sites[ii],]
  sitespan <- max(tmp$Year) - min(tmp$Year) + 1
  pmiss <- 1-length(unique(tmp$Year))/sitespan
  ntax <- length(unique(tmp$harmonized_taxon))
  ntaxa[ii] <- ntax
  
  if(sitespan > minspan & pmiss < maxpmiss & ntax>=mintaxa){ #if data requirements met, calculate synchrony metrics
    loreau[ii] <- codyn::synchrony(tmp, time.var="Year", species.var="harmonized_taxon", 
                                   abundance.var="Value", metric="Loreau")
    gross[ii] <- codyn::synchrony(tmp, time.var="Year", species.var="harmonized_taxon", 
                                  abundance.var="Value", metric="Gross")
  }
}

syncdf <- data.frame(SiteID=sites,
                     ntaxa=ntaxa,
                     loreau=loreau,
                     gross=gross)

syncdf <- left_join(syncdf, site_info)

## Exploratory plotting

hist(syncdf$loreau, main="Loreau synchrony")
hist(syncdf$gross, main="Gross synchrony")

plot(syncdf$loreau, syncdf$gross, ylab="Gross", xlab="Loreau")
abline(lm(syncdf$gross ~ syncdf$loreau), col="red")

plot(syncdf$Lat, syncdf$loreau, ylab="Synchrony", xlab="Latitude", main="Loreau")
abline(lm(syncdf$loreau ~ syncdf$Lat), col="red")

plot(syncdf$Lat, syncdf$gross, ylab="Synchrony", xlab="Latitude", main="Gross")
abline(lm(syncdf$gross ~ syncdf$Lat), col="red")

plot(syncdf$ntaxa, syncdf$loreau, ylab="Synchrony", xlab="N taxa", main="Loreau")
abline(lm(syncdf$loreau ~ syncdf$ntaxa), col="red")

plot(syncdf$ntaxa, syncdf$gross, ylab="Synchrony", xlab="N taxa", main="Gross")
abline(lm(syncdf$gross ~ syncdf$ntaxa), col="red")



## Do analyses of synchrony using time windows ----------------------------------------------------

sites <- unique(as.character(dat_long$SiteID))
SiteID <- NULL
wstart <- NULL
loreau <- NULL
gross <- NULL
ntaxa <- NULL

minspan <- 6 #minimum span of years to run analysis
maxpmiss <- 1/minspan #maximum proportion of missing years
mintaxa <- 2 #minimum number of species
wwidth <- 6 #moving window width, years
wstep <- 1 #steps between moving windows, years

for(ii in 1:length(sites)){
  
  if(ii %in% c(184, 191)){next} #Ignore sites with multiple obs/year under old harmonization script
  
  tmp <- dat_long[dat_long$SiteID==sites[ii],]
  sitespan <- max(tmp$Year) - min(tmp$Year) + 1
  pmiss <- 1-length(unique(tmp$Year))/sitespan
  ntax <- length(unique(tmp$harmonized_taxon)) #TODO put inside windowing

  
  #site has only minimum number of years; compute single values
  if(sitespan==minspan & pmiss < maxpmiss & ntax>=mintaxa){ #if data requirements met, calculate synchrony metrics
    loreau <- c(loreau, codyn::synchrony(tmp, time.var="Year", species.var="harmonized_taxon", 
                                   abundance.var="Value", metric="Loreau"))
    gross <- c(gross, codyn::synchrony(tmp, time.var="Year", species.var="harmonized_taxon", 
                                  abundance.var="Value", metric="Gross"))
    wstart <- c(wstart, min(tmp$Year))
    SiteID <- c(SiteID, sites[ii])
    ntaxa <- c(ntaxa, ntax)
  }
  #site has more years > minimum; compute windowed values
  if(sitespan>minspan & pmiss < maxpmiss & ntax>=mintaxa){ #if data requirements met, calculate synchrony metrics
    wstarts <- seq(from=min(tmp$Year), to=max(tmp$Year)-wwidth+1, by=wstep)
    for(jj in 1:length(wstarts)){
      tmp.ww <- tmp[tmp$Year %in% wstarts[jj]:(wstarts[jj]+(wwidth-1)),]
      pmiss.ww <- 1-length(unique(tmp.ww$Year))/wwidth
      ntax.ww <- length(unique(tmp.ww$harmonized_taxon))
      if(pmiss.ww < maxpmiss & ntax.ww>=mintaxa){ #check window-level data requirements
        loreau <- c(loreau, codyn::synchrony(tmp.ww, time.var="Year", species.var="harmonized_taxon", 
                                             abundance.var="Value", metric="Loreau"))
        gross <- c(gross, codyn::synchrony(tmp.ww, time.var="Year", species.var="harmonized_taxon", 
                                           abundance.var="Value", metric="Gross"))
        wstart <- c(wstart, wstarts[jj])
        SiteID <- c(SiteID, sites[ii])
        ntaxa <- c(ntaxa, ntax.ww)
      }
    }
  }
}

syncdf.ww <- data.frame(SiteID = SiteID,
                        wstart = wstart,
                        ntaxa = ntaxa,
                        loreau = loreau,
                        gross = gross)


syncdf.ww <- left_join(syncdf.ww, site_info)


## exploratory plotting

plot(syncdf.ww$loreau, syncdf.ww$gross)

plot(syncdf.ww$wstart, syncdf.ww$loreau)
abline(lm(syncdf.ww$loreau ~ syncdf.ww$wstart), col="red")

plot(syncdf.ww$wstart, syncdf.ww$gross)
abline(lm(syncdf.ww$gross ~ syncdf.ww$wstart), col="red")
