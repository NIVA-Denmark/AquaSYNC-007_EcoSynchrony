## Explore data and harmonize taxonomy for AquaSync Fish Data

rm(list=ls())

library(arrow)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dat.raw <- read_parquet("../Fish.parquet")

table(dat.raw$Lake.river)
table(dat.raw$Sampled.habitat)
table(dat.raw$Country)
table(dat.raw$SiteName)
table(dat.raw$Taxon)
table(dat.raw$Taxon[dat.raw$Country=="USA"])
table(dat.raw$Taxon[dat.raw$Country=="Norway"])

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



## Sampling through time - USA --------------------------------------------------------
dat.usa <- dat.raw[dat.raw$Country=="USA",]
spp <- unique(as.character(dat.usa$Taxon))

#sampling through time
par(mar=c(4.1,8,2,1))
plot(NA,NA, ylim=c(0.5, length(spp)+0.5), xlim=range(dat.usa$Date),
     xaxt="n", yaxt="n", xlab="Time", ylab="")
axis(2, at=1:length(spp), labels=spp, las=2)
axis(1, at=pretty(dat.usa$Date), labels=as.character(as.Date(pretty(dat.usa$Date))))
abline(h=1:length(spp))
for(ii in 1:length(spp)){
  dat.sub <- dat.usa[dat.usa$Taxon==spp[ii],]
  points(dat.sub$Date, rep(ii, nrow(dat.sub)), pch=16, col="blue")
}
mtext("USA Fishes")

#sampling by location
dat.usa$SiteName <- factor(dat.usa$SiteName)


par(mar=c(12,8,2,1))
plot(NA,NA, ylim=c(0.5, length(spp)+0.5), xlim=c(1,length(unique(dat.usa$SiteName))),
     xaxt="n", yaxt="n", xlab="", ylab="")
axis(2, at=1:length(spp), labels=spp, las=2)
axis(1, at=1:length(unique(dat.usa$SiteName)), labels=unique(dat.usa$SiteName), las=2)
abline(h=1:length(spp))
abline(v=1:length(unique(dat.usa$SiteName)))
for(ii in 1:length(spp)){
  dat.sub <- dat.usa[dat.usa$Taxon==spp[ii],]
  points(as.numeric(dat.sub$SiteName), rep(ii, nrow(dat.sub)), pch=16, col="blue")
}
mtext("USA Fishes")



## Sampling through time - Norway --------------------------------------------------------
dat.nor <- dat.raw[dat.raw$Country=="Norway",]
spp <- unique(as.character(dat.nor$Taxon))

par(mar=c(4.1,8,2,1))
plot(NA,NA, ylim=c(0.5, length(spp)+0.5), xlim=range(dat.nor$Date),
     xaxt="n", yaxt="n", xlab="Time", ylab="")
axis(2, at=1:length(spp), labels=spp, las=2)
axis(1, at=pretty(dat.nor$Date), labels=as.character(as.Date(pretty(dat.nor$Date))))
abline(h=1:length(spp))
for(ii in 1:length(spp)){
  dat.sub <- dat.nor[dat.nor$Taxon==spp[ii],]
  points(dat.sub$Date, rep(ii, nrow(dat.sub)), pch=16, col="blue")
}
mtext("Norway Fishes")


dat.nor$SiteName <- factor(dat.nor$SiteName)

par(mar=c(12,8,2,1))
plot(NA,NA, ylim=c(0.5, length(spp)+0.5), xlim=c(1,length(unique(dat.nor$SiteName))),
     xaxt="n", yaxt="n", xlab="", ylab="")
axis(2, at=1:length(spp), labels=spp, las=2)
axis(1, at=1:length(unique(dat.nor$SiteName)), labels=unique(dat.nor$SiteName), las=2)
abline(h=1:length(spp))
abline(v=1:length(unique(dat.nor$SiteName)))
for(ii in 1:length(spp)){
  dat.sub <- dat.nor[dat.nor$Taxon==spp[ii],]
  points(as.numeric(dat.sub$SiteName), rep(ii, nrow(dat.sub)), pch=16, col="blue")
}
mtext("Norway Fishes")

