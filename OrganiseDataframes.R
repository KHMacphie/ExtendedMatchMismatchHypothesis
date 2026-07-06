#rm(list=ls())
library(tidyverse)

#####################################################
#### Sorting data for mismatch hypothesis models ####
#####################################################

#### Import Data ####
birdphen <- read.csv("Bird_Phenology.csv")
nestlings <- read.csv("Nestlings.csv")
swap <- read.csv("clutchswaps.csv")
adults <- read.csv("Adults.csv")
cater <- read.csv("Branch_Beating.csv")

#Sort bird data link
nestlings$SYB <- paste(nestlings$site, nestlings$box, nestlings$year)
birdphen$SYB <- paste(birdphen$site, birdphen$box, birdphen$year)
swap$SYB <- paste(swap$site, swap$origin.nest, swap$year)
adults$SYB <- paste(adults$site, adults$box, adults$year)

#DF with relevant info
SYB <- data.frame(year=birdphen$year, site=birdphen$site, SYB=birdphen$SYB, species=birdphen$species, fki=birdphen$fki, hd.calc=birdphen$hd_1.45, hd.obs=birdphen$hatching_first_recorded)
SYB$siteyear <- paste(SYB$site,SYB$year)

#remove coal tits 
SYB <- SYB[-which(SYB$species=="coati"),]

# Hatch Date as observed
cater$siteday <- paste(cater$site, cater$date ,cater$year) # Beating days at each site in each year
SYB$calchd <- paste(SYB$site, SYB$hd.calc ,SYB$year) # Days at each site in each year with estimates hatch dates (may or may not be day of observation)
storebeat <- pmatch(SYB$calchd,cater$siteday,duplicates.ok = T)
SYB$beating <- cater$siteday[storebeat] # If the hatch estimate day was a beating day it wont be an NA
SYB$hd.calc <- as.integer(SYB$hd.calc)

SYB$hd.obs.all <- # working out the observed hatch date for those that only have estimated
  ifelse(is.na(SYB$hd.obs), # If there is no record of observed hatch date
         ifelse(is.na(SYB$beating),(SYB$hd.calc+1),SYB$hd.calc), # and there is no beating record for the estimated day then add one to estimated date, if there is a beating record use the estimated date
         SYB$hd.obs) # If there is an observed date use that

SYB$hd.dif <- SYB$hd.calc-SYB$hd.obs # Check the differences between observed and estimated
# table(SYB$hd.dif)
# # Most records with a bigger gap than -1 are from 2020 others I guess are from checks that were forgotten and then hatch date estimated


#remove any that had 2 nests in one season
length(names(which(table(SYB$SYB)==2)))
repnests <- c()
for(i in 1:length(names(which(table(SYB$SYB)==2)))){
  repnests[(i*2-1):(i*2)] <- which(SYB$SYB==c(names(which(table(SYB$SYB)==2)[i])))
}
SYB <- SYB[-c(repnests),]

# remove nests with NAs in fledging
NAnests <- unique(nestlings$SYB[which(is.na(nestlings$fledged)==T)])
nestsNA <- c()
for(i in 1:length(NAnests)){
  nestsNA[i] <- which(SYB$SYB==NAnests[i])
}
SYB <- SYB[-c(nestsNA),]


#remove clutch swap nests
storeswap <- pmatch(SYB$SYB, swap$SYB,duplicates.ok = T)
SYB$swap <- swap$destination.nest[storeswap]
SYB <- subset(SYB, is.na(swap))
SYB$swap <- NULL

#Remove any that dont have either FKI or hatch date (left with all nesting attempts)
SYB$fki <- as.integer(SYB$fki)
SYB <- subset(SYB, !is.na(fki)|!is.na(hd.obs.all)) 

#Match female to box, fill empties with other unique ID (for now remove multifemale nests)
adults <- subset(adults,season=="spring")
adults <- subset(adults,sex=="F")

# identify multi female nests rather than just females caught twice
fem <- adults[,c("ring","SYB")]
fem <- distinct(fem)

which(table(fem$SYB)>1)

multifems <- c()
for(i in 1:length(names(which(table(fem$SYB)>1)))){ #remove nests with multi females
  X <- c(which(fem$SYB==c(names(which(table(fem$SYB)>1)[i]))))  
  multifems <- c(multifems,X)
}

# remove nest(s) with multiple females
fem <- fem[-multifems,]
storefem <- pmatch(SYB$SYB,fem$SYB, duplicates.ok = T)
SYB$femring <- fem$ring[storefem]

SYB$fem <- ifelse(is.na(SYB$femring)==T,SYB$SYB,SYB$femring) #using SYB as unique femID for uncaught females



#Sort fledged data
nestlings$fledged <- as.factor(nestlings$fledged)
nests <- subset(nestlings, fledged=="1", select=c(SYB, fledged))
fledged <- aggregate(as.numeric(as.character(fledged))~SYB, nests, sum)
colnames(fledged)[2] <- "fledged"

storefledge <- pmatch(SYB$SYB,fledged$SYB, duplicates.ok = T)
SYB$fledged <- ifelse(is.na(pmatch(SYB$SYB,fledged$SYB, duplicates.ok = T)==T),0,fledged$fledged[storefledge])

storephen <- pmatch(SYB$SYB,birdphen$SYB)
SYB$fledgedphen <- birdphen$suc[storephen]
SYB <- SYB[-which(SYB$fledgedphen=="-999"),] #-999 signals reason to remove fledging success

# Reduce to just nests with a hatch date
SYBhd <- subset(SYB,!is.na(hd.obs.all)) 


#### Organising Bird and Cater so same site years in each ####
cater$year <- as.factor(cater$year) #year as a factor
cater$siteyear <- paste(cater$site, cater$year)
cater$treeID <- paste(cater$site, cater$tree)

#remove caterpillar data if no bird data
storeSYB <- pmatch(cater$siteyear, SYBhd$siteyear, duplicates.ok = T)
cater$birddat <- SYBhd$siteyear[storeSYB]
cater <- cater[-which(is.na(cater$birddat)==T),]

#set up for siteyear peak estimates within model
SY_id <- data.frame(siteyear_id=as.numeric(as.factor(cater$siteyear)),
                    siteyear=cater$siteyear,
                    site_id=as.numeric(as.factor(cater$site)),
                    site=cater$site,
                    year_id=as.numeric(as.factor(cater$year)),
                    year=cater$year)
SY_id <- distinct(SY_id)
SY_id <- arrange(SY_id,siteyear_id)

#remove birds with no caterpillar peak
storecaterSY <- pmatch(SYBhd$siteyear, cater$siteyear, duplicates.ok = T)
SYBhd$caterSY <- cater$siteyear[storecaterSY]
SYBhd <- subset(SYBhd,!is.na(caterSY))

#remove nests where hatch date was missed
SYBhd <- SYBhd[-c(which(SYBhd$hd.dif==-6),which(SYBhd$hd.dif==-4)),]

#centre cater and bird data on same date
mean(SYBhd$hd.obs.all)
mean(cater$date)
SYBhd$hd_cent <- SYBhd$hd.obs.all-146  
cater$date_cent <- cater$date-146

SYBhd$resid <- 1:nrow(SYBhd)
cater$resid <- 1:nrow(cater) 

#Binomal fledge data 
SYBhd$fledged.binom <- ifelse(SYBhd$fledged>=1,1,0)

#Truncated Poisson fledge data
SYBhd_p <- subset(SYBhd,fledged>0)

#Sorting order for caterpillar peak data for poisson part: length of N_siteyear of Poisson, contain equivalent cater row numbers
SY_id_p <- data.frame(siteyear_id=as.numeric(as.factor(SYBhd_p$siteyear)),
                      siteyear=SYBhd_p$siteyear,
                      site_id=as.numeric(as.factor(SYBhd_p$site)),
                      site=SYBhd_p$site,
                      year_id=as.numeric(as.factor(SYBhd_p$year)),
                      year=SYBhd_p$year)
SY_id_p <- distinct(SY_id_p)
SY_id_p <- arrange(SY_id_p,siteyear_id)
storetrunc <- pmatch(SY_id_p$siteyear,SY_id$siteyear)
SY_id_p$siteyear_id_c <- SY_id$siteyear_id[storetrunc]
SY_id_p$site_id_c <- SY_id$site_id[storetrunc]
SY_id_p$year_id_c <- SY_id$year_id[storetrunc]

rm(list=setdiff(ls(), c("SY_id","SY_id_p","cater","SYBhd","SYBhd_p")))
dev.off()
