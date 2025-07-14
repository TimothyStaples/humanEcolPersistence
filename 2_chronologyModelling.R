# ################################# ####
# Author:    Timothy L Staples      ####
# Collaborators: John Pandolfi      #### 
#                Wolfgang Kiessling ####
# ################################# ####
# Global attributes & working directories ####

rm(list=ls())
setwd('/Users/uqtstapl/Dropbox/Tim/Post-doc/Research projects/novel_comms_legacy/prodCode/')

# Packages & functions ####

# source functions from 'functions' sub-folder
sapply(paste0("./functions/", list.files("./functions", pattern =".R")), source)

package.loader(c("rbacon", "readxl"))

# -------------------------- ####
# IMPORT RAW DATA ####

plant.record.df <- readRDS("./rawdata/processedRecords.rds")
plant.record.df <- droplevels(plant.record.df[plant.record.df$elementtype == "pollen",])

#           Get time-series location ####

# Aggregate sites into continental regions to separate Europe and North America

site.df <- plant.record.df[!duplicated(plant.record.df$siteid) & complete.cases(plant.record.df[,c("long","lat")]),]
#coordinates(site.coords) <- c("long", "lat")

# separate Europe & North America
site.df$REGION = NA

Eu = c(30, 90, -11, 50)
nthAm = c(30, 90, -167, -50)

site.df$REGION[site.df$lat >= Eu[1] & site.df$lat <= Eu[2] & site.df$long >= Eu[3] & site.df$long <= Eu[4]] = "Europe"
site.df$REGION[site.df$lat >= nthAm[1] & site.df$lat <= nthAm[2] & site.df$long >= nthAm[3] & site.df$long <= nthAm[4]] = "North America"

# remove sites not in Europe or Nth America (poor sampling)
site.df <- droplevels(site.df[site.df$REGION %in% c("Europe", "North America"),])

# remove records not in our continental regions, and samples outside of 25000 ybp
plant.record.df <- droplevels(plant.record.df[plant.record.df$siteid %in% site.df$siteid,])
plant.record.df <- droplevels(plant.record.df[plant.record.df$age <= 25100, ])

# Look for duplicate samples
dupeSamps <- paste(plant.record.df$collunitid,
                   plant.record.df$sampleid, 
                   plant.record.df$variablename, 
                   plant.record.df$value, sep=".")

plant.record.df <- droplevels(plant.record.df[!duplicated(dupeSamps),])

plant.record.df = merge(plant.record.df, site.df[,c("siteid", "REGION")],
                         by.x="siteid", by.y="siteid", all.x=TRUE, all.y=FALSE, sort=FALSE)

#           Pollen production adjustments ####

# read in PPE
PPE = as.data.frame(read_excel("./rawdata/RPP_Dataset_v2_Table_5_6.xlsx"))
PPE = PPE[,c(1,2,8,12,16)]
colnames(PPE) = c("type","taxon","ppeAm", "ppeEur", "ppeNthHem")

PPE$taxon[PPE$taxon == "Sambucus nigra-type"] = "Sambucus nigra"
PPE$taxon[PPE$taxon == "Asteraceae"] = "Compositae"
PPE$taxon[PPE$taxon == "Fabaceae"] = "Leguminosae"

# pair as best as possible in the following order - Continent + genus, genus, continent + family, family
taxTable = plant.record.df[,c("genus", "family", "REGION")]
taxTable = taxTable[!duplicated(taxTable[,c("genus", "REGION")]),]

eurTable = taxTable[taxTable$REGION == "Europe",]
eurTable = merge(eurTable, PPE[,c("taxon", "ppeEur", "ppeNthHem")], by.x="genus", by.y = "taxon",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)
colnames(eurTable)[-(1:3)] = c("ppeEurGen", "ppeHemGen")

eurTable = merge(eurTable, PPE[,c("taxon", "ppeEur", "ppeNthHem")], by.x="family", by.y = "taxon",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)
eurPPEraw = eurTable

# sample sizes
firstNona = apply(eurTable[,-(1:3)], 1, function(x){which(!is.na(x))[1]})
table(firstNona)

eurPPE = eurTable[,-(1:3)][cbind(1:nrow(eurTable), apply(eurTable[,-(1:3)], 1, function(x){which(!is.na(x))[1]}))]

eurTable = cbind(eurTable[,1:3],
                 PPE = eurPPE)
eurTable[is.na(eurTable$PPE),]

# Nth Am
amTable = taxTable[taxTable$REGION == "North America",]
amTable = merge(amTable, PPE[,c("taxon", "ppeAm", "ppeNthHem")], by.x="genus", by.y = "taxon",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)
colnames(amTable)[-(1:3)] = c("ppeamGen", "ppeHemGen")

amTable = merge(amTable, PPE[,c("taxon", "ppeAm", "ppeNthHem")], by.x="family", by.y = "taxon",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)
amPPEraw = amTable

# sample sizes

# taxa
firstNona = apply(amTable[,-(1:3)], 1, function(x){which(!is.na(x))[1]})

amPPE = amTable[,-(1:3)][cbind(1:nrow(amTable), apply(amTable[,-(1:3)], 1, function(x){which(!is.na(x))[1]}))]

amTable = cbind(amTable[,1:3],
                 PPE = amPPE)
amTable[is.na(amTable$PPE),]

neoPPE = rbind(eurTable, amTable)
neoPPE = neoPPE[!is.na(neoPPE$family),]

# genSub
genPPE = merge(plant.record.df[!is.na(plant.record.df$genus),], 
               neoPPE[,c("genus","REGION","PPE")], by.x=c("genus", "REGION"), by.y=c("genus", "REGION"),
                all.x=TRUE, all.y=FALSE, sort=FALSE)

neoPPEFam = array2DF(tapply(neoPPE$PPE, list(neoPPE$family, neoPPE$REGION), mean, na.rm=TRUE))
colnames(neoPPEFam) = c("family", "REGION", "PPE")

famPPE = merge(plant.record.df[is.na(plant.record.df$genus),], neoPPEFam[,c("family","REGION","PPE")], by.x=c("family", "REGION"), by.y=c("family", "REGION"),
                        all.x=TRUE, all.y=FALSE, sort=FALSE)
genPPE = genPPE[,match(colnames(famPPE), colnames(genPPE))]

plantPPE = rbind(famPPE, genPPE)

# how many samples do we have PPE for?
ppeStats = tapply(plantPPE$value, is.na(plantPPE$PPE), sum, na.rm=TRUE)
ppeStats / sum(ppeStats)

plantPPE$countPPE = plantPPE$value / plantPPE$PPE

plantPPE = droplevels(plantPPE[!is.na(plantPPE$countPPE),])

# sample sizes by count
continentGenus = as.data.frame(tapply(plant.record.df$value, 
                                      list(plant.record.df$genus,
                                           plant.record.df$REGION),
                        sum))

continentFam = with(plant.record.df[is.na(plant.record.df$genus),],
                      as.data.frame(tapply(value, 
                                      list(family, REGION),
                                      sum)))

# how many have continent specific sums?
contGenCount = sum(continentGenus$Europe[rownames(continentGenus) %in% eurPPEraw$genus[!is.na(eurPPEraw$ppeEurGen)]],
                     continentGenus$`North America`[rownames(continentGenus) %in% amPPEraw$genus[!is.na(amPPEraw$ppeamGen)]])
worldGenCount = sum(continentGenus$Europe[rownames(continentGenus) %in% eurPPEraw$genus[!is.na(eurPPEraw$ppeHemGen) & is.na(eurPPEraw$ppeEurGen)]],
                     continentGenus$`North America`[rownames(continentGenus) %in% amPPEraw$genus[!is.na(amPPEraw$ppeHemGen) & is.na(amPPEraw$ppeamGen)]])

contFamCount = sum(continentGenus$Europe[rownames(continentGenus) %in% eurPPEraw$genus[!is.na(eurPPEraw$ppeEur) & is.na(eurPPEraw$ppeEurGen) & is.na(eurPPEraw$ppeHemGen)]],
                   continentGenus$`North America`[rownames(continentGenus) %in% amPPEraw$genus[!is.na(amPPEraw$ppeAm) & is.na(amPPEraw$ppeHemGen) & is.na(amPPEraw$ppeamGen)]])

contFamCount = sum(continentGenus$Europe[rownames(continentGenus) %in% eurPPEraw$genus[!is.na(eurPPEraw$ppeEur) & is.na(eurPPEraw$ppeEurGen) & is.na(eurPPEraw$ppeHemGen)]],
                   continentGenus$`North America`[rownames(continentGenus) %in% amPPEraw$genus[!is.na(amPPEraw$ppeAm) & is.na(amPPEraw$ppeHemGen) & is.na(amPPEraw$ppeamGen)]])

#           Remove aquatic families ####

# remove aquatic families
plantPPE <- droplevels(plantPPE[!plantPPE$family %in% 
                                                  c("Potamogetonaceae", "Nymphaeaceae", "Typhaceae", "Cabombaceae", "Alismataceae",
                                                    "Haloragaceae"),])

plantPPE = droplevels(plantPPE[!is.na(plantPPE$value),])

# CHRONOLOGY MODELLING ####

# rBacon runs models for individual time series, needs external files and folder
# setup, and can fail. This section queries which datasets need modelling,
# sets up files and runs the chronology model.

sitesToModel <- unique(plant.record.df$datasetid)
sitesDone = list.dirs("./baconRuns")
sitesDone <- substr(sitesDone, regexpr("\\./baconRuns", sitesDone)+12, nchar(sitesDone)) 

sitesToModel <- sitesToModel[!sitesToModel %in% sitesDone]

# assume median error for ages without errors
plant.record.df$ageError = rowMeans(cbind(plant.record.df$ageolder - plant.record.df$age,
                                          plant.record.df$age - plant.record.df$ageyounger))
plant.record.df$ageError[is.na(plant.record.df$ageError)] = median(plant.record.df$ageError, na.rm=TRUE)

plant.record.df = plant.record.df[!is.na(plant.record.df$datasetid),]

# remove any sites in AD/BC years
plant.record.df = droplevels(plant.record.df[plant.record.df$agetype != "Calendar years AD/BC",])

# rBacon process, set up files and run models
loc.list <- lapply(sitesToModel, function(x){
                     
                     print(x)
                     xData = droplevels(plant.record.df[plant.record.df$datasetid == x,])
                     xData = xData[!duplicated(xData$depth),]
                     
                     xSub <- xData[c("sampleid", "age", "ageError", "depth")]

                     if(nrow(xSub) < 2){return(NULL)}
                     
                     # create directory for date files
                     dir.create(paste0("./baconRuns/", x))
                     fileName <- paste0(x, ".csv")
                    
                     xSub = xSub[, c("sampleid", "age", "ageError", "depth")]
                     colnames(xSub) = c("labID", "age", "error", "depth")
                     
                     xSub = xSub[complete.cases(xSub),]
                     
                     # calibrate date
                     needsCal = !grepl("calibrated", tolower(xData$agetype[1]))
                     
                     if(needsCal){
                     xCal = do.call("rbind", lapply(1:nrow(xSub), function(n){
                       print(n)
                       
                       # have to set bombalert to FALSE to allow character 
                       # calibration curves
                       x=try(calibrate(age=xSub$age[n],
                                   error = 10,#xSub$error[n],
                                   cc = ifelse(xSub$age[n] < 0, "nh1", 1),
                                   postbomb = ifelse(xSub$age[n] < 300, TRUE, FALSE),
                                   BCAD = FALSE, 
                                   draw = FALSE,
                                   bombalert = FALSE))
                       
                       if(class(x) == "try-error"){
                         return(data.frame(age = NA,
                                           error = NA))
                       }
                       
                       x=sample(x[[1]][,1], 1e6, replace=TRUE, prob=x[[1]][,2])
                       return(data.frame(age = mean(x),
                                         error = sd(x)))
                     }))
                     
                     xSub$age = xCal[,1]
                     xSub$error = xCal[,2]
                     }
                     
                     xSub <- xSub[!is.na(xSub$age),]
                     xSub <- xSub[order(xSub$depth),]
                     
                     write.csv(xSub,
                               paste0("./baconRuns/", x, "/", fileName),
                               row.names=FALSE)
                     
                     write.table(xSub$depth, paste0("./baconRuns/", x, "/", x, "_depths.txt"), 
                                 row.names=FALSE, col.names=FALSE, quote=FALSE)
                     
                     loc <- try(Bacon(core=x,
                                      coredir=paste0("./baconRuns"), 
                                      cc=0,
                                      depths.file=TRUE,
                                      ask=FALSE,
                                      thick=5,
                                      d.min = min(xSub$depth),
                                      d.max = max(xSub$depth),
                                      acc.mean = 15,
                                      suggest=FALSE,
                                      normal=TRUE,
                                      burnin=500,
                                      ssize=1000))
                     
                     return(loc)
                   })

# rBacon can also look finished but fail invisibly, but we can tell because
# the folders have < 5 files in them.
baconCleanup <- list.files("./baconRuns", full.names = TRUE)
baconLength <- sapply(baconCleanup, function(x){length(list.files(x))})
table(baconLength)
unlink(baconCleanup[baconLength < 5], recursive=TRUE)

# import chronology results and match them to neotoma records

baconAges = list.files("./baconRuns", pattern="_ages.txt", recursive=TRUE, full.names=TRUE)

baconAges = do.call("rbind", lapply(baconAges, function(x){
  tempAge = read.table(x, header=TRUE)
  
  dataName = gsub("\\.\\/baconRuns\\/|_ages.txt", "", x)
  tempAge$dataset = substr(dataName, 1, regexpr("\\/", dataName)-1)
  print(tempAge$dataset[1])
  return(tempAge)
}))
colnames(baconAges) = c("depth", "ageMin", "ageMax", "ageMed", "ageMean", "dataset")

plantPPE = merge(plantPPE, baconAges,
                        by.x=c("datasetid", "depth"), by.y=c("dataset", "depth"),
                        all.x=TRUE, all.y=FALSE, sort=FALSE)

table(is.na(plantPPE$ageMean[!duplicated(plantPPE$datasetid)]))
cor(plantPPE$age, plantPPE$ageMean, use="complete.obs")

# for the sites where calibration and age chronology modelling failed, use
# Neotoma ages
plantPPE$ageMean[is.na(plantPPE$ageMean)] = plantPPE$age[is.na(plantPPE$ageMean)]
plantPPE$neoAge = plantPPE$age
colnames(plantPPE)[colnames(plantPPE) == "ageMean"] = "age"

#           final filtering ####

# remove records with no chronology
plantPPE <- droplevels(plantPPE[!is.na(plantPPE$age), ])

plantPPE <- plantPPE[order(plantPPE$datasetid),]

# remove any NA rows that have creeped in
plantPPE <- plantPPE[complete.cases(plantPPE[,c("datasetid", "sampleid", "variablename", "value")]),]

# WRITE FILES ####

saveRDS(plantPPE, paste0("./rawdata/processed_family_records.rds"))
write.csv(site.df, "./rawdata/siteDf.csv")
