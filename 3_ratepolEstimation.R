# ################################# ####
# Author:    Timothy L Staples      ####
# Collaborators: John Pandolfi      #### 
#                Wolfgang Kiessling ####
# ################################# ####
# SET-UP ####
# Global attributes & working directories ####

rm(list=ls())
setwd('/Users/uqtstapl/Dropbox/Tim/Post-doc/Research projects/novel_comms_legacy/prodCode/')

# Packages & functions ####

# source functions from 'functions' sub-folder
sapply(paste0("./functions/", list.files("./functions", pattern =".R")), source)

package.loader(c("vegan", "hilldiv"))

# load packages (RRatepol is only available via github)
dev = try(library(devtools))
if(class(dev) == "try-error"){
  install.packages("devtools")
  library(devtools)
}

rate = try(library(RRatepol))
if(class(rate) == "try-error"){
devtools::install_github("HOPE-UIB-BIO/R-Ratepol-package")
library(RRatepol)
}

# import data ####
mottlRandN = 999

plantPPE <- readRDS("./rawdata/processed_family_records.rds")
plantHolo = droplevels(plantPPE[plantPPE$age <= 12000 & plantPPE$age > -100, ])

# RUN RATEPOL ESTIMATION (Late glacial - genus) ####

# This approach co-opts Mottl et al. 2021 process from the RRatepol package,
# with additional dissimilarities and statistics calculated on resampled
# compositions and returned as per the ROC values.

trajData = do.call('rbind', lapply(split(plantPPE, f=plantPPE$datasetid), function(x){
  
  print(x$datasetid[1])
  ageData = sort(unique(x$age))
  
  commData = tapply(x$countPPE, list(x$age, x$genus), sum, na.rm=TRUE)
  commData[is.na(commData)] = 0
  commData = cbind(sample_id = as.character(1:nrow(commData)),
                   as.data.frame(commData))
  
  ageData = data.frame(sample_id = as.character(1:nrow(commData)),
                       age = as.numeric(ageData))
  
  if(ncol(commData) < 3){return(NULL)}
  
  ageData = ageData[rowSums(commData[,-1]) > 0,]
  commData = commData[rowSums(commData[,-1]) > 0,]
  rownames(commData) = NULL
  
  # check there are at least 6 500 year bins with observations
  potBins = length(levels(droplevels(cut(ageData$age, breaks=seq(-50,50000,500)))))
  
  if(potBins <= 5 | ncol(commData) < 3){return(NULL)}
  
  suppressMessages({
  commEst = estimate_rocMOD(
    data_source_community = commData,
    data_source_age = ageData,
    age_uncertainty = NULL,
    smooth_method = "shep",
    working_units = "MW",
    number_of_shifts = 1,
    rand = mottlRandN,
    use_parallel = FALSE,
    dissimilarity_coefficient = "bray",
    verbose = FALSE,
    lofK = 3)
  })
  
  commResults = cbind(x[1,c("REGION", "datasetid", "siteid", "depositionalenvironment", "chronologyname", "lat", "long")],
                      commEst[[1]])
  
  #saveRDS(commResults, paste0("./outputs/trajData/", x$datasetid[1], ".rds"))
  return(commResults)
  
  }))
saveRDS(trajData, "./outputs/trajectoryDataGenus.rds")

# RUN RATEPOL ESTIMATION (Holocene - genus) ####
recentTrajData = do.call('rbind', lapply(split(plantHolo, f=plantHolo$datasetid), function(x){
  
  print(x$datasetid[1])
  ageData = sort(unique(x$age))
  
  commData = tapply(x$countPPE, list(x$age, x$genus), sum, na.rm=TRUE)
  commData[is.na(commData)] = 0
  commData = cbind(sample_id = as.character(1:nrow(commData)),
                   as.data.frame(commData))
  
  ageData = data.frame(sample_id = as.character(1:nrow(commData)),
                       age = as.numeric(ageData))
  
  if(ncol(commData) < 3){return(NULL)}
  
  ageData = ageData[rowSums(commData[,-1]) > 0,]
  commData = commData[rowSums(commData[,-1]) > 0,]
  rownames(commData) = NULL
  
  # check there are at least 6 500 year bins with observations
  potBins = length(levels(droplevels(cut(ageData$age, breaks=seq(-50,50000,250)))))
  
  if(potBins <= 5 | ncol(commData) < 3){return(NULL)}
  
  suppressMessages({
    commEst = estimate_rocMOD(
      data_source_community = commData,
      data_source_age = ageData,
      age_uncertainty = NULL,
      smooth_method = "shep",
      working_units = "MW",
      number_of_shifts = 1,
      bin_size = 250,
      rand = mottlRandN,
      use_parallel = FALSE,
      dissimilarity_coefficient = "bray",
      verbose = FALSE,
      lofK = 3)
  })
  
  commResults = cbind(x[1,c("REGION", "datasetid", "siteid", "depositionalenvironment", "chronologyname", "lat", "long")],
                      commEst[[1]])
  
  #saveRDS(commResults, paste0("./outputs/trajData/", x$datasetid[1], ".rds"))
  return(commResults)
  
}))
saveRDS(recentTrajData, "./outputs/recentTrajectoryDataGenus.rds")

# RUN RATEPOL ESTIMATION (Late Glacial - family) ####
trajData = do.call('rbind', lapply(split(plantPPE, f=plantPPE$datasetid), function(x){
  
  print(x$datasetid[1])
  ageData = sort(unique(x$age))
  
  commData = tapply(x$countPPE, list(x$age, x$family), sum, na.rm=TRUE)
  commData[is.na(commData)] = 0
  commData = cbind(sample_id = as.character(1:nrow(commData)),
                   as.data.frame(commData))
  
  ageData = data.frame(sample_id = as.character(1:nrow(commData)),
                       age = as.numeric(ageData))
  
  if(ncol(commData) < 3){return(NULL)}
  
  ageData = ageData[rowSums(commData[,-1]) > 0,]
  commData = commData[rowSums(commData[,-1]) > 0,]
  rownames(commData) = NULL
  
  # check there are at least 6 500 year bins with observations
  potBins = length(levels(droplevels(cut(ageData$age, breaks=seq(-50,50000,500)))))
  
  if(potBins <= 5 | ncol(commData) < 3){return(NULL)}
  
  suppressMessages({
    commEst = estimate_rocMOD(
      data_source_community = commData,
      data_source_age = ageData,
      age_uncertainty = NULL,
      smooth_method = "shep",
      working_units = "MW",
      number_of_shifts = 1,
      rand = mottlRandN,
      use_parallel = FALSE,
      dissimilarity_coefficient = "bray",
      verbose = FALSE,
      lofK = 3)
  })
  
  commResults = cbind(x[1,c("REGION", "datasetid", "siteid", "depositionalenvironment", "chronologyname", "lat", "long")],
                      commEst[[1]])
  
  #saveRDS(commResults, paste0("./outputs/trajData/", x$datasetid[1], ".rds"))
  return(commResults)
  
}))
saveRDS(trajData, "./outputs/trajectoryDataFamily.rds")

# RUN RATEPOL ESTIMATION (Holocene - family) ####
recentTrajData = do.call('rbind', lapply(split(plantHolo, f=plantHolo$datasetid), function(x){
  
  
  print(x$datasetid[1])
  ageData = sort(unique(x$age))
  
  commData = tapply(x$countPPE, list(x$age, x$family), sum, na.rm=TRUE)
  commData[is.na(commData)] = 0
  commData = cbind(sample_id = as.character(1:nrow(commData)),
                   as.data.frame(commData))
  
  ageData = data.frame(sample_id = as.character(1:nrow(commData)),
                       age = as.numeric(ageData))
  
  if(ncol(commData) < 3){return(NULL)}
  
  ageData = ageData[rowSums(commData[,-1]) > 0,]
  commData = commData[rowSums(commData[,-1]) > 0,]
  rownames(commData) = NULL
  
  # check there are at least 6 500 year bins with observations
  potBins = length(levels(droplevels(cut(ageData$age, breaks=seq(-50,50000,250)))))
  
  if(potBins <= 5 | ncol(commData) < 3){return(NULL)}
  
  suppressMessages({
    commEst = estimate_rocMOD(
      data_source_community = commData,
      data_source_age = ageData,
      age_uncertainty = NULL,
      smooth_method = "shep",
      working_units = "MW",
      number_of_shifts = 1,
      bin_size = 250,
      rand = mottlRandN,
      use_parallel = FALSE,
      dissimilarity_coefficient = "bray",
      verbose = FALSE,
      lofK = 3)
  })
  
  commResults = cbind(x[1,c("REGION", "datasetid", "siteid", "depositionalenvironment", "chronologyname", "lat", "long")],
                      commEst[[1]])
  
  #saveRDS(commResults, paste0("./outputs/trajData/", x$datasetid[1], ".rds"))
  return(commResults)
  
}))
saveRDS(recentTrajData, "./outputs/recentTrajectoryDataFamily.rds")