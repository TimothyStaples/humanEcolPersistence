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

package.loader(c("vegan", "DHARMa", "lme4", "multcomp", "plotrix", "shape", "hilldiv",
                 "merTools", "performance", "abind", "sf", "rworldmap", "readxl",
                 "terra", "ncdf4", "tidyr", "raster", "performance", "RRatepol", "rbacon",
                 "piecewiseSEM", "mgcv", "psych", "GPArotation", "FactoMineR", "viridisLite",
                 "glmmTMB", "pBrackets", "gamm4", "grImport"))

# A little custom function to add dates to output files
date.wrap <- function(string, ext){
  paste0(string, " ", Sys.Date(), ext)
}

# converts x and y values into angles on a continuous 360 degree scale
cont.angles <- function(x,y, offset = 0){
  angles <- -atan2(x, y) * (180/pi)
  angles <- ifelse(angles > 0, angles, 360 - abs(angles))
  rot.angles <- angles + 90 + offset
  rot.angles <- ifelse(rot.angles >= 360, rot.angles - 360, rot.angles)
  return(rot.angles)
}

# function to lighten or darken colours using ramps
colorShader <- function(col, scale, direction = "lighten"){
  rgb(colorRamp(c(col, ifelse(direction=="lighten", "white", "black")))(scale)/255)
}

# convert data to scale (default 0-1). For colour ramps for plotting.
unitScale <- function(x, custMin = NULL, custMax=NULL){
  if(is.null(custMax)){custMax = max(x, na.rm=TRUE)}
  if(is.null(custMin)){custMin = min(x, na.rm=TRUE)}
  (x - custMin) / (custMax - custMin)
}

# convert a set of factor levels from the output of cut() into a numeric variable
# with the cut centers
cutToNumeric <- function(x){
  rowMeans(cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", x)),
                 upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", x))))
}

# get objects in environment and rank them by size
envSized = function(){
  obj = ls(envir = globalenv())

  objSize = sapply(obj, function(x){
    object.size(get(x))
    })
  
  obj = obj[order(objSize, decreasing=TRUE)]
  objSize = sort(objSize, decreasing=TRUE)
  return(data.frame(objectSize = objSize / 1000^2))
}

nthAm = c(30, 90, -167, -50)
Eu = c(30, 90, -11, 50)

############################### ####
# LATE GLACIAL ANALYSES - GENUS ####
############################### ####
#           Variable processing & set-up ####

trajData <- readRDS("./outputs/trajectoryDataGenus.rds")

# remove entries that don't have legacy states in 500 year blocks
trajData = droplevels(trajData[trajData$n1gap_1 == -500 &
                                 trajData$n1gap_2 == -1000 &
                                 trajData$n1gap_3 == -1500 &
                                 complete.cases(trajData[,which(grepl("n1gap", colnames(trajData)))[1:3]]),])

testData = trajData

testData$REGION = as.factor(testData$REGION)

# pull out first three post-dissimilarities
testTPMat = testData[, grepl("dn1", colnames(testData))]
testTPMat = testTPMat[, 1:3]

# PCA
tpResPCA = princomp(testTPMat)

# eigenvalues
tpVals = tpResPCA$sdev^2 / sum(tpResPCA$sdev^2)

# reverse first axis if loadings are positive so higher values = more persistent
if(sum(tpResPCA$loadings[,1] > 0) == ncol(testTPMat)){
  tpResPCA$scores[,1] = tpResPCA$scores[,1] * -1
  tpResPCA$loadings[,1] = tpResPCA$loadings[,1] * -1
}

# reverse second axis if late time dissimilarities are positive so
# higher values = more late persistence, reflecting of an "attractor state"
if(tpResPCA$loadings[1,2] < 0 & tpResPCA$loadings[3,2] > 0 ){
  tpResPCA$scores[,2] = tpResPCA$scores[,2] * -1
  tpResPCA$loadings[,2] = tpResPCA$loadings[,2] * -1
}

tpComb = data.frame(tpPC1 = tpResPCA$scores[,1],
                    tpPC2 = tpResPCA$scores[,2])
testData = cbind(testData, tpComb)

saveRDS(testData, "./outputs/modelDataGenus.rds")

#                       tp 1 model ####

# scale continuous predictors
testData$ageS <- as.vector(scale(testData$age))
testData$latS <- as.vector(scale(testData$lat))

# random-effecting raw data (simplified mixed-effect GAM)
timeOnlyMb = gamm4(tpPC1 ~ s(ageS, bs="ps", by=REGION, k=25) + REGION + latS, data=droplevels(testData),
                   random=~(1|datasetid))
summary(timeOnlyMb$gam)
compare_performance(timeOnlyMb$mer, timeOnlyMb$gam)

# make predictions
timeMean = data.frame(age = rep(seq(1500, max(testData$age), len=200), 2),
                      ageS = (rep(seq(1500, max(testData$age), len=200), 2) - mean(testData$age)) / sd(testData$age),
                      latS = 0,
                      REGION = rep(levels(testData$REGION), each=200))
timeMean = cbind(timeMean, as.data.frame(predict(timeOnlyMb$gam, 
                                                 newdata=timeMean, se.fit=TRUE)))
timeMean$lci = timeMean$fit - 1.96 * timeMean$se.fit  
timeMean$uci = timeMean$fit + 1.96 * timeMean$se.fit  

#                       tp 2 model ####

# random-effecting raw data (simplified mixed-effect GAM)
timeOnlyMb2 = gamm4(tpPC2 ~ s(ageS, bs="ps", by=REGION, k=20) + REGION + latS, data=droplevels(testData),
                   random=~(1|datasetid))
summary(timeOnlyMb2$gam)
compare_performance(timeOnlyMb2$mer, timeOnlyMb2$gam)

timeMean2 = data.frame(age = rep(seq(1500, max(testData$age), len=200), 2),
                       ageS = (rep(seq(1500, max(testData$age), len=200), 2) - mean(testData$age)) / sd(testData$age),
                       latS = 0,
                      REGION = rep(levels(testData$REGION), each=200))
timeMean2 = cbind(timeMean2, as.data.frame(predict(timeOnlyMb2$gam, 
                                                 newdata=timeMean, se.fit=TRUE)))
timeMean2$lci = timeMean2$fit - 1.96 * timeMean2$se.fit  
timeMean2$uci = timeMean2$fit + 1.96 * timeMean2$se.fit  

#           continental temperature models ####

# use TRACE cells with random effects to estimate continental temperature change over time
regTempRast<-stack("./rawdata/climate/trace.01-36.22000BP.clm2.TSA.22000BP_decavg_400BCE.nc",
                   varname="TSA")

extent(regTempRast) = extent(c(0,360,-90,90))
# convert longitude to -180:180
regTempRast = rotate(regTempRast)

# convert to celcius
regTempRast = regTempRast - 273.15

nthAm = c(30, 90, -167, -50)
Eu = c(30, 90, -11, 50)
europeTempRast = raster::crop(x=regTempRast, y=extent(Eu[3], Eu[4], Eu[1], Eu[2]))
europeTemp = values(europeTempRast)
europeTemp = europeTemp[rowSums(is.na(europeTemp)) < ncol(europeTemp),]
europeTemp = t(apply(europeTemp, 1, function(x){x - x[1]}))
europeTime = as.numeric(gsub("X\\.|X", "", colnames(europeTemp)))

europeTempDF = data.frame(temp = as.vector(europeTemp),
                          cellID = paste0("Europe", rep(1:nrow(europeTemp), ncol(europeTemp))),
                          age = rep(europeTime, each=nrow(europeTemp)),
                          REGION = "Europe")

nthAmTempRast = raster::crop(x=regTempRast, y=extent(nthAm[3], nthAm[4], nthAm[1], nthAm[2]))
nthAmTemp = values(nthAmTempRast)
nthAmTemp = nthAmTemp[rowSums(is.na(nthAmTemp)) < ncol(nthAmTemp),]
nthAmTemp = t(apply(nthAmTemp, 1, function(x){x - x[1]}))
nthAmTime = as.numeric(gsub("X\\.|X", "", colnames(nthAmTemp)))

nthAmTempDF = data.frame(temp = as.vector(nthAmTemp),
                         cellID = paste0("nthAm", rep(1:nrow(nthAmTemp), ncol(nthAmTemp))),
                         age = rep(nthAmTime, each=nrow(nthAmTemp)),
                         REGION = "North America")

tempDF = rbind(europeTempDF, nthAmTempDF)
tempDF$REGION = as.factor(tempDF$REGION)

tempM <- gam(temp ~ s(age, by=REGION) + REGION, data=tempDF)

#  FIGURE 2 - GENUS ####

tempPred = data.frame(age = rep(seq(min(tempDF$age),
                                            max(tempDF$age),
                                            len=1000), 2),
                              REGION = rep(levels(tempDF$REGION), each = 1000))
tempPred = cbind(tempPred, as.data.frame(predict(tempM, newdata=tempPred, se.fit=TRUE)))
tempPred = split(tempPred, f=tempPred$REGION)

plot(tempPred[[1]]$fit ~ tempPred[[1]]$age, type="l", ylim=c(-1,14))
lines(tempPred[[2]]$fit ~ tempPred[[2]]$age, col="red")

pdf("./plots/trendsOverTimeRawGenus.pdf", height=5.5, width=7)

xlims=c(28500,-1750)
ylim=rbind(c(-0.35,0.25),
           c(-0.1,0.1))
braceHeight = 1
contCol = c("blue","darkgreen")

tempCol = c(rgb(colorRamp(c("blue", "red"))(1)/255),
            rgb(colorRamp(c('darkgreen', "red"))(1)/255))
contAlpha = col2rgb(contCol)/255
contAlpha = apply(contAlpha, 2, function(x){rgb(x[1],x[2],x[3],0.25)})

split.screen(rbind(t(replicate(2, c(0.12,0.9,0.565,0.99))), # TP1 Europe
                   t(replicate(2, c(0.12,0.9,0.1,0.525))),
                   c(0.12,0.9,0.525,0.565))) # TP2 Europe
                   
sapply(1:2, function(n){

tempAxis = n
yAxis = ylim[tempAxis,]
  
print(n)

# each odd screen sets up the global temp underlay
screen(n + 1 * (n-1))
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

plot(NULL, xlim=xlims, ylim=c(-1,13), xlab="", ylab="", axes=FALSE, xaxs="i")
 
rect(xright=10000,
      xleft=0,
      ybottom=par("usr")[3], ytop=relative.axis.point(braceHeight, "y"),
      col="grey90", border="black", lty="31")
if(n==1){
  text(x=mean(c(0,10000)), y=relative.axis.point(0.065, "y"), lheight=0.75,
       labels="Ellis et al. 2021\nanthrome coverage", col="black", font=1, cex=0.8)
}

lines(y=tempPred[[1]]$fit, x=(tempPred[[1]]$age * 1000), col=tempCol[1], lwd=1.5, lty="31")
text(y=tempPred[[1]]$fit[1], x=tempPred[[1]]$age[1]*1000, labels="Europe\ntemp",
     col=tempCol[1], pos=4)
lines(y=tempPred[[2]]$fit, x=(tempPred[[2]]$age * 1000), col=tempCol[2], lwd=1.5, lty="31")
text(y=tempPred[[2]]$fit[1], x=tempPred[[2]]$age[1]*1000, labels="Nth Am\ntemp",
     col=tempCol[2], pos=4)

box()

# only put axis if it's a right-hand plot
axis(side=4, col="red", col.ticks="red", col.axis="red")
if(n==1){
    mtext(side=4, line=1.5, at=par("usr")[3], text="Global temperature anomaly", col="red", las=0)
    mtext(side=4, line=2.5, at=par("usr")[3], text=expression("("*Delta*degree*"C from 1961-1990 mean)"), col="red", las=0)
  }

close.screen(n + 1 * (n-1))

screen(n + 1 * (n-1) + 1)

par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(NULL, xlim=xlims, ylim=yAxis, xlab="", ylab="", axes=FALSE, xaxs="i")

yHeight = diff(par("usr")[3:4])

abline(h=0, lty="31", col="grey")

if(n == 2){
  axis(side=1, at=seq(5000,25000,5000), labels=format(seq(5000,25000,5000), big.mark=","), mgp=c(3,0.2,0))
  axis(side=1, at=0, labels=0, mgp=c(3,0.2,0))
  axis(side=3, at=seq(0,25000,5000), labels=NA, mgp=c(3,0.2,0), tcl=0.25)
  axis(side=3, at=seq(0,25000,1000), labels=NA, tcl=0.125)
  mtext(side=1, line=1.5, text="Years before 1950AD")
  } else {axis(side=1, at=seq(0,25000,5000), labels=NA)}
axis(side=1, at=seq(0,25000,1000), labels=NA, tcl=-0.125)

axis(side=2)

if(n==1){
mtext(side=2, line=3, at=par("usr")[3], text="Assemblage persistence axes", las=0, font=2)
}

if(n == 1){
  mtext(side=2, line=2, text="PC1: Persistence", las=0)
} else {
  mtext(side=2, line=2, text="PC2: State attractor strength", las=0)  
}

# plot observed
if(n == 1){obs = timeMean} else {obs = timeMean2}
with(obs[obs$REGION == "Europe",], {
  lines(fit ~ age, lwd=2, col=contCol[1])
  polygon(x=c(age, rev(age)),
          y=c(lci, rev(uci)), border=NA, col=contAlpha[1])
  text(y=fit[1], x=age[1], labels="Europe", col=contCol[1], pos=4)  
})

with(obs[obs$REGION == "North America",], {
  lines(fit ~ age, lwd=2, col=contCol[2])
  polygon(x=c(age, rev(age)),
          y=c(lci, rev(uci)), border=NA, col=contAlpha[2])
  text(y=fit[1], x=age[1], labels="North\nAmerica", col=contCol[2], pos=4)  
})

text(x=relative.axis.point(0.01, "x"),
     y=relative.axis.point(0.935, "y"),
     labels=paste0("(", LETTERS[c(1:2)][n], ")"), font=2, adj=0)

arrowLabs = list(c("More\ntransient", "More\npersistent"),
                 c("Repulsor\nstate", "Attractor\nstate"))[[n]]
Arrows(y0=0, y1= 0 + (c(-0.15, 0.15) * c(1,0.35)[n]),
       x0=relative.axis.point(0.06, "x"), x1=relative.axis.point(0.06, "x"),
       arr.type="triangle", col="grey80", arr.width=0.1, arr.length=0.1, lwd=2)
text(y=0 + (c(-0.15, 0.15) * c(1,0.35)[n]), x=relative.axis.point(0.06, "x"),
     pos=c(1,3), adj=0,
     labels=arrowLabs,
     font=3, col="grey80", cex=0.8)

close.screen(n + 1 * (n-1) + 1)
})

screen(5)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(NULL, xlim=xlims, ylim=c(0,1), axes=FALSE, xlab="", ylab="", xaxs="i")

stad <- read.csv("./rawdata/stadial.csv")
stad$age <- as.numeric(gsub(",", "", stad$age))
stadAge = c(11700, stad$age, 25000)
statCat = as.factor(c(stad$cat, "stadial"))

rect(xleft=stadAge[-1], xright=stadAge[-length(stadAge)],
     ybottom=par('usr')[3], ytop=par("usr")[4], 
     col=c("white", "grey")[statCat], border=NA)
rect(xleft=11700, xright=0, ybottom=par('usr')[3], ytop=par("usr")[4],
     col="white")
text(x=c(mean(rev(stadAge)[c(1,2)]),
         mean(rev(stadAge)[c(5,8)]),
         mean(rev(stadAge)[c(8,14)]),
         mean(c(0, 11700))),
     y=0.5, labels=c("Stad", "Stad", "Inter", "Holocene"), cex=0.75, col=rgb(0.5,0.5,0.5,0.75))

box()
close.screen(all.screens=TRUE)
dev.off()

################################ ####
# LATE GLACIAL ANALYSES - FAMILY ####
################################ ####
#           Variable processing & set-up ####

trajData <- readRDS("./outputs/trajectoryDataFamily.rds")

# remove entries that don't have legacy states in 500 year blocks
trajData = droplevels(trajData[trajData$n1gap_1 == -500 &
                                 trajData$n1gap_2 == -1000 &
                                 trajData$n1gap_3 == -1500 &
                                 complete.cases(trajData[,which(grepl("n1gap", colnames(trajData)))[1:3]]),])

testData = trajData

testData$REGION = as.factor(testData$REGION)

# pull out first three post-dissimilarities
testTPMat = testData[, grepl("dn1", colnames(testData))]
testTPMat = testTPMat[, 1:3]

# PCA
tpResPCA = princomp(testTPMat)

# eigenvalues
tpVals = tpResPCA$sdev^2 / sum(tpResPCA$sdev^2)

# reverse first axis if loadings are positive so higher values = more persistent
if(sum(tpResPCA$loadings[,1] > 0) == ncol(testTPMat)){
  tpResPCA$scores[,1] = tpResPCA$scores[,1] * -1
  tpResPCA$loadings[,1] = tpResPCA$loadings[,1] * -1
}

# reverse second axis if late time dissimilarities are positive so
# higher values = more late persistence, reflecting of an "attractor state"
if(tpResPCA$loadings[1,2] < 0 & tpResPCA$loadings[3,2] > 0 ){
  tpResPCA$scores[,2] = tpResPCA$scores[,2] * -1
  tpResPCA$loadings[,2] = tpResPCA$loadings[,2] * -1
}

tpComb = data.frame(tpPC1 = tpResPCA$scores[,1],
                    tpPC2 = tpResPCA$scores[,2])
testData = cbind(testData, tpComb)

saveRDS(testData, "./outputs/modelDataFamily.rds")

#                       tp 1 model ####

# scale continuous predictors
testData$ageS <- as.vector(scale(testData$age))
testData$latS <- as.vector(scale(testData$lat))

# random-effecting raw data (simplified mixed-effect GAM)
timeOnlyMb = gamm4(tpPC1 ~ s(ageS, bs="ps", by=REGION, k=25) + REGION + latS, data=droplevels(testData),
                   random=~(1|datasetid))
summary(timeOnlyMb$gam)
compare_performance(timeOnlyMb$mer, timeOnlyMb$gam)

# make predictions
timeMean = data.frame(age = rep(seq(1500, max(testData$age), len=200), 2),
                      ageS = (rep(seq(1500, max(testData$age), len=200), 2) - mean(testData$age)) / sd(testData$age),
                      latS = 0,
                      REGION = rep(levels(testData$REGION), each=200))
timeMean = cbind(timeMean, as.data.frame(predict(timeOnlyMb$gam, 
                                                 newdata=timeMean, se.fit=TRUE)))
timeMean$lci = timeMean$fit - 1.96 * timeMean$se.fit  
timeMean$uci = timeMean$fit + 1.96 * timeMean$se.fit  

#                       tp 2 model ####

# random-effecting raw data (simplified mixed-effect GAM)
timeOnlyMb2 = gamm4(tpPC2 ~ s(ageS, bs="ps", by=REGION, k=20) + REGION + latS, data=droplevels(testData),
                    random=~(1|datasetid))
summary(timeOnlyMb2$gam)
compare_performance(timeOnlyMb2$mer, timeOnlyMb2$gam)

timeMean2 = data.frame(age = rep(seq(1500, max(testData$age), len=200), 2),
                       ageS = (rep(seq(1500, max(testData$age), len=200), 2) - mean(testData$age)) / sd(testData$age),
                       latS = 0,
                       REGION = rep(levels(testData$REGION), each=200))
timeMean2 = cbind(timeMean2, as.data.frame(predict(timeOnlyMb2$gam, 
                                                   newdata=timeMean, se.fit=TRUE)))
timeMean2$lci = timeMean2$fit - 1.96 * timeMean2$se.fit  
timeMean2$uci = timeMean2$fit + 1.96 * timeMean2$se.fit  

#           continental temperature models ####

# use TRACE cells with random effects to estimate continental temperature change over time
regTempRast<-stack("./rawdata/climate/trace.01-36.22000BP.clm2.TSA.22000BP_decavg_400BCE.nc",
                   varname="TSA")

extent(regTempRast) = extent(c(0,360,-90,90))
# convert longitude to -180:180
regTempRast = rotate(regTempRast)

# convert to celcius
regTempRast = regTempRast - 273.15

europeTempRast = raster::crop(x=regTempRast, y=extent(Eu[3], Eu[4], Eu[1], Eu[2]))
europeTemp = values(europeTempRast)
europeTemp = europeTemp[rowSums(is.na(europeTemp)) < ncol(europeTemp),]
europeTemp = t(apply(europeTemp, 1, function(x){x - x[1]}))
europeTime = as.numeric(gsub("X\\.|X", "", colnames(europeTemp)))

europeTempDF = data.frame(temp = as.vector(europeTemp),
                          cellID = paste0("Europe", rep(1:nrow(europeTemp), ncol(europeTemp))),
                          age = rep(europeTime, each=nrow(europeTemp)),
                          REGION = "Europe")

nthAmTempRast = raster::crop(x=regTempRast, y=extent(nthAm[3], nthAm[4], nthAm[1], nthAm[2]))
nthAmTemp = values(nthAmTempRast)
nthAmTemp = nthAmTemp[rowSums(is.na(nthAmTemp)) < ncol(nthAmTemp),]
nthAmTemp = t(apply(nthAmTemp, 1, function(x){x - x[1]}))
nthAmTime = as.numeric(gsub("X\\.|X", "", colnames(nthAmTemp)))

nthAmTempDF = data.frame(temp = as.vector(nthAmTemp),
                         cellID = paste0("nthAm", rep(1:nrow(nthAmTemp), ncol(nthAmTemp))),
                         age = rep(nthAmTime, each=nrow(nthAmTemp)),
                         REGION = "North America")

tempDF = rbind(europeTempDF, nthAmTempDF)
tempDF$REGION = as.factor(tempDF$REGION)

tempM <- gam(temp ~ s(age, by=REGION) + REGION, data=tempDF)

#  FIGURE 2 - Family ####

tempPred = data.frame(age = rep(seq(min(tempDF$age),
                                    max(tempDF$age),
                                    len=1000), 2),
                      REGION = rep(levels(tempDF$REGION), each = 1000))
tempPred = cbind(tempPred, as.data.frame(predict(tempM, newdata=tempPred, se.fit=TRUE)))
tempPred = split(tempPred, f=tempPred$REGION)

plot(tempPred[[1]]$fit ~ tempPred[[1]]$age, type="l", ylim=c(-1,14))
lines(tempPred[[2]]$fit ~ tempPred[[2]]$age, col="red")

pdf("./plots/trendsOverTimeRawFamily.pdf", height=5.5, width=7)

xlims=c(28500,-1750)
ylim=rbind(c(-0.25,0.2),
           c(-0.05,0.05))
braceHeight = 1
contCol = c("blue","darkgreen")

tempCol = c(rgb(colorRamp(c("blue", "red"))(1)/255),
            rgb(colorRamp(c('darkgreen', "red"))(1)/255))
contAlpha = col2rgb(contCol)/255
contAlpha = apply(contAlpha, 2, function(x){rgb(x[1],x[2],x[3],0.25)})

split.screen(rbind(t(replicate(2, c(0.12,0.9,0.565,0.99))), # TP1 Europe
                   t(replicate(2, c(0.12,0.9,0.1,0.525))),
                   c(0.12,0.9,0.525,0.565))) # TP2 Europe

sapply(1:2, function(n){
  
  tempAxis = n
  yAxis = ylim[tempAxis,]
  
  print(n)
  
  # each odd screen sets up the global temp underlay
  screen(n + 1 * (n-1))
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  
  plot(NULL, xlim=xlims, ylim=c(-1,13), xlab="", ylab="", axes=FALSE, xaxs="i")
  
  rect(xright=10000,
       xleft=0,
       ybottom=par("usr")[3], ytop=relative.axis.point(braceHeight, "y"),
       col="grey90", border="black", lty="31")
  if(n==1){
    text(x=mean(c(0,10000)), y=relative.axis.point(0.065, "y"), lheight=0.75,
         labels="Ellis et al. 2021\nanthrome coverage", col="black", font=1, cex=0.8)
  }
  
  lines(y=tempPred[[1]]$fit, x=(tempPred[[1]]$age * 1000), col=tempCol[1], lwd=1.5, lty="31")
  text(y=tempPred[[1]]$fit[1], x=tempPred[[1]]$age[1]*1000, labels="Europe\ntemp",
       col=tempCol[1], pos=4)
  lines(y=tempPred[[2]]$fit, x=(tempPred[[2]]$age * 1000), col=tempCol[2], lwd=1.5, lty="31")
  text(y=tempPred[[2]]$fit[1], x=tempPred[[2]]$age[1]*1000, labels="Nth Am\ntemp",
       col=tempCol[2], pos=4)
  
  box()
  
  # only put axis if it's a right-hand plot
  axis(side=4, col="red", col.ticks="red", col.axis="red")
  if(n==1){
    mtext(side=4, line=1.5, at=par("usr")[3], text="Global temperature anomaly", col="red", las=0)
    mtext(side=4, line=2.5, at=par("usr")[3], text=expression("("*Delta*degree*"C from 1961-1990 mean)"), col="red", las=0)
  }
  
  close.screen(n + 1 * (n-1))
  
  screen(n + 1 * (n-1) + 1)
  
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  plot(NULL, xlim=xlims, ylim=yAxis, xlab="", ylab="", axes=FALSE, xaxs="i")
  
  yHeight = diff(par("usr")[3:4])
  
  abline(h=0, lty="31", col="grey")
  
  if(n == 2){
    axis(side=1, at=seq(5000,25000,5000), labels=format(seq(5000,25000,5000), big.mark=","), mgp=c(3,0.2,0))
    axis(side=1, at=0, labels=0, mgp=c(3,0.2,0))
    axis(side=3, at=seq(0,25000,5000), labels=NA, mgp=c(3,0.2,0), tcl=0.25)
    axis(side=3, at=seq(0,25000,1000), labels=NA, tcl=0.125)
    mtext(side=1, line=1.5, text="Years before 1950AD")
  } else {axis(side=1, at=seq(0,25000,5000), labels=NA)}
  axis(side=1, at=seq(0,25000,1000), labels=NA, tcl=-0.125)
  
  axis(side=2)
  
  if(n==1){
    mtext(side=2, line=3, at=par("usr")[3], text="Assemblage persistence axes", las=0, font=2)
  }
  
  if(n == 1){
    mtext(side=2, line=2, text="PC1: Persistence", las=0)
  } else {
    mtext(side=2, line=2, text="PC2: State attractor strength", las=0)  
  }
  
  # plot observed
  if(n == 1){obs = timeMean} else {obs = timeMean2}
  with(obs[obs$REGION == "Europe",], {
    lines(fit ~ age, lwd=2, col=contCol[1])
    polygon(x=c(age, rev(age)),
            y=c(lci, rev(uci)), border=NA, col=contAlpha[1])
    text(y=fit[1], x=age[1], labels="Europe", col=contCol[1], pos=4)  
  })
  
  with(obs[obs$REGION == "North America",], {
    lines(fit ~ age, lwd=2, col=contCol[2])
    polygon(x=c(age, rev(age)),
            y=c(lci, rev(uci)), border=NA, col=contAlpha[2])
    text(y=fit[1], x=age[1], labels="North\nAmerica", col=contCol[2], pos=4)  
  })
  
  text(x=relative.axis.point(0.01, "x"),
       y=relative.axis.point(0.935, "y"),
       labels=paste0("(", LETTERS[c(1:2)][n], ")"), font=2, adj=0)
  
  arrowLabs = list(c("More\ntransient", "More\npersistent"),
                   c("Repulsor\nstate", "Attractor\nstate"))[[n]]
  Arrows(y0=0, y1= 0 + (c(-0.15, 0.15) * c(1,0.35)[n]),
         x0=relative.axis.point(0.06, "x"), x1=relative.axis.point(0.06, "x"),
         arr.type="triangle", col="grey80", arr.width=0.1, arr.length=0.1, lwd=2)
  text(y=0 + (c(-0.15, 0.15) * c(1,0.35)[n]), x=relative.axis.point(0.06, "x"),
       pos=c(1,3), adj=0,
       labels=arrowLabs,
       font=3, col="grey80", cex=0.8)
  
  close.screen(n + 1 * (n-1) + 1)
})

screen(5)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(NULL, xlim=xlims, ylim=c(0,1), axes=FALSE, xlab="", ylab="", xaxs="i")

stad <- read.csv("./rawdata/stadial.csv")
stad$age <- as.numeric(gsub(",", "", stad$age))
stadAge = c(11700, stad$age, 25000)
statCat = as.factor(c(stad$cat, "stadial"))

rect(xleft=stadAge[-1], xright=stadAge[-length(stadAge)],
     ybottom=par('usr')[3], ytop=par("usr")[4], 
     col=c("white", "grey")[statCat], border=NA)
rect(xleft=11700, xright=0, ybottom=par('usr')[3], ytop=par("usr")[4],
     col="white")
text(x=c(mean(rev(stadAge)[c(1,2)]),
         mean(rev(stadAge)[c(5,8)]),
         mean(rev(stadAge)[c(8,14)]),
         mean(c(0, 11700))),
     y=0.5, labels=c("Stad", "Stad", "Inter", "Holocene"), cex=0.75, col=rgb(0.5,0.5,0.5,0.75))

box()
close.screen(all.screens=TRUE)
dev.off()

########################### ####
# HOLOCENE ANALYSES - GENUS ####
########################### ####

rData <- readRDS("./outputs/recentTrajectoryDataGenus.rds")

rDatasetID = rData[!duplicated(rData$datasetid),]

# run persistence PCA
testTPMat = rData[, grepl("dn1", colnames(rData))]
testTPMat = testTPMat[, 1:3]

rData = rData[complete.cases(testTPMat),]
testTPMat = testTPMat[complete.cases(testTPMat),]

#tpResPCA = princomp(sapply(testTPMat, function(x){logit(beta.tr(x))}))
tpResPCA = princomp(sapply(testTPMat, function(x){x}))

tpVals = tpResPCA$sdev^2 / sum(tpResPCA$sdev^2)

# reverse first axis if loadings are positive so higher values = more persistent
if(sum(tpResPCA$loadings[,1] > 0) == ncol(testTPMat)){
  tpResPCA$scores[,1] = tpResPCA$scores[,1] * -1
  tpResPCA$loadings[,1] = tpResPCA$loadings[,1] * -1
}

# reverse second axis if late time dissimilarities are positive so
# higher values = more late persistence, reflecting of an "attractor state"
if(tpResPCA$loadings[1,2] < 0 & tpResPCA$loadings[3,2] > 0 ){
  tpResPCA$scores[,2] = tpResPCA$scores[,2] * -1
  tpResPCA$loadings[,2] = tpResPCA$loadings[,2] * -1
}

print(tpResPCA$loadings)

tpComb = data.frame(tpPC1 = tpResPCA$scores[,1],
                    tpPC2 = tpResPCA$scores[,2])

rData = cbind(rData, tpComb)

#           Add modern population density ####

library(terra)
popDens = rast("/Users/uqtstapl/Downloads/gpw-v4-population-density-rev11_totpop_2pt5_min_nc/gpw_v4_population_density_rev11_2pt5_min.nc")

rData$popDens = extract(popDens, rData[,c("long","lat")])[,2]

#           Ellis Anthrome pop dens estimates ####
ellis <- read_sf("/Users/uqtstapl/Downloads/Anthromes-12k-DGG/an12_dgg_inputs.shp")
ellisData = read.csv("~/Downloads/dataverse_files/an12_dgg_baseline.csv")

ellisCat = data.frame(ellisID = c(11, 12, 21, 22, 23, 24, 31, 32, 33, 34,
                                  41, 42, 43, 51, 52, 53, 54, 61, 62, 63, 70),
                      ellisName = c("Urban", "Mixed settlements", "Rice villages", 
                                    "Irrigated villages", "Rainfed villages", "Pastoral villages",
                                    "Residential irrigated croplands", "Residential rainfed croplands",
                                    "Populated croplands", "Remote croplands", "Residental rangelands",
                                    "Populated rangelands", "Remote rangelands", "Residential woodlands",
                                    " Populated woodlands", "Remote woodlands", "Inhabited drylands",
                                    "Wild woodlands", "Wild drylands", "Ice, uninhabited", "No land"),
                      ellisClass = c("urban", rep("village", 5), rep("cropland", 4),
                                     rep("rangeland", 3), rep("woodland", 4),
                                     rep("woodland", 2), "ice", "no land"),
                      ellisPopDens = c(5000, 
                                       rep(2500, 5), 
                                       100, 100, 10, 1,
                                       100, 10, 1,
                                       100, 10, 1, 1,
                                       0, 0, 0, 0))

rDataCoords = rData[!duplicated(rData$datasetid),c("long", "lat")]
c = st_as_sf(rDataCoords, coords=c("long", "lat"))
c = st_set_crs(c, st_crs(ellis))
rDataEllisID = st_intersects(c, ellis$geometry)
rDataEllisID = as.data.frame(rDataEllisID)

rDataEllisID$datasetid = unique(rData$datasetid)[rDataEllisID$row.id]
colnames(rDataEllisID) = c("id", "ellisID", "datasetid")

rData = merge(rData, rDataEllisID[,c("ellisID", "datasetid")],
              by.x="datasetid", by.y="datasetid", all.x=TRUE, all.y=FALSE, sort=FALSE)

# now get the closest age.

ellisAges = gsub("X", "", colnames(ellisData)[-1])
ellisAges = cbind(substr(ellisAges, 1, nchar(ellisAges)-2),
                  substr(ellisAges, nchar(ellisAges)-1, nchar(ellisAges)))
ellisAges = ifelse(ellisAges[,2] == "BC", 
                   as.numeric(ellisAges[,1]) + 1950,
                   1950 - as.numeric(ellisAges[,1]))

#Adding 750 years so we get max density across the "persistence window"
rDataEllisCols = sapply(rData$age, function(x){
  which.min(abs((x - 750) - ellisAges))
})

rData$ellisCat = ellisData[cbind(rData$ellisID,rDataEllisCols)]

# now grab metadata from table with pop density etc, and look for trends
# in tp1 based on these cats.
rData = merge(rData, ellisCat, by.x = "ellisCat", by.y="ellisID",
              all.x=TRUE, all.y=FALSE, sort=FALSE)

rData$ellisDensRegion = as.factor(paste0(rData$REGION, ":", rData$ellisPopDens))
rData$ellisPopNum = rData$ellisPopDens
rData$ellisPopDens = as.factor(rData$ellisPopDens)

rData = rData[complete.cases(rData[,c("age", "lat", "popDens", "ellisPopDens")]),]
rData = droplevels(rData[rData$lat > 25, ])

rData$logDensS = as.vector(scale(log(rData$popDens+1)))
rData$ageS = as.vector(scale(rData$age))
rData$latS = as.vector(scale(rData$lat))

tsPopChange = tapply(as.numeric(as.character(rData$ellisPopDens)), rData$datasetid, range)
tsDS = names(tsPopChange)
tsPopChange = do.call("rbind", tsPopChange)
rownames(tsPopChange) = tsDS

tsPopSplit = split(as.data.frame(tsPopChange), f=rData$REGION[match(rownames(tsPopChange), rData$datasetid)])
table(tsPopSplit[[1]][,1], tsPopSplit[[1]][,2])
table(tsPopSplit[[2]][,1], tsPopSplit[[2]][,2])

sum(tsPopSplit[[2]][,1] < tsPopSplit[[2]][,2]) / nrow(tsPopSplit[[2]])

#           ellis exploration ####

contCol = c("blue","darkgreen")
contAlpha = col2rgb(contCol)/255
contAlpha = apply(contAlpha, 2, function(x){rgb(x[1],x[2],x[3],0.25)})

# consistent vs changing time series
a = table(cut(rData$age, breaks=seq(0,8000,1000)), rData$ellisPopDens, rData$REGION)
#a = tapply(rData$ellisPopDens, rData$datasetid, function(x){length(unique(x))})
barplot(t(prop.table(a[,,1], 1)))
barplot(t(prop.table(a[,,2], 1)))

# modern pop versus historical pop categories
boxplot(log(rData$popDens+1) ~ rData$ellisPopDen * rData$REGION,
        col=c(rep("blue", 5), rep("darkgreen", 5)), yaxt="n", xlab="", ylab="")
axis(side=2, at=log(c(0, 1,10,100,1000,10000)+1), labels=c(0,1,10,100,1000,10000), las=1)
axis(side=2, at=log(c(seq(1,10,1), seq(10,100,10), seq(100,1000,100), seq(1000,10000,1000))),
     labels=NA, tcl=-0.125)
mtext(side=2, line=2, text=expression("Population density at 2000AD (persons km"^-2*")"))

# test for modern pop ability to explain historical prop
rData$ellisPopNum <- as.numeric(as.character(rData$ellisPopDens))
histPopM <- lmer(log(ellisPopNum+1) ~ log(popDens+1) * REGION + (1|datasetid), data=rData)
performance(histPopM)
summary(histPopM)

saveRDS(rData, "./outputs/holoDataGenus.rds")

# ellis pop dens as proportions across Holocene
ellisEurope = st_crop(ellis, st_bbox(c(xmin=Eu[3], xmax=Eu[4], ymin=Eu[1], ymax=Eu[2]), crs=st_crs(4326)))
ellisNthAm = st_crop(ellis, st_bbox(c(xmin=nthAm[3], xmax=nthAm[4], ymin=nthAm[1], ymax=nthAm[2]), crs=st_crs(4326)))

# pull out and aggregate ellis data for each 1000 years
ellisDataEurope = ellisData[match(ellisEurope$id, ellisData$id),]
ellisDataNthAm = ellisData[match(ellisNthAm$id, ellisData$id),]

# now aggregate data for each time aggregate (thousand years?)
ellisColsCut = as.factor(round(ellisAges, -2))

ellisFullDens = do.call("rbind", lapply(unique(ellisColsCut), function(n){
  
  a = ellisDataEurope[,-1][,ellisColsCut == n & !is.na(ellisColsCut)]
  if(class(a) == "data.frame"){
    a = unlist(a)
  }
  
  aDens = factor(ellisCat$ellisPopDens[match(a, ellisCat$ellisID)], levels=unique(ellisCat$ellisPopDens))
  aDens = table(aDens) / length(aDens)
  
  b = ellisDataNthAm[,-1][,ellisColsCut == n & !is.na(ellisColsCut)]
  if(class(b) == "data.frame"){
    b = unlist(b)
  }
  bDens = factor(ellisCat$ellisPopDens[match(b, ellisCat$ellisID)], levels=unique(ellisCat$ellisPopDens))
  bDens = table(bDens) / length(bDens)
  
  return(data.frame(REGION = rep(c("Europe", "North America"), each = length(aDens)),
                    time = n,
                    popCat = rep(names(aDens), 2),
                    popDens = c(aDens, bDens)))
  
}))

#           Model ####

# remove missing pop dens & very low sampled high pop dens categories
rModel = droplevels(rData[!is.na(rData$ellisPopDens) & rData$ellisPopDens != "2500",])

tempA = gamm4(tpPC1 ~ s(ageS, k = 12, by=ellisDensRegion) + ellisPopDens + REGION * latS, data=rModel,
              random = ~(1|datasetid))

predDF = cbind(expand.grid(ageS = seq(min(rData$ageS, na.rm=TRUE), max(rData$ageS, na.rm=TRUE), len=1000),
                           ellisDensRegion = sort(unique(rModel$ellisDensRegion))),
               latS = 0,
               tempPred = 0)
regCats = do.call('rbind', strsplit(as.character(predDF$ellisDensRegion), ":"))
colnames(regCats) = c("REGION", "ellisPopDens")
predDF = cbind(predDF, regCats)

predDF = cbind(predDF, as.data.frame(predict(tempA$gam, newdata=predDF, se.fit=TRUE)))
predDF$age = predDF$ageS * sd(rData$age, na.rm=TRUE) + mean(rData$age, na.rm=TRUE)
predDF = predDF[predDF$age <= 12000,]

gradientCutoff = floor(0.005 * length(unique(rData$datasetid)))

europePopSamplingLims = sapply(split(rData[rData$REGION == "Europe",], 
                               f=rData$ellisPopDens[rData$REGION == "Europe"]),function(x){
                                 sort(tapply(x$age, x$datasetid, max), decreasing=TRUE)[gradientCutoff]
                               })
AmpopSamplingLims = sapply(split(rData[rData$REGION == "North America",], 
                               f=rData$ellisPopDens[rData$REGION == "North America"]),function(x){
                                 sort(tapply(x$age, x$datasetid, max), decreasing=TRUE)[gradientCutoff]
                               })

#           Plot ####

library(pBrackets)

pdf("./plots/holoceneLineGenus.pdf", height=6, width=11, useDingbats=FALSE)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

split.screen(rbind(c(0.12,0.535,0.3,0.95),
                   c(0.575,0.99,0.3,0.95),
                   
                   c(0.12,0.535,0.1,0.3), # population densities
                   c(0.575,0.99,0.1,0.3),
                   
                   c(0.88, 0.92, 0.325, 0.5))) # pop dens gradient

# try a linear plot instead of a surface plot
screen(1)
plot(NULL, xlim=c(12200,0), xlab="", ylab="", xaxt="n", xaxs="i", ylim=c(-0.3,0.25))  #ylim=c(-1.65,1)
abline(h=0, lty="31")
axis(side=1, at=seq(0,12000,2000), labels=NA)
mtext(side=2, line=4, text="Assemblage persistence PC1", las=0)
mtext(side=3, line=0, at=par('usr')[1], text="Europe", font=2,col="blue", adj=0)

par(xpd=NA)
arrowLabs = c("More\ntransient", "More\npersistent")
Arrows(y0=0, y1= 0 + c(-0.15, 0.15),
       x0=relative.axis.point(-0.12, "x"), x1=relative.axis.point(-0.12, "x"),
       arr.type="triangle", col="black", arr.width=0.1, arr.length=0.1, lwd=2)
text(y=0 + c(-0.15, 0.15), x=relative.axis.point(-0.12, "x"),
     pos=c(1,3), adj=0,
     labels=arrowLabs,
     font=3, col="black", cex=0.8)
par(xpd=FALSE)

europeCols = c(colorRampPalette(c("grey70", "blue"))(5)[-1], "darkblue")
x = split(predDF[grepl("Europe", predDF$ellisDensRegion),], f=predDF$ellisDensRegion[grepl("Europe", predDF$ellisDensRegion)])
x = x[sapply(x, nrow) > 0]
mapply(y=x, n=seq_along(x), function(y,n){

  y = y[y$age <= europePopSamplingLims[n],]
   
  colrgb = col2rgb(europeCols[n])/255

  polygon(x=c(y$age, rev(y$age)),
           y=c(y$fit + 1.96 * y$se.fit, rev(y$fit - 1.96 * y$se.fit)),
          border=NA, col=rgb(colrgb[1],colrgb[2],colrgb[3],0.25))
  lines(y$fit ~ y$age, col=europeCols[n], lwd=2, lty=ifelse(n==1, "31", "solid"))
  text(x=y$age[1], y=y$fit[1], labels=c("0","<1","<10","<100", ">100")[n],
       col=europeCols[n], pos=4, offset=0.2) 
 })

# add in 2500 cat as individual non-modelled point
points(x=750, y=mean(rData$tpPC1[rData$ellisDensRegion == "Europe:2500"]),
       pch=16, col="darkblue")
text(x=750, y=mean(rData$tpPC1[rData$ellisDensRegion == "Europe:2500"]), 
     labels=">100*", col="darkblue", pos=2, offset=0.2) 

#brackets(x2=max(x[[1]]$age), x1=8800, y1=relative.axis.point(0.75, "y"), y2 = relative.axis.point(0.75, "y"), h=0.15, type=1)
#segments(x0=8800, x1=8800, y0=relative.axis.point(0.75, "y"), y1=x[[2]]$fit[which.min(abs(8800-x[[2]]$age))], lty="dotted")
#par(lheight=0.8)
#text(x=mean(c(max(x[[1]]$age), 8800)), y=relative.axis.point(0.87, "y"), labels="(1) Increased\npersistence with\nhuman presence", font=3, adj=0.5)

brackets(x2=max(x[[1]]$age), x1=2900, y1=relative.axis.point(0.75, "y"), y2 = relative.axis.point(0.75, "y"), h=0.03, type=1)
segments(x0=max(x[[1]]$age), x1=max(x[[1]]$age), y0=relative.axis.point(0.75, "y"), y1=x[[2]]$fit[which.max(x[[2]]$age)], lty="dotted")
segments(x0=2900, x1=2900, y0=relative.axis.point(0.75, "y"), y1=x[[4]]$fit[which.min(abs(europePopSamplingLims[4]-x[[4]]$age))], lty="dotted")
text(x=mean(c(max(x[[1]]$age), 2900)), y=relative.axis.point(0.87, "y"), labels="(0) Human presence\nhas no detectable\npersistence signal", font=3, adj=0.5)

brackets(x2=2900, x1=min(x[[1]]$age), y1=relative.axis.point(0.75, "y"), y2 = relative.axis.point(0.75, "y"), h=0.03, type=1)
segments(x0=min(x[[1]]$age), x1=min(x[[1]]$age), y0=relative.axis.point(0.75, "y"), y1=x[[4]]$fit[which.min(x[[4]]$age)], lty="dotted")
text(x=mean(c(min(x[[1]]$age), 2900)), y=relative.axis.point(0.87, "y"), labels="(-) Decreased\npersistence with\nhuman presence", font=3, adj=0.5)
text(x=relative.axis.point(0.01, "x"), y=relative.axis.point(0.975, "y"),
     adj=0, labels="(A)", font=2)
close.screen(1)

screen(2)
plot(NULL, xlim=c(12200,0), xlab="", ylab="", xaxt="n", xaxs="i", ylim=c(-0.3,0.25)) #ylim=c(-1.25,1.75)
abline(h=0, lty="31")
axis(side=1, at=seq(0,12000,2000), labels=NA)
mtext(side=3, line=0, at=par('usr')[1], text="North America", font=2,col="darkgreen", adj=0)

# lines(tempPred[,2] ~ tempPred[,1], col="red")
amCols = c(colorRampPalette(c("grey70", "darkgreen"))(5)[-1], "#003000FF")
x = split(predDF[grepl("North America", predDF$ellisDensRegion),], f=predDF$ellisDensRegion[grepl("North America", predDF$ellisDensRegion)])
x = x[sapply(x, nrow) > 0]

mapply(y=x, n=seq_along(x), function(y,n){
  
  y = y[y$age <= AmpopSamplingLims[n],]
  
  colrgb = col2rgb(amCols[n])/255
  
  polygon(x=c(y$age, rev(y$age)),
          y=c(y$fit + 1.96 * y$se.fit, rev(y$fit - 1.96 * y$se.fit)),
          border=NA, col=rgb(colrgb[1],colrgb[2],colrgb[3],0.25), lwd=0.5)
  lines(y$fit ~ y$age, col=amCols[n], lwd=2, lty=ifelse(n==1, "31", "solid"))
  text(x=y$age[1], y=y$fit[1], labels=c("0","<1","<10","<100", ">100")[n],
       col=amCols[n], pos=4, offset=0.2) 
})

brackets(x2=max(x[[1]]$age), x1=9700, y1=relative.axis.point(0.75, "y"), y2 = relative.axis.point(0.75, "y"), h=0.03, type=1)
segments(x0=max(x[[1]]$age), x1=max(x[[1]]$age), y0=relative.axis.point(0.75, "y"), y1=x[[2]]$fit[which.max(x[[2]]$age)], lty="dotted")
segments(x0=9700, x1=9700, y0=relative.axis.point(0.75, "y"), y1=x[[2]]$fit[which.min(abs(9700-x[[2]]$age))], lty="dotted")
par(lheight=0.8)
text(x=mean(c(max(x[[1]]$age), 9700)), y=relative.axis.point(0.87, "y"), labels="(-) Decreased\npersistence with\nhuman presence", adj=0.5, font=3)

brackets(x2=9700, x1=6000, y1=relative.axis.point(0.75, "y"), y2 = relative.axis.point(0.75, "y"), h=0.03, type=1)
segments(x0=6000, x1=6000, y0=relative.axis.point(0.75, "y"), y1=x[[2]]$fit[which.min(abs(6000-x[[2]]$age))], lty="dotted")
par(lheight=0.8)
text(x=mean(c(9700, 6000)), y=relative.axis.point(0.87, "y"), labels="(+) Increased\npersistence with\nhuman presence", adj=0.5, font=3)

brackets(x2=6000, x1=1200, y1=relative.axis.point(0.75, "y"), y2 = relative.axis.point(0.75, "y"), h=0.03, type=1)
segments(x0=1200, x1=1200, y0=relative.axis.point(0.75, "y"), y1=x[[3]]$fit[which.min(abs(1200-x[[2]]$age))], lty="dotted")
text(x=mean(c(min(x[[1]]$age), 6000)), y=relative.axis.point(0.87, "y"), labels="(0) Human presence\nhas no detectable\npersistence signal", font=3, adj=0.5)

brackets(x2=1200, x1=min(x[[1]]$age), y1=relative.axis.point(0.75, "y"), y2 = relative.axis.point(0.75, "y"), h=0.03, type=1)
segments(x0=min(x[[1]]$age), x1=min(x[[1]]$age), y0=relative.axis.point(0.75, "y"), y1=x[[4]]$fit[which.min(x[[4]]$age)], lty="dotted")
text(x=mean(c(min(x[[1]]$age), 1200)), y=relative.axis.point(0.87, "y"), labels="(-)", font=3, adj=0.5)

text(x=relative.axis.point(0.01, "x"), y=relative.axis.point(0.975, "y"),
     adj=0, labels="(B)", font=2)
close.screen(2)

screen(3)
# unsampled windows of higher density human populations
plot(NULL, xlim=c(12200, 0), ylim=c(0,1), axes=FALSE, xlab="", ylab="", yaxs="i", xaxs="i")
axis(side=1, at=seq(2000,11500,2000), tcl=-0.25, 
     labels=format(seq(2000,11500,2000), big.mark=","), mgp=c(3,0.2,0))
axis(side=1, at=0, labels=0, mgp=c(3,0.2,0))
axis(side=1, at=seq(0,12000,1000), tcl=-0.125, labels=NA)

mtext(side=1, line=1.25, text="Years before present")
axis(side=2, at=seq(0,1,0.2))
mtext(side=2, line=4, text="Human population density", las=0, lheight=0.75)
mtext(side=2, line=2.5, text="(Proportion of\ncontinental area)", las=0, lheight=0.75, font=3)

# europe pop dens
ellisEurWide = ellisFullDens[ellisFullDens$REGION=="Europe",]
ellisEurWide = as.data.frame(pivot_wider(ellisEurWide[,-1], values_from="popDens", names_from="popCat"))
ellisEurWide[,"2500"] = rowSums(ellisEurWide[,c("2500","5000")])
rownames(ellisEurWide) = ellisEurWide[,1]
ellisEurWide = ellisEurWide[,!colnames(ellisEurWide) %in% c("time", "5000")]
ellisEurWide = ellisEurWide[,ncol(ellisEurWide):1]

europeEllisPlot = rbind(0, apply(ellisEurWide, 1, cumsum))

europeCols = c(colorRampPalette(c("grey70", "blue"))(5)[-1], "darkblue")

sapply(2:(nrow(europeEllisPlot)), function(n){
  print(n)
  polygon(x=c(as.numeric(colnames(europeEllisPlot)),
              rev(as.numeric(colnames(europeEllisPlot)))),
          y=c(europeEllisPlot[n-1,], rev(europeEllisPlot[n,])),
          col=europeCols[n-1])
})

sapply(1:ncol(ellisEurWide), function(n){
  labX = weighted.mean(as.numeric(rownames(ellisEurWide)), ellisEurWide[,n])
  labY = mean(europeEllisPlot[,which.min(abs(as.numeric(rownames(ellisEurWide)) - labX))][n:(n+1)])
  text(x=labX, y=labY, label=c("0", "<1", "<10", "<100")[n],
       col=c("black","black","white","white")[n], cex=0.75)
})
box()
text(x=relative.axis.point(0.01, "x"), y=relative.axis.point(0.925, "y"),
     adj=0, labels="(C)", font=2)

axis(side=3, at=seq(0,12000,2000), tcl=0.25, labels=NA)
axis(side=3, at=seq(0,12000,1000), tcl=0.125, labels=NA)

close.screen(3)

screen(4)
# unsampled windows of higher density human populations
plot(NULL, xlim=c(12200, 0), ylim=c(0,1), axes=FALSE, xlab="", ylab="", yaxs="i", xaxs="i")
axis(side=1, at=seq(2000,11500,2000), tcl=-0.25, 
     labels=format(seq(2000,11500,2000), big.mark=","), mgp=c(3,0.2,0))
axis(side=1, at=seq(0,12000,1000), tcl=-0.125, labels=NA)
axis(side=1, at=0, labels=0, mgp=c(3,0.2,0))
mtext(side=1, line=1.25, text="Years before present")
axis(side=2, at=seq(0,1,0.2))

# europe pop dens
ellisEurWide = ellisFullDens[ellisFullDens$REGION=="North America",]
ellisEurWide = as.data.frame(pivot_wider(ellisEurWide[,-1], values_from="popDens", names_from="popCat"))
ellisEurWide[,"2500"] = rowSums(ellisEurWide[,c("2500","5000")])
rownames(ellisEurWide) = ellisEurWide[,1]
ellisEurWide = ellisEurWide[,!colnames(ellisEurWide) %in% c("time", "5000")]
ellisEurWide = ellisEurWide[,ncol(ellisEurWide):1]

europeEllisPlot = rbind(0, apply(ellisEurWide, 1, cumsum))

europeCols = c(colorRampPalette(c("grey70", "darkgreen"))(5)[-1], "#003000FF")

sapply(2:(nrow(europeEllisPlot)), function(n){
  print(n)
  polygon(x=c(as.numeric(colnames(europeEllisPlot)),
              rev(as.numeric(colnames(europeEllisPlot)))),
          y=c(europeEllisPlot[n-1,], rev(europeEllisPlot[n,])),
          col=europeCols[n-1])
})

sapply(1:ncol(ellisEurWide), function(n){
  labX = weighted.mean(as.numeric(rownames(ellisEurWide)), ellisEurWide[,n])
  labY = mean(europeEllisPlot[,which.min(abs(as.numeric(rownames(ellisEurWide)) - labX))][n:(n+1)])
  text(x=labX, y=labY, label=c("0", "<1", "<10", "<100")[n],
       col=c("black","black","white","white")[n], cex=0.75)
})
box()
text(x=relative.axis.point(0.01, "x"), y=relative.axis.point(0.925, "y"),
     adj=0, labels="(D)", font=2)
axis(side=3, at=seq(0,12000,2000), tcl=0.25, labels=NA)
axis(side=3, at=seq(0,12000,1000), tcl=0.125, labels=NA)
close.screen(4)

screen(5)
image(y=seq(0,1,len=5),
      x=c(0,1),
      z=matrix(1:10, byrow=TRUE, nrow=2),
      col=c(c(colorRampPalette(c("grey70", "blue"))(5)[-1], "darkblue"),
            c(colorRampPalette(c("grey70", "darkgreen"))(5)[-1], "#003000FF")),
      useRaster=TRUE, axes=FALSE, xlab="", ylab="")
axis(side=2, mgp=c(3,0.5,0), at=seq(0,1,len=5),
     labels=c("0", "<1",  "<10", "<100", ">100"))
abline(v=0.5,)
box()
mtext(side=2, line=3.25, text="Pop density", las=0)
mtext(side=2, line=2.25, text=expression("(persons km"^-2*")"), las=0)

close.screen(5)

close.screen(all.screens=TRUE)
dev.off()

############################ ####
# HOLOCENE ANALYSES - FAMILY ####
############################ ####

rData <- readRDS("./outputs/recentTrajectoryDataFamily.rds")

rDatasetID = rData[!duplicated(rData$datasetid),]

# run persistence PCA
testTPMat = rData[, grepl("dn1", colnames(rData))]
testTPMat = testTPMat[, 1:3]

rData = rData[complete.cases(testTPMat),]
testTPMat = testTPMat[complete.cases(testTPMat),]

#tpResPCA = princomp(sapply(testTPMat, function(x){logit(beta.tr(x))}))
tpResPCA = princomp(sapply(testTPMat, function(x){x}))

tpVals = tpResPCA$sdev^2 / sum(tpResPCA$sdev^2)

# reverse first axis if loadings are positive so higher values = more persistent
if(sum(tpResPCA$loadings[,1] > 0) == ncol(testTPMat)){
  tpResPCA$scores[,1] = tpResPCA$scores[,1] * -1
  tpResPCA$loadings[,1] = tpResPCA$loadings[,1] * -1
}

# reverse second axis if late time dissimilarities are positive so
# higher values = more late persistence, reflecting of an "attractor state"
if(tpResPCA$loadings[1,2] < 0 & tpResPCA$loadings[3,2] > 0 ){
  tpResPCA$scores[,2] = tpResPCA$scores[,2] * -1
  tpResPCA$loadings[,2] = tpResPCA$loadings[,2] * -1
}

print(tpResPCA$loadings)

tpComb = data.frame(tpPC1 = tpResPCA$scores[,1],
                    tpPC2 = tpResPCA$scores[,2])

rData = cbind(rData, tpComb)

#           Add modern population density ####

library(terra)
popDens = rast("/Users/uqtstapl/Downloads/gpw-v4-population-density-rev11_totpop_2pt5_min_nc/gpw_v4_population_density_rev11_2pt5_min.nc")

rData$popDens = extract(popDens, rData[,c("long","lat")])[,2]

#           Ellis Anthrome pop dens estimates ####
ellis <- read_sf("/Users/uqtstapl/Downloads/Anthromes-12k-DGG/an12_dgg_inputs.shp")
ellisData = read.csv("~/Downloads/dataverse_files/an12_dgg_baseline.csv")

ellisCat = data.frame(ellisID = c(11, 12, 21, 22, 23, 24, 31, 32, 33, 34,
                                  41, 42, 43, 51, 52, 53, 54, 61, 62, 63, 70),
                      ellisName = c("Urban", "Mixed settlements", "Rice villages", 
                                    "Irrigated villages", "Rainfed villages", "Pastoral villages",
                                    "Residential irrigated croplands", "Residential rainfed croplands",
                                    "Populated croplands", "Remote croplands", "Residental rangelands",
                                    "Populated rangelands", "Remote rangelands", "Residential woodlands",
                                    " Populated woodlands", "Remote woodlands", "Inhabited drylands",
                                    "Wild woodlands", "Wild drylands", "Ice, uninhabited", "No land"),
                      ellisClass = c("urban", rep("village", 5), rep("cropland", 4),
                                     rep("rangeland", 3), rep("woodland", 4),
                                     rep("woodland", 2), "ice", "no land"),
                      ellisPopDens = c(5000, 
                                       rep(2500, 5), 
                                       100, 100, 10, 1,
                                       100, 10, 1,
                                       100, 10, 1, 1,
                                       0, 0, 0, 0))

rDataCoords = rData[!duplicated(rData$datasetid),c("long", "lat")]
c = st_as_sf(rDataCoords, coords=c("long", "lat"))
c = st_set_crs(c, st_crs(ellis))
rDataEllisID = st_intersects(c, ellis$geometry)
rDataEllisID = as.data.frame(rDataEllisID)

rDataEllisID$datasetid = unique(rData$datasetid)[rDataEllisID$row.id]
colnames(rDataEllisID) = c("id", "ellisID", "datasetid")

rData = merge(rData, rDataEllisID[,c("ellisID", "datasetid")],
              by.x="datasetid", by.y="datasetid", all.x=TRUE, all.y=FALSE, sort=FALSE)

# now get the closest age.

ellisAges = gsub("X", "", colnames(ellisData)[-1])
ellisAges = cbind(substr(ellisAges, 1, nchar(ellisAges)-2),
                  substr(ellisAges, nchar(ellisAges)-1, nchar(ellisAges)))
ellisAges = ifelse(ellisAges[,2] == "BC", 
                   as.numeric(ellisAges[,1]) + 1950,
                   1950 - as.numeric(ellisAges[,1]))

#Adding 750 years so we get max density across the "persistence window"
rDataEllisCols = sapply(rData$age, function(x){
  which.min(abs((x - 750) - ellisAges))
})

rData$ellisCat = ellisData[cbind(rData$ellisID,rDataEllisCols)]

# now grab metadata from table with pop density etc, and look for trends
# in tp1 based on these cats.
rData = merge(rData, ellisCat, by.x = "ellisCat", by.y="ellisID",
              all.x=TRUE, all.y=FALSE, sort=FALSE)

rData$ellisDensRegion = as.factor(paste0(rData$REGION, ":", rData$ellisPopDens))
rData$ellisPopNum = rData$ellisPopDens
rData$ellisPopDens = as.factor(rData$ellisPopDens)

rData = rData[complete.cases(rData[,c("age", "lat", "popDens", "ellisPopDens")]),]
rData = droplevels(rData[rData$lat > 25, ])

rData$logDensS = as.vector(scale(log(rData$popDens+1)))
rData$ageS = as.vector(scale(rData$age))
rData$latS = as.vector(scale(rData$lat))

tsPopChange = tapply(as.numeric(as.character(rData$ellisPopDens)), rData$datasetid, range)
tsDS = names(tsPopChange)
tsPopChange = do.call("rbind", tsPopChange)
rownames(tsPopChange) = tsDS

tsPopSplit = split(as.data.frame(tsPopChange), f=rData$REGION[match(rownames(tsPopChange), rData$datasetid)])
table(tsPopSplit[[1]][,1], tsPopSplit[[1]][,2])
table(tsPopSplit[[2]][,1], tsPopSplit[[2]][,2])

sum(tsPopSplit[[2]][,1] < tsPopSplit[[2]][,2]) / nrow(tsPopSplit[[2]])

#           ellis exploration ####

contCol = c("blue","darkgreen")
contAlpha = col2rgb(contCol)/255
contAlpha = apply(contAlpha, 2, function(x){rgb(x[1],x[2],x[3],0.25)})

# consistent vs changing time series
a = table(cut(rData$age, breaks=seq(0,8000,1000)), rData$ellisPopDens, rData$REGION)
#a = tapply(rData$ellisPopDens, rData$datasetid, function(x){length(unique(x))})
barplot(t(prop.table(a[,,1], 1)))
barplot(t(prop.table(a[,,2], 1)))

# modern pop versus historical pop categories
boxplot(log(rData$popDens+1) ~ rData$ellisPopDen * rData$REGION,
        col=c(rep("blue", 5), rep("darkgreen", 5)), yaxt="n", xlab="", ylab="")
axis(side=2, at=log(c(0, 1,10,100,1000,10000)+1), labels=c(0,1,10,100,1000,10000), las=1)
axis(side=2, at=log(c(seq(1,10,1), seq(10,100,10), seq(100,1000,100), seq(1000,10000,1000))),
     labels=NA, tcl=-0.125)
mtext(side=2, line=2, text=expression("Population density at 2000AD (persons km"^-2*")"))

# test for modern pop ability to explain historical prop
rData$ellisPopNum <- as.numeric(as.character(rData$ellisPopDens))
histPopM <- lmer(log(ellisPopNum+1) ~ log(popDens+1) * REGION + (1|datasetid), data=rData)
performance(histPopM)
summary(histPopM)

saveRDS(rData, "./outputs/holoDataGenus.rds")

# ellis pop dens as proportions across Holocene
ellisEurope = st_crop(ellis, st_bbox(c(xmin=Eu[3], xmax=Eu[4], ymin=Eu[1], ymax=Eu[2]), crs=st_crs(4326)))
ellisNthAm = st_crop(ellis, st_bbox(c(xmin=nthAm[3], xmax=nthAm[4], ymin=nthAm[1], ymax=nthAm[2]), crs=st_crs(4326)))

# pull out and aggregate ellis data for each 1000 years
ellisDataEurope = ellisData[match(ellisEurope$id, ellisData$id),]
ellisDataNthAm = ellisData[match(ellisNthAm$id, ellisData$id),]

# now aggregate data for each time aggregate (thousand years?)
ellisColsCut = as.factor(round(ellisAges, -2))

ellisFullDens = do.call("rbind", lapply(unique(ellisColsCut), function(n){
  
  a = ellisDataEurope[,-1][,ellisColsCut == n & !is.na(ellisColsCut)]
  if(class(a) == "data.frame"){
    a = unlist(a)
  }
  
  aDens = factor(ellisCat$ellisPopDens[match(a, ellisCat$ellisID)], levels=unique(ellisCat$ellisPopDens))
  aDens = table(aDens) / length(aDens)
  
  b = ellisDataNthAm[,-1][,ellisColsCut == n & !is.na(ellisColsCut)]
  if(class(b) == "data.frame"){
    b = unlist(b)
  }
  bDens = factor(ellisCat$ellisPopDens[match(b, ellisCat$ellisID)], levels=unique(ellisCat$ellisPopDens))
  bDens = table(bDens) / length(bDens)
  
  return(data.frame(REGION = rep(c("Europe", "North America"), each = length(aDens)),
                    time = n,
                    popCat = rep(names(aDens), 2),
                    popDens = c(aDens, bDens)))
  
}))

#           Model ####

# remove missing pop dens & very low sampled high pop dens categories
rModel = droplevels(rData[!is.na(rData$ellisPopDens) & rData$ellisPopDens != "2500",])

tempA = gamm4(tpPC1 ~ s(ageS, k = 12, by=ellisDensRegion) + ellisPopDens + REGION * latS, data=rModel,
              random = ~(1|datasetid))

predDF = cbind(expand.grid(ageS = seq(min(rData$ageS, na.rm=TRUE), max(rData$ageS, na.rm=TRUE), len=1000),
                           ellisDensRegion = sort(unique(rModel$ellisDensRegion))),
               latS = 0,
               tempPred = 0)
regCats = do.call('rbind', strsplit(as.character(predDF$ellisDensRegion), ":"))
colnames(regCats) = c("REGION", "ellisPopDens")
predDF = cbind(predDF, regCats)

predDF = cbind(predDF, as.data.frame(predict(tempA$gam, newdata=predDF, se.fit=TRUE)))
predDF$age = predDF$ageS * sd(rData$age, na.rm=TRUE) + mean(rData$age, na.rm=TRUE)
predDF = predDF[predDF$age <= 12000,]

gradientCutoff = floor(0.005 * length(unique(rData$datasetid)))

europePopSamplingLims = sapply(split(rData[rData$REGION == "Europe",], 
                                     f=rData$ellisPopDens[rData$REGION == "Europe"]),function(x){
                                       sort(tapply(x$age, x$datasetid, max), decreasing=TRUE)[gradientCutoff]
                                     })
AmpopSamplingLims = sapply(split(rData[rData$REGION == "North America",], 
                                 f=rData$ellisPopDens[rData$REGION == "North America"]),function(x){
                                   sort(tapply(x$age, x$datasetid, max), decreasing=TRUE)[gradientCutoff]
                                 })

#           Plot ####

library(pBrackets)

pdf("./plots/holoceneLineFamily.pdf", height=6, width=11, useDingbats=FALSE)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

split.screen(rbind(c(0.12,0.535,0.3,0.95),
                   c(0.575,0.99,0.3,0.95),
                   
                   c(0.12,0.535,0.1,0.3), # population densities
                   c(0.575,0.99,0.1,0.3),
                   
                   c(0.88, 0.92, 0.325, 0.5))) # pop dens gradient

# try a linear plot instead of a surface plot
screen(1)
plot(NULL, xlim=c(12200,0), xlab="", ylab="", xaxt="n", xaxs="i", ylim=c(-0.3,0.2))  #ylim=c(-1.65,1)
abline(h=0, lty="31")
axis(side=1, at=seq(0,12000,2000), labels=NA)
mtext(side=2, line=4, text="Assemblage persistence PC1", las=0)
mtext(side=3, line=0, at=par('usr')[1], text="Europe", font=2,col="blue", adj=0)

par(xpd=NA)
arrowLabs = c("More\ntransient", "More\npersistent")
Arrows(y0=0, y1= 0 + c(-0.15, 0.15),
       x0=relative.axis.point(-0.12, "x"), x1=relative.axis.point(-0.12, "x"),
       arr.type="triangle", col="black", arr.width=0.1, arr.length=0.1, lwd=2)
text(y=0 + c(-0.15, 0.15), x=relative.axis.point(-0.12, "x"),
     pos=c(1,3), adj=0,
     labels=arrowLabs,
     font=3, col="black", cex=0.8)
par(xpd=FALSE)

europeCols = c(colorRampPalette(c("grey70", "blue"))(5)[-1], "darkblue")
x = split(predDF[grepl("Europe", predDF$ellisDensRegion),], f=predDF$ellisDensRegion[grepl("Europe", predDF$ellisDensRegion)])
x = x[sapply(x, nrow) > 0]
mapply(y=x, n=seq_along(x), function(y,n){
  
  y = y[y$age <= europePopSamplingLims[n],]
  
  colrgb = col2rgb(europeCols[n])/255
  
  polygon(x=c(y$age, rev(y$age)),
          y=c(y$fit + 1.96 * y$se.fit, rev(y$fit - 1.96 * y$se.fit)),
          border=NA, col=rgb(colrgb[1],colrgb[2],colrgb[3],0.25))
  lines(y$fit ~ y$age, col=europeCols[n], lwd=2, lty=ifelse(n==1, "31", "solid"))
  text(x=y$age[1], y=y$fit[1], labels=c("0","<1","<10","<100", ">100")[n],
       col=europeCols[n], pos=4, offset=0.2) 
})

# add in 2500 cat as individual non-modelled point
points(x=750, y=mean(rData$tpPC1[rData$ellisDensRegion == "Europe:2500"]),
       pch=16, col="darkblue")
text(x=750, y=mean(rData$tpPC1[rData$ellisDensRegion == "Europe:2500"]), 
     labels=">100*", col="darkblue", pos=2, offset=0.2) 

#brackets(x2=max(x[[1]]$age), x1=8800, y1=relative.axis.point(0.75, "y"), y2 = relative.axis.point(0.75, "y"), h=0.15, type=1)
#segments(x0=8800, x1=8800, y0=relative.axis.point(0.75, "y"), y1=x[[2]]$fit[which.min(abs(8800-x[[2]]$age))], lty="dotted")
#par(lheight=0.8)
#text(x=mean(c(max(x[[1]]$age), 8800)), y=relative.axis.point(0.87, "y"), labels="(1) Increased\npersistence with\nhuman presence", font=3, adj=0.5)

brackets(x2=max(x[[1]]$age), x1=europePopSamplingLims[4], y1=relative.axis.point(0.75, "y"), y2 = relative.axis.point(0.75, "y"), h=0.03, type=1)
segments(x0=max(x[[1]]$age), x1=max(x[[1]]$age), y0=relative.axis.point(0.75, "y"), y1=x[[2]]$fit[which.max(x[[2]]$age)], lty="dotted")
segments(x0=europePopSamplingLims[4], x1=europePopSamplingLims[4], y0=relative.axis.point(0.75, "y"), y1=x[[4]]$fit[which.min(abs(europePopSamplingLims[4]-x[[4]]$age))], lty="dotted")
text(x=mean(c(max(x[[1]]$age), 2500)), y=relative.axis.point(0.87, "y"), labels="(0) Human presence\nhas no detectable\npersistence signal", font=3, adj=0.5)

brackets(x2=europePopSamplingLims[4], x1=min(x[[1]]$age), y1=relative.axis.point(0.75, "y"), y2 = relative.axis.point(0.75, "y"), h=0.03, type=1)
segments(x0=min(x[[1]]$age), x1=min(x[[1]]$age), y0=relative.axis.point(0.75, "y"), y1=x[[4]]$fit[which.min(x[[4]]$age)], lty="dotted")
text(x=mean(c(min(x[[1]]$age), europePopSamplingLims[4])), y=relative.axis.point(0.87, "y"), labels="(-) Decreased\npersistence with\nhuman presence", font=3, adj=0.5)
text(x=relative.axis.point(0.01, "x"), y=relative.axis.point(0.975, "y"),
     adj=0, labels="(A)", font=2)
close.screen(1)

screen(2)
plot(NULL, xlim=c(12200,0), xlab="", ylab="", xaxt="n", xaxs="i", ylim=c(-0.3,0.3)) #ylim=c(-1.25,1.75)
abline(h=0, lty="31")
axis(side=1, at=seq(0,12000,2000), labels=NA)
mtext(side=3, line=0, at=par('usr')[1], text="North America", font=2,col="darkgreen", adj=0)

# lines(tempPred[,2] ~ tempPred[,1], col="red")
amCols = c(colorRampPalette(c("grey70", "darkgreen"))(5)[-1], "#003000FF")
x = split(predDF[grepl("North America", predDF$ellisDensRegion),], f=predDF$ellisDensRegion[grepl("North America", predDF$ellisDensRegion)])
x = x[sapply(x, nrow) > 0]

mapply(y=x, n=seq_along(x), function(y,n){
  
  y = y[y$age <= AmpopSamplingLims[n],]
  
  colrgb = col2rgb(amCols[n])/255
  
  polygon(x=c(y$age, rev(y$age)),
          y=c(y$fit + 1.96 * y$se.fit, rev(y$fit - 1.96 * y$se.fit)),
          border=NA, col=rgb(colrgb[1],colrgb[2],colrgb[3],0.25), lwd=0.5)
  lines(y$fit ~ y$age, col=amCols[n], lwd=2, lty=ifelse(n==1, "31", "solid"))
  text(x=y$age[1], y=y$fit[1], labels=c("0","<1","<10","<100", ">100")[n],
       col=amCols[n], pos=4, offset=0.2) 
})

brackets(x2=max(x[[1]]$age), x1=9400, y1=relative.axis.point(0.75, "y"), y2 = relative.axis.point(0.75, "y"), h=0.03, type=1)
segments(x0=max(x[[1]]$age), x1=max(x[[1]]$age), y0=relative.axis.point(0.75, "y"), y1=x[[2]]$fit[which.max(x[[2]]$age)], lty="dotted")
segments(x0=9400, x1=9400, y0=relative.axis.point(0.75, "y"), y1=x[[2]]$fit[which.min(abs(9400-x[[2]]$age))], lty="dotted")
par(lheight=0.8)
text(x=mean(c(max(x[[1]]$age), 9400)), y=relative.axis.point(0.87, "y"), labels="(0) Human presence\nhas no detectable\npersistence signal", adj=0.5, font=3)

brackets(x2=9400, x1=5400, y1=relative.axis.point(0.75, "y"), y2 = relative.axis.point(0.75, "y"), h=0.03, type=1)
segments(x0=5400, x1=5400, y0=relative.axis.point(0.75, "y"), y1=x[[2]]$fit[which.min(abs(5400-x[[2]]$age))], lty="dotted")
par(lheight=0.8)
text(x=mean(c(max(x[[1]]$age), 5400)), y=relative.axis.point(0.87, "y"), labels="(+) Increased\npersistence with\nhuman presence", adj=0.5, font=3)


brackets(x2=5400, x1=min(x[[1]]$age), y1=relative.axis.point(0.75, "y"), y2 = relative.axis.point(0.75, "y"), h=0.03, type=1)
segments(x0=min(x[[1]]$age), x1=min(x[[1]]$age), y0=relative.axis.point(0.75, "y"), y1=x[[4]]$fit[which.min(x[[4]]$age)], lty="dotted")
text(x=mean(c(min(x[[1]]$age), 5400)), y=relative.axis.point(0.87, "y"), labels="(0) Human presence\nhas no detectable\npersistence signal", font=3, adj=0.5)
text(x=relative.axis.point(0.01, "x"), y=relative.axis.point(0.975, "y"),
     adj=0, labels="(B)", font=2)
close.screen(2)

screen(3)
# unsampled windows of higher density human populations
plot(NULL, xlim=c(12200, 0), ylim=c(0,1), axes=FALSE, xlab="", ylab="", yaxs="i", xaxs="i")
axis(side=1, at=seq(2000,11500,2000), tcl=-0.25, 
     labels=format(seq(2000,11500,2000), big.mark=","), mgp=c(3,0.2,0))
axis(side=1, at=0, labels=0, mgp=c(3,0.2,0))
axis(side=1, at=seq(0,12000,1000), tcl=-0.125, labels=NA)

mtext(side=1, line=1.25, text="Years before present")
axis(side=2, at=seq(0,1,0.2))
mtext(side=2, line=4, text="Human population density", las=0, lheight=0.75)
mtext(side=2, line=2.5, text="(Proportion of\ncontinental area)", las=0, lheight=0.75, font=3)

# europe pop dens
ellisEurWide = ellisFullDens[ellisFullDens$REGION=="Europe",]
ellisEurWide = as.data.frame(pivot_wider(ellisEurWide[,-1], values_from="popDens", names_from="popCat"))
ellisEurWide[,"2500"] = rowSums(ellisEurWide[,c("2500","5000")])
rownames(ellisEurWide) = ellisEurWide[,1]
ellisEurWide = ellisEurWide[,!colnames(ellisEurWide) %in% c("time", "5000")]
ellisEurWide = ellisEurWide[,ncol(ellisEurWide):1]

europeEllisPlot = rbind(0, apply(ellisEurWide, 1, cumsum))

europeCols = c(colorRampPalette(c("grey70", "blue"))(5)[-1], "darkblue")

sapply(2:(nrow(europeEllisPlot)), function(n){
  print(n)
  polygon(x=c(as.numeric(colnames(europeEllisPlot)),
              rev(as.numeric(colnames(europeEllisPlot)))),
          y=c(europeEllisPlot[n-1,], rev(europeEllisPlot[n,])),
          col=europeCols[n-1])
})

sapply(1:ncol(ellisEurWide), function(n){
  labX = weighted.mean(as.numeric(rownames(ellisEurWide)), ellisEurWide[,n])
  labY = mean(europeEllisPlot[,which.min(abs(as.numeric(rownames(ellisEurWide)) - labX))][n:(n+1)])
  text(x=labX, y=labY, label=c("0", "<1", "<10", "<100")[n],
       col=c("black","black","white","white")[n], cex=0.75)
})
box()
text(x=relative.axis.point(0.01, "x"), y=relative.axis.point(0.925, "y"),
     adj=0, labels="(C)", font=2)

axis(side=3, at=seq(0,12000,2000), tcl=0.25, labels=NA)
axis(side=3, at=seq(0,12000,1000), tcl=0.125, labels=NA)

close.screen(3)

screen(4)
# unsampled windows of higher density human populations
plot(NULL, xlim=c(12200, 0), ylim=c(0,1), axes=FALSE, xlab="", ylab="", yaxs="i", xaxs="i")
axis(side=1, at=seq(2000,11500,2000), tcl=-0.25, 
     labels=format(seq(2000,11500,2000), big.mark=","), mgp=c(3,0.2,0))
axis(side=1, at=seq(0,12000,1000), tcl=-0.125, labels=NA)
axis(side=1, at=0, labels=0, mgp=c(3,0.2,0))
mtext(side=1, line=1.25, text="Years before present")
axis(side=2, at=seq(0,1,0.2))

# europe pop dens
ellisEurWide = ellisFullDens[ellisFullDens$REGION=="North America",]
ellisEurWide = as.data.frame(pivot_wider(ellisEurWide[,-1], values_from="popDens", names_from="popCat"))
ellisEurWide[,"2500"] = rowSums(ellisEurWide[,c("2500","5000")])
rownames(ellisEurWide) = ellisEurWide[,1]
ellisEurWide = ellisEurWide[,!colnames(ellisEurWide) %in% c("time", "5000")]
ellisEurWide = ellisEurWide[,ncol(ellisEurWide):1]

europeEllisPlot = rbind(0, apply(ellisEurWide, 1, cumsum))

europeCols = c(colorRampPalette(c("grey70", "darkgreen"))(5)[-1], "#003000FF")

sapply(2:(nrow(europeEllisPlot)), function(n){
  print(n)
  polygon(x=c(as.numeric(colnames(europeEllisPlot)),
              rev(as.numeric(colnames(europeEllisPlot)))),
          y=c(europeEllisPlot[n-1,], rev(europeEllisPlot[n,])),
          col=europeCols[n-1])
})

sapply(1:ncol(ellisEurWide), function(n){
  labX = weighted.mean(as.numeric(rownames(ellisEurWide)), ellisEurWide[,n])
  labY = mean(europeEllisPlot[,which.min(abs(as.numeric(rownames(ellisEurWide)) - labX))][n:(n+1)])
  text(x=labX, y=labY, label=c("0", "<1", "<10", "<100")[n],
       col=c("black","black","white","white")[n], cex=0.75)
})
box()
text(x=relative.axis.point(0.01, "x"), y=relative.axis.point(0.925, "y"),
     adj=0, labels="(D)", font=2)
axis(side=3, at=seq(0,12000,2000), tcl=0.25, labels=NA)
axis(side=3, at=seq(0,12000,1000), tcl=0.125, labels=NA)
close.screen(4)

screen(5)
image(y=seq(0,1,len=5),
      x=c(0,1),
      z=matrix(1:10, byrow=TRUE, nrow=2),
      col=c(c(colorRampPalette(c("grey70", "blue"))(5)[-1], "darkblue"),
            c(colorRampPalette(c("grey70", "darkgreen"))(5)[-1], "#003000FF")),
      useRaster=TRUE, axes=FALSE, xlab="", ylab="")
axis(side=2, mgp=c(3,0.5,0), at=seq(0,1,len=5),
     labels=c("0", "<1",  "<10", "<100", ">100"))
abline(v=0.5,)
box()
mtext(side=2, line=3.25, text="Pop density", las=0)
mtext(side=2, line=2.25, text=expression("(persons km"^-2*")"), las=0)

close.screen(5)

close.screen(all.screens=TRUE)
dev.off()

####################################### ####
# FIG 4: CONCEPTUAL HUMAN DRIVER FIGURE ####
####################################### ####

pdf("./plots/persConcept.pdf", height=5, width=4.5, useDingbats = FALSE)
par(mar=c(3,3,0.5,0.5), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

plot(NULL, xlim=c(0.09,0.52), ylim=c(-1,1), xaxt="n", ylab="", xlab="")
mtext(side=2, line=2, las=0, text="Persistence")

xDivs = c(par("usr")[1], 0.225, 0.375, par("usr")[2])
rect(xleft=xDivs[-length(xDivs)], xright=xDivs[-1],
     ybottom=par("usr")[3], ytop=par("usr")[4], col=c("grey80", "white", "grey80"),
     border = NA)
abline(h=0, lty="31")

axis(side=1, at=seq(0.1,0.5,0.05), labels=NA)

sapply(1:3, function(x){
  xLims = list(c(0.09, 0.21),
               c(0.24,0.36),
               c(0.39,0.51))[[x]]
  
  axis(side=1, line=2.5, labels=NA, at=xLims, tcl=0.25)
  mtext(side=1, line=2.65,
        at=mean(xLims),
        text=c("Early Holocene", "Mid Holocene", "Late Holocene")[x])
})
box()
# EARlY

earlyNo = c(-0.3,0,0,0,0)
Arrows(x0=seq(0.09,0.11,len=5)[3], 
       x1=seq(0.09,0.11,len=5)[3], 
       y0=sum(earlyNo), 
       y1=sum(earlyNo) + earlyNo[earlyNo != 0],
       arr.type = "triangle", arr.length=0.15, arr.width = 0.15, col="red",
       lwd=2)
points(x=0.1, y = sum(earlyNo), pch=16)
#segments(x0=0.08, x1=0.12, y0=-0.2, y1=-0.2, lwd=2)

earlyLight = c(-0.3,0.3,0,0,0)
Arrows(x0=seq(0.135,0.165,len=5)[c(3,3)], 
       x1=seq(0.135,0.165,len=5)[c(3,3)], 
       y0=sum(earlyLight), 
       y1=sum(earlyLight) + earlyLight[earlyLight != 0],
       col=c("red", "darkgreen"),
       arr.type = "triangle", arr.length=0.15, arr.width = 0.15, 
       lwd=2)
points(x=0.15, y = sum(earlyLight), pch=16)

earlyHeavy = c(-0.3,0.6,-0.1,0,0)
Arrows(x0=seq(0.185,0.215,len=5)[c(1,3,5)], 
       x1=seq(0.185,0.215,len=5)[c(1,3,5)], 
       y0=sum(earlyHeavy), 
       y1=sum(earlyHeavy) + earlyHeavy[earlyHeavy != 0],
       col=c("red", "darkgreen", "orange"),
       arr.type = "triangle", arr.length=0.15, arr.width = 0.15, 
       lwd=2)
points(x=0.2, y = sum(earlyHeavy), pch=16)
segments(x0=0.185, x1=0.215, y0=sum(earlyHeavy), y1=sum(earlyHeavy), lwd=2)

midNo = c(-0.05,0,0,0,0)
Arrows(x0=seq(0.235,0.265,len=5)[3], 
       x1=seq(0.235,0.265,len=5)[3], 
       y0=sum(midNo), 
       y1=sum(midNo) + midNo[midNo != 0],
       col=c("red"),
       arr.type = "triangle", arr.length=0.15, arr.width = 0.15, 
       lwd=2)
points(x=0.25, y = sum(midNo), pch=16)

midLight = c(-0.05,0.3,-0.1,-0.1)
Arrows(x0=seq(0.285,0.315,len=4), 
       x1=seq(0.285,0.315,len=4), 
       y0=sum(midLight), 
       y1=sum(midLight) + midLight[midLight != 0],
       col=c("red", "darkgreen", "orange", "blue"),
       arr.type = "triangle", arr.length=0.15, arr.width = 0.15, 
       lwd=2)
points(x=0.3, y = sum(midLight), pch=16)
segments(x0=0.285, x1=0.315, y0=sum(midLight), y1=sum(midLight), lwd=2)

midHeavy = c(-0.05,0.6,-0.3,-0.3,0)
Arrows(x0=seq(0.335,0.365,len=4), 
       x1=seq(0.335,0.365,len=4), 
       y0=sum(midHeavy), 
       y1=sum(midHeavy) + midHeavy[midHeavy != 0],
       col=c("red", "darkgreen", "orange", "blue"),
       arr.type = "triangle", arr.length=0.15, arr.width = 0.15, 
       lwd=2)
text(x=seq(0.335,0.365,len=4), 
     y=sum(midHeavy) + midHeavy[midHeavy != 0],
     pos=c(1,3,1,1),
     col=c("red", "darkgreen", "orange", "blue"),
     labels=c("Clim", "Use", "Coll", expression(Delta*"Use")), cex=0.75)
points(x=0.35, y = sum(midHeavy), pch=16)
segments(x0=0.335, x1=0.365, y0=sum(midHeavy), y1=sum(midHeavy), lwd=2)


lateNo = c(0.1,0,0,0,0)
Arrows(x0=seq(0.39,0.41,len=5)[3], 
       x1=seq(0.39,0.41,len=5)[3], 
       y0=sum(lateNo), 
       y1=sum(lateNo) + lateNo[lateNo != 0],
       col=c("red"),
       arr.type = "triangle", arr.length=0.15, arr.width = 0.15, 
       lwd=2)
points(x=0.4, y = sum(lateNo), pch=16)

lateMid = c(0.1,0.3,-0.2,-0.2)  
Arrows(x0=seq(0.435,0.465,len=4), 
       x1=seq(0.435,0.465,len=4), 
       y0=sum(lateMid), 
       y1=sum(lateMid) + lateMid[lateMid != 0],
       col=c("red", "darkgreen", "orange", "blue"),
       arr.type = "triangle", arr.length=0.15, arr.width = 0.15, 
       lwd=2)
points(x=0.45, y = sum(lateMid), pch=16)
segments(x0=0.435, x1=0.465, y0=sum(lateMid), y1=sum(lateMid), lwd=2)


lateHeavy = c(0.1,0.6,-0.5,-0.5)  
Arrows(x0=seq(0.485,0.515,len=4), 
       x1=seq(0.485,0.515,len=4), 
       y0=sum(lateHeavy), 
       y1=sum(lateHeavy) + lateHeavy[lateHeavy != 0],
       col=c("red", "darkgreen", "orange", "blue"),
       arr.type = "triangle", arr.length=0.15, arr.width = 0.15, 
       lwd=2)
points(x=0.5, y = sum(lateHeavy), pch=16)
segments(x0=0.485, x1=0.515, y0=sum(lateHeavy), y1=sum(lateHeavy), lwd=2)

dev.off()

############################### ####
# FIG 1: PCA PERSISTENCE FIGURE ####
############################### ####

set.seed(110037)

fakeMat = t(matrix(rpois(5, 2:6), nrow=5, ncol=10)) + rnorm(50, 0, 1)
fakeMat[fakeMat<0] = 0
fakeMat = cbind(fakeMat,
                seq(1,10, len=10) + rnorm(10, 0, 1))
fakeMat = cbind(fakeMat,
                rev(seq(1,10, len=10)) + rnorm(10, 0, 1))
fakeMat = prop.table(fakeMat, 1)

testOrd = metaMDS(fakeMat)

examplePos = as.matrix(vegdist(fakeMat))[5,][6:8]
examplePC = logit(examplePos) %*% tpResPCA$loadings
exampleMean = mean(examplePos)
exampleRatio = examplePos[1] / examplePos[3]

plotBG <- function(col, plotBorder=TRUE){
  rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop = par("usr")[4],
       col=col, border=ifelse(plotBorder, "black", NA))
}

pdf("./plots/pcaDidactic2Genus.pdf", height=4.25, width=6.25)

arrowScale = 1

split.screen(rbind(c(0.35,0.6,0.7,0.99),
                   c(0.09,0.6,0.15,0.99),
                   c(0.7,0.99,0.62,0.99),
                   c(0.7,0.99,0.15,0.52)))

pointSub = sample(1:nrow(testData),10000)
xlims=c(-1.5,1)
ylims=c(-0.6,1.3)

screen(2)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(testData$tpPC2[pointSub] ~ testData$tpPC1[pointSub],  xlab="", ylab="", axes=FALSE,type="n",
     xlim=xlims, ylim=ylims)
points(testData$tpPC2[pointSub] ~ testData$tpPC1[pointSub],pch=16, col="grey40", cex=0.6)
segments(x0=0, x1=0, y0=par("usr")[3], y1=relative.axis.point(0.9, "y"), lty="31")
abline(h=0,lty="31")

tpPCval = tpResStore$sdev^2 / sum(tpResStore$sdev^2)

Arrows(x0=0, y0=0, x1=tpResStore$loadings[,1]*arrowScale, y1=tpResStore$loadings[,2]*arrowScale, arr.type="triangle",
       arr.width=0.15, arr.length=0.15, lwd=2, col="blue")
text(x=tpResStore$loadings[,1]*arrowScale, y=tpResStore$loadings[,2]*arrowScale,
     labels=c(expression("D"["+3"]), expression("D"["+2"]), expression("D"["+1"])), 
     pos=2, col="blue")
axis(side=1, mgp=c(3,0.1,0))
axis(side=1, at=-10:10, tcl=-0.125, labels=NA)
axis(side=2)
axis(side=2, at=-10:10, tcl=-0.125, labels=NA)
mtext(side=1, line=0.75, text=paste0("PC1 (", sprintf("%.2f", tpPCval[1]*100), "%)"))
mtext(side=1, line=1.5, text='"State persistence"')
mtext(side=2, line=2, text=paste0("PC2 (", sprintf("%.2f", tpPCval[2]*100), "%)"), las=0)
mtext(side=2, line=1.25, text='"State attractor strength"', las=0)
box()
text(x=relative.axis.point(0.02, "x"), y=relative.axis.point(0.965, "y"),
     labels=expression(bold("(A)")*" Assemblage persistence"), adj=0)
text(x=relative.axis.point(0.02, "x"), y=relative.axis.point(0.93, "y"),
     labels="principal component analysis", adj=0)

close.screen(2)

tpPCval = tpResPCA$sdev^2 / sum(tpResPCA$sdev^2)

screen(1)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)

plot(testOrd$points, axes=FALSE, xlab="", ylab="", 
     xlim=range(testOrd$points[,1])*1.1, ylim=range(testOrd$points[,2])*1.1 + c(0,0.05),
     type="n")

plotBG("white")

sapply(2:nrow(testOrd$points), function(n){
  lines(testOrd$points[(n-1):n,], col = c(rep("black",6), rep("black", 4))[n])
})

box()
mtext(side=1, line=0.1, text="nMDS1")
mtext(side=2, line=0.1, text="nMDS2", las=0)

backOrd = testOrd$points[nrow(testOrd$points):1,]
Arrows(x0 = backOrd[-nrow(backOrd),1],
       y0 = backOrd[-nrow(backOrd),2],
       x1 = backOrd[-nrow(backOrd),1] + 0.5 * diff(backOrd[,1]),
       y1 = backOrd[-nrow(backOrd),2] + 0.5 * diff(backOrd[,2]),
       arr.type="triangle", arr.length=0.1, arr.width=0.1)

points(testOrd$points, axes=FALSE,
       pch=c(rep(21,6), rep(16, 4)), 
       cex=0.75,
       col = c(rep("black",6), rep("black", 4)),
       bg="white")

segments(x0=testOrd$points[6,1], y0=testOrd$points[6,2],
         x1=testOrd$points[3:5,1], y1=testOrd$points[3:5,2], lty="31", col="blue", lwd=2)
#text(testOrd$points, labels=format(5000 + 500 * 1:nrow(testOrd$points), big.mark=","), pos=4)

points(y=testOrd$points[3:5,2], x=testOrd$points[3:5,1], pch=21, col="blue", bg="white", cex=1, lwd=2)
text(y=testOrd$points[3:5,2], x=testOrd$points[3:5,1], 
     labels=c(expression("D"["+3"]), expression("D"["+2"]), expression("D"["+1"])), pos=c(1,2,4), col="blue")
box()
points(y=testOrd$points[6,2], x=testOrd$points[6,1], pch=21, bg="red", cex=1.5)
text(x=relative.axis.point(0.02, "x"), y=relative.axis.point(0.935, "y"),
     labels=expression(bold("(B)")*" Persistence window"), adj=0)
close.screen(1)

screen(3)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)

plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(-1.5,1), xlab="", ylab="", xaxt="n")
dissMean = rowMeans(testData[,c("dn1_1", "dn1_2", "dn1_3")])

points(testData$tpPC1[pointSub] ~ dissMean[pointSub], pch=16, col="grey40", cex=0.6)

axis(side=1, mgp=c(3,0.2,0))
mtext(side=1, line=1, text = expression("Mean D"["+1-3"]))
mtext(side=2, line=2, las=0, text="PC1 (persistence)")
text(x=relative.axis.point(0.02, "x"), y=relative.axis.point(0.05, "y"),
     adj=0, labels = paste0("R = ", sprintf("%.3f", cor(testData$tpPC1, dissMean))))
text(x=relative.axis.point(0.02, "x"), y=relative.axis.point(0.935, "y"),
     labels=expression(bold("(C)")), adj=0)
close.screen(3)

screen(4)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)

plot(x=NULL, y=NULL, xlim=c(-3,4), ylim=c(-0.5,0.65), xlab="", ylab="", xaxt="n")

dissRat = log(beta.tr(testData$dn1_3) / beta.tr(testData$dn1_1))

points(testData$tpPC2[pointSub] ~ dissRat[pointSub], pch=16, col="grey40", cex=0.6)

axis(side=1, mgp=c(3,0.2,0))
mtext(side=1, line=1, text = expression("Log ratio (D"["+3"]*":D"["+1"]*")"))
mtext(side=2, line=2, las=0, text="PC2 (attractor strength)")
text(x=relative.axis.point(0.02, "x"), y=relative.axis.point(0.05, "y"),
     adj=0, labels = paste0("R = ", sprintf("%.3f", cor(testData$tpPC2, dissRat))))
text(x=relative.axis.point(0.02, "x"), y=relative.axis.point(0.935, "y"),
     labels=expression(bold("(D)")), adj=0)
close.screen(4)

close.screen(all.screens=TRUE)
dev.off()
######################## ####
# SUPPLEMENTARY ANALYSES ####
######################## ####
# random intercepts comparison ####

# check both persistence models (late glacial vs holocene) ####

lgRI = ranef(timeOnlyMb$mer)$datasetid
lgRI$datasetid = rownames(lgRI)
colnames(lgRI)[1] = "lateGlac"

holoRI = ranef(tempA$mer)$datasetid
holoRI$datasetid = rownames(holoRI)
colnames(holoRI)[1] = "holo"

combRI = merge(lgRI, holoRI, all.x=TRUE, all.y=TRUE, sort=FALSE)
plot(combRI$lateGlac ~ combRI$holo)
abline(0,1)

combRI$diff = combRI$lateGlac - combRI$holo

cor(combRI$holo, combRI$diff, use="complete.obs")
# what's up with large deviations?

rData[rData$datasetid %in% combRI$datasetid[combRI$diff < -0.2 & !is.na(combRI$diff)],]





# where are the high density sites? ####

highSites = rData[rData$ellisPopNum >= 10 & !duplicated(rData$datasetid), ]

plot(europeMap)
points(highSites[highSites$REGION == "Europe",c("long","lat")], pch=16, col="blue")

plot(nthAmMap)
points(highSites[highSites$REGION == "North America",c("long","lat")], pch=16, col="darkgreen")


# old trends model ####

testPred = split(predDFB, f=grepl("Europe", predDFB$ellisDensRegion))
dev.off()
plot(NULL, xlim=c(12000,0), ylim=c(-0.5,0.5))
points(testPred[[1]]$fit ~ testPred[[1]]$age, pch=16,
       col=testPred[[1]]$ellisDensRegion)

plot(NULL, xlim=c(12000,0), ylim=c(-0.5,0.5))
points(testPred[[2]]$fit ~ testPred[[2]]$age, pch=16,
       col=testPred[[2]]$ellisDensRegion)


pdf("./plots/ellisTrendsGenus.pdf", height=5, width=10, useDingbats=FALSE)

par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

split.screen(rbind(c(0.1,0.475,0.70,0.97),
                   c(0.525,0.9,0.70,0.97),
                   c(0.1,0.475,0.40,0.67),
                   c(0.525,0.9,0.40,0.67), # time trends
                   
                   c(0.91,0.93,0.735,0.935), # color gradients
                   c(0.91,0.93,0.435,0.635),
                   
                   c(0.1,0.475,0.1,0.37), # population densities
                   c(0.525,0.9,0.1,0.37),
                   
                   c(0.91, 0.93, 0.135, 0.335))) # pop dens gradient


tp1Lims = c(-1.4,1.4)
tp2Lims = c(-0.16,0.16)
contLevelsTP1 = seq(tp1Lims[1],tp1Lims[2],0.1)
contLevelsTP2 = seq(tp2Lims[1],tp2Lims[2],0.02)

tp1Ramp = colorRampPalette(c("red","white","blue"))(length(contLevelsTP1)-1)
tp2Ramp = colorRampPalette(c("orange","white","purple"))(length(contLevelsTP2)-1)

# what ages do we stop reporting predictions because there aren't sufficient
# observations for them to be reliable?
gradientCutoff = floor(0.01 * length(unique(rData$datasetid)))

# TP ramps
popCutoffs = sapply(1:4, function(n){
print(n)
screen(n)

if(n < 3){
fitMat = matrix(testPred[[n]]$fit, nrow=200, ncol=4)
subLims = tp1Lims
subContour = contLevelsTP1
subRamp = tp1Ramp
  
} else {
  fitMat = matrix(testPred2[[n-2]]$fit, nrow=200, ncol=4)
  subLims = tp2Lims
  subContour = contLevelsTP2
  subRamp = tp2Ramp
}
  
region = ifelse(n %% 2 == 1, "Europe", "North America")

# hide predictions we don't have observations for
popSamplingLims = sapply(split(rData[rData$REGION == region,], 
                               f=rData$ellisPopDens[rData$REGION == region]),function(x){
                                 sort(tapply(x$age, x$datasetid, max), decreasing=TRUE)[gradientCutoff]
                               })
popSamplingLims[is.na(popSamplingLims)] = min(rData$age, na.rm=TRUE)
fitMat = sapply(1:length(popSamplingLims), function(n){
  tempMat = fitMat[,n]
  tempMat[sort(unique(testPred[[1]]$age)) > popSamplingLims[n]] = NA
  return(tempMat)
}) 
#europeMat = europeMat[nrow(europeMat):1,]

image(x=sort(unique(testPred[[1]]$age)),
      y=sort(unique(testPred[[1]]$ellisMag)),
      z=fitMat,
      col=subRamp,
      zlim=subLims, xlim=c(11500,0),
      axes=FALSE, xlab="", ylab="", useRaster=TRUE)

sapply(1:length(popSamplingLims), function(n){

  if(which.max(is.na(fitMat[,n])) == 1){return(NULL)}
  
  contour(x=sort(unique(testPred[[1]]$age)),
          y=sort(unique(testPred[[1]]$ellisMag))[n] + c(-0.5,0.5),
          z=matrix(fitMat[,n], nrow=200, ncol=2),
          add=TRUE, levels=subContour, method="edge")
  
  rect(xleft=sort(unique(testPred[[1]]$age))[which.max(is.na(fitMat[,n]))-1],
       xright=min(sort(unique(testPred[[1]]$age))),
       ybottom=sort(unique(testPred[[1]]$ellisMag))[n] - 0.5,
       ytop=sort(unique(testPred[[1]]$ellisMag))[n] + 0.5)

  # add in a turning point as bold dashed line
  fitDiff = diff(fitMat[,n])
  segments(x0=testPred[[1]]$age[which.min(fitDiff > 0) + 1],
           x1=testPred[[1]]$age[which.min(fitDiff > 0) + 1],
           y0 = sort(unique(testPred[[1]]$ellisMag))[n] - 0.5,
           y1 = sort(unique(testPred[[1]]$ellisMag))[n] + 0.5,
           lty="31", lwd=2, lend="butt")
         
})

box()
axis(side=1, at=seq(0,11500,2000), tcl=-0.25, labels=NA)

axis(side=2, at=0:3, labels=c("0", "<1", "<10", "<100"), las=1)

if(n == 1){
mtext(side=2, line=3.25, at=par("usr")[3], text = "Population density estimate", las=0)
mtext(side=2, line=2.25, at=par("usr")[3], text = expression("(persons km"^-2*")"), las=0)
mtext(side=3, at=par("usr")[1], text="Europe", col="blue", font=2, adj=0)
}

if(n==2){
  mtext(side=3, at=par("usr")[1], text="North America", col="darkgreen", font=2, adj=0)
}

text(x=relative.axis.point(0.01, "x"), y=relative.axis.point(0.925, "y"),
     labels=paste0("(", LETTERS[n], ")"), font=2, adj=0)

close.screen(n)
return(popSamplingLims)
})

# color gradient legends
screen(5)
image(y=seq(tp1Lims[1], tp1Lims[2], len=length(contLevelsTP1)-1),
      x=c(0,1),
      z=matrix(1:(length(contLevelsTP1)-1), nrow=1),
      col=tp1Ramp, useRaster=TRUE, axes=FALSE, xlab="", ylab="")
axis(side=4, mgp=c(3,0.5,0))
box()
mtext(side=4, line=2, text="Persistence", las=0)
close.screen(5)

screen(6)
image(y=seq(tp2Lims[1], tp2Lims[2], len=length(contLevelsTP2)-1),
      x=c(0,1),
      z=matrix(1:(length(contLevelsTP2)-1), nrow=1),
      col=tp2Ramp, useRaster=TRUE, axes=FALSE, xlab="", ylab="")
axis(side=4, mgp=c(3,0.5,0))
box()
mtext(side=4, line=2.25, text="Attractor strength", las=0)
close.screen(6)

screen(7)
plot(NULL, xlim=c(11500, 0), ylim=c(0,1), axes=FALSE, xlab="", ylab="", yaxs="i", xaxs="i")
axis(side=1, at=seq(2000,11500,2000), tcl=-0.25, 
     labels=format(seq(2000,11500,2000), big.mark=","), mgp=c(3,0.2,0))
axis(side=1, at=0, labels=0, mgp=c(3,0.2,0))
mtext(side=1, line=1.25, text="Years before present")
axis(side=2, at=seq(0,1,0.25))
mtext(side=2, line=2, text="Proportion of\ntime series", las=0, lheight=0.75)

# europe pop dens
europeEllis = with(rData[rData$REGION == "Europe",], table(cut(age, breaks=seq(0,11500,250)), ellisPopNum))
europeEllis = europeEllis[rowSums(europeEllis) > 0 ,]
europeEllis = prop.table(europeEllis, 1)
rownames(europeEllis) = cutToNumeric(rownames(europeEllis))
europeEllisPlot = rbind(0, apply(europeEllis, 1, cumsum))
ageSave = as.numeric(colnames(europeEllisPlot))
europeEllisPlot = do.call("cbind", apply(europeEllisPlot, 2, function(x){cbind(x,x)}, simplify=FALSE))
colnames(europeEllisPlot) = rep(ageSave, each=2) + c(-125,125)

europeCols = colorRampPalette(c("white", "blue"))(4)

sapply(2:(nrow(europeEllisPlot)), function(n){
  print(n)
  polygon(x=c(as.numeric(colnames(europeEllisPlot)),
              rev(as.numeric(colnames(europeEllisPlot)))),
          y=c(europeEllisPlot[n-1,], rev(europeEllisPlot[n,])),
          col=europeCols[n-1])
})

# modern population
europeModCut = table(cut(rData$popDens[rData$REGION == "Europe"], breaks=c(-1,0,1,10,100,Inf)))
europeModCut = europeModCut / sum(europeModCut)
europeModCut = c(0, cumsum(europeModCut))
rect(xleft=0+350, xright=0,
     ybottom=europeModCut[-length(europeModCut)],
     ytop = europeModCut[-1],
     col=c(europeCols, "black"))

sapply(1:ncol(europeEllis), function(n){
labX = weighted.mean(as.numeric(rownames(europeEllis)), europeEllis[,n])
labY = mean(europeEllisPlot[,which.min(abs(as.numeric(rownames(europeEllis)) - labX))][n:(n+1)])
text(x=labX, y=labY, label=c("0", "<1", "<10", "<100")[n],
     col=c("black","black","white","white")[n], cex=0.75)
})
box()
text(x=relative.axis.point(0.01, "x"), y=relative.axis.point(0.925, "y"),
     adj=0, labels="(E)", font=2)
close.screen(7)

screen(8)
plot(NULL, xlim=c(11500, 0), ylim=c(0,1), axes=FALSE, xlab="", ylab="", yaxs="i", xaxs="i")
axis(side=1, at=seq(2000,11500,2000), tcl=-0.25, 
     labels=format(seq(2000,11500,2000), big.mark=","), mgp=c(3,0.2,0))
axis(side=1, at=0, labels=0, mgp=c(3,0.2,0))
mtext(side=1, line=1.25, text="Years before present")
axis(side=2, at=seq(0,1,0.25))
# europe pop dens
europeEllis = with(rData[rData$REGION == "North America",], table(cut(age, breaks=seq(0,11500,250)), ellisPopNum))
europeEllis = europeEllis[rowSums(europeEllis) > 0 ,]
europeEllis = prop.table(europeEllis, 1)
rownames(europeEllis) = cutToNumeric(rownames(europeEllis))
europeEllisPlot = rbind(0, apply(europeEllis, 1, cumsum))
ageSave = as.numeric(colnames(europeEllisPlot))
europeEllisPlot = do.call("cbind", apply(europeEllisPlot, 2, function(x){cbind(x,x)}, simplify=FALSE))
colnames(europeEllisPlot) = rep(ageSave, each=2) + c(-125,125)

europeCols = colorRampPalette(c("white", "darkgreen"))(4)

sapply(2:(nrow(europeEllisPlot)), function(n){
  print(n)
  polygon(x=c(as.numeric(colnames(europeEllisPlot)),
              rev(as.numeric(colnames(europeEllisPlot)))),
          y=c(europeEllisPlot[n-1,], rev(europeEllisPlot[n,])),
          col=europeCols[n-1])
})

# modern population
europeModCut = table(cut(rData$popDens[rData$REGION == "North America"], breaks=c(-1,0,1,10,100,Inf)))
europeModCut = europeModCut / sum(europeModCut)
europeModCut = c(0, cumsum(europeModCut))
rect(xleft=0+350, xright=0,
     ybottom=europeModCut[-length(europeModCut)],
     ytop = europeModCut[-1],
     col=c(europeCols, "black"))

sapply(1:ncol(europeEllis), function(n){
  labX = weighted.mean(as.numeric(rownames(europeEllis)), europeEllis[,n])
  labY = mean(europeEllisPlot[,which.min(abs(as.numeric(rownames(europeEllis)) - labX))][n:(n+1)])
  text(x=labX, y=labY, label=c("0", "<1", "<10", "<100")[n],
       col=c("black","black","white","white")[n], cex=0.75)
})

box()
text(x=relative.axis.point(0.01, "x"), y=relative.axis.point(0.925, "y"),
     adj=0, labels="(F)", font=2)
close.screen(8)

screen(9)
image(y=seq(0,1,len=5),
      x=c(0,1),
      z=matrix(1:10, byrow=TRUE, nrow=2),
      col=c(colorRampPalette(c("white", "blue"))(4), "black",
            colorRampPalette(c("white", "darkgreen"))(4), "black"), 
      useRaster=TRUE, axes=FALSE, xlab="", ylab="")
axis(side=4, mgp=c(3,0.5,0), at=seq(0,1,len=5),
     labels=c("0", "<1",  "<10", "<100", ">100"))
abline(v=0.5,)
box()
mtext(side=4, line=2.25, text="Population density", las=0)
close.screen(9)

close.screen(all.screens=TRUE)
dev.off()

# persistence and attractor strength correlations ####

plot(NULL, xlim=c(0,1), ylim=c(0,1))

persPredsEur = matrix(testPred[[1]]$fit, nrow=200)
persPredsEur = sapply(1:4, function(n){
  x = persPredsEur[,n]
  xAges = sort(unique(testPred[[1]]$age))
  x[xAges > popCutoffs[n,1]] = NA
  return(x)
})

persPredsAm = matrix(testPred[[2]]$fit, nrow=200)
persPredsAm = sapply(1:4, function(n){
  x = persPredsAm[,n]
  xAges = sort(unique(testPred[[1]]$age))
  x[xAges > popCutoffs[n,2]] = NA
  return(x)
})
persPredsAm = persPredsAm[,colSums(is.na(persPredsAm)) < nrow(persPredsAm)]

persPredsCor = cor(cbind(persPredsEur, persPredsAm), use="complete.obs")
dimnames(persPredsCor) = paste0(rep(0,1,10,100), )

plot(y=(persPredsEur[,1] - persPredsEur[1,1]), x=(persPredsAm[,1] - persPredsAm[1,1]), type="l", col="red",
     xlim=c(-1,1), ylim=c(-1,1))
lines(y=(persPredsEur[,1] - persPredsEur[1,1]), x=(persPredsAm[,2] - persPredsAm[1,2]), type="l", col="blue")
lines(y=(persPredsEur[,1] - persPredsEur[1,1]), x=(persPredsAm[,3] - persPredsAm[1,3]), type="l")

abline(0,1, lty='31')

# individual plots? ####

plot(NULL, xlim=c(12000,0), ylim=c(-3,3))
abline(h=0, lty="31")
sapply(split(rData, f=rData$datasetid)[1:10], function(x){
  x = x[order(x$age),]
  if(nrow(x)>10){
    lines(x$tpPC1 ~ x$age, col="grey70")
  }
})

# north america final persistence values ####

lateRdata = rData[rData$age < 3000,]

lateRdata$ageCut = as.factor(cutToNumeric(cut(lateRdata$age, breaks=seq(0,3000,250))))

latePC1 = lmer(tpPC1 ~ ellisDensRegion * ageCut + latS + (1|datasetid) + (1|depositionalenvironment),
               data=droplevels(lateRdata))
latePC2 = lmer(tpPC2 ~ ellisDensRegion * ageCut + latS + (1|datasetid) + (1|depositionalenvironment),
               data=droplevels(lateRdata))

plot(simulateResiduals(latePC1))
summary(latePC1)
performance(latePC1)

predDF = expand.grid(ellisDensRegion = levels(lateRdata$ellisDensRegion),
                     ageCut = levels(lateRdata$ageCut),
                     latS = 0)
predKeep = table(lateRdata$ellisDensRegion, lateRdata$ageCut) > 0
predDF = predDF[predKeep,]        
predDF = cbind(predDF, as.data.frame(predict(latePC1, newdata=predDF, re.form=NA, se.fit=TRUE)))

pdf("./plots/lateHoloTrends.pdf", height=4, width=5.5, useDingbats = FALSE)
par(mar=c(3,3.5,0.5,0.5), las=1, tcl=-0.25, mgp=c(3,0.5,0), ps=10)
plot(predDF$fit ~ as.numeric(as.character(predDF$ageCut)), type="n", xlim=c(3000,250), xlab="", ylab="", xaxt="n")
axis(side=1, mgp=c(3,0.2,0))
mtext(side=1, line=1.5, text="Years before present")
mtext(side=2, line=2.5, las=0, text="Assemblage persistence")

abline(h=0, lty="31")
sapply(split(predDF, f=predDF$ellisDensRegion), function(x){
  
  colVect = list(europeCols, amCols)[[ifelse(grepl("Europe", x$ellisDensRegion[1]), 1, 2)]]
  
  densCat = strsplit(as.character(x$ellisDensRegion[1]), ":")[[1]][2]
  col = colVect[which(levels(lateRdata$ellisPopDens) == densCat)]
  
  lines(x$fit ~ as.numeric(as.character(x$ageCut)), lwd=2, col=col,
        lty=ifelse(densCat=="0", "31", "solid"))
  points(x$fit ~ as.numeric(as.character(x$ageCut)), pch=21, bg=col)
  text(y=x$fit[1], x=as.numeric(as.character(x$ageCut))[1],
       pos=4, labels=paste0(ifelse(densCat=="0", densCat, paste0("<",densCat)),
                            ifelse(grepl("Europe", x$ellisDensRegion[1]), " (Eu)", " (Am)")),
       col=col)
})
dev.off()
# SUMMARY STATS ####

# what does persistence mean in raw dissimilarity terms?
rData$logitPers = beta.tr(rowMeans(rData[,c("dn1_1", "dn1_2", "dn1_3")]))

library(betareg)
tp1LogitM <- betareg(logitPers ~ tpPC1, data=rData)

# Nth Am human occupied vs unoccupied
persMean = tapply(predDF$fit[predDF$age > 7000 & predDF$age < 7500],
                  predDF$ellisDensRegion[predDF$age > 7000 & predDF$age < 7500],
                  mean)
tp1Preds = predict(tp1LogitM, newdata=data.frame(tpPC1 = persMean))
tp1Preds[6] / tp1Preds[5]

# europe start of Holo vs stabilization at 8,000 years
persMean = c(predDF$fit[predDF$ellisDensRegion == "Europe:0" & predDF$age == max(predDF$age)],
             predDF$fit[predDF$ellisDensRegion == "Europe:0" & predDF$age == predDF$age[which.min(abs(predDF$age - 8000))]])

tp1Preds = plogis(predict(tp1LogitM, newdata=data.frame(tpPC1 = persMean)))
tp1Preds[2] / tp1Preds[1]

plot(rData$tpPC1 ~ logit(rowMeans(rData[,c("dn1_1", "dn1_2", "dn1_3")])))
abline(h=-0.027, col="red")
abline(h=0.197, col="red")


