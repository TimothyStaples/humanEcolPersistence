modelDiagTests <- function(model, time, coords, iter){

  require(ape)
  require(lmtest)

  # great circle distance function
  gc.dist<-function(lat1, lon1, lat2, lon2){
    
    deg2rad <- function(deg) {(deg * pi) / (180)} # custom function to turn degrees to radians
    
    R<-6371 # Radius of the earth in km
    dLat<-deg2rad(lat2-lat1) # deg2rad below
    dLon<-deg2rad(lon2-lon1); 
    a<-sin(dLat/2) * sin(dLat/2) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dLon/2) * sin(dLon/2)
    
    c = 2 * atan2(as.numeric(sqrt(a)), as.numeric(sqrt(1-a)))
    d = R * c # Distance in km
    return(d)
  }
    
  modelSim <- simulateResiduals(model)
  
  unif.test = testUniformity(modelSim)
  print(ifelse(abs(unif.test$statistic) > 0.1, 
               "Uniformity not okay :(",
               "Uniformity okay :)"))
  
  disp.test = testDispersion(modelSim)
  print(ifelse(disp.test$p.value <=0.05, 
               "Dispersal not okay :(",
               "Dispersal okay :)"))
  
  time.data <- time
  
  tauto.test <- lapply(1:iter, function(n){
    
  subRes <- sapply(split(residuals(model), f=time), function(x){
    
    x[sample(1:length(x),1)]
    
  })

  names(subRes) = substr(names(subRes), 1, regexpr("\\.", names(subRes))-1)
  
  a <- dwtest(subRes[match(names(subRes), sort(unique(time)))] ~ sort(unique(time)))
  
  b <- acf(subRes[match(names(subRes), sort(unique(time)))])
  
  return(list(data.frame(statistic = a$statistic,
                    method = a$method,
                    alternative = a$alternative,
                    p = a$p.value,
                    data.name = a$data.name),
              b))
  })
  
  tauto.df <- do.call("rbind", lapply(tauto.test, function(x){x[[1]]}))
  tauto.acf <- do.call("rbind", lapply(tauto.test, function(x){x[[2]]$acf}))
  
  # sauto tests
  coordID <- paste0(coords[,1],":",coords[,2])
    sauto.test <- do.call('rbind', lapply(1:iter, function(n){
    
    subRes <- sapply(split(residuals(model), f=coordID), function(x){
      
      names(x) <- NULL
      
      x[sample(1:length(x),1)]
      
    })
    
    #names(subRes) = substr(names(subRes), 1, regexpr("\\.", names(subRes))-1)
    
    subCoords <- apply(do.call("rbind",strsplit(names(subRes), ":")),2, as.numeric)
    
    subRes = subRes[complete.cases(subCoords)]
    subCoords = subCoords[complete.cases(subCoords),]
    
    distmat <- do.call("rbind", apply(subCoords, 1, function(x){
      gc.dist(x[1],x[2], subCoords[,1], subCoords[,2])
    }, simplify=FALSE))
    
    iDist <- 1/distmat
    diag(iDist) = 0 #as per Moran.I help page
    
    a <- Moran.I(subRes, weight=iDist)
    
    return(data.frame(observed = a$observed,
                           expected = a$expected,
                           sd = a$sd,
                           p = a$p.value))
    
    }))

    return(list(unif = unif.test,
                disp = disp.test,
                tauto = tauto.df,
                sauto = sauto.test,
                tauto.acf=tauto.acf))
      
}