ecoLOF <- function(x, k, method, ...){

  require(vegan)

  if(!is.numeric(k)){stop("k must be numeric")}
  if(k > nrow(x)){stop("k must be less than the dataset size")}

  if(class(x)[1] != "matrix"){

    x <- tryCatch(as.matrix(x), error=function(e) e)
    if(inherits(x, "error")){
      stop("Reference set must be a matrix or matrix-like object")
    }}

  if(!mode(x) %in% c("numeric", "double", "logical")){
    stop("Dissimilarity can only be calculated using binary or numeric values")
  }

  distMat <- as.matrix(noveldist(x, method=method, ...))

  distSort <- t(apply(distMat, 1, sort))

  kNeigh <- lapply(1:nrow(distMat), function(n){

    kDist <- distSort[n, k+1]

    nK <- colnames(distMat)[which(distMat[n,] <= kDist)]

    nVect <- distMat[n,nK]
    return(sort(nVect[nVect > 0]))

  })

  LRD <- sapply(kNeigh, function(y){

    kN <- length(y)

    neighKDist <- distSort[names(y), kN+1]
    temp <- rbind(neighKDist, y)
    numneigh <- ncol(temp)

    return(1/(sum(apply(temp, 2, max))/numneigh))

  })
  names(LRD) <- rownames(distMat)

  # this returns the mean neighbour reachability distance (LOF denominator),
  # focal reach distance (LOF numerator) and LOF. 
  LOF <- do.call('rbind', mapply(N=kNeigh, L=LRD, function(N, L){
    return(data.frame(meanNRD = mean(LRD[names(N)]),
                      LRD = L,
                      LOF = mean(LRD[names(N)]) / L))
    }, SIMPLIFY = FALSE))
    
  rownames(LOF) = rownames(distMat)

  return(LOF)

}