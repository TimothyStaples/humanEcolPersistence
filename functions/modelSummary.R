modelSummary <- function(model, savePath,
                         diagTest = FALSE, time = NULL, coords = NULL, diagIter=99, digits=3){

  modelClass = class(model)
  digMod = paste0("%.", digits, "f")
  
  if(modelClass %in% c("lmerMod", "glmerMod")){
    mFix = summary(model)$coefficients
    mFix = apply(mFix, 2, function(x){sprintf(digMod, x)})
    rownames(mFix) = names(fixef(model))
    
    # SEPARATE RANDOM SLOPES IF NEEDED?
    
    mRan = as.data.frame(VarCorr(model))
    mRan = mRan[,colSums(is.na(mRan)) == 0]
    colnames(mRan) = c("Group", "Variance", "Std.Dev")
    mRan[,-1] = sapply(mRan[,-1], function(x){sprintf(digMod, x)})
  }
  
  mDiags = modelDiagTests(model, time, coords, diagIter)
  
  mDiagTable = data.frame(Test = c("Uniformity", "Dispersion", "TempAuto", "SpatAuto"),
                         iterationN = c(NA, NA, diagIter, diagIter),
                         obsStatistic = sprintf(digMod, c(mDiags$unif$statistic,
                                                          mDiags$disp$statistic,
                                                          mean(mDiags$tauto$statistic),
                                                          mean(mDiags$sauto$observed))),
                         expStatistic = c(NA, NA, NA, sprintf(digMod, mean(mDiags$sauto$expected))),
                         p.value = sprintf(digMod, c(mDiags$unif$p.value,
                                                     mDiags$disp$p.value,
                                                     sum(mDiags$tauto$p <= 0.05) / diagIter,
                                                     sum(mDiags$sauto$p <= 0.05) / diagIter)))
  mDiagTable$iterationN[is.na(mDiagTable$iterationN)] = "-"
  mDiagTable$expStatistic[is.na(mDiagTable$expStatistic)] = "-"
  
  # structure model summary
  tabCols = 5
  mTable = rbind(c("Fixed effects", rep("", tabCols-1)),
                 c("Variable", colnames(mFix), rep("", tabCols - ncol(mFix) - 1)),
                 cbind(rownames(mFix), mFix, rep("", tabCols - ncol(mFix) - 1)),
                 rep("", tabCols),
                 
                 c("Random effects", rep("", tabCols-1)),
                 c(colnames(mRan), rep("", tabCols - ncol(mRan))),
                   cbind(as.matrix(mRan), matrix("", nrow=nrow(mRan), ncol=tabCols - ncol(mRan))),
                 
                 c("Model diagnostic tests", rep("", tabCols-1)),
                 colnames(mDiagTable),
                 as.matrix(mDiagTable))
  
  rownames(mTable) = NULL
  colnames(mTable) = NULL
  
  write.csv(mTable, savePath)
  
}
