#bayesMCClust
install.packages("bayesMCClust")
library(bayesMCClust)


#rm(list=ls(all=TRUE))

# ==================================================================================
if ( TRUE ) {
  # ==================================================================================
  
  # set working directory
  oldDir <- getwd()
  curDir <- tempdir()
  setwd(curDir)
  
  if ( !file.exists("bayesMCClust-wd") ) dir.create("bayesMCClust-wd")
  setwd("bayesMCClust-wd")
  myOutfilesDir <- "mcClust-Example-Outfiles" 
  
  # load data 
  data(MCCExampleData)
  
  # ==================================================================================
  
  # function call 
  system.time(
    outList <- mcClust(    # parameter lists (every four) must be complete!
      Data=list(dataFile=MCCExampleData$Njk.i, 
                storeDir=myOutfilesDir,
                priorFile= NULL), 
      Prior=list(H=2, # sample(2:6, 1), # 4 
                 e0=4, 
                 c=1,
                 cOff=1,
                 usePriorFile=FALSE,
                 xiPooled=FALSE,
                 N0=5), 
      Initial=list(xi.start.ind=3, 
                   pers=0.7), 
      Mcmc=list(M=100, 
                M0=20, 
                mOut=5, 
                mSave=50, 
                seed=sample(1:100000, 1) # 123 
      ) 
    )
  )
  
  str(outList)
  
  #outFileName
  #results <- load(outFileName)
  #results
  #estTransProb
  
  allocList <- calcAllocationsMCC(outList, thin=1, maxi=50) # , plotPathsForEta=TRUE
  str(allocList)
  
  myTransProbs <- calcTransProbs(outList, estGroupSize=allocList$estGroupSize, thin=1, 
                                 printXtable=FALSE, printSd=FALSE, printTogether=TRUE ) 
  # , plotPaths=TRUE, grLabels=paste("Group", 1:Prior$H)
  str(myTransProbs)
  
  myTransList <- plotTransProbs(outList, estTransProb=myTransProbs$estTransProb, 
                                estGroupSize=allocList$estGroupSize, class=allocList$class, plotPooled=TRUE, 
                                plotContTable=TRUE, printContTable=TRUE, plotContPooled=TRUE) 
  # , grLabels=paste("Group", 1:Prior$H)
  str(myTransList)
  
  (equiDist <- calcEquiDist(outList, thin=1, maxi=50)) 
  #, printEquiDist=TRUE, plotEquiDist=TRUE , grLabels=paste("Group", 1:Prior$H)
  
  myLongRunDistList <- calcLongRunDist(outList, 
                                       initialStateData=MCCExampleData$initialState, 
                                       class=allocList$class, equiDist=equiDist, maxi=50) 
  # , printLongRunDist=TRUE, grLabels=paste("Group", 1:Prior$H)
  str(myLongRunDistList)
  
  myTypicalMembs <- plotTypicalMembers(outList, moreTypMemb=c(10,25,40,55,70,85,100), 
                                       myObsList=MCCExampleData$obsList, classProbs=allocList$classProbs) # noTypMemb=7
  str(myTypicalMembs)
  
  plotScatter(outList, thin=1, xi11=c(1,1), xi12=c(2,2), xi21=c(2,2), xi22=c(3,3), 
              xi31=c(1,1), xi32=c(3,3) )
  
  mySegPower <- calcSegmentationPower(outList, classProbs=allocList$classProbs, 
                                      class=allocList$class, printXtable=TRUE, calcSharp=TRUE, printSharpXtable=TRUE ) 
  # , grLabels=paste("Group", 1:Prior$H)
  str(mySegPower)
  
  myEntropy <- calcEntropy(outList, classProbs=allocList$classProbs, 
                           class=allocList$class, printXtable=TRUE ) 
  # , grLabels=paste("Group", 1:Prior$H)
  myEntropy
  
  plotLikeliPaths(outList, from=10, by=1 )
  
  myNumEffTables <- calcNumEff( outList, thin=1, printXi=TRUE, printE=TRUE, 
                                printBeta=TRUE, grLabels=paste("Group", 1:outList$Prior$H) ) 
  str(myNumEffTables)
  
  myMSCrits <- calcMSCritMCC(workDir=myOutfilesDir, myLabel="mcClust-Example", H0=4, 
                             whatToDoList=c("approxML", "approxMCL", "postMode") ) 
  str(myMSCrits)
  
  setwd(oldDir)
  
} # end if

# ==================================================================================
# ==================================================================================
# ==================================================================================

# ==================================================================================
if ( FALSE ) {
  # ==================================================================================
  
  rm(list=ls(all=TRUE))
  
  # set working directory
  oldDir <- getwd()
  curDir <- tempdir()
  setwd(curDir)
  
  if ( !file.exists("bayesMCClust-wd") ) dir.create("bayesMCClust-wd")
  setwd("bayesMCClust-wd") 
  myOutfilesDir <- "mcClustExtended-Example-Outfiles"      
  
  # load data 
  data(MCCExtExampleData)
  if (!is.element("MCCExtExampleData$covariates", search())) { 
    attach(MCCExtExampleData$covariates)
  }
  
  # ==================================================================================
  
  groupNr <- 2 # sample(2:6, 1) # 3
  
  # ==================================================================================
  
  results <- kmeans( log( MCCExtExampleData$NjkiMat + 0.5 ) , groupNr, nstart=2)
  
  # ==================================================================================
  
  require(nnet, quietly = TRUE)
  H <- groupNr
  X = cbind( intercept=1, alrateBezNew, unskilled, skilled, angStart ) 
  
  N <- dim(X)[1]
  mX <- data.frame( cbind(group=as.factor( results$cluster ), X[,-1], 
                          matrix(sample(1:H,H*N,replace=TRUE),N,H)) )
  
  colnames(mX)[6:(6+groupNr-1)] <- 
    c( "as.1", "as.2", "as.3", "as.4", "as.5", "as.6" )[1:groupNr] 
  
  tempMNom <- multinom(group ~ alrateBezNew+ unskilled+ skilled+ angStart, 
                       data=as.data.frame(mX)) 
  
  toStartBeta <- t(rbind(0,coef( tempMNom )))
  
  # ==================================================================================
  # function call 
  outList <- mcClustExtended(      
    Data=list(dataFile=MCCExtExampleData$Njk.i, # parameter lists must be complete!!!
              storeDir=myOutfilesDir,
              priorFile= NULL,
              X = cbind( intercept=1, alrateBezNew, unskilled, skilled, angStart ) ), 
    Prior=list(H=groupNr, 
               c=1,
               cOff=1,
               usePriorFile=FALSE,
               xiPooled=FALSE,
               N0=5,
               betaPrior = "informative", # N(0,1)
               betaPriorMean = 0,
               betaPriorVar = 1),
    Initial=list(xi.start.ind=3, 
                 pers=0.7,
                 S.i.start = results$cluster,
                 Beta.start = toStartBeta ), 
    Mcmc=list(M=100, 
              M0=50, 
              mOut=10, 
              mSave=50, 
              seed=sample(1:100000, 1) # 69814651 
    ) 
  )
  
  str(outList)
  
  #outFileName <- outList$workspaceFile
  #results <- load(outFileName)
  #results
  #estTransProb
  
  allocList <- calcAllocationsMCCExt(outList, thin=1, maxi=50) 
  str(allocList)
  
  myTransProbs <- calcTransProbs(outList, estGroupSize=allocList$estGroupSize, thin=1, 
                                 printXtable=FALSE, printSd=FALSE, printTogether=TRUE ) 
  # plotPaths=TRUE, grLabels=paste("Group", 1:Prior$H)
  str(myTransProbs)
  
  myTransList <- plotTransProbs(outList, estTransProb=myTransProbs$estTransProb, 
                                estGroupSize=allocList$estGroupSize, class=allocList$class, plotPooled=TRUE, 
                                plotContTable=TRUE, printContTable=TRUE, plotContPooled=TRUE) 
  # , grLabels=paste("Group", 1:Prior$H)
  str(myTransList)
  
  (equiDist <- calcEquiDist(outList, thin=1, maxi=50)) 
  # , printEquiDist=TRUE, plotEquiDist=TRUE, grLabels=paste("Group", 1:Prior$H)
  
  myRegCoeffs <- calcRegCoeffs(outList, hBase=2, thin=1) 
  #, M0=Mcmc$M0, grLabels=paste("Group", 1:Prior$H), 
  # printHPD=TRUE, plotPaths=TRUE, plotACFs=TRUE
  str(myRegCoeffs)
  
  myLongRunDistList <- calcLongRunDist(outList, initialStateData=initialState, 
                                       class=allocList$class, equiDist=equiDist, maxi=50) 
  # , printLongRunDist=TRUE
  str(myLongRunDistList)
  
  myTypicalMembs <- plotTypicalMembers(outList, myObsList=MCCExtExampleData$obsList, 
                                       classProbs=allocList$classProbs) 
  # , noTypMemb=7, moreTypMemb=c(10,25,50,100,200,500,1000)
  str(myTypicalMembs)
  
  plotScatter(outList, thin=1, xi11=c(1,1), xi12=c(2,2), xi21=c(2,2), xi22=c(3,3), 
              xi31=c(1,1), xi32=c(3,3) )
  
  mySegPower <- calcSegmentationPower(outList, classProbs=allocList$classProbs, 
                                      class=allocList$class, printXtable=TRUE, calcSharp=TRUE, printSharpXtable=TRUE ) 
  # , grLabels=paste("Group", 1:Prior$H)
  str(mySegPower)
  
  myEntropy <- calcEntropy(outList, classProbs=allocList$classProbs, 
                           class=allocList$class, printXtable=TRUE ) 
  # , grLabels=paste("Group", 1:Prior$H)
  myEntropy
  
  plotLikeliPaths(outList, from=10, by=1 )
  
  myNumEffTables <- calcNumEff( outList, thin=1, printXi=TRUE, printE=TRUE, 
                                printBeta=TRUE, grLabels=paste("Group", 1:outList$Prior$H) ) 
  str(myNumEffTables)
  
  myMSCrits <- calcMSCritMCCExt(workDir=myOutfilesDir, NN=outList$N, 
                                myLabel="mcClustExtended-Example", ISdraws=100, H0=3, 
                                whatToDoList=c("approxML", "approxMCL", "postMode" ) ) 
  str(myMSCrits)
  
  setwd(oldDir)
  
  # ==================================================================================
  
  if (is.element("MCCExtExampleData$covariates", search())) { 
    detach(MCCExtExampleData$covariates)
  }
  
  # ==================================================================================
} # end if
# ==================================================================================

# ==================================================================================
# ==================================================================================

