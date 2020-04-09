library(RSiena) # or RSienaTest
library("mice")
source("./simulation/siena07ToConvergence.R")

getNet <- function(observedNet,edgeList) {
  # observedNet = observed network as adjacency with missing data
  # edgeList = edgelist that is returned by siena07(...)$sims
  observedNet[is.na(observedNet)] <- 0
  for (i in 1:nrow(edgeList)) {
    observedNet[edgeList[i,1],edgeList[i,2]] <- 1
  }
  return(observedNet)
}

load("./data/simulated/Data30.RData")


################################################################################
#########                                                             ##########
#########                   Imputing behavior with mice               ##########
#########                                                             ##########
################################################################################

miceImpAlco30.10.n <- array(rep(NA, 30*50*100), c(30, 50, 100))
impNets.1.30.10.n <- list()
impAlco.1.30.10.n <- list()

for (i in 1:2) {
  indegree1 <- colSums(fr.30.1.mis.10.n, na.rm = TRUE)
  indegree2 <- colSums(fr.30.2.sim.mis.10.n[,,i], na.rm = TRUE)
  indegree3 <- colSums(fr.30.3.sim.mis.10.n[,,i], na.rm = TRUE)
  
  avgAltA1 <- rowSums(sweep(t(fr.30.1.mis.10.n),
                            MARGIN = 2, alco.30.1.mis.10.n,'*'), na.rm = TRUE) /
                            rowSums(t(fr.30.1.mis.10.n), na.rm = TRUE)
  
  avgAltA2 <- rowSums(sweep(t(fr.30.2.sim.mis.10.n[,,i]),
                    MARGIN = 2, alco.30.sim.mis.10.n[,1,i],'*'), na.rm = TRUE) /
                            rowSums(t(fr.30.2.sim.mis.10.n[,,i]), na.rm = TRUE)
  
  avgAltA3 <- rowSums(sweep(t(fr.30.3.sim.mis.10.n[,,i]),
                     MARGIN = 2, alco.30.sim.mis.10.n[,2,i],'*'),na.rm = TRUE) /
                            rowSums(t(fr.30.3.sim.mis.10.n[,,i]), na.rm = TRUE)
  
  avgAltA1[is.nan(avgAltA1)] <- NA
  avgAltA2[is.nan(avgAltA2)] <- NA
  avgAltA3[is.nan(avgAltA3)] <- NA
  
  miceData <- cbind(alco.30.1.mis.10.n, alco.30.sim.mis.10.n[,,i], indegree1,
                    indegree2, indegree3, avgAltA1, avgAltA2, avgAltA3)
  
  set.seed(11019)
  miceImp <- mice(miceData, m = 50, defaultMethod = "pmm", maxit = 20)
  for (d in 1:50) {
    miceImpAlco30.10.n[,d,i] <- complete(miceImp, d)$alco.30.1.mis.10.n
  }
  
  friendship <- sienaDependent(array(c(fr.30.1.mis.10.n, fr.30.1.mis.10.n),
                               dim = c(30, 30, 2)) ,
                               allowOnly = FALSE)
  
  w2 <- coDyadCovar(fr.30.2.sim.mis.10.n[,,i]) # the 2nd wave incomplete
  # network as covariate
  a2 <- coCovar(alco.30.sim.mis.10.n[,1,i]) # the 2nd wave incomplete
  # behavior as covariate
  
  
  stationaryDataList <- list()
  
  for (d in 1:50) {
    drinkingbeh <- sienaDependent(cbind(miceImpAlco30.10.n[,d,i],
                                        miceImpAlco30.10.n[,d,i]), 
                                  type = "behavior", allowOnly = FALSE)
    
    stationaryDataList[[d]] <- sienaDataCreate(friendship,drinkingbeh, 
                                               w2,a2)
  }
  
  Data.stationary <- sienaGroupCreate(stationaryDataList)
  
  effects.stationary <- getEffects(Data.stationary)
  effects.stationary <- includeEffects(effects.stationary, outActSqrt, inPopSqrt,
                                       gwespFF, gwespBB)
  
  # 2nd wave as covariate
  effects.stationary <- includeEffects(effects.stationary, X, name ="friendship",
                                       interaction1 = "w2")
  effects.stationary <- includeEffects(effects.stationary, effFrom, 
                                       name = "drinkingbeh", interaction1 ="a2")
  
  # influence
  effects.stationary <- includeEffects(effects.stationary, name = "drinkingbeh",
                                       avAlt, interaction1 = "friendship")
  
  #selection
  effects.stationary <- includeEffects(effects.stationary,egoX,altX,
                                       name = "friendship",
                                       interaction1 = "drinkingbeh")
  
  for (d in 1:50) {
    effects.stationary <- setEffect(effects.stationary, Rate, initialValue = 5,
                                    name = "friendship",fix = TRUE, 
                                    group = d,type = "rate",test = FALSE)
    
    effects.stationary <- setEffect(effects.stationary, Rate, initialValue = 3,
                                    name = "drinkingbeh",fix = TRUE,
                                    group = d,type = "rate",test = FALSE)
  }
  
  estimation.options.st <- sienaAlgorithmCreate(useStdInits = FALSE,
                                                seed = 214,
                                                n3 = 3000, maxlike = FALSE,
                                                cond = FALSE, diagonalize = 0.6,
                                                firstg = 0.02,
                                                behModelType = c(drinkingbeh=2),
                                                lessMem = TRUE)
  
  period0saom <- siena07ToConvergence(alg = estimation.options.st,
                                      dat = Data.stationary,nodes = 7,
                                      eff = effects.stationary, threshold=0.25)
  
  imputation.options <- sienaAlgorithmCreate(useStdInits = FALSE,
                                             seed = 214,
                                             cond = FALSE, 
                                             behModelType = c(drinkingbeh = 2),
                                             maxlike = TRUE,
                                             nsub = 0,
                                             simOnly = TRUE,
                                             n3 = 10)
  
  set.seed(142)
  # adding obeserved as a covariate
  
  stationaryImpDataList <- list()
  
  for (d in 1:50) {
    n1 <- fr.30.1.mis.10.n
    n1 <- n1 + 10
    diag(n1) <- 0
    n2 <- n1
    tieList <- c(1:(nrow(n1)**2))[c(n1 == 11)]
    tieList <- tieList[!is.na(tieList)]
    
    changedTie <- sample(tieList,1)
    
    n1[changedTie] <- 0
    n2[changedTie] <- 1
    
    friendship <- sienaDependent(array(c(n1,n2), dim = c(N,N, 2)),
                                 allowOnly = FALSE )
    
    
    a1 <- alco.30.1.mis.10.n
    a1.3s <- c(1:N)[a1 == 3 & !is.na(a1)]
    a1c <- sample(a1.3s,1)
    a1change <- miceImpAlco30.10.n[,d,i]
    a1change[a1c] <- sample(c(2,4),1)
    
    
    
    drinkingbeh <- sienaDependent(cbind(a1change,a1), type = "behavior",
                                  allowOnly = FALSE)
    
    stationaryImpDataList[[d]] <- sienaDataCreate(friendship, drinkingbeh, 
                                                  w2,a2)
  }
  
  Data.stationary.imp <- sienaGroupCreate(stationaryImpDataList)
  
  sims <- siena07(imputation.options, data = Data.stationary.imp,
                  effects = effects.stationary,
                  prevAns = period0saom,
                  returnDeps = TRUE)$sims[[10]]
  
  
  net1imp <- list()
  alc1imp <- matrix(NA,30,50)
  
  for (d in 1:50) {
    net1imp[[d]] = getNet(fr.30.1.mis.10.n, sims[[d]][[1]][[1]]) 
    alc1imp[,d] = sims[[d]][[1]][[2]]
  }
  
  impNets.1.30.10.n[[i]] = net1imp
  impAlco.1.30.10.n[[i]] = alc1imp
  
  }
