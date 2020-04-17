library(RSiena) # or RSienaTest
library("mice")
source("./simulation/siena07ToConvergence.R")

Nnodes = 31
M = 2 # number of waves
N = 30 # number of nodes
D = 3 # number of imputations
S = 10 # number of dataSet

getNet <- function(observedNet,edgeList) {
  # observedNet = observed network as adjacency with missing data
  # edgeList = edgelist that is returned by siena07(...)$sims
  observedNet[is.na(observedNet)] <- 0
  for (i in 1:nrow(edgeList)) {
    observedNet[edgeList[i,1],edgeList[i,2]] <- 1
  }
  return(observedNet)
}

load("./data/simulated/Data30_2waves_v2.RData")


################################################################################
#########                                                             ##########
#########                   Imputing behavior with mice               ##########
#########                                                             ##########
################################################################################

miceImpAlco.30.1.10.n <- array(rep(NA, N*D*S), c(N, D, S))

miceImpAlco.30.2.10.n <- array(rep(NA, N*D*S), c(N, D, S))

impNets.30.1.10.n <- list()
impAlco.30.1.10.n <- list()

impNets.30.2.10.n <- list()
impAlco.30.2.10.n <- list()

for (i in 1:S) {
  indegree1 <- colSums(fr.30.1.sim.mis.10.n[,,i], na.rm = TRUE)
  indegree2 <- colSums(fr.30.2.sim.mis.10.n[,,i], na.rm = TRUE)
  
  avgAltA1 <- rowSums(sweep(t(fr.30.1.sim.mis.10.n[,,i]),
    MARGIN = 2, alco.30.1.sim.mis.10.n[,i],'*'), na.rm = TRUE) /
    rowSums(t(fr.30.1.sim.mis.10.n[,,i]), na.rm = TRUE)
  
  avgAltA2 <- rowSums(sweep(t(fr.30.2.sim.mis.10.n[,,i]),
    MARGIN = 2, alco.30.2.sim.mis.10.n[,i],'*'), na.rm = TRUE) /
    rowSums(t(fr.30.2.sim.mis.10.n[,,i]), na.rm = TRUE)

  
  
  avgAltA1[is.nan(avgAltA1)] <- NA
  avgAltA2[is.nan(avgAltA2)] <- NA
  
  
  miceData <- cbind(alco.30.1.sim.mis.10.n[,i], alco.30.2.sim.mis.10.n[,i],
                    sex.F.30,
                    indegree1, indegree2,
                    avgAltA1, avgAltA2)
  
  set.seed(11019)
  miceImp <- mice(miceData, m = D, defaultMethod = "pmm", maxit = 20)
  for (d in 1:D) {
    miceImpAlco.30.1.10.n[,d,i] <- complete(miceImp, d)$V1
    
    miceImpAlco.30.2.10.n[,d,i] <- complete(miceImp, d)$V2

  }
  
  friendship <- sienaDependent(array(c(fr.30.1.sim.mis.10.n[,,i],
                                        fr.30.1.sim.mis.10.n[,,i]),
                                     dim = c(N, N, M)) ,
                               allowOnly = FALSE)
  
  w2 <- coDyadCovar(fr.30.2.sim.mis.10.n[,,i]) # the 2nd wave incomplete
  # network as covariate
  a2 <- coCovar(alco.30.2.sim.mis.10.n[,i]) # the 2nd wave incomplete
  
  gender <- coCovar(sex.F.30)
  
  stationaryDataList <- list()
  
  for (d in 1:D) {
    drinkingbeh <- sienaDependent(cbind(miceImpAlco.30.1.10.n[,d,i],
                                        miceImpAlco.30.1.10.n[,d,i]), 
                                  type = "behavior", allowOnly = FALSE)
    
    stationaryDataList[[d]] <- sienaDataCreate(friendship,
                                               drinkingbeh,
                                               w2, a2, gender)
  }
  
  Data.stationary <- sienaGroupCreate(stationaryDataList)
  
  effects.stationary <- getEffects(Data.stationary)
  effects.stationary <- includeEffects(effects.stationary,
                                       outActSqrt, inPopSqrt,
                                      gwespFF, gwespBB)
  
  # 2nd wave as covariate
  effects.stationary <- includeEffects(effects.stationary, X,
                                       name ="friendship",interaction1 = "w2")
  
  effects.stationary <- includeEffects(effects.stationary, effFrom, 
                                       name = "drinkingbeh", interaction1 ="a2")
  
  
  effects.stationary <- includeEffects(effects.stationary, sameX, 
                                       interaction1 ="gender") 
  # influence
  effects.stationary <- includeEffects(effects.stationary, name = "drinkingbeh",
                                       avAlt,
                                       interaction1 = "friendship")
  

  #selection
  effects.stationary <- includeEffects(effects.stationary, egoX, altX, egoXaltX, 
                                       name = "friendship",
                                       interaction1 = "drinkingbeh")
  
  effects.stationary <- includeEffects(effects.stationary, effFrom,
                                       name = "drinkingbeh",
                                       interaction1 = "gender")

  
  for (d in 1:D) {
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
                                                behModelType =
                                                c(drinkingbeh=2),
                                                lessMem = TRUE)

  period0saom <- siena07ToConvergence(alg = estimation.options.st,
                                      dat = Data.stationary, cluster = TRUE,
                                      nodes = Nnodes,
                                      eff = effects.stationary, threshold=0.25)

  
  # 3 tconv  max: 0.108 
  # Time difference of 37.08413 mins
  imputation.options <- sienaAlgorithmCreate(useStdInits = FALSE,
                                             seed = 214,
                                             cond = FALSE, 
                                             behModelType =
                                               c(drinkingbeh=2),
                                             maxlike = TRUE,
                                             nsub = 0,
                                             simOnly = TRUE,
                                             n3 = 10)
  
  set.seed(142)
  # adding obeserved as a covariate
  
  stationaryImpDataList <- list()
  
  for (d in 1:D) {
    n1 <- fr.30.1.sim.mis.10.n[,,i]
    n1 <- n1 + 10
    diag(n1) <- 0
    n2 <- n1
    tieList <- c(1:(nrow(n1)**2))[c(n1 == 11)]
    tieList <- tieList[!is.na(tieList)]
    
    changedTie <- sample(tieList,1)
    
    n1[changedTie] <- 0
    n2[changedTie] <- 1
    
    friendship <- sienaDependent(array(c(n1,n2), dim = c(N, N, M)),
                                 allowOnly = FALSE )
    
    
    a1 <- alco.30.1.sim.mis.10.n[,i]
    a1.3s <- c(1:N)[a1 == 3 & !is.na(a1)]
    a1c <- sample(a1.3s,1)
    a1change <- miceImpAlco.30.1.10.n[,d,i]
    a1change[a1c] <- sample(c(4,5),1)
    
    drinkingbeh <- sienaDependent(cbind(a1change,a1), type = "behavior",
                                  allowOnly = FALSE)
    

    gender <- coCovar(sex.F.30)
    
    stationaryImpDataList[[d]] <- sienaDataCreate(friendship,
                                                  drinkingbeh,
                                                  w2,a2, gender)
  }
  
  Data.stationary.imp <- sienaGroupCreate(stationaryImpDataList)
  
  sims <- siena07(imputation.options, data = Data.stationary.imp,
                  effects = effects.stationary,
                  prevAns = period0saom,
                  returnDeps = TRUE)$sims[[10]]
  
  
  net1imp <- list()
  alc1imp <- matrix(NA,N,D)
  
  for (d in 1:D) {
    net1imp[[d]] = getNet(fr.30.1.sim.mis.10.n[,,i], sims[[d]][[1]][[1]]) 
    alc1imp[,d] = sims[[d]][[1]][[2]]
  }
  
  impNets.30.1.10.n[[i]] = net1imp
  impAlco.30.1.10.n[[i]] = alc1imp
  
  ########################### later waves imputation ###########################
  
  alc2imp <- matrix(NA,N,D)
  net2imp <- list()
  
  set.seed(1402)
  
  estimation.options <- sienaAlgorithmCreate(useStdInits = FALSE, seed = 214,
                                             n3 = 3000, maxlike = FALSE,
                                             cond = FALSE, diagonalize = 0.3,
                                             firstg = 0.02, lessMem = TRUE,
                                             behModelType =
                                             c(drinkingbeh=2))

  
  for (d in 1:D) {
    
    cat('imputation',d,'\n')
    
    # now impute wave2
    
    friendship <- sienaDependent(array(c(impNets.30.1.10.n[[i]][[d]],
                                         fr.30.2.sim.mis.10.n[,,i]),
                                       dim = c(N,N,M)))
    
    drinkingbeh <- sienaDependent(cbind(impAlco.30.1.10.n[[i]][,d],
                                        alco.30.2.sim.mis.10.n[,i]),
                                  type = "behavior")
  
    gender <- coCovar(sex.F.30)
    
    Data.w2  <- sienaDataCreate(friendship, drinkingbeh, gender)
    
    effects.twoWaves <- getEffects(Data.w2)

    effects.twoWaves <- includeEffects(effects.twoWaves, inPopSqrt, cycle3,
                                       transTrip, transRecTrip, inActSqrt)
    
    effects.twoWaves <- includeEffects(effects.twoWaves, sameX,
                                       interaction1 = "gender")  
    
    effects.twoWaves <- includeEffects(effects.twoWaves, egoX,  altX, egoXaltX,
                                       interaction1 =  "drinkingbeh")
    
    effects.twoWaves <- includeEffects(effects.twoWaves, avAlt, avSim,
                                       name = 'drinkingbeh',
                                       interaction1 =  "friendship")
    
    effects.twoWaves <- includeEffects(effects.twoWaves, effFrom,
                                        name = 'drinkingbeh',
                                        interaction1 =  "gender")
  
    
    
    if (d == 1) {
      period1saom <- siena07ToConvergence(alg = estimation.options,
                                          dat = Data.w2, nodes = Nnodes,
                                          eff = effects.twoWaves,
                                          threshold = 0.25)
    } else {
      period1saom <- siena07ToConvergence(alg = estimation.options,
                                          dat = Data.w2,
                                          eff = effects.twoWaves,
                                          threshold = 0.25, nodes = Nnodes,
                                          ans0 = period1saom)
    }
    
    sims <- siena07(imputation.options, data = Data.w2,
                    effects = effects.twoWaves, prevAns = period1saom,
                    returnDeps = TRUE)$sims[[10]]
    
    net2imp[[d]] <- getNet(impNets.30.1.10.n[[i]][[d]], sims[[1]])
    #a1changelc2imp[,d] <- sims[[2]]
    alc2imp[,d] <- sims[[2]]
    
  }
  save.image('mi.RData')
  impNets.30.2.10.n[[i]] = list(net2imp)
  impAlco.30.2.10.n[[i]] = list(alc2imp) 
}

