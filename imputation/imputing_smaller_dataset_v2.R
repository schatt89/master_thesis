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

load("./data/simulated/Data30_2waves.RData")


################################################################################
#########                                                             ##########
#########                   Imputing behavior with mice               ##########
#########                                                             ##########
################################################################################

miceImpAlco30.10.n <- array(rep(NA, 30*50*100), c(30, 50, 100))
miceImpToba30.10.n <- array(rep(NA, 30*50*100), c(30, 50, 100))

miceImpAlco2.30.10.n <- array(rep(NA, 30*50*100), c(30, 50, 100))
miceImpToba2.30.10.n <- array(rep(NA, 30*50*100), c(30, 50, 100))

impNets.1.30.10.n <- list()
impAlco.1.30.10.n <- list()
impToba.1.30.10.n <- list()
impNets.2.30.10.n <- list()
impAlco.2.30.10.n <- list()
impToba.2.30.10.n <- list()
impNets.3.30.10.n <- list()
impAlco.3.30.10.n <- list()
impToba.3.30.10.n <- list()

for (i in 1:2) {
  indegree1 <- colSums(fr.30.1.mis.10.n, na.rm = TRUE)
  indegree2 <- colSums(fr.30.2.sim.mis.10.n[,,i], na.rm = TRUE)
  
  avgAltA1 <- rowSums(sweep(t(fr.30.1.mis.10.n),
                            MARGIN = 2, alco.30.1.mis.10.n,'*'), na.rm = TRUE) /
    rowSums(t(fr.30.1.mis.10.n), na.rm = TRUE)
  
  avgAltA2 <- rowSums(sweep(t(fr.30.2.sim.mis.10.n[,,i]),
                            MARGIN = 2, alco.30.2.sim.mis.10.n[,i],'*'), na.rm = TRUE) /
    rowSums(t(fr.30.2.sim.mis.10.n[,,i]), na.rm = TRUE)
  
  avgAltT1 <- rowSums(sweep(t(fr.30.1.mis.10.n),
                            MARGIN = 2, toba.30.1.mis.10.n,'*'), na.rm = TRUE) /
    rowSums(t(fr.30.1.mis.10.n), na.rm = TRUE)
  
  avgAltT2 <- rowSums(sweep(t(fr.30.2.sim.mis.10.n[,,i]),
                            MARGIN = 2, toba.30.2.sim.mis.10.n[,i],'*'), na.rm = TRUE) /
    rowSums(t(fr.30.2.sim.mis.10.n[,,i]), na.rm = TRUE)
  
  
  avgAltA1[is.nan(avgAltA1)] <- NA
  avgAltA2[is.nan(avgAltA2)] <- NA
  
  avgAltT1[is.nan(avgAltT1)] <- NA
  avgAltT2[is.nan(avgAltT2)] <- NA
  
  miceData <- cbind(alco.30.1.mis.10.n, alco.30.2.sim.mis.10.n[,i],
                    toba.30.1.mis.10.n, toba.30.2.sim.mis.10.n[,i],
                    sex.F.30,
                    indegree1, indegree2,
                    avgAltA1, avgAltA2, avgAltT1, avgAltT2)
  
  set.seed(11019)
  miceImp <- mice(miceData, m = 50, defaultMethod = "pmm", maxit = 20)
  for (d in 1:50) {
    miceImpAlco30.10.n[,d,i] <- complete(miceImp, d)$alco.30.1.mis.10.n
    miceImpToba30.10.n[,d,i] <- complete(miceImp, d)$toba.30.1.mis.10.n
    
    miceImpAlco2.30.10.n[,d,i] <- complete(miceImp, d)$V2
    miceImpToba2.30.10.n[,d,i] <- complete(miceImp, d)$V4
  }
  
  friendship <- sienaDependent(array(c(fr.30.1.mis.10.n, fr.30.1.mis.10.n),
                                     dim = c(30, 30, 2)) ,
                               allowOnly = FALSE)
  
  w2 <- coDyadCovar(fr.30.2.sim.mis.10.n[,,i]) # the 2nd wave incomplete
  # network as covariate
  a2 <- coCovar(alco.30.2.sim.mis.10.n[,i]) # the 2nd wave incomplete
  t2 <- coCovar(toba.30.2.sim.mis.10.n[,i]) # 2nd wave incomlete
  
  gender <- coCovar(sex.F.30)
  
  stationaryDataList <- list()
  
  for (d in 1:50) {
    drinkingbeh <- sienaDependent(cbind(miceImpAlco30.10.n[,d,i],
                                        miceImpAlco30.10.n[,d,i]), 
                                  type = "behavior", allowOnly = FALSE)
    
    smokingbeh <- sienaDependent(cbind(miceImpToba30.10.n[,d,i],
                                        miceImpToba30.10.n[,d,i]), 
                                  type = "behavior", allowOnly = FALSE)
    
    stationaryDataList[[d]] <- sienaDataCreate(friendship,
                                               drinkingbeh, smokingbeh,
                                               w2, a2, t2, gender)
  }
  
  Data.stationary <- sienaGroupCreate(stationaryDataList)
  
  effects.stationary <- getEffects(Data.stationary)
  effects.stationary <- includeEffects(effects.stationary,
                                       gwespFF, gwespBB)
  
  # 2nd wave as covariate
  effects.stationary <- includeEffects(effects.stationary, X,
                                       name ="friendship",interaction1 = "w2")
  
  effects.stationary <- includeEffects(effects.stationary, effFrom, 
                                       name = "drinkingbeh", interaction1 ="a2")
  
  effects.stationary <- includeEffects(effects.stationary, effFrom, 
                                       name = "smokingbeh", interaction1 ="t2")
  
  effects.stationary <- includeEffects(effects.stationary, sameX, 
                                       interaction1 ="gender") 
  # influence
  effects.stationary <- includeEffects(effects.stationary, name = "drinkingbeh",
                                       avAlt,
                                       interaction1 = "friendship")
  
  effects.stationary <- includeEffects(effects.stationary, name = "smokingbeh",
                                       avAlt,
                                       interaction1 = "friendship")
  
  #selection
  effects.stationary <- includeEffects(effects.stationary, egoX, altX, egoXaltX, 
                                       name = "friendship",
                                       interaction1 = "drinkingbeh")
  
  effects.stationary <- includeEffects(effects.stationary, effFrom,
                                       name = "drinkingbeh",
                                       interaction1 = "gender")
  
  effects.stationary <- includeEffects(effects.stationary, effFrom,
                                       name = "drinkingbeh",
                                       interaction1 = "smokingbeh")
  
  effects.stationary <- includeEffects(effects.stationary, egoX, altX, egoXaltX, 
                                       name = "friendship",
                                       interaction1 = "smokingbeh")
  
  effects.stationary <- includeEffects(effects.stationary, effFrom,
                                       name = "smokingbeh",
                                       interaction1 = "gender")
  
  effects.stationary <- includeEffects(effects.stationary, effFrom,
                                       name = "smokingbeh",
                                       interaction1 = "drinkingbeh")
  
  for (d in 1:50) {
    effects.stationary <- setEffect(effects.stationary, Rate, initialValue = 5,
                                    name = "friendship",fix = TRUE, 
                                    group = d,type = "rate",test = FALSE)
    
    effects.stationary <- setEffect(effects.stationary, Rate, initialValue = 3,
                                    name = "drinkingbeh",fix = TRUE,
                                    group = d,type = "rate",test = FALSE)
    
    effects.stationary <- setEffect(effects.stationary, Rate, initialValue = 3,
                                    name = "smokingbeh",fix = TRUE,
                                    group = d,type = "rate",test = FALSE)
  }
  
  estimation.options.st <- sienaAlgorithmCreate(useStdInits = FALSE,
                                                seed = 214,
                                                n3 = 3000, maxlike = FALSE,
                                                cond = FALSE, diagonalize = 0.6,
                                                firstg = 0.02,
                                                behModelType =
                                                c(drinkingbeh=2, smokingbeh=2),
                                                lessMem = TRUE)
  
  period0saom <- siena07ToConvergence(alg = estimation.options.st,
                                      dat = Data.stationary, nodes = 7,
                                      eff = effects.stationary, threshold=0.25)
  
  imputation.options <- sienaAlgorithmCreate(useStdInits = FALSE,
                                             seed = 214,
                                             cond = FALSE, 
                                             behModelType =
                                               c(drinkingbeh=2, smokingbeh=2),
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
    
    friendship <- sienaDependent(array(c(n1,n2), dim = c(30,30, 2)),
                                 allowOnly = FALSE )
    
    
    a1 <- alco.30.1.mis.10.n
    a1.3s <- c(1:30)[a1 == 3 & !is.na(a1)]
    a1c <- sample(a1.3s,1)
    a1change <- miceImpAlco30.10.n[,d,i]
    a1change[a1c] <- sample(c(4,5),1)
    
    t1 <- toba.30.1.mis.10.n
    t1.3s <- c(1:30)[t1 == 3 & !is.na(t1)]
    t1c <- sample(t1.3s,1)
    t1change <- miceImpToba30.10.n[,d,i]
    t1change[a1c] <- sample(c(4,5),1)
    
    
    drinkingbeh <- sienaDependent(cbind(a1change,a1), type = "behavior",
                                  allowOnly = FALSE)
    
    smokingbeh <- sienaDependent(cbind(t1change,t1), type = "behavior",
                                  allowOnly = FALSE)
    
    gender <- coCovar(sex.F.30)
    
    stationaryImpDataList[[d]] <- sienaDataCreate(friendship,
                                                  drinkingbeh, smokingbeh,
                                                  w2,a2,t2, gender)
  }
  
  Data.stationary.imp <- sienaGroupCreate(stationaryImpDataList)
  
  sims <- siena07(imputation.options, data = Data.stationary.imp,
                  effects = effects.stationary,
                  prevAns = period0saom,
                  returnDeps = TRUE)$sims[[10]]
  
  
  net1imp <- list()
  alc1imp <- matrix(NA,30,50)
  toba1imp <- matrix(NA,30,50)
  
  for (d in 1:50) {
    net1imp[[d]] = getNet(fr.30.1.mis.10.n, sims[[d]][[1]][[1]]) 
    alc1imp[,d] = sims[[d]][[1]][[2]]
    toba1imp[,d] = sims[[d]][[1]][[3]]
  }
  
  impNets.1.30.10.n[[i]] = net1imp
  impAlco.1.30.10.n[[i]] = alc1imp
  impToba.1.30.10.n[[i]] = toba1imp
  
  ########################### later waves imputation ###########################
  
  alc2imp <- matrix(NA,30,50)
  toba2imp <- matrix(NA,30,50)
  net2imp <- list()
  
  set.seed(1402)
  
  estimation.options <- sienaAlgorithmCreate(useStdInits = FALSE, seed = 214,
                                             n3 = 3000, maxlike = FALSE,
                                             cond = FALSE, diagonalize = 0.3,
                                             firstg = 0.02, lessMem = TRUE
                                             #behModelType =
                                             # c(drinkingbeh=2, smokingbeh=2)
                                             )
  
  # estimation.options <- sienaAlgorithmCreate(useStdInits = FALSE,
  #                                            seed = 2214,
  #                                            n3 = 1000, maxlike = FALSE,
  #                                            cond = FALSE)
  
  for (d in 1:50) {
    
    cat('imputation',d,'\n')
    
    # now impute wave2
    
    friendship <- sienaDependent(array(c(impNets.1.30.10.n[[i]][[d]],
                                         fr.30.2.sim.mis.10.n[,,d]),
                                       dim = c(30,30,2)))
    
    drinkingbeh <- sienaDependent(cbind(impAlco.1.30.10.n[[i]][,d],
                                        alco.30.2.sim.mis.10.n[,i]),
                                  type = "behavior")
    
    impToba.1.30.10.n[[i]][,d][sample(1:30, 1)] <- 4 # essential
    
    smokingbeh <- sienaDependent(cbind(impToba.1.30.10.n[[i]][,d],
                                       toba.30.2.sim.mis.10.n[,i]),
                                  type = "behavior")
    gender <- coCovar(sex.F.30)
    
    Data.w2  <- sienaDataCreate(friendship, drinkingbeh, smokingbeh, gender)
    
    effects.twoWaves <- getEffects(Data.w2)
    effects.twoWaves <- includeEffects(effects.twoWaves, inPopSqrt, cycle3,
                                       transTrip, transRecTrip, inActSqrt)
    
    effects.twoWaves <- includeEffects(effects.twoWaves, sameX,
                                       interaction1 = "gender")  
    
    effects.twoWaves <- includeEffects(effects.twoWaves, egoX,  altX, egoXaltX,
                                       interaction1 =  "drinkingbeh")
    
    effects.twoWaves <- includeEffects(effects.twoWaves, avAlt,
                                       name = 'drinkingbeh',
                                       interaction1 =  "friendship")
    
    effects.twoWaves <- includeEffects(effects.twoWaves, effFrom,
                                       name = 'drinkingbeh',
                                       interaction1 =  "gender")
    effects.twoWaves <- includeEffects(effects.twoWaves, effFrom,
                                       name = 'drinkingbeh',
                                       interaction1 =  "smokingbeh")
    
    effects.twoWaves <- includeEffects(effects.twoWaves, egoX,  altX, egoXaltX,
                                       interaction1 =  "smokingbeh")
    
    effects.twoWaves <- includeEffects(effects.twoWaves, avAlt,
                                       name = 'smokingbeh',
                                       interaction1 =  "friendship")
    
    effects.twoWaves <- includeEffects(effects.twoWaves, effFrom,
                                       name = 'smokingbeh',
                                       interaction1 =  "gender")
    effects.twoWaves <- includeEffects(effects.twoWaves, effFrom,
                                       name = 'smokingbeh',
                                       interaction1 =  "drinkingbeh")
    
    
    if (d == 1) {
      period1saom <- siena07ToConvergence(alg = estimation.options,
                                          dat = Data.w2,nodes = 7,
                                          eff = effects.twoWaves,
                                          threshold = 0.25)
    } else {
      period1saom <- siena07ToConvergence(alg = estimation.options,
                                          dat = Data.w2,
                                          eff = effects.twoWaves,
                                          threshold = 0.25,nodes = 7,
                                          ans0 = period1saom)
    }
    
    sims <- siena07(imputation.options, data = Data.w2,
                    effects = effects.twoWaves, prevAns = period1saom,
                    returnDeps = TRUE)$sims[[10]]
    
    net2imp[[d]] <- getNet(impNets.1.30.10.n[[i]][[d]], sims[[1]])
    #a1changelc2imp[,d] <- sims[[2]]
    alc2imp[,d] <- sims[[2]]
    # impute wave 3
    
    friendship <- sienaDependent(array( c(net2imp[[d]],
                                          fr.30.3.sim.mis.10.n[,,d]),
                                        dim = c(30,30, 2)))
    drinkingbeh <- sienaDependent(cbind(alc2imp[,d],
                                        alco.30.sim.mis.10.n[,2,i]),
                                  type = "behavior")
    
    Data.w3  <- sienaDataCreate(friendship, drinkingbeh, smoke1, gender)
    
    if (d == 1) {
      period2saom <- siena07ToConvergence(alg = estimation.options,
                                          dat = Data.w3,nodes = 7,
                                          eff = effects.twoWaves,
                                          threshold = 0.25)
    } else {
      period2saom <- siena07ToConvergence(alg = estimation.options,
                                          dat = Data.w3,
                                          eff = effects.twoWaves,
                                          threshold = 0.25,nodes = 7,
                                          ans0 = period2saom)
    }
    
    
    sims <- siena07(imputation.options, data = Data.w3,
                    effects = effects.twoWaves,
                    prevAns = period2saom, returnDeps = TRUE)$sims[[10]]
    
    net3imp[[d]] <- getNet(fr.30.3.sim.mis.10.n[,,d], sims[[1]])
    alc3imp[,d] <- sims[[2]]
    
  }
  save.image('mi.RData')
  impNets.2.30.10.n[[i]] = list(net2imp)
  impAlco.2.30.10.n[[i]] = list(alc2imp)
  impNets.3.30.10.n[[i]] = list(net3imp)
  impAlco.3.30.10.n[[i]] = list(alc3imp)
  
}
