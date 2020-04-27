# rm(list = ls())

library(RSiena) # or RSienaTest
source('./simulation/siena07ToConvergence.R')
source('./smaller_experiment/simulateNetworkBehavior.R')

load("./data/simulated/smaller_exp.RData")

S = 100 # number of datasets (can be set to 10 for speed)
N = 60 # number of nodes
M = 2 # number of waves

################################################################################
########                                                               #########
########           Missing data generation                             #########
########                                                               #########
################################################################################

################## Missings depend on the network ##############################

### 30% missing data ####

N_miss = 18

fr.1.mis.30.n <- fr.1
alco.1.mis.30.n <- alco[,1]

probs <- ((rowSums(fr.1)-max(rowSums(fr.1)))*-1)/N
to_remove.1 <- sample(1:N, N_miss, prob = probs, replace=F)

missing.30 <- rep(1, N)
missing.30[to_remove.1] <- 0

fr.1.mis.30.n[to_remove.1,] <- NA
alco.1.mis.30.n[to_remove.1] <- NA

fr.60.2.sim.mis.30.n <- fr.60.2.sim
alco.60.2.sim.mis.30.n <- alco.60.2.sim

for (i in 1:S) {
  probs <- ((rowSums(fr.60.2.sim[,,i])-max(rowSums(fr.60.2.sim[,,i])))*-1)/N
  to_remove.2 <- sample(1:N, N_miss, prob = probs, replace=F)
  
  fr.60.2.sim.mis.30.n[,,i][to_remove.2,] <- NA
  alco.60.2.sim.mis.30.n[,i][to_remove.2] <- NA
}


################################################################################
#########                                                             ##########
#########                   Imputing datasets                         ##########
#########                                                             ##########
################################################################################

library(mice)

getNet <- function(observedNet,edgeList) {
  # observedNet = observed network as adjacency with missing data
  # edgeList = edgelist that is returned by siena07(...)$sims
  observedNet[is.na(observedNet)] <- 0
  for (i in 1:nrow(edgeList)) {
    observedNet[edgeList[i,1],edgeList[i,2]] <- 1
  }
  return(observedNet)
}

Jaccard = function (x, y) {
  M.11 = sum(x == 1 & y == 1)
  M.10 = sum(x == 1 & y == 0)
  M.01 = sum(x == 0 & y == 1)
  return (M.11 / (M.11 + M.10 + M.01))
}


################################################################################
#########                                                             ##########
#########                   Imputing behavior with mice               ##########
#########                                                             ##########
################################################################################

########   fr.60.2.sim.mis.20.n, alco.60.2.sim.mis.20.n

D = 50

miceImpAlco.60.1.20.n <- array(rep(NA, N*D), c(N, D))
miceImpAlco.60.2.20.n <- array(rep(NA, N*D*S), c(N, D, S))

impNets.60.1.20.n <- list()
impAlco.60.1.20.n <- list()

impNets.60.2.20.n <- list()
impAlco.60.2.20.n <- list()

for (i in 1:S) {
  indegree1 <- colSums(fr.1.mis.20.n, na.rm = TRUE)
  indegree2 <- colSums(fr.60.2.sim.mis.20.n[,,i], na.rm = TRUE)
  
  avgAltA1 <- rowSums(sweep(t(fr.1.mis.20.n),
                            MARGIN = 2, alco.1.mis.20.n,'*'), na.rm = TRUE) /
    rowSums(t(fr.1.mis.20.n), na.rm = TRUE)
  
  avgAltA2 <- rowSums(sweep(t(fr.60.2.sim.mis.20.n[,,i]),
                            MARGIN = 2, alco.60.2.sim.mis.20.n[,i],'*'), na.rm = TRUE) /
    rowSums(t(fr.60.2.sim.mis.20.n[,,i]), na.rm = TRUE)
  
  
  
  avgAltA1[is.nan(avgAltA1)] <- NA
  avgAltA2[is.nan(avgAltA2)] <- NA
  
  
  miceData <- cbind(alco.1.mis.20.n, alco.60.2.sim.mis.20.n[,i],
                    indegree1, indegree2,
                    avgAltA1, avgAltA2)
  
  set.seed(11019)
  miceImp <- mice(miceData, m = D, defaultMethod = "pmm", maxit = 20)
  for (d in 1:D) {
    miceImpAlco.60.1.20.n[,d] <- complete(miceImp, d)$alco.1.mis.20.n
    
    miceImpAlco.60.2.20.n[,d,i] <- complete(miceImp, d)$V2
    
  }
  
  friendship <- sienaDependent(array(c(fr.1.mis.20.n,
                                       fr.1.mis.20.n),
                                     dim = c(N, N, M)) ,
                               allowOnly = FALSE)
  
  w2 <- coDyadCovar(fr.60.2.sim.mis.20.n[,,i]) # the 2nd wave incomplete
  # network as covariate
  a2 <- coCovar(alco.60.2.sim.mis.20.n[,i]) # the 2nd wave incomplete
  # missing data covariate
  m1 <- coCovar(missing.20, center = F)
  
  
  stationaryDataList <- list()
  
  for (d in 1:D) {
    drinkingbeh <- sienaDependent(cbind(miceImpAlco.60.1.20.n[,d],
                                        miceImpAlco.60.1.20.n[,d]), 
                                  type = "behavior", allowOnly = FALSE)
    
    stationaryDataList[[d]] <- sienaDataCreate(friendship,
                                               drinkingbeh,
                                               w2, a2, m1)
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
  
  # influence
  effects.stationary <- includeEffects(effects.stationary, name = "drinkingbeh",
                                       avAlt,
                                       interaction1 = "friendship")
  
  
  #selection
  effects.stationary <- includeEffects(effects.stationary, egoX, altX, egoXaltX,
                                       name = "friendship",
                                       interaction1 = "drinkingbeh")
  
  
  for (d in 1:D) {
    effects.stationary <- setEffect(effects.stationary, Rate, initialValue = 5,
                                    name = "friendship",fix = TRUE, 
                                    group = d,type = "rate",test = FALSE)
    
    effects.stationary <- setEffect(effects.stationary, Rate, initialValue = 3,
                                    name = "drinkingbeh",fix = TRUE,
                                    group = d,type = "rate",test = FALSE)
  }
  
  estimation.options.st <- sienaAlgorithmCreate(useStdInits = FALSE,
                                                seed = 218,
                                                n3 = 3000, maxlike = FALSE,
                                                cond = FALSE, diagonalize = 0.6,
                                                firstg = 0.02,
                                                behModelType =
                                                  c(drinkingbeh=2),
                                                lessMem = TRUE)
  source('./simulation/siena07ToConvergence.R')
  
  period0saom <- siena07ToConvergence(alg = estimation.options.st,
                                      dat = Data.stationary, cluster = TRUE,
                                      nodes = Nnodes - 1,
                                      eff = effects.stationary, threshold=0.25)
  
  
  # 3 tconv  max: 0.108 
  # Time difference of 37.08413 mins
  imputation.options <- sienaAlgorithmCreate(useStdInits = FALSE, seed = 214,
                                             cond = FALSE, maxlike = FALSE,
                                             behModelType = c(drinkingbeh = 2),
                                             nsub = 0, simOnly = TRUE, n3 = 10)
  
  effects.stationary <- includeEffects(effects.stationary, RateX,
                                       name = "drinkingbeh",
                                       type = "rate", interaction1 = "m1",
                                       fix = TRUE,
                                       test = FALSE)
  effects.stationary <- setEffect(effects.stationary, RateX,
                                  name = "drinkingbeh",
                                  type = "rate", interaction1 = "m1",
                                  fix = TRUE,
                                  initialValue = -1000)
  
  set.seed(142)
  # adding obeserved as a covariate
  
  stationaryImpDataList <- list()
  
  for (d in 1:D) {
    n1 <- fr.1.mis.20.n
    n1 <- n1 + 10
    diag(n1) <- 0
    n2 <- n1
    
    friendship <- sienaDependent(array(c(n1,n2), dim = c(N,N, 2)),
                                 allowOnly = FALSE )
    
    
    a1 <- alco.1.mis.20.n
    a1.3s <- c(1:N)[a1 == 3 & !is.na(a1)]
    a1c <- sample(a1.3s,1)
    a1change <- miceImpAlco.60.1.20.n[,d]
    a1change[a1c] <- sample(c(2,4),1)
    
    
    
    drinkingbeh <- sienaDependent(cbind(a1change,a1), type = "behavior",
                                  allowOnly = FALSE)
    
    stationaryImpDataList[[d]] <- sienaDataCreate(friendship, drinkingbeh,
                                                  w2,a2,m1)
  }
  
  
  Data.stationary.imp <- sienaGroupCreate(stationaryImpDataList)
  
  # Now we can proceed and impute the data for the fist wave:
  
  imp.ans <- siena07(imputation.options, data = Data.stationary.imp,
                     effects = effects.stationary, prevAns = period0saom,
                     returnChains = TRUE,
                     returnDeps = FALSE)
  
  sims <- imp.ans$sims[[10]]
  
  net1imp <- list()
  alc1imp <- array(rep(NA,N*D), c(N,D))
  
  for (d in 1:D) {
    net1imp[[d]] = getNet(fr.1.mis.20.n, sims[[d]][[1]][[1]])
    alc1imp[,d] = sims[[d]][[2]][[1]]
  }
  
  for (d in 1:D) {
    a <- sum(((alco.1.mis.20.n - alc1imp[,d]) != 0), na.rm = T)
    changed <- which(((alco.1.mis.20.n - alc1imp[,d]) != 0))[[1]]
    alc1imp[,d][changed] <- alco.1.mis.20.n[changed]
    print(changed)
    print(a)
  }
  
  impNets.60.1.20.n[[i]] <- net1imp
  impAlco.60.1.20.n[[i]] <- alc1imp
  
  ########################### later waves imputation ###########################
  
  alc2imp <- matrix(NA,N,D)
  net2imp <- list()
  
  set.seed(1402)
  
  estimation.options <- sienaAlgorithmCreate(useStdInits = FALSE, seed = 214,
                                             n3 = 3000, maxlike = FALSE,
                                             cond = FALSE, diagonalize = 0.6,
                                             firstg = 0.02, lessMem = TRUE,
                                             behModelType =
                                               c(drinkingbeh=2))
  
  source('./simulation/siena07ToConvergence_v2.R')
  
  for (d in 1:D) {
    
    cat('imputation',d,'\n')
    
    # now impute wave2
    
    friendship <- sienaDependent(array(c(net1imp[[d]],
                                         fr.60.2.sim.mis.20.n[,,i]),
                                       dim = c(N,N,M)))
    
    drinkingbeh <- sienaDependent(cbind(alc1imp[,d],
                                        alco.60.2.sim.mis.20.n[,i]),
                                  type = "behavior")
    
    
    Data.w2  <- sienaDataCreate(friendship, drinkingbeh)
    
    effects.twoWaves <- getEffects(Data.w2)
    
    effects.twoWaves <- includeEffects(effects.twoWaves,
                                       outActSqrt, inPopSqrt,
                                       gwespFF, gwespBB)
    
    effects.twoWaves <- includeEffects(effects.twoWaves, egoX,  altX, egoXaltX,
                                       interaction1 =  "drinkingbeh")
    
    effects.twoWaves <- includeEffects(effects.twoWaves, avAlt,
                                       name = 'drinkingbeh',
                                       interaction1 =  "friendship")
    effects.twoWaves
    
    tryCatch({
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
      net2imp[[d]] <- getNet(fr.60.2.sim.mis.20.n[,,i], sims[[1]][[1]]$`1`)
      #a1changelc2imp[,d] <- sims[[2]]
      alc2imp[,d] <- sims[[1]][[2]]$`1`
      
    }, error = function(e) {
      print(e)
    })
    
  }
  impNets.60.2.20.n[[i]] <- net2imp
  impAlco.60.2.20.n[[i]] <- alc2imp
}

save.image('mi.RData')