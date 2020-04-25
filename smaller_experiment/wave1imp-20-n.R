# rm(list = ls())

library(RSiena) # or RSienaTest
source('./simulation/siena07ToConvergence.R')
source('./smaller_experiment/simulateNetworkBehavior.R')

load("./data/simulated/smaller_exp.RData")

Nnodes = 8 # n of cores

S = 100 # number of datasets (can be set to 10 for speed)
N = 60 # number of nodes
M = 2 # number of waves

D = 50 # number of imputations

################################################################################
########                                                               #########
########           Missing data generation                             #########
########                                                               #########
################################################################################

################## Missings depend on the network ##############################

### 20% missing data ####

N_miss = 12

fr.1.mis.20.n <- fr.1
alco.1.mis.20.n <- alco[,1]

probs <- ((rowSums(fr.1)-max(rowSums(fr.1)))*-1)/N
to_remove.1 <- sample(1:N, N_miss, prob = probs, replace=F)

missing.20 <- rep(1, N)
missing.20[to_remove.1] <- 0

missing.20.inv <- rep(0, N)
missing.20.inv[to_remove.1] <- 1

fr.1.mis.20.n[to_remove.1,] <- NA
alco.1.mis.20.n[to_remove.1] <- NA

fr.60.2.sim.mis.20.n <- fr.60.2.sim
alco.60.2.sim.mis.20.n <- alco.60.2.sim

missing.20.2.inv <- array(rep(NA, N*S), c(N, S))
missing.20.2 <- array(rep(NA, N*S), c(N, S))

for (i in 1:S) {
  probs <- ((rowSums(fr.60.2.sim[,,i])-max(rowSums(fr.60.2.sim[,,i])))*-1)/N
  to_remove.2 <- sample(1:N, N_miss, prob = probs, replace=F)
  
  fr.60.2.sim.mis.20.n[,,i][to_remove.2,] <- NA
  alco.60.2.sim.mis.20.n[,i][to_remove.2] <- NA
  
  missing.2 <- rep(1, N)
  missing.2[to_remove.2] <- 0
  
  missing.20.2[,i] <- missing.2 
  
  missing.2 <- rep(0, N)
  missing.2[to_remove.2] <- 1
  
  missing.20.2.inv[,i] <- missing.2
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

miceImpAlco.60.1.20.n <- array(rep(NA, N*D), c(N, D))

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
}

friendship <- sienaDependent(array(c(fr.1.mis.20.n,
                                     fr.1.mis.20.n),
                                   dim = c(N, N, M)) ,
                             allowOnly = FALSE)

w2 <- coDyadCovar(fr.60.2.sim.mis.20.n[,,i]) # the 2nd wave incomplete
# network as covariate
a2 <- coCovar(alco.60.2.sim.mis.20.n[,i]) # the 2nd wave incomplete
# missing data covariates 
m1 <- coCovar(missing.20, center = F)
m1.inv <- coCovar(missing.20.inv, center = F)


stationaryDataList <- list()

for (d in 1:D) {
  drinkingbeh <- sienaDependent(cbind(miceImpAlco.60.1.20.n[,d],
                                      miceImpAlco.60.1.20.n[,d]), 
                                type = "behavior", allowOnly = FALSE)
  
  stationaryDataList[[d]] <- sienaDataCreate(friendship,
                                             drinkingbeh,
                                             w2, a2, m1, m1.inv)
}

Data.stationary <- sienaGroupCreate(stationaryDataList)

effects.stationary <- getEffects(Data.stationary)
effects.stationary <- includeEffects(effects.stationary,
                                     gwespFF,
                                     outActSqrt, inPopSqrt)

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

effects.stationary <- includeEffects(effects.stationary, egoX,
                                     name = "friendship",
                                     interaction1 = "m1.inv",
                                     fix = TRUE,
                                     test = FALSE)

thetas <- c(-2, -1)
for(t in 1:2) {
  effects.stationary <- setEffect(effects.stationary, egoX,
                                  name = "friendship",
                                  interaction1 = "m1.inv",
                                  fix = TRUE,
                                  initialValue = thetas[t])
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
                                                  w2,a2,m1,m1.inv)
  }
  
  
  Data.stationary.imp <- sienaGroupCreate(stationaryImpDataList)
  
  # Now we can proceed and impute the data for the fist wave:
  
  imp.ans <- siena07(imputation.options, data = Data.stationary.imp,
                     effects = effects.stationary, prevAns = period0saom,
                     returnChains = TRUE,
                     returnDeps = TRUE)
  
  sims <- imp.ans$sims[[10]]
  
  net1imp.t1 <- list()
  alc1imp.t1 <- array(rep(NA,N*D), c(N,D))
  
  net1imp.t2 <- list()
  alc1imp.t2 <- array(rep(NA,N*D), c(N,D))
  
  if (t == 1) {
    for (d in 1:D) {
      net1imp.t1[[d]] = getNet(fr.1.mis.20.n, sims[[d]][[1]][[1]])
      alc1imp.t1[,d] = sims[[d]][[2]][[1]]
      
      changed <- which(((alco.1.mis.20.n - alc1imp.t1[,d]) != 0))[[1]]
      alc1imp.t1[,d][changed] <- alco.1.mis.20.n[changed]
    }
  } else {
    for (d in 1:D) {
      net1imp.t2[[d]] = getNet(fr.1.mis.20.n, sims[[d]][[1]][[1]])
      alc1imp.t2[,d] = sims[[d]][[2]][[1]]
      
      changed <- which(((alco.1.mis.20.n - alc1imp.t2[,d]) != 0))[[1]]
      alc1imp.t2[,d][changed] <- alco.1.mis.20.n[changed]
    }
  }
}

save(net1imp.t1, alco1imp.t1, net1imp.t2, alco1imp.t2,
     fr.60.2.sim.mis.20.n,
     alco.60.2.sim.mis.20.n,
     missing.20.2,
     missing.20.2.inv,
     file = "./data/results/wave1imp-20-n.RData")

