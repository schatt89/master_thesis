# rm(list = ls())

library(RSiena) # or RSienaTest

source('./smaller_experiment/simulateNetworkBehavior.R')

load("./data/results/wave1imp-20-b.RData")

Nnodes = 12 # n of cores

S = 100 # number of datasets (can be set to 10 for speed)
N = 60 # number of nodes
M = 2 # number of waves

D = 50 # number of imputations

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

########   fr.60.2.sim.mis.20.b, alco.60.2.sim.mis.20.b


impNets.60.2.20.b.t1 <- list()
impAlco.60.2.20.b.t1 <- list()

impNets.60.2.20.b.t2 <- list()
impAlco.60.2.20.b.t2 <- list()

saom.results.20.b.t1 <- list()
saom.results.20.b.t2 <- list()

fr.60.2.sim.mis.20.b[,,6] <- fr.60.2.sim.mis.20.b[,,70]
alco.60.2.sim.mis.20.b[,6] <- alco.60.2.sim.mis.20.b[,70]

thetas <- c(1, 2)

for (i in 6:S) {

  ########################### later waves imputation ###########################
  
  alc2imp.t1 <- matrix(NA,N,D)
  net2imp.t1 <- list()
  
  alc2imp.t2 <- matrix(NA,N,D)
  net2imp.t2 <- list()

  set.seed(1402)
  
  estimation.options <- sienaAlgorithmCreate(useStdInits = FALSE, seed = 2149,
                                             n3 = 3000, maxlike = FALSE,
                                             cond = FALSE, diagonalize = 0.6,
                                             firstg = 0.02, lessMem = TRUE,
                                             MaxDegree = c(friendship = 6),
                                             behModelType =
                                               c(drinkingbeh=2))
  

  
  saom.results.t1 <- list()
  saom.results.t2 <- list()

### here was for each theta
for (t in 1:2) { 
  for (d in 1:D) {

    cat('dataset', i, "\n")
    cat('theta',t,'\n')
    cat('imputation',d,'\n')
    
    # now impute wave2
    if (t == 1) {
    friendship <- sienaDependent(array(c(net1imp.t1[[d]],
                                        fr.60.2.sim.mis.20.b[,,i]),
                                      dim = c(N,N,M)))
    
    drinkingbeh <- sienaDependent(cbind(alc1imp.t1[,d],
                                        alco.60.2.sim.mis.20.b[,i]),
                                                  type = "behavior")
    } else {
    friendship <- sienaDependent(array(c(net1imp.t2[[d]],
                                        fr.60.2.sim.mis.20.b[,,i]),
                                      dim = c(N,N,M)))
    
    drinkingbeh <- sienaDependent(cbind(alc1imp.t2[,d],
                                        alco.60.2.sim.mis.20.b[,i]),
                                                  type = "behavior")
    }

    m2.inv <- coCovar(missing.20.2.inv[,i], center = FALSE)
    
    Data.w2  <- sienaDataCreate(friendship, drinkingbeh, m2.inv)
    # , m2 removed
    
    effects.twoWaves <- getEffects(Data.w2)
    
    effects.twoWaves <- includeEffects(effects.twoWaves,
                                      outActSqrt, inPopSqrt,
                                      gwespFF)
    
    effects.twoWaves <- includeEffects(effects.twoWaves, egoX,  altX,
                                      egoXaltX,
                                      interaction1 =  "drinkingbeh")
    
    effects.twoWaves <- includeEffects(effects.twoWaves, avAlt,
                                      name = 'drinkingbeh',
                                      interaction1 =  "friendship")
    effects.twoWaves
    
    if (d == 1) {
      source('./simulation/siena07ToConvergence.R')
      period1saom <- siena07ToConvergence(alg = estimation.options,
                                          dat = Data.w2, nodes = Nnodes,
                                          eff = effects.twoWaves,
                                          threshold = 0.25)
      
    } else {
      source('./simulation/siena07ToConvergence_v2.R')
      period1saom <- tryCatch({
        siena07ToConvergence(alg = estimation.options,
                             dat = Data.w2,
                             eff = effects.twoWaves,
                             threshold = 0.25,nodes = Nnodes,
                             ans0 = period1saom)       
      }, error = function(e) {
        period1saom
      })
    }

    effects.twoWaves <- includeEffects(effects.twoWaves, effFrom,
                                      name = "drinkingbeh",
                                      interaction1 = "m2.inv",
                                      fix = TRUE,
                                      test = FALSE)
    if (t == 1) {
      effects.twoWaves <- setEffect(effects.twoWaves, effFrom,
                                    name = "drinkingbeh",
                                    interaction1 = "m2.inv",
                                    fix = TRUE,
                                    initialValue = thetas[t])
      
      imputation.options <- sienaAlgorithmCreate(useStdInits = FALSE,
                                             seed = 214,
                                             cond = FALSE, maxlike = TRUE,
                                             behModelType = c(drinkingbeh = 2),
                                             nsub = 0, simOnly = TRUE, n3 = 10)
      
      sims <- siena07(imputation.options, data = Data.w2,
                      effects = effects.twoWaves, prevAns = period1saom,
                      returnDeps = TRUE)$sims[[10]]
      
      net2imp.t1[[d]] <- getNet(fr.60.2.sim.mis.20.b[,,i], sims[[1]])
      alc2imp.t1[,d] <- sims[[2]]
    } else {
      effects.twoWaves <- setEffect(effects.twoWaves, effFrom,
                                    name = "drinkingbeh",
                                    interaction1 = "m2.inv",
                                    fix = TRUE,
                                    initialValue = thetas[t])
      
      imputation.options <- sienaAlgorithmCreate(useStdInits = FALSE,
                                             seed = 214,
                                             cond = FALSE, maxlike = TRUE,
                                             behModelType = c(drinkingbeh = 2),
                                             nsub = 0, simOnly = TRUE, n3 = 10)
      
      sims <- siena07(imputation.options, data = Data.w2,
                      effects = effects.twoWaves, prevAns = period1saom,
                      returnDeps = TRUE)$sims[[10]]
      
      net2imp.t2[[d]] <- getNet(fr.60.2.sim.mis.20.b[,,i], sims[[1]])
      alc2imp.t2[,d] <- sims[[2]]
    }
    
  }
}

save(net1imp.t1, net2imp.t1, net1imp.t2, net2imp.t2,
alc1imp.t1, alc2imp.t1, alc1imp.t2, alc2imp.t2,
file = "./data/results/20misB_before_estimation.RData")

  ###################### completed models estimation part ######################
for (t in 1:2) {
  for (d in 1:D) {
    cat('dataset', i, "\n")
    cat('theta', t, '\n')
    cat('estimation',d,'\n')
    if (t == 1) {
      friendship <- sienaDependent(array(c(net1imp.t1[[d]],
                                           net2imp.t1[[d]]),
                                         dim = c(N,N,M)))
      drinkingbeh <- sienaDependent(cbind(alc1imp.t1[,d],
                                          alc2imp.t1[,d]),
                                    type = "behavior")
    } else {
      friendship <- sienaDependent(array(c(net1imp.t2[[d]],
                                           net2imp.t2[[d]]),
                                         dim = c(N,N,M)))
      drinkingbeh <- sienaDependent(cbind(alc1imp.t2[,d],
                                          alc2imp.t2[,d]),
                                    type = "behavior")
    }
    
    Data.imputed  <- sienaDataCreate(friendship, drinkingbeh)
    
    effects.imputed <- getEffects(Data.imputed)
    
    effects.imputed <- includeEffects(effects.imputed, gwespFF,
                                       outActSqrt, inPopSqrt)
    
    effects.imputed <- includeEffects(effects.imputed,
                                       egoX,  altX, egoXaltX,
                                       interaction1 =  "drinkingbeh")
    
    effects.imputed <- includeEffects(effects.imputed, avAlt,
                                       name = 'drinkingbeh',
                                       interaction1 =  "friendship")
    
    options.imputed <- sienaAlgorithmCreate(projname = "model", seed = d+209)
    
    if (t == 1) {
      if (d == 1) {
        source('./simulation/siena07ToConvergence.R')
        saom.results.t1[[d]] <- siena07ToConvergence(
                                            alg = options.imputed,
                                            dat = Data.imputed, nodes = Nnodes,
                                            eff = effects.imputed,
                                            threshold = 0.25)
      } else {
        source('./simulation/siena07ToConvergence_v2.R')
        saom.results.t1[[d]] <- tryCatch({siena07ToConvergence(
                                        alg = options.imputed,
                                        dat = Data.imputed, nodes = Nnodes,
                                        eff = effects.imputed,
                                        ans0 = saom.results.t1[[d - 1]],
                                        threshold = 0.25) 
        }, error = function(e) {
          tryCatch({
              source('./simulation/siena07ToConvergence_v4.R')
              siena07ToConvergence(alg = options.imputed,
                                            dat = Data.imputed, nodes = Nnodes,
                                            eff = effects.imputed,
                                            threshold = 0.25)
          }, error = function(e) {
              saom.results.t1[[d-1]]
          })
        })
       
      }

    } else {
      if (d == 1) {
        source('./simulation/siena07ToConvergence.R')
        saom.results.t2[[d]] <- siena07ToConvergence(
                                            alg = options.imputed,
                                            dat = Data.imputed, nodes = Nnodes,
                                            eff = effects.imputed,
                                            threshold = 0.25)
      } else {
        source('./simulation/siena07ToConvergence_v2.R')
        saom.results.t2[[d]] <- tryCatch({siena07ToConvergence(
                                            alg = options.imputed,
                                            dat = Data.imputed, nodes = Nnodes,
                                            eff = effects.imputed,
                                            ans0 = saom.results.t2[[d - 1]],
                                            threshold = 0.25)    
        }, error = function(e) {
          tryCatch({
              source('./simulation/siena07ToConvergence_v4.R')
              siena07ToConvergence(alg = options.imputed,
                                            dat = Data.imputed, nodes = Nnodes,
                                            eff = effects.imputed,
                                            threshold = 0.25)
          }, error = function(e) {
              saom.results.t2[[d-1]]
          })
        })
    
      }
    }

  }

  }

  thetas.t1 <- list()
  covthetas.t1 <- list()
  
  thetas.t2 <- list()
  covthetas.t2 <- list()
  
    for (d in 1:D) {
      thetas.t1[[d]] <- saom.results.t1[[d]]$theta
      covthetas.t1[[d]] <- saom.results.t1[[d]]$covtheta
      
      thetas.t2[[d]] <- saom.results.t2[[d]]$theta
      covthetas.t2[[d]] <- saom.results.t2[[d]]$covtheta
    }
  
  load('./data/results/result-20-b.RData')
  
  saom.results.20.b.t1[[i]] <- list(thetas.t1, covthetas.t1)
  saom.results.20.b.t2[[i]] <- list(thetas.t2, covthetas.t2)
  
  save(effects.imputed,
       saom.results.20.b.t1, saom.results.20.b.t2, # save after each dataset imp
       file = './data/results/result-20-b.RData')

}