# rm(list = ls())

library(RSiena) # or RSienaTest
source('./simulation/siena07ToConvergence.R')
source('./smaller_experiment/simulateNetworkBehavior.R')

load("./data/simulated/smaller_exp.RData")
# load("./data/results/wave1imp-20-n.RData")

Nnodes = 8 # n of cores

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

########   fr.60.2.sim.mis.20.n, alco.60.2.sim.mis.20.n


impNets.60.2.20.n.t1 <- list()
impAlco.60.2.20.n.t1 <- list()

impNets.60.2.20.n.t2 <- list()
impAlco.60.2.20.n.t2 <- list()

saom.results.20.n.t1 <- list()
saom.results.20.n.t2 <- list()

for (i in 1:S) {

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

  for (t in 1:2) {
  
      for (d in 1:D) {
        cat('theta', t, '\n')
        cat('imputation',d,'\n')
        
        # now impute wave2
        if (t == 1) {
        friendship <- sienaDependent(array(c(net1imp.t1[[d]],
                                            fr.60.2.sim.mis.20.n[,,i]),
                                          dim = c(N,N,M)))
        
        drinkingbeh <- sienaDependent(cbind(alco1imp.t1[,d],
                                            alco.60.2.sim.mis.20.n[,i]),
                                                      type = "behavior")
        } else {
        friendship <- sienaDependent(array(c(net1imp.t2[[d]],
                                            fr.60.2.sim.mis.20.n[,,i]),
                                          dim = c(N,N,M)))
        
        drinkingbeh <- sienaDependent(cbind(alco1imp.t2[,d],
                                            alco.60.2.sim.mis.20.n[,i]),
                                                      type = "behavior")
        }
        missing.20.2.inv
        
        m2 <- coCovar(missing.20.2[,i], center = FALSE)
        
        m2.inv <- coCovar(missing.20.2.inv[,i], center = FALSE)
        
        Data.w2  <- sienaDataCreate(friendship, drinkingbeh, m2, m2.inv)
        
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
                                                threshold = 0.25,nodes = Nnodes,
                                                ans0 = period1saom)
          }
        }, error = function(e) {
          print(e)
        })
        effects.twoWaves <- includeEffects(effects.twoWaves, RateX,
                                             name = "drinkingbeh",
                                             type = "rate", interaction1 = "m2",
                                             fix = TRUE,
                                             test = FALSE)
        effects.twoWaves <- setEffect(effects.twoWaves, RateX,
                                        name = "drinkingbeh",
                                        type = "rate", interaction1 = "m2",
                                        fix = TRUE,
                                        initialValue = -1000)     
        
        
        effects.twoWaves <- includeEffects(effects.twoWaves, egoX,
                                          name = "friendship",
                                          interaction1 = "m2.inv",
                                          fix = TRUE,
                                          test = FALSE)

        effects.twoWaves <- setEffect(effects.twoWaves, egoX,
                                          name = "friendship",
                                          interaction1 = "m2.inv",
                                          fix = TRUE,
                                          initialValue = thetas[t])
        
        sims <- siena07(imputation.options, data = Data.w2,
                      effects = effects.twoWaves, prevAns = period1saom,
                      returnDeps = TRUE)$sims[[10]]
            
        net2imp[[d]] <- getNet(fr.60.2.sim.mis.20.n[,,i], sims[[1]][[1]]$`1`)
        alc2imp[,d] <- sims[[1]][[2]]$`1`
        
        changed <- which(((alco.60.2.sim.mis.20.n[,i] - alc2imp[,d]) != 0))[[1]]
        alc2imp[,d][changed] <- alco.60.2.sim.mis.20.n[,i][changed]
        
      }
  if (t == 1) {
      impNets.60.2.20.n.t1[[i]] <- net2imp
      impAlco.60.2.20.n.t1[[i]] <- alc2imp
  } else {
      impNets.60.2.20.n.t2[[i]] <- net2imp
      impAlco.60.2.20.n.t2[[i]] <- alc2imp
  }
    
  ###################### completed models estimation part ######################
  source('./simulation/siena07ToConvergence.R')
  saom.results.t1 <- list()
  saom.results.t2 <- list()
  for (d in 1:D) {
    cat('estimation',d,'\n')
    if (t == 1) {
      friendship <- sienaDependent(array(c(net1imp.t1[[d]],
                                           impNets.60.2.20.n.t1[[i]][[d]]),
                                         dim = c(N,N,M)))
      drinkingbeh <- sienaDependent(cbind(alco1imp.t1[,d],
                                          impAlco.60.2.20.n.t1[[i]][,d]),
                                    type = "behavior")
    } else {
      friendship <- sienaDependent(array(c(net1imp.t2[[d]],
                                           impNets.60.2.20.n.t2[[i]][[d]]),
                                         dim = c(N,N,M)))
      drinkingbeh <- sienaDependent(cbind(alco1imp.t2[,d],
                                          impAlco.60.2.20.n.t2[[i]][,d]),
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
    
    options.imputed <- sienaAlgorithmCreate(projname = "model",seed = d+239)
    
    if (t == 1) {
      if (d == 1) {
        saom.results.t1[[d]] <- siena07ToConvergence(
                                            alg = options.imputed,
                                            dat = Data.imputed, nodes = Nnodes,
                                            eff = effects.imputed,
                                            threshold = 0.25)
      } else {
        saom.results.t1[[d]] <- siena07ToConvergence(
                                      alg = options.imputed,
                                      dat = Data.imputed, nodes = Nnodes,
                                      eff = effects.imputed,
                                      ans0 = saom.results.t1[[d - 1]],
                                      threshold = 0.25)        
      }

    } else {
      if (d == 1) {
        saom.results.t2[[d]] <- siena07ToConvergence(
                                            alg = options.imputed,
                                            dat = Data.imputed, nodes = Nnodes,
                                            eff = effects.imputed,
                                            threshold = 0.25)
      } else {
        saom.results.t2[[d]] <- siena07ToConvergence(
                                      alg = options.imputed,
                                      dat = Data.imputed, nodes = Nnodes,
                                      eff = effects.imputed,
                                      ans0 = saom.results.t2[[d - 1]],
                                      threshold = 0.25)        
      }
    }

  }
  
  saom.results.20.n.t1[[i]] <- saom.results.t1
  saom.results.20.n.t2[[i]] <- saom.results.t2

  }
  
  save(saom.results.20.n.t1, saom.results.20.n.t2, # save after each dataset imp
       file = './data/results/result-20-n.RData')

}

