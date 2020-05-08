# rm(list = ls())

library(RSiena) # or RSienaTest
source('./simulation/siena07ToConvergence.R')
source('./smaller_experiment/simulateNetworkBehavior.R')

load("./data/results/wave1imp-30-mcar.RData")

Nnodes = 16 # n of cores

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

########   fr.60.2.sim.mis.30.mcar, alco.60.2.sim.mis.30.mcar


impNets.60.2.30.mcar <- list()
impAlco.60.2.30.mcar <- list()

saom.results.30.mcar <- list()

fr.60.2.sim.mis.30.mcar[,,4] <- fr.60.2.sim.mis.30.mcar[,,48]
alco.60.2.sim.mis.30.mcar[,4] <- alco.60.2.sim.mis.30.mcar[,48]

for (i in 4:S) {

  ########################### later waves imputation ###########################
  
  alc2imp <- matrix(NA,N,D)
  net2imp <- list()
  
  alc2imp <- matrix(NA,N,D)
  net2imp <- list()

  set.seed(1402)
  
  estimation.options <- sienaAlgorithmCreate(useStdInits = FALSE, seed = 49,
                                             n3 = 3000, maxlike = FALSE,
                                             cond = FALSE, diagonalize = 0.6,
                                             firstg = 0.02, lessMem = TRUE,
                                             MaxDegree = c(friendship = 6),
                                             behModelType =
                                               c(drinkingbeh=2))

  
  saom.results <- list()
  saom.results <- list()

### here was for each theta

for (d in 1:D) {

cat('dataset', i, "\n")
cat('imputation',d,'\n')

# now impute wave2

friendship <- sienaDependent(array(c(net1imp[[d]],
                                    fr.60.2.sim.mis.30.mcar[,,i]),
                                    dim = c(N,N,M)))

drinkingbeh <- sienaDependent(cbind(alc1imp[,d],
                                    alco.60.2.sim.mis.30.mcar[,i]),
                                                type = "behavior")

Data.w2  <- sienaDataCreate(friendship, drinkingbeh)
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
    period1saom <- tryCatch({siena07ToConvergence(alg = estimation.options,
                            dat = Data.w2,
                            eff = effects.twoWaves,
                            threshold = 0.25,nodes = Nnodes,
                            ans0 = period1saom)       
    }, error = function(e) {
      tryCatch({
        source('./simulation/siena07ToConvergence_v4.R')
        siena07ToConvergence(alg = estimation.options,
                                        dat = Data.w2, nodes = Nnodes,
                                        eff = effects.twoWaves,
                                        threshold = 0.25)
      }, error = function(e) {
           period1saom
      })
    })
}

      
imputation.options <- sienaAlgorithmCreate(useStdInits = FALSE,
                                        seed = 214,
                                        cond = FALSE, maxlike = TRUE,
                                        behModelType = c(drinkingbeh = 2),
                                        nsub = 0, simOnly = TRUE, n3 = 10)

sims <- siena07(imputation.options, data = Data.w2,
                effects = effects.twoWaves, prevAns = period1saom,
                returnDeps = TRUE)$sims[[10]]

net2imp[[d]] <- getNet(fr.60.2.sim.mis.30.mcar[,,i], sims[[1]])
alc2imp[,d] <- sims[[2]]
    
}


  ###################### completed models estimation part ######################
for (d in 1:D) {
    cat('dataset', i, "\n")
    cat('estimation',d,'\n')

    friendship <- sienaDependent(array(c(net1imp[[d]],
                                        net2imp[[d]]),
                                        dim = c(N,N,M)))
    drinkingbeh <- sienaDependent(cbind(alc1imp[,d],
                                        alc2imp[,d]),
                                type = "behavior")
    
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
    
    options.imputed <- sienaAlgorithmCreate(projname = "model", seed = d+309)
    

    if (d == 1) {
    source('./simulation/siena07ToConvergence.R')
    saom.results[[d]] <- siena07ToConvergence(
                                        alg = options.imputed,
                                        dat = Data.imputed, nodes = Nnodes,
                                        eff = effects.imputed,
                                        threshold = 0.25)
    } else {
    source('./simulation/siena07ToConvergence_v2.R')
    saom.results[[d]] <- tryCatch({siena07ToConvergence(
                                    alg = options.imputed,
                                    dat = Data.imputed, nodes = Nnodes,
                                    eff = effects.imputed,
                                    ans0 = saom.results[[d - 1]],
                                    threshold = 0.25) 
    }, error = function(e) {
      tryCatch({
        source('./simulation/siena07ToConvergence_v4.R')
        siena07ToConvergence(alg = options.imputed,
                                        dat = Data.imputed, nodes = Nnodes,
                                        eff = effects.imputed,
                                        threshold = 0.25)
      }, error = function(e) {
        saom.results[[d-1]]
      })
    })
    }



}


thetas <- list()
covthetas <- list()

  
for (d in 1:D) {
    thetas[[d]] <- saom.results[[d]]$theta
    covthetas[[d]] <- saom.results[[d]]$covtheta
}
   
  load('./data/results/result-30-mcar.RData')
  
  saom.results.30.mcar[[i]] <- list(thetas, covthetas)
  
  save(effects.imputed,
       saom.results.30.mcar, # save after each dataset imp
       file = './data/results/result-30-mcar.RData')

}
