# install.packages('RSiena', repos='http://cran.us.r-project.org')
# conda install -c conda-forge r-mice

library(RSiena) # or RSienaTest
source('./simulation/SimulateNetworksBehavior.R')
source('./simulation/siena07ToConvergence.R')

################################################################################
#######                                                                   ######
#######                      Data preparation                             ######
#######                                                                   ######
################################################################################

load("./data/Glasgow_data/Glasgow-friendship.RData")
load("./data/Glasgow_data/Glasgow-substances.RData")
load("./data/Glasgow_data/Glasgow-demographic.RData")

sex.F <- sex.F - 1

friendship.1[friendship.1 == 2] <- 1
friendship.2[friendship.2 == 2] <- 1
friendship.3[friendship.3 == 2] <- 1

friendship.1[friendship.1 == 10] <- NA
friendship.2[friendship.2 == 10] <- NA
friendship.3[friendship.3 == 10] <- NA

not_missing <- rowSums(is.na(friendship.1) | is.na(friendship.2) | 
                      is.na(friendship.3)) != ncol(friendship.1)

friendship.1 <- friendship.1[not_missing, not_missing]
friendship.2 <- friendship.2[not_missing, not_missing]
friendship.3 <- friendship.3[not_missing, not_missing]

alcohol <- alcohol[not_missing, ]

sex.F <- sex.F[not_missing]

#################################################################################
##################   Select only boys ###########################################

friendship.1 <- friendship.1[sex.F == 0, sex.F == 0]
friendship.2 <- friendship.2[sex.F == 0, sex.F == 0]
friendship.3 <- friendship.3[sex.F == 0, sex.F == 0]

alcohol <- alcohol[sex.F == 0, ]

# library(igraph)
# 
# g1 <- graph_from_adjacency_matrix(friendship.1)
# lcc1 <- decompose.graph(g1)[[1]]
# lcc1_nodes <- V(lcc1)$name
# 
# g2 <- graph_from_adjacency_matrix(friendship.2)
# lcc2 <- decompose.graph(g2)[[1]]
# lcc2_nodes <- V(lcc2)$name
# 
# names70 <- intersect(lcc1_nodes, lcc2_nodes)
# set.seed(151) # seed that creates one connected component
# names <- sample(names70, 60, replace = F)
# 
# save(names, file = "./data/names60.Rdata")

load("./data/names60.Rdata")

fr.1 <- subset(friendship.1, rownames(friendship.1) %in% names)
fr.1 <- fr.1[, (colnames(fr.1) %in% names)]

fr.2 <- subset(friendship.2, rownames(friendship.2) %in% names)
fr.2 <- fr.2[, (colnames(fr.2) %in% names)]

alco <- subset(alcohol, rownames(alcohol) %in% names)

################################################################################
#######                                                                   ######
#######                   Model Estimation                                ######
#######                                                                   ######
################################################################################

N = 60
M = 2
Nnodes = 32

friendship <- sienaDependent(array(c(fr.1, fr.2),
                                   dim = c(N, N, M)))

drinkingbeh <- sienaDependent(alco[,1:M], type = "behavior" )

myCoEvolutionData <- sienaDataCreate(friendship,
                                     drinkingbeh)

myCoEvolutionEff <- getEffects(myCoEvolutionData)

print01Report(myCoEvolutionData, modelname = "complete_data_model_60")

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, dens, recip,
                                   transTrip, transRecTrip, cycle3,
                                   outActSqrt, inPopSqrt)

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, altX, egoX, egoXaltX,
                                   interaction1 = "drinkingbeh" )

# If we want to parse out whether there is a selection or influence (or both)
# effect for drinking behaviour,
# we need to also include sender, receiver and homophily effects
# of drinking for friendship formation:
myCoEvolutionEff <- includeEffects(myCoEvolutionEff,
                                   name = "drinkingbeh", avAlt,
                                   interaction1 = "friendship")

# Check what effects you have decided to include:

myCoEvolutionEff

# Now we have to define the algorithm settings.
# The defaults are adequate. You only have to specify the filename
# that will receive the results in text format.

myCoEvAlgorithm <- sienaAlgorithmCreate(projname = "model", seed = 300)

# Finally, estimate the model; the whole command is put in parentheses
# to have the results printed directly to the screen.

(ans60 <- siena07ToConvergence(alg = myCoEvAlgorithm, dat = myCoEvolutionData,
                             eff = myCoEvolutionEff, threshold = 0.25,
                             nodes = Nnodes))

SimulateNetworksBehavior <- function(net.w1, b1.w1, n, c1,
                                     rate, dens, rec, tt, tRt, c3, inPopSq,
                                     outActSq,
                                     altX.b1, egoX.b1, egoXaltX.b1,
                                     rate.b1, lin.b1, qu.b1, avalt.b1){
  b1.w1[is.na(b1.w1)] <- 1

  # Create initial 2-wave data to get a suitable data structure.
  # arbitrarily, this initial network has an expected average degree of 3
  X0 <- matrix(rbinom(n * n, 1, 3 / (n - 1)), n, n)
  diag(X0) <- 0
  X1 <- X0
  # but X0 and X1 should not be identical for use in sienaDependent
  X0[1, 2] <- 0
  X0[2, 1] <- 1
  X1[1, 2] <- 1
  X1[2, 1] <- 0
  XX <- array(NA, c(n, n, 2))
  XX[, , 1] <- X0
  XX[, , 2] <- X1
  # Create behavior variable; initial distribution uniform on {1, ..., c}.
  ZZ1 <- pmin(matrix(trunc(c1 * runif(n * 2)) + 1, n, 2), c1)
  # With this data structure, we now can create the data.
  X  <- sienaDependent(XX, allowOnly = FALSE)
  Z1  <- sienaDependent(ZZ1, type = "behavior", allowOnly = FALSE)
  InitData <- sienaDataCreate(X, Z1)
  InitEff0 <- getEffects(InitData)
  # sink to avoid printing to the screen
  sink("eff.txt")
  # Specify the parameters.
  # The rate parameters are first multiplied by 10,
  # which will be used only to get from the totally random network XX[,,1] = X0
  # to the network that will be the simulated first wave.
  InitEff0 <- setEffect(InitEff0, Rate, type = "rate", initialValue = 10 * rate)
  InitEff0 <- setEffect(InitEff0, density, initialValue = dens)
  InitEff0 <- setEffect(InitEff0, recip, initialValue = rec)
  InitEff0 <- setEffect(InitEff0, transTrip, initialValue = tt)
  InitEff0 <- setEffect(InitEff0, cycle3, initialValue = c3)
  
  # altX, egoX, egoXaltX
  
  InitEff0 <- setEffect(InitEff0, altX, initialValue = altX.b1,
                        interaction1 = "Z1")
  InitEff0 <- setEffect(InitEff0, egoX, initialValue = egoX.b1,
                        interaction1 = "Z1")
  InitEff0 <- setEffect(InitEff0, egoXaltX, initialValue = egoXaltX.b1,
                        interaction1 = "Z1")
  
  InitEff0 <- setEffect(InitEff0, name = "Z1", Rate, type = "rate",
                        initialValue = 10 * rate.b1)
  InitEff0 <- setEffect(InitEff0, name = "Z1", linear, initialValue = lin.b1)
  InitEff0 <- setEffect(InitEff0, name = "Z1", quad, initialValue = qu.b1)
  InitEff0 <- setEffect(InitEff0, name = "Z1", avAlt, interaction1 = "X",
                        initialValue = avalt.b1)
  
  
  ## The parameter given for n3 should be larger than sum(InitEff0$include)
  nthree <- sum(InitEff0$include)	+ 5
  InitAlg <- sienaAlgorithmCreate(projname = "Init", useStdInits = FALSE,
                                  cond = FALSE, nsub = 0, n3 = nthree,
                                  simOnly = TRUE)
  # Simulate the first wave.
  # InitSim   <- siena07(InitAlg, data = InitData, eff = InitEff0,
  #                      returnDeps = TRUE, batch = TRUE, silent = TRUE)
  # Now prepare for simulating waves 2 to M.
  # Create empty result network and behavior matrices
  Xs <- array(NA, dim = c(n, n, M))
  Z1s <- array(NA, dim = c(n, M))
  # The rate parameter values from the function call are reinstated in InitEff.
  InitEff <- InitEff0
  InitEff <- setEffect(InitEff, Rate, type = "rate", initialValue = rate)
  InitEff <- setEffect(InitEff, name = "Z1", Rate, type = "rate",
                       initialValue = rate.b1)

  sink()
  for (m in 1:M) {
    # Note that we start this loop with a previously simulated network.
    # Transform the previously simulated network
    # from edge list into adjacency matrix
    if (m == 1) {
      XXsim <- net.w1
      Z1sim <- b1.w1
      # Put simulated network and behavior into the result matrix.
      Xs[, , 1] <- XXsim
      Z1s[, 1] <- Z1sim
      
      # Put simulated network in desired places for the next simulation
      XX[, , 2] <- XX[, , 1] # used only to get the data structure
      XX[, , 1] <- XXsim
      ZZ1[, 2] <- ZZ1[, 1]
      ZZ1[, 1] <- Z1sim
    } else {
      XXsim <- matrix(0, n, n)
      nsim  <- InitAlg$n3
      XXsim[InitSim$sims[[nsim]][[1]]$X[[1]][, 1:2]]  <-
        InitSim$sims[[nsim]][[1]]$X[[1]][, 3]
      Z1sim <- InitSim$sims[[nsim]][[1]]$Z1[[1]]
      # Put simulated network and behavior into the result matrix.
      Xs[, , m] <- XXsim
      Z1s[, m] <- Z1sim
      # Put simulated network in desired places for the next simulation
      XX[, , 2] <- XX[, , 1] # used only to get the data structure
      XX[, , 1] <- XXsim
      ZZ1[, 2] <- ZZ1[, 1]
      ZZ1[, 1] <- Z1sim
    }
    
    if (m < M) {
      # The following is only to prevent the error that would occur
      # in the very unlikely event XX[,,1] == XX[,,2] or ZZ[,1] == ZZ[,2].
      if (identical(XX[, , 1], XX[, , 2])) {
        XX[1, 2, 1] <- 1 - XX[1, 2, 2]
      }
      if (identical(ZZ1[, 1], ZZ1[, 2])) {
        ZZ1[1, 1] <- ifelse((ZZ1[1, 1] == 1), 2, 1)
      }
      # Specify the two-wave network data set starting with XX[,,1].
      X <- sienaDependent(XX, allowOnly = FALSE)
      Z1 <- sienaDependent(ZZ1, type = "behavior", allowOnly = FALSE)
      # Simulate wave m+1 starting at XX[,,1] which is the previous XXsim
      InitData  <- sienaDataCreate(X, Z1)
      InitSim <- siena07(InitAlg, data = InitData, eff = InitEff,
                         returnDeps = TRUE, batch = TRUE, silent = TRUE)
    }
  }
  # Present the average degrees to facilitate tuning the outdegree parameter
  # to achieve a desired average value for the average degrees.
  closeAllConnections()
  cat("Average degrees ", round(colSums(Xs, dims = 2) / n, digits = 2), "\n")
  cat("Average behavior ", round(colSums(Z1s) / n, digits = 2), "\n")
  # Result: simulated data set; covara and covarb are vectors of length n;
  # networks is an array of dimension nxnxM;
  # behaviors is a matrix of dimension nxM
  list(networks = Xs[,,2:M], behavior = Z1s[,2])
}

net.w1 = fr.1
b1.w1 = alco[,1]
n = 60
c1 = 5
rate = 2
dens = ans60$theta[2] - 2
rec = ans60$theta[3]
tt = ans60$theta[4]
tRt = ans60$theta[5]
c3 = ans60$theta[6]
inPopSq = ans60$theta[7]
outActSq = ans60$theta[8]
altX.b1 = ans60$theta[9] 
egoX.b1 = ans60$theta[10] + 1
egoXaltX.b1 = ans60$theta[11] + 1
rate.b1 = ans60$theta[12]
lin.b1 = ans60$theta[13]
qu.b1 = ans60$theta[14]
avalt.b1 = ans60$theta[15]

S = 100

fr.60.2.sim <- array(rep(0, n*n*S), c(n, n, S))
alco.60.2.sim <- array(rep(0, n*S), c(n, S))

for (i in 1:S) { 
  SN <- SimulateNetworksBehavior(net.w1, b1.w1, n, c1,
                                 rate, dens, rec, tt, tRt, c3, inPopSq, outActSq,
                                 altX.b1, egoX.b1, egoXaltX.b1,
                                 rate.b1, lin.b1, qu.b1, avalt.b1)
  
  fr.60.2.sim[,,i] <- SN$networks
  alco.60.2.sim[,i] <- SN$behavior
  print(i)
  print(which(rowSums(fr.60.2.sim[,,i]) == 0 & colSums(fr.60.2.sim[,,i]) == 0)) 
}


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

probs <- ((rowSums(fr.1)-max(rowSums(fr.1)))*-1)/n
to_remove.1 <- sample(1:n, N_miss, prob = probs, replace=F)

missing.20 <- rep(1, N)
missing.20[to_remove.1] <- 0

fr.1.mis.20.n[to_remove.1,] <- NA
alco.1.mis.20.n[to_remove.1] <- NA

fr.60.2.sim.mis.20.n <- fr.60.2.sim
alco.60.2.sim.mis.20.n <- alco.60.2.sim



for (i in 1:S) {
  probs <- ((rowSums(fr.60.2.sim[,,i])-max(rowSums(fr.60.2.sim[,,i])))*-1)/n
  to_remove.2 <- sample(1:n, N_miss, prob = probs, replace=F)

  fr.60.2.sim.mis.20.n[,,i][to_remove.2,] <- NA
  alco.60.2.sim.mis.20.n[,i][to_remove.2] <- NA
}

### 30% missing data ####

N_miss = 18

fr.1.mis.30.n <- fr.1
alco.1.mis.30.n <- alco[,1]

probs <- ((rowSums(fr.1)-max(rowSums(fr.1)))*-1)/n
to_remove.1 <- sample(1:n, N_miss, prob = probs, replace=F)

missing.30 <- rep(1, N)
missing.30[to_remove.1] <- 0

fr.1.mis.30.n[to_remove.1,] <- NA
alco.1.mis.30.n[to_remove.1] <- NA

fr.60.2.sim.mis.30.n <- fr.60.2.sim
alco.60.2.sim.mis.30.n <- alco.60.2.sim

for (i in 1:S) {
  probs <- ((rowSums(fr.60.2.sim[,,i])-max(rowSums(fr.60.2.sim[,,i])))*-1)/n
  to_remove.2 <- sample(1:n, N_miss, prob = probs, replace=F)
  
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

for (i in 2:S) {
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
                     returnDeps = TRUE)
  
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
    
    failed.imps = list()
    
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