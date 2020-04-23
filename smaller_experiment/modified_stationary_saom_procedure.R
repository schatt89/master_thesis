################################################################################
## Multiple Imputation of Missing Data in RSiena - Network and Behavior        #
## Script prepared by Robert Krause and Anna Iashina                           #
## Date: September 2019                                                        #
## Version: 2                                                                  #
################################################################################


library("RSiena")
library("mice")
packageVersion("RSiena")

# (1) Preparatory steps in R

# Before we start, we create some missing data in the s50 data set. The s50 data
# set is included in the RSiena package. See  ?s50 .

# Select the data internal to RSiena, set the number of imputations, and set a
# variable for the number of actors.

s501miss <- s501
s502miss <- s502
s503miss <- s503
alc <- s50a

D <- 5
N <- nrow(s501)

# Arbitrarily create some missing data

set.seed(1)
wave1miss <- sample(1:50,10)
wave2miss <- sample(1:50,10)
wave3miss <- sample(1:50,10)

s501miss[wave1miss,] <- NA
s502miss[wave2miss,] <- NA
s503miss[wave3miss,] <- NA

alc[wave1miss,1] <- NA
alc[wave2miss,2] <- NA
alc[wave3miss,3] <- NA

# Covariate for missing nodes

missing1 <- rep(1, 50)
missing2 <- rep(1, 50)
missing3 <- rep(1, 50)

missing1[wave1miss] <- 0
missing2[wave2miss] <- 0
missing3[wave3miss] <- 0

# We will define the following  siena07ToConvergence  function to (hopefully)
# ensure convergence. This function will run a given SAOM until the algorithm is
# converged. A circuit breaker is applied after 20 runs, or when the divergence
# is seen as too large. You may use any other specifications of  siena07()  to
# get converged estimates. See also the manual section 6.2.
# (https://www.stats.ox.ac.uk/~snijders/siena/RSiena_Manual.pdf).

# For the code below please take care to use the number of nodes 
# (number of processes for parallel processing) that is appropriate for your machine.

siena07ToConvergence <- function(alg, dat, eff, ans0=NULL, threshold, nodes = 3,
                                 cluster = TRUE, n2startPrev = 1000, ...) {
  # parameters are:
  # alg, dat, eff: Arguments for siena07: algorithm, data, effects object.
  # ans0: previous answer, if available; used as prevAns in siena07.
  # threshold: largest satisfactory value
  #            for overall maximum convergence ratio (indicating convergence).
  # nodes: number of processes for parallel processing.
  numr <- 0
  if (is.null(ans0)) {
    ans <- siena07(alg, data = dat, effects = eff, prevAns = ans0,nbrNodes = nodes,
                   returnDeps = TRUE, useCluster = cluster, ...) # the first run
  } else {
    alg$nsub <- 1
    alg$n2start <- n2startPrev
    ans <- siena07(alg, data = dat, effects = eff, prevAns = ans0,nbrNodes = nodes,
                   returnDeps = TRUE, useCluster = cluster, ...)
  }
  repeat {
    #save(ans, file = paste("ans",numr,".RData",sep = "")) # to be safe
    numr <- numr + 1           # count number of repeated runs
    tm <- ans$tconv.max      # convergence indicator
    cat(numr,"tconv  max:", round(tm,3),"\n")       # report how far we are
    if (tm < threshold) {break}   # success
    if (tm > 10) {stop()}     # divergence without much hope
    # of returning to good parameter values
    if (numr > 100) {stop()}  # now it has lasted too long
    alg$nsub <- 1
    alg$n2start <- 1000 + numr * 1000
    alg$n3 <- 2000 + numr * 1000
    ans <- siena07(alg, data = dat,effects = eff,prevAns = ans,nbrNodes = nodes,
                   returnDeps = TRUE, useCluster = cluster, ...)
  }
  if (tm > threshold) {
    stop("Warning: convergence inadequate.\n")
  }
  ans
}

# Since version 1.2-12, Maximum Likelihood (ML) estimation by RSiena with
# returnDeps = TRUE returns an edgelist of the final network at the end of the
# phase 3 simulation. The following function,  getNet() , uses this edgelist to
# impute the data. (We need a version higher than 1.2-12 to also get coevolving
# behaviors returned.)

getNet <- function(observedNet,edgeList) {
  # observedNet = observed network as adjacency with missing data
  # edgeList = edgelist that is returned by siena07(...)$sims
  observedNet[is.na(observedNet)] <- 0
  for (i in 1:nrow(edgeList)) {
    observedNet[edgeList[i,1],edgeList[i,2]] <- 1
  }
  return(observedNet)
}



#Indegree:

indegree1 <- colSums(s501miss, na.rm = TRUE)
indegree2 <- colSums(s502miss, na.rm = TRUE)
indegree3 <- colSums(s503miss, na.rm = TRUE)


#Average in-alter behavior:

avgAltA1 <- rowSums(sweep(t(s501miss), MARGIN = 2, alc[,1],'*'),na.rm = TRUE) /
            rowSums(t(s501miss), na.rm = TRUE)
avgAltA2 <- rowSums(sweep(t(s502miss), MARGIN = 2, alc[,2],'*'),na.rm = TRUE) /
            rowSums(t(s502miss), na.rm = TRUE)
avgAltA3 <- rowSums(sweep(t(s503miss), MARGIN = 2, alc[,3],'*'),na.rm = TRUE) /
            rowSums(t(s503miss), na.rm = TRUE)

avgAltA1[is.nan(avgAltA1)] <- NA
avgAltA2[is.nan(avgAltA2)] <- NA
avgAltA3[is.nan(avgAltA3)] <- NA


# Now we combine these into one data.frame:

miceData <- cbind(alc,indegree1, indegree2, indegree3, avgAltA1, avgAltA2,
                  avgAltA3)



#We can pass the obtained data to mice(). Do not forget to use $\textsf{pmm}$...
set.seed(11019)
miceImp <- mice(miceData, m = D ,defaultMethod = "pmm", maxit = 20)

#... and check convergence:

plot(miceImp)


# The algorithm is converged, if there is no clear pattern in the plots.




## Stationary SAOM imputation of the first wave



friendship <- sienaDependent(array(c(s501miss, s501miss), dim = c(N,N, 2)) ,
allowOnly = FALSE)

w2 <- coDyadCovar(s502miss) # the 2nd wave incomplete network as covariate
a2 <- coCovar(alc[,2]) # the 2nd wave incomplete behavior as covariate
m1 <- coCovar(missing1, center = FALSE) # missing node indicator for wst wave


stationaryDataList <- list()

for (d in 1:D) {
  drinkingbeh <- sienaDependent(cbind(complete(miceImp,d)$V1,
                                      complete(miceImp,d)$V1),
                                type = "behavior", allowOnly = FALSE)

  stationaryDataList[[d]] <- sienaDataCreate(friendship,drinkingbeh,w2,a2,m1)
}

Data.stationary <- sienaGroupCreate(stationaryDataList)


# After creating the data, we can choose the effects for the stationary model.
# The suggestion is to include all effects also chosen for the longitudinal
# model. Note that creation and maintenance cannot be separated in a stationary
# model. Do not forget to include effects of the second wave.

# In the stationary model it is recommended to use the geometrically weighted
# triangle effects (gwespFF, gwespBB, gwespFB, or gwespBF), instead of the
# linear additie terms (transTrip).



effects.stationary <- getEffects(Data.stationary)
effects.stationary <- includeEffects(effects.stationary, outActSqrt, inPopSqrt,
                                      gwespFF, gwespBB)

# 2nd wave as covariate
effects.stationary <- includeEffects(effects.stationary, X, name = "friendship",
                                      interaction1 = "w2")
effects.stationary <- includeEffects(effects.stationary, effFrom,
                                      name = "drinkingbeh", interaction1 = "a2")

# influence
effects.stationary <- includeEffects(effects.stationary, name = "drinkingbeh",
                                      avAlt, interaction1 = "friendship")

#selection
effects.stationary <- includeEffects(effects.stationary,egoX,altX,
                                     name = "friendship",
                                     interaction1 = "drinkingbeh")


# An important step for estimating stationary SAOMs with RSiena is to fix the
# rate function. For the imputation rather small values should be sufficient:


for (d in 1:D) {
  effects.stationary <- setEffect(effects.stationary, Rate, initialValue = 5,
                                  name = "friendship",fix = TRUE,
                                  group = d,type = "rate",test = FALSE)

  effects.stationary <- setEffect(effects.stationary, Rate, initialValue = 3,
                                  name = "drinkingbeh",fix = TRUE,
                                  group = d,type = "rate",test = FALSE)
}


# Before we can start the estimation of the imputation model we need to set the
# estimation options. It is important to set useStdInits = FALSE and
# cond = FALSE. We further recommend to set lessMem = TRUE to reduce the memory
# load during the estimation, which is especially useful when many imputations
# are done.

# Setting firstg = 0.02 and diagonalize = 0.6 may help with obtaining
# convergence faster. This will depend on the data set.
# See the RSiena manual for more information on these parameters.

# Further, we use the behavior model type 2 (behModelType = c(drinkingbeh = 2)),
# see the manual for the different behavior model types.


estimation.options.st <- sienaAlgorithmCreate(useStdInits = FALSE, seed = 214,
                                              n3 = 3000, maxlike = FALSE,
                                              cond = FALSE, diagonalize = 0.6,
                                              behModelType = c(drinkingbeh = 2),
                                              firstg = 0.02, lessMem = TRUE)


# Now we can start the estimation of the stationary SAOM. This is likely to take
# several estimation runs until it is converged.

period0saom <- siena07ToConvergence(alg = estimation.options.st,
                                    dat = Data.stationary,nodes = 7,
                                    eff = effects.stationary, threshold = 0.25)


# After obtaining one converged model we change the RSiena algorithm to
# imputation. This is achieved by the following settings:


imputation.options <- sienaAlgorithmCreate(useStdInits = FALSE, seed = 214,
                                           cond = FALSE, maxlike = FALSE,
                                           behModelType = c(drinkingbeh = 2),
                                           nsub = 0, simOnly = TRUE, n3 = 10)


# All that is left to do now is to obtain the D imputations. However, one minor
# caveat is the problem that Maximum Likelihood (ML) simulation will only start
# if there is at least one (1) change in each dependent variable. We obtain that
# by setting one random tie in the "starting network" from 1 to 0. The ML
# algorithm will guarantee that the tie will have changed back to 1 correctly by
# the end of the simulation. We do this similarly for the behavior, here a
# random value of the mean value of the scale (in this example the scale goes
# from $1$ to $5$, with a scale mean of $3$) is randomly either increased or
# decreased by one (i.e., in this example set to $2$ or $4$).

# Further, to reduce random changes in our imputation we will set all observed
# ties (except the changed one) to structurally fixed values (1 -> 11, 0 -> 10),
# with exception of the previously changed tie. This way no redundant changes
# will be simulated on the observed values. Such an elegant option is not
# available for the behavior.

# These steps are only done for the stationary imputation and will not be
# repeated for the imputation of later waves.


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
  n1 <- s501miss
  n1 <- n1 + 10
  diag(n1) <- 0
  n2 <- n1

  friendship <- sienaDependent(array(c(n1,n2), dim = c(N,N, 2)),
                                allowOnly = FALSE )


  a1 <- alc[,1]
  a1.3s <- c(1:N)[a1 == 3 & !is.na(a1)]
  a1c <- sample(a1.3s,1)
  a1change <- complete(miceImp,d)$V1
  a1change[a1c] <- sample(c(2,4),1)



  drinkingbeh <- sienaDependent(cbind(a1change,a1), type = "behavior",
                               allowOnly = FALSE)

  stationaryImpDataList[[d]] <- sienaDataCreate(friendship, drinkingbeh,w2,a2,m1)
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
  net1imp[[d]] = getNet(s501miss, sims[[d]][[1]][[1]])
  alc1imp[,d] = sims[[d]][[2]][[1]]
}

for (d in 1:D) {
    a <- sum(((alc[,1] - alc1imp[,d]) != 0), na.rm = T)
    changed <- which(((alc[,1] - alc1imp[,d]) != 0))[[1]]
    alc1imp[,d][changed] <- alc[,1][changed]
    print(changed)
    print(a)
}




#Now we have D joint imputations of the network (stored in net1imp) and behavior
# (stored in alc1imp).



# Multiple Imputation - Imputing later waves

# We first start by preparing some lists to save the networks in:


alc2imp <- matrix(NA,N,D)
alc3imp <- matrix(NA,N,D)
net2imp <- list()
net3imp <- list()



# The procedure goes wave by wave and is straightforward for RSiena veterans,
# like you. We first define the data and effects objects as usual, using one
# imputed network for wave 1 and the observed wave 2, with the missing data, as
# the end network. Then we estimate using MoM. After we have a converged model
# we use the estimated parameters in a ML simulation to retain one imputation of
# wave 2. This imputation is passed on to the next period as the starting
# network. We define a new data object with the imputed wave 2 as the starting
# network and the observed wave 3, with missing data, as the target. We, again,
# estimate with MoM and then impute with ML and retain one imputation. This
# process is repeated D times, starting once with each of the previously
# obtained 1st wave imputations.

set.seed(1402)

estimation.options <- sienaAlgorithmCreate(useStdInits = FALSE, seed = 214,
                                           n3 = 3000, maxlike = FALSE,
                                           cond = FALSE, diagonalize = 0.3,
                                           firstg = 0.02, lessMem = TRUE,
                                           behModelType = c(drinkingbeh = 2))

for (d in 1:D) {

  cat('imputation',d,'\n')

# now impute wave2

  friendship <- sienaDependent(array(c(net1imp[[d]],s502miss), dim = c(N,N,2)))
  drinkingbeh <- sienaDependent(cbind(alc1imp[,d], alc[,2]), type = "behavior")
  Data.w2  <- sienaDataCreate(friendship, drinkingbeh)


  effects.twoWaves <- getEffects(Data.w2)
  effects.twoWaves <- includeEffects(effects.twoWaves, density, recip,
                                      outActSqrt, inPopSqrt, gwespFF, gwespBB)
  effects.twoWaves <- includeEffects(effects.twoWaves, egoX,  altX, egoXaltX,
                                      interaction1 =  "drinkingbeh")
  effects.twoWaves <- includeEffects(effects.twoWaves, avAlt,
                                     name = 'drinkingbeh',
                                     interaction1 =  "friendship")

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

  net2imp[[d]] <- getNet(s502miss, sims[[1]])
  a1changelc2imp[,d] <- sims[[2]]
# impute wave 3

  friendship <- sienaDependent(array( c(net2imp[[d]], s503miss),
                                      dim = c(N,N, 2)))
  drinkingbeh <- sienaDependent(cbind(alc2imp[,d],alc[,3]), type = "behavior")
  Data.w3  <- sienaDataCreate(friendship, drinkingbeh)

  if (d == 1) {
    period2saom <- siena07ToConvergence(alg = estimation.options,
                                        dat = Data.w3,nodes = 7,
                                        eff = effects.twoWaves,
                                        threshold = 0.25)
  } else {
    period2saom <- siena07ToConvergence(alg = estimation.options, dat = Data.w3,
                                        eff = effects.twoWaves,
                                        threshold = 0.25,nodes = 7,
                                        ans0 = period2saom)
  }


  sims <- siena07(imputation.options, data = Data.w3,effects = effects.twoWaves,
                  prevAns = period2saom, returnDeps = TRUE)$sims[[10]]

  net3imp[[d]] <- getNet(s503miss, sims[[1]])
  alc3imp[,d] <- sims[[2]]
  save.image('mi.RData')
}



# The  save.image('mi.RData')  is inside the loop to save the data after every
# completed imputation. This way we are able to restart the imputation process
# in case of any crash at the current d  imputation. If we do so, we should set
# the random number seed again  set.seed(####)  and note down the new seed and
# at which d we restarted. This way we can ensure that we will later be able to
# obtain exactly the same results again.

# Keeping the  save.image()  outside the loop is more time efficient,
# especially in larger data sets where  save.image()  will take a lot of time to
# back up your data. Personally, I prefer the security that I do not have to
# restart from scratch in case of a crash to some minutes (or a maybe an hour?)
# longer computation time.


# Estimating the analysis models and combining results

# Now that we have obtained our imputed data sets we need to run our model on
# the  complete data sets. The procedure works normally as any RSiena analysis,
# except that we are running it in a loop, once for each of the D imputed data
# sets and save the results in a list.
# Use the value of n3 in the algorithm object that, in your experience,
# is sufficient to obtain stable standard error estimates for this data set.


saomResults <- list()

for (d in 1:D) {
  cat('Imputation',d,'\n')

  friendship <- sienaDependent(array(c(net1imp[[d]], net2imp[[d]],
                                net3imp[[d]]), dim = c(N,N, 3)))
  drinkingbeh <- sienaDependent(cbind(alc1imp[,d],alc2imp[,d],alc3imp[,d]),
                                type = "behavior")

  Data <- sienaDataCreate(friendship, drinkingbeh)
  effectsData <- getEffects(Data)
  effectsData <- includeEffects(effectsData, density, recip, outActSqrt,
                                inPopSqrt, gwespFF, gwespBB)
  effectsData <- includeEffects(effectsData, egoX,  altX, egoXaltX,
                                interaction1 =  "drinkingbeh")
  effectsData <- includeEffects(effectsData, avAlt, name = 'drinkingbeh',
                                interaction1 =  "friendship")

  estimation.options <- sienaAlgorithmCreate(useStdInits = FALSE, seed = 312,
                                              n3 = 3000, maxlike = FALSE,
                                              lessMem = TRUE)
  if (d == 1) {
    saomResults[[d]] <- siena07ToConvergence(alg = estimation.options,
                                     dat = Data, eff = effectsData,
                                      threshold = 0.25,nodes = 7)
  } else {
    saomResults[[d]] <- siena07ToConvergence(alg = estimation.options,
                                              dat = Data, eff = effectsData,
                                            threshold = 0.25,nodes = 7,
                                            ans0 = saomResults[[d - 1]])
  }
  save.image('mi.RData')
}


# Now we have D RSiena results and all that is left is to combine them. We need
# to extract all parameter and standard error estimates from the models and
# combine them using Rubin's Rules.

# Rubin's rules for combining results include combining parameter estimates and
# covariances. Let $\hat{\gamma}_d$ denote the th estimate of the parameter
# $\gamma$ and $W_d= \mbox{cov}(\hat{\gamma}_d \, | \, x_d)$ the
# (within-imputation) covariance matrix of the parameters of data set $x_d$. The
# combined estimate for the parameters is the average of the estimates of the
# analyses:

# $\bar{\gamma}_D = \frac{1}{D} \sum_{d=1}^D \hat{\gamma}_d.$

# Obtaining the proper standard errors is a bit less straightforward. The
# combined estimate for the standard error needs to take into account the
# variance within and between imputations. It requires the average
# within-imputation covariance matrix $\bar{W}_D$ and the between-imputation
# covariance matrix $B_D$.The average within-imputation covariance matrix is
# given by

# $\bar{W}_D = \frac{1}{D} \sum_{d=1}^D W_d$

# and the between covariance matrix by

#$B_D = \frac{1}{D-1} \sum_{d=1}^D (\hat{\gamma}_d - \bar{\gamma}_D)
# (\hat{\gamma}_d - \bar{\gamma}_D)'.$

# The total variability for $\bar{\gamma}_D$ is estimated by

# $T_D = \hat{\mbox{cov}}(\bar{\gamma}_D) = \bar{W}_D + \left(1 + \frac{1}{D}
# \right) B_D.$

# The standard errors for the parameters are given by the square roots of the
# diagonal elements of $T_D$.

# To obtain this in R, a function for the row variance will be helpful:


rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}


# How many parameters do we have?

npar <- sum(effectsData$include)


# We start by creating a dataframe for our results and fill it with the D
# estimated parameters and standard errors.


MIResults <- as.data.frame(matrix(,npar,(2 * D)))

for (d in 1:D) {
  names(MIResults)[d * 2 - 1] <- paste("imp" , "mean", sep = as.character(d))
  names(MIResults)[d * 2] <- paste("imp" , "se", sep = as.character(d))
  MIResults[,d * 2 - 1] <- saomResults[[d]]$theta
  MIResults[,d * 2] <-  sqrt(diag(saomResults[[d]]$covtheta))
}


#Now we get the average covariance structure between the parameters.


WDMIs <- matrix(0,npar,npar)

for (d in 1:D) {
  WDMIs <- WDMIs + saomResults[[d]]$covtheta
}

WDMIs <- (1/D) * WDMIs


# Using Rubin's Rules we combine the parameters and standard errors and complete
# the procedure.


finalResults <- as.data.frame(matrix(,npar,2))
names(finalResults) <- c("combinedEstimate", "combinedSE")
rownames(finalResults) <- effectsData$effectName[effectsData$include]
finalResults$combinedEstimate <- rowMeans(MIResults[,seq(1,2*D,2)])
finalResults$combinedSE <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                                  rowVar(MIResults[,seq(1,2*D,2)]))
kable(round(finalResults, 3))


#And let us not forget to save...

save.image('mi.RData')