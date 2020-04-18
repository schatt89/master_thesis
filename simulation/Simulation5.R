library(RSiena) # or RSienaTest
source("./simulation/SimulateNetworksBehavior_v5.R")
source("./simulation/siena07ToConvergence.R")

Nnodes = 31
M = 2 # number of waves
D = 50 # number of imputations
S = 100 # number of dataSet

################################################################################
#######                                                                   ######
#######                      Data preparation                             ######
#######                                                                   ######
################################################################################

load("./data/Glasgow_data/Glasgow-friendship.RData")
load("./data/Glasgow_data/Glasgow-substances.RData")
load("./data/Glasgow_data/Glasgow-demographic.RData")

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

alcohol <- alcohol[not_missing, ]
tobacco <- tobacco[not_missing, ]
age <- age[not_missing]
sex.F <- sex.F[not_missing] - 1 # to 0s and 1s 

sex.F.30 <- sex.F[1:30]


################################################################################
#######                                                                   ######
#######             Full data model estimation                            ######
#######                                                                   ######
################################################################################

friendship <- sienaDependent(array(c(friendship.1, friendship.2),
                                   dim = c(129, 129, M)))

drinkingbeh <- sienaDependent(alcohol[,1:2], type = "behavior" )

gender <- coCovar(sex.F)

myCoEvolutionData <- sienaDataCreate(friendship, gender,
                                     drinkingbeh)

myCoEvolutionEff <- getEffects(myCoEvolutionData)

print01Report(myCoEvolutionData, modelname = "complete_data_model_30")

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, dens, recip,
                                   transTrip, transRecTrip, cycle3,
                                   outActSqrt, inPopSqrt)

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, sameX,
                                   interaction1 = "gender" )

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, altX, egoX, egoXaltX,
                                   interaction1 = "drinkingbeh" )

# If we want to parse out whether there is a selection or influence (or both)
# effect for drinking behaviour,
# we need to also include sender, receiver and homophily effects
# of drinking for friendship formation:
myCoEvolutionEff <- includeEffects(myCoEvolutionEff,
                                   name = "drinkingbeh", avAlt,
                                   interaction1 = "friendship")

# myCoEvolutionEff <- includeEffects(myCoEvolutionEff, effFrom,
                                  #  name = "drinkingbeh",
                                  #  interaction1 = "smokingbeh")

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, effFrom,
                                   name = "drinkingbeh",
                                   interaction1 = "gender")

# Check what effects you have decided to include:

myCoEvolutionEff

# Now we have to define the algorithm settings.
# The defaults are adequate. You only have to specify the filename
# that will receive the results in text format.

myCoEvAlgorithm <- sienaAlgorithmCreate(projname = "model", seed = 300)

# Finally, estimate the model; the whole command is put in parentheses
# to have the results printed directly to the screen.

(ans <- siena07ToConvergence(alg = myCoEvAlgorithm, dat = myCoEvolutionData,
                               eff = myCoEvolutionEff, threshold = 0.25,
                                nodes = Nnodes))

################################################################################
#######                                                                   ######
#######             Smaller dataset simulation                            ######
#######                                                                   ######
################################################################################


covar = sex.F.30
m = 2
n = 30
c1 = 5

rate = ans$theta[1] - 7.63
dens = ans$theta[2] - 1.1
rec = ans$theta[3] - 0.37
tt = ans$theta[4] - 0.49
tRt = ans$theta[5] 
c3 = ans$theta[6] - 0.5
inPopSq = ans$theta[7]
outActSq = ans$theta[8]

Vsame = ans$theta[9]

altX.b1 = ans$theta[10]
egoX.b1 = ans$theta[11]
egoXaltX.b1 = ans$theta[12]

rate.b1 = ans$theta[13] - 1
lin.b1 = ans$theta[14] - 0.4
qu.b1 = ans$theta[15]
avalt.b1 = ans$theta[16] - 0.57
effF.b1.V = ans$theta[17]


fr.30.1.sim <- array(rep(0, n*n*S), c(n, n, S))
fr.30.2.sim <- array(rep(0, n*n*S), c(n, n, S))

alco.30.1.sim <- array(rep(0, n*S), c(n, S))
alco.30.2.sim <- array(rep(0, n*S), c(n, S))


for (i in 1:S) { 
  SN <- SimulateNetworksBehavior(covar, M = m, n, c1,
                       rate, dens, rec, tt, tRt, c3, inPopSq, outActSq,
                       Vsame,
                       altX.b1, egoX.b1, egoXaltX.b1,
                       rate.b1, lin.b1, qu.b1, avalt.b1, effF.b1.V)

  fr.30.1.sim[,,i] <- SN$networks[,,1]
  fr.30.2.sim[,,i] <- SN$networks[,,2]

  alco.30.1.sim[,i] <- SN$behavior1[,1]
  alco.30.2.sim[,i] <- SN$behavior1[,2]

}

################################################################################
#######                                                                   ######
#######             Bigger dataset simulation                             ######
#######                                                                   ######
################################################################################

sex.F.60 <- sex.F[sample(1:129, 60, replace=F)]

covar = sex.F.60
m = 2
n = 60
c1 = 5

rate = ans$theta[1] - 7.63
dens = ans$theta[2] - 1.25

fr.60.1.sim <- array(rep(0, n*n*S), c(n, n, S))
fr.60.2.sim <- array(rep(0, n*n*S), c(n, n, S))

alco.60.1.sim <- array(rep(0, n*S), c(n, S))
alco.60.2.sim <- array(rep(0, n*S), c(n, S))


for (i in 1:S) {
  SN <- SimulateNetworksBehavior(covar, M = m, n, c1,
                       rate, dens, rec, tt, tRt, c3, inPopSq, outActSq,
                       Vsame,
                       altX.b1, egoX.b1, egoXaltX.b1,
                       rate.b1, lin.b1, qu.b1, avalt.b1, effF.b1.V)
  
  fr.60.1.sim[,,i] <- SN$networks[,,1]
  fr.60.2.sim[,,i] <- SN$networks[,,2]

  alco.60.1.sim[,i] <- SN$behavior1[,1]
  alco.60.2.sim[,i] <- SN$behavior1[,2]
  
}

################################################################################
########                                                               #########
########           Missing data generation - smaller network           #########
########                                                               #########
################################################################################

################## Missings depend on the network ##############################

### 10% missing data ####

fr.30.1.sim.mis.10.n <- fr.30.1.sim
fr.30.2.sim.mis.10.n <- fr.30.2.sim

alco.30.1.sim.mis.10.n <- alco.30.1.sim

alco.30.2.sim.mis.10.n <- alco.30.2.sim

for (i in 1:S) {
  to_remove.1 <- sample(which(rowSums(fr.30.1.sim[,,i]) < 5), 3, replace=F)
  to_remove.2 <- sample(which(rowSums(fr.30.2.sim[,,i]) < 5), 3, replace=F)
  for (j in 1:3) {
    fr.30.1.sim.mis.10.n[,,i][to_remove.1[j],] <- NA
    alco.30.1.sim.mis.10.n[,i][to_remove.1[j]] <- NA

    fr.30.2.sim.mis.10.n[,,i][to_remove.2[j],] <- NA
    alco.30.2.sim.mis.10.n[,i][to_remove.2[j]] <- NA
  }
}

### 20% missing data ####

fr.30.1.sim.mis.20.n <- fr.30.1.sim
fr.30.2.sim.mis.20.n <- fr.30.2.sim

alco.30.1.sim.mis.20.n <- alco.30.1.sim

alco.30.2.sim.mis.20.n <- alco.30.2.sim

for (i in 1:S) {
  to_remove.1 <- sample(which(rowSums(fr.30.1.sim[,,i]) < 5), 6, replace=F)
  to_remove.2 <- sample(which(rowSums(fr.30.2.sim[,,i]) < 5), 6, replace=F)
  for (j in 1:6) {
    fr.30.1.sim.mis.20.n[,,i][to_remove.1[j],] <- NA
    alco.30.1.sim.mis.20.n[,i][to_remove.1[j]] <- NA

    fr.30.2.sim.mis.20.n[,,i][to_remove.2[j],] <- NA
    alco.30.2.sim.mis.20.n[,i][to_remove.2[j]] <- NA
  }
}

### 30% missing data ####

fr.30.1.sim.mis.30.n <- fr.30.1.sim
fr.30.2.sim.mis.30.n <- fr.30.2.sim

alco.30.1.sim.mis.30.n <- alco.30.1.sim

alco.30.2.sim.mis.30.n <- alco.30.2.sim

for (i in 1:S) {
  to_remove.1 <- sample(which(rowSums(fr.30.1.sim[,,i]) < 5), 9, replace=F)
  to_remove.2 <- sample(which(rowSums(fr.30.2.sim[,,i]) < 5), 9, replace=F)
  for (j in 1:9) {
    fr.30.1.sim.mis.30.n[,,i][to_remove.1[j],] <- NA
    alco.30.1.sim.mis.30.n[,i][to_remove.1[j]] <- NA

    fr.30.2.sim.mis.30.n[,,i][to_remove.2[j],] <- NA
    alco.30.2.sim.mis.30.n[,i][to_remove.2[j]] <- NA
  }
}

################## Missings depend on the behavior #############################

### 10% missing data ####

fr.30.1.sim.mis.10.b <- fr.30.1.sim
fr.30.2.sim.mis.10.b <- fr.30.2.sim

alco.30.1.sim.mis.10.b <- alco.30.1.sim

alco.30.2.sim.mis.10.b <- alco.30.2.sim

for (i in 1:S) {
  to_remove.1 <- sample(which(alco.30.1.sim[,i] > 2), 3, replace=F)
  to_remove.2 <- sample(which(alco.30.2.sim[,i] > 2), 3, replace=F)
  for (j in 1:3) {
    fr.30.1.sim.mis.10.b[,,i][to_remove.1[j],] <- NA
    alco.30.1.sim.mis.10.b[,i][to_remove.1[j]] <- NA

    fr.30.2.sim.mis.10.b[,,i][to_remove.2[j],] <- NA
    alco.30.2.sim.mis.10.b[,i][to_remove.2[j]] <- NA
  }
}

### 20% missing data ####

fr.30.1.sim.mis.20.b <- fr.30.1.sim
fr.30.2.sim.mis.20.b <- fr.30.2.sim

alco.30.1.sim.mis.20.b <- alco.30.1.sim

alco.30.2.sim.mis.20.b <- alco.30.2.sim

for (i in 1:S) {
  to_remove.1 <- sample(which(alco.30.1.sim[,i] > 2), 6, replace=F)
  to_remove.2 <- sample(which(alco.30.2.sim[,i] > 2), 6, replace=F)
  for (j in 1:6) {
    fr.30.1.sim.mis.20.b[,,i][to_remove.1[j],] <- NA
    alco.30.1.sim.mis.20.b[,i][to_remove.1[j]] <- NA

    fr.30.2.sim.mis.20.b[,,i][to_remove.2[j],] <- NA
    alco.30.2.sim.mis.20.b[,i][to_remove.2[j]] <- NA
  }
}

### 30% missing data ####

fr.30.1.sim.mis.30.b <- fr.30.1.sim
fr.30.2.sim.mis.30.b <- fr.30.2.sim

alco.30.1.sim.mis.30.b <- alco.30.1.sim

alco.30.2.sim.mis.30.b <- alco.30.2.sim

for (i in 1:S) {
  to_remove.1 <- sample(which(alco.30.1.sim[,i] > 2), 9, replace=F)
  to_remove.2 <- sample(which(alco.30.2.sim[,i] > 2), 9, replace=F)
  for (j in 1:9) {
    fr.30.1.sim.mis.30.b[,,i][to_remove.1[j],] <- NA
    alco.30.1.sim.mis.30.b[,i][to_remove.1[j]] <- NA

    fr.30.2.sim.mis.30.b[,,i][to_remove.2[j],] <- NA
    alco.30.2.sim.mis.30.b[,i][to_remove.2[j]] <- NA
  }
}

################## Missings depend on both net and beh #########################

### 10% data ####

fr.30.1.sim.mis.10.nb <- fr.30.1.sim
fr.30.2.sim.mis.10.nb <- fr.30.2.sim

alco.30.1.sim.mis.10.nb <- alco.30.1.sim

alco.30.2.sim.mis.10.nb <- alco.30.2.sim


for (i in 1:S) {
  to_remove.1 <- sample(intersect(which(alco.30.1.sim[,i] > 2),
                                  which(rowSums(fr.30.1.sim[,,i]) < 5)),
                        3, replace = F)
  to_remove.2 <- sample(intersect(which(alco.30.2.sim[,i] > 2),
                                  which(rowSums(fr.30.2.sim[,,i]) < 5)),
                        3, replace = F)
  for (j in 1:3) {
    fr.30.1.sim.mis.10.nb[,,i][to_remove.1[j],] <- NA
    alco.30.1.sim.mis.10.nb[,i][to_remove.1[j]] <- NA

    fr.30.2.sim.mis.10.nb[,,i][to_remove.2[j],] <- NA
    alco.30.2.sim.mis.10.nb[,i][to_remove.2[j]] <- NA
  }
}

### 20% data ####

fr.30.1.sim.mis.20.nb <- fr.30.1.sim
fr.30.2.sim.mis.20.nb <- fr.30.2.sim

alco.30.1.sim.mis.20.nb <- alco.30.1.sim

alco.30.2.sim.mis.20.nb <- alco.30.2.sim

for (i in 1:S) {
  to_remove.1 <- sample(intersect(which(alco.30.1.sim[,i] > 1),
                                  which(rowSums(fr.30.1.sim[,,i]) < 5)),
                        6, replace = F)
  to_remove.2 <- sample(intersect(which(alco.30.2.sim[,i] > 1),
                                  which(rowSums(fr.30.2.sim[,,i]) < 5)),
                        6, replace = F)
  for (j in 1:6) {
    fr.30.1.sim.mis.20.nb[,,i][to_remove.1[j],] <- NA
    alco.30.1.sim.mis.20.nb[,i][to_remove.1[j]] <- NA

    fr.30.2.sim.mis.20.nb[,,i][to_remove.2[j],] <- NA
    alco.30.2.sim.mis.20.nb[,i][to_remove.2[j]] <- NA
  }
}

### 30% data ####

fr.30.1.sim.mis.30.nb <- fr.30.1.sim
fr.30.2.sim.mis.30.nb <- fr.30.2.sim

alco.30.1.sim.mis.30.nb <- alco.30.1.sim

alco.30.2.sim.mis.30.nb <- alco.30.2.sim

for (i in 1:S) {
  to_remove.1 <- sample(intersect(which(alco.30.1.sim[,i] > 1),
                                  which(rowSums(fr.30.1.sim[,,i]) < 5)),
                        9, replace = F)
  to_remove.2 <- sample(intersect(which(alco.30.2.sim[,i] > 1),
                                  which(rowSums(fr.30.2.sim[,,i]) < 5)),
                        9, replace = F)
  for (j in 1:9) {
    fr.30.1.sim.mis.30.nb[,,i][to_remove.1[j],] <- NA
    alco.30.1.sim.mis.30.nb[,i][to_remove.1[j]] <- NA

    fr.30.2.sim.mis.30.nb[,,i][to_remove.2[j],] <- NA
    alco.30.2.sim.mis.30.nb[,i][to_remove.2[j]] <- NA
  }
}


######## SAVING SMALLER DATASET ################################################

save(alco.30.1.sim, alco.30.2.sim,
     alco.30.1.sim.mis.10.b, alco.30.1.sim.mis.10.n,
     alco.30.1.sim.mis.10.nb, alco.30.1.sim.mis.20.b,
     alco.30.1.sim.mis.20.n, alco.30.1.sim.mis.20.nb,
     alco.30.1.sim.mis.30.b, alco.30.1.sim.mis.30.n,
     alco.30.1.sim.mis.30.nb, alco.30.2.sim.mis.10.b,
     alco.30.2.sim.mis.10.n, alco.30.2.sim.mis.10.nb,
     alco.30.2.sim.mis.20.b, alco.30.2.sim.mis.20.n,
     alco.30.2.sim.mis.20.nb, alco.30.2.sim.mis.30.b,
     alco.30.2.sim.mis.30.n, alco.30.2.sim.mis.30.nb,
     ans,
     fr.30.1.sim,
     fr.30.1.sim.mis.10.b, fr.30.1.sim.mis.10.n,
     fr.30.1.sim.mis.10.nb, fr.30.1.sim.mis.20.b,
     fr.30.1.sim.mis.20.n, fr.30.1.sim.mis.20.nb,
     fr.30.1.sim.mis.30.b, fr.30.1.sim.mis.30.n,
     fr.30.1.sim.mis.30.nb,
     fr.30.2.sim,
     fr.30.2.sim.mis.10.b,
     fr.30.2.sim.mis.10.n, fr.30.2.sim.mis.10.nb,
     fr.30.2.sim.mis.20.b, fr.30.2.sim.mis.20.n,
     fr.30.2.sim.mis.20.nb, fr.30.2.sim.mis.30.b,
     fr.30.2.sim.mis.30.n, fr.30.2.sim.mis.30.nb,
     sex.F.30,
     file = "./data/simulated/Data30_2waves_v2.RData")


################################################################################
########                                                               #########
########           Missing data generation - bigger network            #########
########                                                               #########
################################################################################

################## Missings depend on the network ##############################

### 10% missing data ####

fr.60.1.sim.mis.10.n <- fr.60.1.sim
fr.60.2.sim.mis.10.n <- fr.60.2.sim

alco.60.1.sim.mis.10.n <- alco.60.1.sim

alco.60.2.sim.mis.10.n <- alco.60.2.sim

for (i in 1:S) {
  to_remove.1 <- sample(which(rowSums(fr.60.1.sim[,,i]) < 5), 3, replace=F)
  to_remove.2 <- sample(which(rowSums(fr.60.2.sim[,,i]) < 5), 3, replace=F)
  for (j in 1:3) {
    fr.60.1.sim.mis.10.n[,,i][to_remove.1[j],] <- NA
    alco.60.1.sim.mis.10.n[,i][to_remove.1[j]] <- NA

    fr.60.2.sim.mis.10.n[,,i][to_remove.2[j],] <- NA
    alco.60.2.sim.mis.10.n[,i][to_remove.2[j]] <- NA
  }
}

### 20% missing data ####

fr.60.1.sim.mis.20.n <- fr.60.1.sim
fr.60.2.sim.mis.20.n <- fr.60.2.sim

alco.60.1.sim.mis.20.n <- alco.60.1.sim

alco.60.2.sim.mis.20.n <- alco.60.2.sim

for (i in 1:S) {
  to_remove.1 <- sample(which(rowSums(fr.60.1.sim[,,i]) < 5), 6, replace=F)
  to_remove.2 <- sample(which(rowSums(fr.60.2.sim[,,i]) < 5), 6, replace=F)
  for (j in 1:6) {
    fr.60.1.sim.mis.20.n[,,i][to_remove.1[j],] <- NA
    alco.60.1.sim.mis.20.n[,i][to_remove.1[j]] <- NA

    fr.60.2.sim.mis.20.n[,,i][to_remove.2[j],] <- NA
    alco.60.2.sim.mis.20.n[,i][to_remove.2[j]] <- NA
  }
}

### 30% missing data ####

fr.60.1.sim.mis.30.n <- fr.60.1.sim
fr.60.2.sim.mis.30.n <- fr.60.2.sim

alco.60.1.sim.mis.30.n <- alco.60.1.sim

alco.60.2.sim.mis.30.n <- alco.60.2.sim

for (i in 1:S) {
  to_remove.1 <- sample(which(rowSums(fr.60.1.sim[,,i]) < 5), 9, replace=F)
  to_remove.2 <- sample(which(rowSums(fr.60.2.sim[,,i]) < 5), 9, replace=F)
  for (j in 1:9) {
    fr.60.1.sim.mis.30.n[,,i][to_remove.1[j],] <- NA
    alco.60.1.sim.mis.30.n[,i][to_remove.1[j]] <- NA

    fr.60.2.sim.mis.30.n[,,i][to_remove.2[j],] <- NA
    alco.60.2.sim.mis.30.n[,i][to_remove.2[j]] <- NA
  }
}

################## Missings depend on the behavior #############################

### 10% missing data ####

fr.60.1.sim.mis.10.b <- fr.60.1.sim
fr.60.2.sim.mis.10.b <- fr.60.2.sim

alco.60.1.sim.mis.10.b <- alco.60.1.sim

alco.60.2.sim.mis.10.b <- alco.60.2.sim

for (i in 1:S) {
  to_remove.1 <- sample(which(alco.60.1.sim[,i] > 2), 3, replace=F)
  to_remove.2 <- sample(which(alco.60.2.sim[,i] > 2), 3, replace=F)
  for (j in 1:3) {
    fr.60.1.sim.mis.10.b[,,i][to_remove.1[j],] <- NA
    alco.60.1.sim.mis.10.b[,i][to_remove.1[j]] <- NA

    fr.60.2.sim.mis.10.b[,,i][to_remove.2[j],] <- NA
    alco.60.2.sim.mis.10.b[,i][to_remove.2[j]] <- NA
  }
}

### 20% missing data ####

fr.60.1.sim.mis.20.b <- fr.60.1.sim
fr.60.2.sim.mis.20.b <- fr.60.2.sim

alco.60.1.sim.mis.20.b <- alco.60.1.sim

alco.60.2.sim.mis.20.b <- alco.60.2.sim

for (i in 1:S) {
  to_remove.1 <- sample(which(alco.60.1.sim[,i] > 2), 6, replace=F)
  to_remove.2 <- sample(which(alco.60.2.sim[,i] > 2), 6, replace=F)
  for (j in 1:6) {
    fr.60.1.sim.mis.20.b[,,i][to_remove.1[j],] <- NA
    alco.60.1.sim.mis.20.b[,i][to_remove.1[j]] <- NA

    fr.60.2.sim.mis.20.b[,,i][to_remove.2[j],] <- NA
    alco.60.2.sim.mis.20.b[,i][to_remove.2[j]] <- NA
  }
}

### 30% missing data ####

fr.60.1.sim.mis.30.b <- fr.60.1.sim
fr.60.2.sim.mis.30.b <- fr.60.2.sim

alco.60.1.sim.mis.30.b <- alco.60.1.sim

alco.60.2.sim.mis.30.b <- alco.60.2.sim

for (i in 1:S) {
  to_remove.1 <- sample(which(alco.60.1.sim[,i] > 2), 9, replace=F)
  to_remove.2 <- sample(which(alco.60.2.sim[,i] > 2), 9, replace=F)
  for (j in 1:9) {
    fr.60.1.sim.mis.30.b[,,i][to_remove.1[j],] <- NA
    alco.60.1.sim.mis.30.b[,i][to_remove.1[j]] <- NA

    fr.60.2.sim.mis.30.b[,,i][to_remove.2[j],] <- NA
    alco.60.2.sim.mis.30.b[,i][to_remove.2[j]] <- NA
  }
}

################## Missings depend on both net and beh #########################

### 10% data ####

fr.60.1.sim.mis.10.nb <- fr.60.1.sim
fr.60.2.sim.mis.10.nb <- fr.60.2.sim

alco.60.1.sim.mis.10.nb <- alco.60.1.sim

alco.60.2.sim.mis.10.nb <- alco.60.2.sim


for (i in 1:S) {
  to_remove.1 <- sample(intersect(which(alco.60.1.sim[,i] > 2),
                                  which(rowSums(fr.60.1.sim[,,i]) < 5)),
                        3, replace = F)
  to_remove.2 <- sample(intersect(which(alco.60.2.sim[,i] > 2),
                                  which(rowSums(fr.60.2.sim[,,i]) < 5)),
                        3, replace = F)
  for (j in 1:3) {
    fr.60.1.sim.mis.10.nb[,,i][to_remove.1[j],] <- NA
    alco.60.1.sim.mis.10.nb[,i][to_remove.1[j]] <- NA

    fr.60.2.sim.mis.10.nb[,,i][to_remove.2[j],] <- NA
    alco.60.2.sim.mis.10.nb[,i][to_remove.2[j]] <- NA
  }
}

### 20% data ####

fr.60.1.sim.mis.20.nb <- fr.60.1.sim
fr.60.2.sim.mis.20.nb <- fr.60.2.sim

alco.60.1.sim.mis.20.nb <- alco.60.1.sim

alco.60.2.sim.mis.20.nb <- alco.60.2.sim

for (i in 1:S) {
  to_remove.1 <- sample(intersect(which(alco.60.1.sim[,i] > 1),
                                  which(rowSums(fr.60.1.sim[,,i]) < 5)),
                        6, replace = F)
  to_remove.2 <- sample(intersect(which(alco.60.2.sim[,i] > 1),
                                  which(rowSums(fr.60.2.sim[,,i]) < 5)),
                        6, replace = F)
  for (j in 1:6) {
    fr.60.1.sim.mis.20.nb[,,i][to_remove.1[j],] <- NA
    alco.60.1.sim.mis.20.nb[,i][to_remove.1[j]] <- NA

    fr.60.2.sim.mis.20.nb[,,i][to_remove.2[j],] <- NA
    alco.60.2.sim.mis.20.nb[,i][to_remove.2[j]] <- NA
  }
}

### 30% data ####

fr.60.1.sim.mis.30.nb <- fr.60.1.sim
fr.60.2.sim.mis.30.nb <- fr.60.2.sim

alco.60.1.sim.mis.30.nb <- alco.60.1.sim

alco.60.2.sim.mis.30.nb <- alco.60.2.sim

for (i in 1:S) {
  to_remove.1 <- sample(intersect(which(alco.60.1.sim[,i] > 1),
                                  which(rowSums(fr.60.1.sim[,,i]) < 5)),
                        9, replace = F)
  to_remove.2 <- sample(intersect(which(alco.60.2.sim[,i] > 1),
                                  which(rowSums(fr.60.2.sim[,,i]) < 5)),
                        9, replace = F)
  for (j in 1:9) {
    fr.60.1.sim.mis.30.nb[,,i][to_remove.1[j],] <- NA
    alco.60.1.sim.mis.30.nb[,i][to_remove.1[j]] <- NA

    fr.60.2.sim.mis.30.nb[,,i][to_remove.2[j],] <- NA
    alco.60.2.sim.mis.30.nb[,i][to_remove.2[j]] <- NA
  }
}


######## SAVING SMALLER DATASET ################################################

save(alco.60.1.sim, alco.60.2.sim,
     alco.60.1.sim.mis.10.b, alco.60.1.sim.mis.10.n,
     alco.60.1.sim.mis.10.nb, alco.60.1.sim.mis.20.b,
     alco.60.1.sim.mis.20.n, alco.60.1.sim.mis.20.nb,
     alco.60.1.sim.mis.30.b, alco.60.1.sim.mis.30.n,
     alco.60.1.sim.mis.30.nb, alco.60.2.sim.mis.10.b,
     alco.60.2.sim.mis.10.n, alco.60.2.sim.mis.10.nb,
     alco.60.2.sim.mis.20.b, alco.60.2.sim.mis.20.n,
     alco.60.2.sim.mis.20.nb, alco.60.2.sim.mis.30.b,
     alco.60.2.sim.mis.30.n, alco.60.2.sim.mis.30.nb,
     ans,
     fr.60.1.sim,
     fr.60.1.sim.mis.10.b, fr.60.1.sim.mis.10.n,
     fr.60.1.sim.mis.10.nb, fr.60.1.sim.mis.20.b,
     fr.60.1.sim.mis.20.n, fr.60.1.sim.mis.20.nb,
     fr.60.1.sim.mis.30.b, fr.60.1.sim.mis.30.n,
     fr.60.1.sim.mis.30.nb,
     fr.60.2.sim,
     fr.60.2.sim.mis.10.b,
     fr.60.2.sim.mis.10.n, fr.60.2.sim.mis.10.nb,
     fr.60.2.sim.mis.20.b, fr.60.2.sim.mis.20.n,
     fr.60.2.sim.mis.20.nb, fr.60.2.sim.mis.30.b,
     fr.60.2.sim.mis.30.n, fr.60.2.sim.mis.30.nb,
     sex.F.60,
     file = "./data/simulated/Data60_2waves_v2.RData")

