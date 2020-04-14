library(RSiena) # or RSienaTest
source("./simulation/SimulateNetworksBehavior_v4.R")
source("./simulation/siena07ToConvergence.R")

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

############################# Smaller dataset ##################################

fr.30.1 <- friendship.1[1:30, 1:30]
fr.30.2 <- friendship.2[1:30, 1:30]

alco.30 <- alcohol[1:30,1:2]
toba.30 <- tobacco[1:30,1:2]
sex.F.30 <- sex.F[1:30]

############################## Bigger dataset ##################################
set.seed(234)
sample <- sample((1:129), 60, replace = F)
fr.60.1 <- friendship.1[sample, sample]
fr.60.2 <- friendship.2[sample, sample]

alco.60 <- alcohol[sample,1:2]
toba.60 <- tobacco[sample,1:2]
sex.F.60 <- sex.F[sample]

################################################################################
#######                                                                   ######
#######             Smaller dataset model estimation                      ######
#######                                                                   ######
################################################################################

friendship <- sienaDependent(array(c(fr.30.1, fr.30.2),
                                   dim = c(30, 30, 2)))

drinkingbeh <- sienaDependent(alco.30, type = "behavior" )
smokingbeh <- sienaDependent(toba.30, type = "behavior" )
gender <- coCovar(sex.F.30)

myCoEvolutionData <- sienaDataCreate(friendship, gender, smokingbeh,
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

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, altX, egoX, egoXaltX,
                                   interaction1 = "smokingbeh" )
# If we want to parse out whether there is a selection or influence (or both)
# effect for drinking behaviour,
# we need to also include sender, receiver and homophily effects
# of drinking for friendship formation:
myCoEvolutionEff <- includeEffects(myCoEvolutionEff,
                                   name = "drinkingbeh", avAlt,
                                   interaction1 = "friendship")

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, effFrom,
                                   name = "drinkingbeh",
                                   interaction1 = "smokingbeh")

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, effFrom,
                                   name = "drinkingbeh",
                                   interaction1 = "gender")

myCoEvolutionEff <- includeEffects(myCoEvolutionEff,
                                   name = "smokingbeh", avAlt,
                                   interaction1 = "friendship")


myCoEvolutionEff <- includeEffects(myCoEvolutionEff, effFrom,
                                   name = "smokingbeh",
                                   interaction1 = "gender")

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, effFrom,
                                   name = "smokingbeh",
                                   interaction1 = "drinkingbeh")

# Check what effects you have decided to include:

myCoEvolutionEff

# Now we have to define the algorithm settings.
# The defaults are adequate. You only have to specify the filename
# that will receive the results in text format.

myCoEvAlgorithm <- sienaAlgorithmCreate(projname = "model_30", seed = 500)

# Finally, estimate the model; the whole command is put in parentheses
# to have the results printed directly to the screen.

(ans30 <- siena07ToConvergence(alg = myCoEvAlgorithm, dat = myCoEvolutionData,
                               eff = myCoEvolutionEff, threshold = 0.25,
                                nodes = 16))

################################################################################
#######                                                                   ######
#######             Bigger dataset model estimation                       ######
#######                                                                   ######
################################################################################

friendship <- sienaDependent(array(c(fr.60.1, fr.60.2),
                                   dim = c(60, 60, 2)))

drinkingbeh <- sienaDependent(alco.60, type = "behavior")
smokingbeh <- sienaDependent(toba.60, type = "behavior")
gender <- coCovar(sex.F.60)

myCoEvolutionData <- sienaDataCreate(friendship, gender, smokingbeh, drinkingbeh)
myCoEvolutionEff <- getEffects(myCoEvolutionData)

print01Report(myCoEvolutionData, modelname = "complete_data_model_60")

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, dens, recip,
                                   transTrip, transRecTrip, cycle3,
                                   outActSqrt, inPopSqrt)

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, sameX,
                                   interaction1 = "gender" )

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, altX, egoX, egoXaltX,
                                   interaction1 = "drinkingbeh" )

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, altX, egoX, egoXaltX,
                                   interaction1 = "smokingbeh" )
# If we want to parse out whether there is a selection or influence (or both)
# effect for drinking behaviour,
# we need to also include sender, receiver and homophily effects
# of drinking for friendship formation:

myCoEvolutionEff <- includeEffects(myCoEvolutionEff,
                                   name = "drinkingbeh", avAlt,
                                   interaction1 = "friendship")

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, effFrom,
                                   name = "drinkingbeh",
                                   interaction1 = "smokingbeh")

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, effFrom,
                                   name = "drinkingbeh",
                                   interaction1 = "gender")

myCoEvolutionEff <- includeEffects(myCoEvolutionEff,
                                   name = "smokingbeh", avAlt,
                                   interaction1 = "friendship")


myCoEvolutionEff <- includeEffects(myCoEvolutionEff, effFrom,
                                   name = "smokingbeh",
                                   interaction1 = "gender")

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, effFrom,
                                   name = "smokingbeh",
                                   interaction1 = "drinkingbeh")

# Check what effects you have decided to include:

myCoEvolutionEff


# Now we have to define the algorithm settings.
# The defaults are adequate. You only have to specify the filename
# that will receive the results in text format.

myCoEvAlgorithm <- sienaAlgorithmCreate(projname = "model_60", seed = 300)

(ans60 <- siena07ToConvergence(alg = myCoEvAlgorithm, dat = myCoEvolutionData,
                               eff = myCoEvolutionEff, threshold = 0.25,
                               nodes = 16))

################################################################################
#######                                                                   ######
#######             Smaller dataset simulation                            ######
#######                                                                   ######
################################################################################

M <- 2

net.w1 = fr.30.1
covar = sex.F.30
b1.w1 = toba.30[,1]
b2.w1 = alco.30[,1]
n = 30
c1 = 5
c2 = 5

rate = ans30$theta[1]
dens = ans30$theta[2]
rec = ans30$theta[3]
tt = ans30$theta[4]
tRt = ans30$theta[5]
c3 = ans30$theta[6]
inPopSq = ans30$theta[7]
outActSq = ans30$theta[8]

Vsame = ans30$theta[9]

altX.b1 = ans30$theta[10]
egoX.b1 = ans30$theta[11]
egoXaltX.b1 = ans30$theta[12]

rate.b1 = ans30$theta[16]
lin.b1 = ans30$theta[17]
qu.b1 = ans30$theta[18]
avalt.b1 = ans30$theta[19]
effF.b1.V = ans30$theta[20]
effF.b1.b2 = ans30$theta[21]

altX.b2 = ans30$theta[13]
egoX.b2 = ans30$theta[14]
egoXaltX.b2 = ans30$theta[15]

rate.b2 = ans30$theta[22]
lin.b2 = ans30$theta[23]
qu.b2 = ans30$theta[24]
avalt.b2 = ans30$theta[25]
effF.b2.V = ans30$theta[26]
effF.b2.b1 = ans30$theta[27]

fr.30.2.sim <- array(rep(0, n*n*100), c(n, n, 100))
fr.30.3.sim <- array(rep(0, n*n*100), c(n, n, 100))

alco.30.sim <- array(rep(0, n*100), c(n, 100))
toba.30.sim <- array(rep(0, n*100), c(n, 100))

for (i in 1:100) { 
  SN <- SimulateNetworksBehavior(net.w1, covar, b1.w1, b2.w1, n, c1, c2,
                       rate, dens, rec, tt, tRt, c3, inPopSq, outActSq,
                       Vsame,
                       altX.b1, egoX.b1, egoXaltX.b1,
                       rate.b1, lin.b1, qu.b1, avalt.b1, effF.b1.V, effF.b1.b2,
                       altX.b2, egoX.b2, egoXaltX.b2,
                       rate.b2, lin.b2, qu.b2, avalt.b2, effF.b2.V, effF.b2.b1)

  fr.30.2.sim[,,i] <- SN$networks

  toba.30.sim[,i] <- SN$behavior1
  alco.30.sim[,i] <- SN$behavior2
}

################################################################################
#######                                                                   ######
#######             Bigger dataset simulation                             ######
#######                                                                   ######
################################################################################

M <- 2

net.w1 = fr.60.1
covar = sex.F.60
b1.w1 = toba.60[,1]
b2.w1 = alco.60[,1]
n = 60
c1 = 5
c2 = 5

rate = ans60$theta[1] - 3.2
dens = ans60$theta[2]
rec = ans60$theta[3]
tt = ans60$theta[4]
tRt = ans60$theta[5]
c3 = ans60$theta[6]
inPopSq = ans60$theta[7]
outActSq = ans60$theta[8]

Vsame = ans60$theta[9]

altX.b1 = ans60$theta[10]
egoX.b1 = ans60$theta[11]
egoXaltX.b1 = ans60$theta[12]

rate.b1 = ans60$theta[16]
lin.b1 = ans60$theta[17]
qu.b1 = ans60$theta[18]
avalt.b1 = ans60$theta[19]
effF.b1.V = ans60$theta[20]
effF.b1.b2 = ans60$theta[21]

altX.b2 = ans60$theta[13]
egoX.b2 = ans60$theta[14]
egoXaltX.b2 = ans60$theta[15]

rate.b2 = ans60$theta[22]
lin.b2 = ans60$theta[23]
qu.b2 = ans60$theta[24]
avalt.b2 = ans60$theta[25]
effF.b2.V = ans60$theta[26]
effF.b2.b1 = ans60$theta[27]


fr.60.2.sim = array(rep(0, n*n*100), c(n, n, 100))
fr.60.3.sim = array(rep(0, n*n*100), c(n, n, 100))

alco.60.sim = array(rep(0, n*100), c(n, 100))
toba.60.sim = array(rep(0, n*100), c(n, 100))

for (i in 1:100) {
  SN <- SimulateNetworksBehavior(net.w1, covar, b1.w1, b2.w1, n, c1, c2,
                                 rate, dens, rec, tt, tRt, c3, inPopSq, outActSq,
                                 Vsame,
                                 altX.b1, egoX.b1, egoXaltX.b1,
                                 rate.b1, lin.b1, qu.b1, avalt.b1, effF.b1.V, effF.b1.b2,
                                 altX.b2, egoX.b2, egoXaltX.b2,
                                 rate.b2, lin.b2, qu.b2, avalt.b2, effF.b2.V, effF.b2.b1)
  
  fr.60.2.sim[,,i] <- SN$networks
  
  toba.60.sim[,i] <- SN$behavior1
  alco.60.sim[,i] <- SN$behavior2
  
}

################################################################################
########                                                               #########
########           Missing data generation - smaller network           #########
########                                                               #########
################################################################################

################## Missings depend on the network ##############################

### 10% missing data ####

fr.30.1.mis.10.n <- fr.30.1
fr.30.1.mis.10.n[rowSums(fr.30.1) < 1,] <- NA # 3 removed

alco.30.1.mis.10.n <- alco.30[,1]
alco.30.1.mis.10.n[rowSums(fr.30.1) < 1] <- NA

toba.30.1.mis.10.n <- toba.30[,1]
toba.30.1.mis.10.n[rowSums(fr.30.1) < 1] <- NA

fr.30.2.sim.mis.10.n <- fr.30.2.sim

alco.30.2.sim.mis.10.n <- alco.30.sim
toba.30.2.sim.mis.10.n <- toba.30.sim

for (i in 1:100) {
  to_remove.2 <- sample(which(rowSums(fr.30.2.sim[,,i]) < 4), 3, replace=F)
  for (j in 1:3) {
    fr.30.2.sim.mis.10.n[,,i][to_remove.2[j],] <- NA
    alco.30.2.sim.mis.10.n[,i][to_remove.2[j]] <- NA
    toba.30.2.sim.mis.10.n[,i][to_remove.2[j]] <- NA
  }
}

### 20% missing data ####

to_remove.1 <- sample(which(rowSums(fr.30.1) < 2), 6, replace=F)

fr.30.1.mis.20.n <- fr.30.1
fr.30.1.mis.20.n[rowSums(fr.30.1) < 1,] <- NA # 6 removed

alco.30.1.mis.20.n <- alco.30[,1]
alco.30.1.mis.20.n[rowSums(fr.30.1) < 1] <- NA

toba.30.1.mis.20.n <- toba.30[,1]
toba.30.1.mis.20.n[rowSums(fr.30.1) < 1] <- NA

fr.30.2.sim.mis.20.n <- fr.30.2.sim

alco.30.2.sim.mis.20.n <- alco.30.sim
toba.30.2.sim.mis.20.n <- toba.30.sim

for (i in 1:100) {
  to_remove.2 <- sample(which(rowSums(fr.30.2.sim[,,i]) < 5), 6, replace=F)
  for (j in 1:6) {
    fr.30.2.sim.mis.20.n[,,i][to_remove.2[j],] <- NA
    alco.30.2.sim.mis.20.n[,i][to_remove.2[j]] <- NA
    toba.30.2.sim.mis.20.n[,i][to_remove.2[j]] <- NA
  }
}

### 30% missing data ####

to_remove.1 <- sample(which(rowSums(fr.30.1) < 2), 9, replace=F)

fr.30.1.mis.30.n <- fr.30.1
fr.30.1.mis.30.n[rowSums(fr.30.1) < 1,] <- NA # 9 removed

alco.30.1.mis.30.n <- alco.30[,1]
alco.30.1.mis.30.n[rowSums(fr.30.1) < 1] <- NA

toba.30.1.mis.30.n <- toba.30[,1]
toba.30.1.mis.30.n[rowSums(fr.30.1) < 1] <- NA

fr.30.2.sim.mis.30.n <- fr.30.2.sim

alco.30.2.sim.mis.30.n <- alco.30.sim
toba.30.2.sim.mis.30.n <- toba.30.sim

for (i in 1:100) {
  to_remove.2 <- sample(which(rowSums(fr.30.2.sim[,,i]) < 6), 9, replace=F)
  for (j in 1:9) {
    fr.30.2.sim.mis.30.n[,,i][to_remove.2[j],] <- NA
    alco.30.2.sim.mis.30.n[,i][to_remove.2[j]] <- NA
    toba.30.2.sim.mis.30.n[,i][to_remove.2[j]] <- NA
  }
}

################## Missings depend on the behavior #############################

### 10% missing data ####

to_remove.1 <- sample(which(alco.30[,1] > 3), 3, replace = F)

fr.30.1.mis.10.b <- fr.30.1
fr.30.1.mis.10.b[to_remove.1,] <- NA # 3 removed

alco.30.1.mis.10.b <- alco.30[,1]
alco.30.1.mis.10.b[to_remove.1] <- NA

toba.30.1.mis.10.b <- toba.30[,1]
toba.30.1.mis.10.b[to_remove.1] <- NA

fr.30.2.sim.mis.10.b <- fr.30.2.sim

alco.30.2.sim.mis.10.b <- alco.30.sim
toba.30.2.sim.mis.10.b <- toba.30.sim

for (i in 1:100) {
  to_remove.2 <- sample(which(alco.30.sim[,i] > 3), 3, replace=F)
  for (j in 1:3) {
    fr.30.2.sim.mis.10.b[,,i][to_remove.2[j],] <- NA
    alco.30.2.sim.mis.10.b[,i][to_remove.2[j]] <- NA
    toba.30.2.sim.mis.10.b[,i][to_remove.2[j]] <- NA
  }
}

### 20% missing data ####

to_remove.1 <- sample(which(alco.30[,1] > 3), 6, replace = F)

fr.30.1.mis.20.b <- fr.30.1
fr.30.1.mis.20.b[to_remove.1,] <- NA # 6 removed

alco.30.1.mis.20.b <- alco.30[,1]
alco.30.1.mis.20.b[to_remove.1] <- NA

toba.30.1.mis.20.b <- toba.30[,1]
toba.30.1.mis.20.b[to_remove.1] <- NA

fr.30.2.sim.mis.20.b <- fr.30.2.sim

alco.30.2.sim.mis.20.b <- alco.30.sim
toba.30.2.sim.mis.20.b <- toba.30.sim

for (i in 1:100) {
  to_remove.2 <- sample(which(alco.30.sim[,i] > 2), 6, replace=F)
  for (j in 1:6) {
    fr.30.2.sim.mis.20.b[,,i][to_remove.2[j],] <- NA
    alco.30.2.sim.mis.20.b[,i][to_remove.2[j]] <- NA
    toba.30.2.sim.mis.20.b[,i][to_remove.2[j]] <- NA
  }
}

### 30% missing data ####

to_remove.1 <- sample(which(alco.30[,1] > 2), 9, replace = F)

fr.30.1.mis.30.b <- fr.30.1
fr.30.1.mis.30.b[to_remove.1,] <- NA # 9 removed

alco.30.1.mis.30.b <- alco.30[,1]
alco.30.1.mis.30.b[to_remove.1] <- NA

toba.30.1.mis.30.b <- toba.30[,1]
toba.30.1.mis.30.b[to_remove.1] <- NA

fr.30.2.sim.mis.30.b <- fr.30.2.sim

alco.30.2.sim.mis.30.b <- alco.30.sim
toba.30.2.sim.mis.30.b <- toba.30.sim

for (i in 1:100) {
  to_remove.2 <- sample(which(alco.30.sim[,i] > 2), 9, replace=F)
  for (j in 1:9) {
    fr.30.2.sim.mis.30.b[,,i][to_remove.2[j],] <- NA
    alco.30.2.sim.mis.30.b[,i][to_remove.2[j]] <- NA
    toba.30.2.sim.mis.30.b[,i][to_remove.2[j]] <- NA
  }
}

################## Missings depend on both net and beh #########################

### 10% data ####

to_remove.1 <- sample(intersect(which(alco.30[,1] > 3),
                                which(rowSums(fr.30.1) < 2)),
                      3, replace = F)

fr.30.1.mis.10.nb <- fr.30.1
fr.30.1.mis.10.nb[to_remove.1,] <- NA # 3 removed

alco.30.1.mis.10.nb <- alco.30[,1]
alco.30.1.mis.10.nb[to_remove.1] <- NA

toba.30.1.mis.10.nb <- toba.30[,1]
toba.30.1.mis.10.nb[to_remove.1] <- NA

fr.30.2.sim.mis.10.nb <- fr.30.2.sim

alco.30.2.sim.mis.10.nb <- alco.30.sim
toba.30.2.sim.mis.10.nb <- toba.30.sim


for (i in 1:100) {
  to_remove.2 <- sample(intersect(which(alco.30.sim[,i] > 1),
                                  which(rowSums(fr.30.2.sim[,,i]) < 4)),
                        3, replace = F)
  for (j in 1:3) {
    fr.30.2.sim.mis.10.nb[,,i][to_remove.2[j],] <- NA
    alco.30.2.sim.mis.10.nb[,i][to_remove.2[j]] <- NA
    toba.30.2.sim.mis.10.nb[,i][to_remove.2[j]] <- NA
  }
}

### 20% data ####

to_remove.1 <- sample(intersect(which(alco.30[,1] > 1),
                                which(rowSums(fr.30.1) < 2)),
                      6, replace = F)

fr.30.1.mis.20.nb <- fr.30.1
fr.30.1.mis.20.nb[to_remove.1,] <- NA # 6 removed

alco.30.1.mis.20.nb <- alco.30[,1]
alco.30.1.mis.20.nb[to_remove.1] <- NA

toba.30.1.mis.20.nb <- toba.30[,1]
toba.30.1.mis.20.nb[to_remove.1] <- NA

fr.30.2.sim.mis.20.nb <- fr.30.2.sim

alco.30.2.sim.mis.20.nb <- alco.30.sim
toba.30.2.sim.mis.20.nb <- toba.30.sim

for (i in 1:100) {
  to_remove.2 <- sample(intersect(which(alco.30.sim[,i] > 1),
                                  which(rowSums(fr.30.2.sim[,,i]) < 6)),
                        6, replace = F)
  for (j in 1:6) {
    fr.30.2.sim.mis.20.nb[,,i][to_remove.2[j],] <- NA
    alco.30.2.sim.mis.20.nb[,i][to_remove.2[j]] <- NA
    toba.30.2.sim.mis.20.nb[,i][to_remove.2[j]] <- NA
  }
}

### 30% data ####

to_remove.1 <- sample(intersect(which(alco.30[,1] > 1),
                                which(rowSums(fr.30.1) < 2)),
                      9, replace = F)

fr.30.1.mis.30.nb <- fr.30.1
fr.30.1.mis.30.nb[to_remove.1,] <- NA # 9 removed

alco.30.1.mis.30.nb <- alco.30[,1]
alco.30.1.mis.30.nb[to_remove.1] <- NA

toba.30.1.mis.30.nb <- toba.30[,1]
toba.30.1.mis.30.nb[to_remove.1] <- NA

fr.30.2.sim.mis.30.nb <- fr.30.2.sim

alco.30.2.sim.mis.30.nb <- alco.30.sim
toba.30.2.sim.mis.30.nb <- toba.30.sim

for (i in 1:100) {
  to_remove.2 <- sample(intersect(which(alco.30.sim[,i] > 1),
                                  which(rowSums(fr.30.2.sim[,,i]) < 7)),
                        9, replace = F)
  for (j in 1:9) {
    fr.30.2.sim.mis.20.nb[,,i][to_remove.2[j],] <- NA
    alco.30.2.sim.mis.20.nb[,i][to_remove.2[j]] <- NA
    toba.30.2.sim.mis.20.nb[,i][to_remove.2[j]] <- NA
  }
}

a <- ls()

######## SAVING SMALLER DATASET ################################################

save(alco.30,
     alco.30.1.mis.10.b,alco.30.1.mis.10.n,
     alco.30.1.mis.10.nb,alco.30.1.mis.20.b,
     alco.30.1.mis.20.n,alco.30.1.mis.20.nb,
     alco.30.1.mis.30.b,alco.30.1.mis.30.n,
     alco.30.1.mis.30.nb,alco.30.2.sim.mis.10.b,
     alco.30.2.sim.mis.10.n,alco.30.2.sim.mis.10.nb,
     alco.30.2.sim.mis.20.b,alco.30.2.sim.mis.20.n,
     alco.30.2.sim.mis.20.nb,alco.30.2.sim.mis.30.b,
     alco.30.2.sim.mis.30.n,alco.30.2.sim.mis.30.nb,
     ans30,
     fr.30.1,
     fr.30.1.mis.10.b,fr.30.1.mis.10.n,
     fr.30.1.mis.10.nb,fr.30.1.mis.20.b,
     fr.30.1.mis.20.n,fr.30.1.mis.20.nb,
     fr.30.1.mis.30.b,fr.30.1.mis.30.n,
     fr.30.1.mis.30.nb,fr.30.2,
     fr.30.2.sim.mis.10.b,
     fr.30.2,
     fr.30.2.sim.mis.10.n,fr.30.2.sim.mis.10.nb,
     fr.30.2.sim.mis.20.b,fr.30.2.sim.mis.20.n,
     fr.30.2.sim.mis.20.nb,fr.30.2.sim.mis.30.b,
     fr.30.2.sim.mis.30.n,fr.30.2.sim.mis.30.nb,
     toba.30,
     toba.30.1.mis.10.b,toba.30.1.mis.10.n,
     toba.30.1.mis.10.nb,toba.30.1.mis.20.b,
     toba.30.1.mis.20.n,toba.30.1.mis.20.nb,
     toba.30.1.mis.30.b,toba.30.1.mis.30.n,
     toba.30.1.mis.30.nb,toba.30.2.sim.mis.10.b,
     toba.30.2.sim.mis.10.n,toba.30.2.sim.mis.10.nb,
     toba.30.2.sim.mis.20.b,toba.30.2.sim.mis.20.n,
     toba.30.2.sim.mis.20.nb,toba.30.2.sim.mis.30.b,
     toba.30.2.sim.mis.30.n,toba.30.2.sim.mis.30.nb,
     sex.F.30,
     file = "./data/simulated/Data30_2waves.RData")


################################################################################
########                                                               #########
########           Missing data generation - bigger network            #########
########                                                               #########
################################################################################

################## Missings depend on the network ##############################

### 10% missing data ####

fr.60.1.mis.10.n <- fr.60.1
fr.60.1.mis.10.n[rowSums(fr.60.1) < 1,] <- NA # 6 removed

alco.60.1.mis.10.n <- alco.60[,1]
alco.60.1.mis.10.n[rowSums(fr.60.1) < 1] <- NA

toba.60.1.mis.10.n <- toba.60[,1]
toba.60.1.mis.10.n[rowSums(fr.60.1) < 1] <- NA

fr.60.2.sim.mis.10.n <- fr.60.2.sim

alco.60.2.sim.mis.10.n <- alco.60.sim
toba.60.2.sim.mis.10.n <- toba.60.sim

for (i in 1:100) {
  to_remove.2 <- sample(which(rowSums(fr.60.2.sim[,,i]) < 4), 6, replace=F)
  for (j in 1:6) {
    fr.60.2.sim.mis.10.n[,,i][to_remove.2[j],] <- NA
    alco.60.2.sim.mis.10.n[,i][to_remove.2[j]] <- NA
    toba.60.2.sim.mis.10.n[,i][to_remove.2[j]] <- NA
  }
}

### 20% missing data ####

to_remove.1 <- sample(which(rowSums(fr.60.1) < 2), 12, replace=F)

fr.60.1.mis.20.n <- fr.60.1
fr.60.1.mis.20.n[rowSums(fr.60.1) < 1,] <- NA # 12 removed

alco.60.1.mis.20.n <- alco.60[,1]
alco.60.1.mis.20.n[rowSums(fr.60.1) < 1] <- NA

toba.60.1.mis.20.n <- toba.60[,1]
toba.60.1.mis.20.n[rowSums(fr.60.1) < 1] <- NA

fr.60.2.sim.mis.20.n <- fr.60.2.sim

alco.60.2.sim.mis.20.n <- alco.60.sim
toba.60.2.sim.mis.20.n <- toba.60.sim

for (i in 1:100) {
  to_remove.2 <- sample(which(rowSums(fr.60.2.sim[,,i]) < 5), 12, replace=F)
  for (j in 1:12) {
    fr.60.2.sim.mis.20.n[,,i][to_remove.2[j],] <- NA
    alco.60.2.sim.mis.20.n[,i][to_remove.2[j]] <- NA
    toba.60.2.sim.mis.20.n[,i][to_remove.2[j]] <- NA
  }
}

### 30% missing data ####

to_remove.1 <- sample(which(rowSums(fr.60.1) < 2), 18, replace=F)

fr.60.1.mis.30.n <- fr.60.1
fr.60.1.mis.30.n[rowSums(fr.60.1) < 1,] <- NA # 18 removed

alco.60.1.mis.30.n <- alco.60[,1]
alco.60.1.mis.30.n[rowSums(fr.60.1) < 1] <- NA

toba.60.1.mis.30.n <- toba.60[,1]
toba.60.1.mis.30.n[rowSums(fr.60.1) < 1] <- NA

fr.60.2.sim.mis.30.n <- fr.60.2.sim

alco.60.2.sim.mis.30.n <- alco.60.sim
toba.60.2.sim.mis.30.n <- toba.60.sim

for (i in 1:100) {
  to_remove.2 <- sample(which(rowSums(fr.60.2.sim[,,i]) < 6), 18, replace=F)
  for (j in 1:18) {
    fr.60.2.sim.mis.30.n[,,i][to_remove.2[j],] <- NA
    alco.60.2.sim.mis.30.n[,i][to_remove.2[j]] <- NA
    toba.60.2.sim.mis.30.n[,i][to_remove.2[j]] <- NA
  }
}

################## Missings depend on the behavior #############################

### 10% missing data ####

to_remove.1 <- sample(which(alco.60[,1] > 3), 6, replace = F)

fr.60.1.mis.10.b <- fr.60.1
fr.60.1.mis.10.b[to_remove.1,] <- NA # 6 removed

alco.60.1.mis.10.b <- alco.60[,1]
alco.60.1.mis.10.b[to_remove.1] <- NA

toba.60.1.mis.10.b <- toba.60[,1]
toba.60.1.mis.10.b[to_remove.1] <- NA

fr.60.2.sim.mis.10.b <- fr.60.2.sim

alco.60.2.sim.mis.10.b <- alco.60.sim
toba.60.2.sim.mis.10.b <- toba.60.sim

for (i in 1:100) {
  to_remove.2 <- sample(which(alco.60.sim[,i] > 3), 6, replace=F)
  for (j in 1:6) {
    fr.60.2.sim.mis.10.b[,,i][to_remove.2[j],] <- NA
    alco.60.2.sim.mis.10.b[,i][to_remove.2[j]] <- NA
    toba.60.2.sim.mis.10.b[,i][to_remove.2[j]] <- NA
  }
}

### 20% missing data ####

to_remove.1 <- sample(which(alco.60[,1] > 2), 12, replace = F)

fr.60.1.mis.20.b <- fr.60.1
fr.60.1.mis.20.b[to_remove.1,] <- NA # 12 removed

alco.60.1.mis.20.b <- alco.60[,1]
alco.60.1.mis.20.b[to_remove.1] <- NA

toba.60.1.mis.20.b <- toba.60[,1]
toba.60.1.mis.20.b[to_remove.1] <- NA

fr.60.2.sim.mis.20.b <- fr.60.2.sim

alco.60.2.sim.mis.20.b <- alco.60.sim
toba.60.2.sim.mis.20.b <- toba.60.sim

for (i in 1:100) {
  to_remove.2 <- sample(which(alco.60.sim[,i] > 2), 12, replace=F)
  for (j in 1:12) {
    fr.60.2.sim.mis.20.b[,,i][to_remove.2[j],] <- NA
    alco.60.2.sim.mis.20.b[,i][to_remove.2[j]] <- NA
    toba.60.2.sim.mis.20.b[,i][to_remove.2[j]] <- NA
  }
}

### 30% missing data ####

to_remove.1 <- sample(which(alco.60[,1] > 2), 18, replace = F)

fr.60.1.mis.30.b <- fr.60.1
fr.60.1.mis.30.b[to_remove.1,] <- NA # 18 removed

alco.60.1.mis.30.b <- alco.60[,1]
alco.60.1.mis.30.b[to_remove.1] <- NA

toba.60.1.mis.30.b <- toba.60[,1]
toba.60.1.mis.30.b[to_remove.1] <- NA

fr.60.2.sim.mis.30.b <- fr.60.2.sim

alco.60.2.sim.mis.30.b <- alco.60.sim
toba.60.2.sim.mis.30.b <- toba.60.sim

for (i in 1:100) {
  to_remove.2 <- sample(which(alco.60.sim[,i] > 2), 18, replace=F)
  for (j in 1:18) {
    fr.60.2.sim.mis.30.b[,,i][to_remove.2[j],] <- NA
    alco.60.2.sim.mis.30.b[,i][to_remove.2[j]] <- NA
    toba.60.2.sim.mis.30.b[,i][to_remove.2[j]] <- NA
  }
}

################## Missings depend on both net and beh #########################

### 10% data ####

to_remove.1 <- sample(intersect(which(alco.60[,1] > 3),
                                which(rowSums(fr.60.1) < 2)),
                      6, replace = F)

fr.60.1.mis.10.nb <- fr.60.1
fr.60.1.mis.10.nb[to_remove.1,] <- NA # 6 removed

alco.60.1.mis.10.nb <- alco.60[,1]
alco.60.1.mis.10.nb[to_remove.1] <- NA

toba.60.1.mis.10.nb <- toba.60[,1]
toba.60.1.mis.10.nb[to_remove.1] <- NA

fr.60.2.sim.mis.10.nb <- fr.60.2.sim

alco.60.2.sim.mis.10.nb <- alco.60.sim
toba.60.2.sim.mis.10.nb <- toba.60.sim


for (i in 1:100) {
  to_remove.2 <- sample(intersect(which(alco.60.sim[,i] > 1),
                                  which(rowSums(fr.60.2.sim[,,i]) < 4)),
                        6, replace = F)
  for (j in 1:6) {
    fr.60.2.sim.mis.10.nb[,,i][to_remove.2[j],] <- NA
    alco.60.2.sim.mis.10.nb[,i][to_remove.2[j]] <- NA
    toba.60.2.sim.mis.10.nb[,i][to_remove.2[j]] <- NA
  }
}

### 20% data ####

to_remove.1 <- sample(intersect(which(alco.60[,1] > 1),
                                which(rowSums(fr.60.1) < 2)),
                      12, replace = F)

fr.60.1.mis.20.nb <- fr.60.1
fr.60.1.mis.20.nb[to_remove.1,] <- NA # 12 removed

alco.60.1.mis.20.nb <- alco.60[,1]
alco.60.1.mis.20.nb[to_remove.1] <- NA

toba.60.1.mis.20.nb <- toba.60[,1]
toba.60.1.mis.20.nb[to_remove.1] <- NA

fr.60.2.sim.mis.20.nb <- fr.60.2.sim

alco.60.2.sim.mis.20.nb <- alco.60.sim
toba.60.2.sim.mis.20.nb <- toba.60.sim

for (i in 1:100) {
  to_remove.2 <- sample(intersect(which(alco.60.sim[,i] > 1),
                                  which(rowSums(fr.60.2.sim[,,i]) < 6)),
                        12, replace = F)
  for (j in 1:12) {
    fr.60.2.sim.mis.20.nb[,,i][to_remove.2[j],] <- NA
    alco.60.2.sim.mis.20.nb[,i][to_remove.2[j]] <- NA
    toba.60.2.sim.mis.20.nb[,i][to_remove.2[j]] <- NA
  }
}

### 30% data ####

to_remove.1 <- sample(intersect(which(alco.60[,1] > 1),
                                which(rowSums(fr.60.1) < 2)),
                      18, replace = F)

fr.60.1.mis.30.nb <- fr.60.1
fr.60.1.mis.30.nb[to_remove.1,] <- NA # 18 removed

alco.60.1.mis.30.nb <- alco.60[,1]
alco.60.1.mis.30.nb[to_remove.1] <- NA

toba.60.1.mis.30.nb <- toba.60[,1]
toba.60.1.mis.30.nb[to_remove.1] <- NA

fr.60.2.sim.mis.30.nb <- fr.60.2.sim

alco.60.2.sim.mis.30.nb <- alco.60.sim
toba.60.2.sim.mis.30.nb <- toba.60.sim

for (i in 1:100) {
  to_remove.2 <- sample(intersect(which(alco.60.sim[,i] > 1),
                                  which(rowSums(fr.60.2.sim[,,i]) < 7)),
                        18, replace = F)
  for (j in 1:18) {
    fr.60.2.sim.mis.20.nb[,,i][to_remove.2[j],] <- NA
    alco.60.2.sim.mis.20.nb[,i][to_remove.2[j]] <- NA
    toba.60.2.sim.mis.20.nb[,i][to_remove.2[j]] <- NA
  }
}


######## SAVING SMALLER DATASET ################################################

save(alco.60,
     alco.60.1.mis.10.b,alco.60.1.mis.10.n,
     alco.60.1.mis.10.nb,alco.60.1.mis.20.b,
     alco.60.1.mis.20.n,alco.60.1.mis.20.nb,
     alco.60.1.mis.30.b,alco.60.1.mis.30.n,
     alco.60.1.mis.30.nb,alco.60.2.sim.mis.10.b,
     alco.60.2.sim.mis.10.n,alco.60.2.sim.mis.10.nb,
     alco.60.2.sim.mis.20.b,alco.60.2.sim.mis.20.n,
     alco.60.2.sim.mis.20.nb,alco.60.2.sim.mis.30.b,
     alco.60.2.sim.mis.30.n,alco.60.2.sim.mis.30.nb,
     ans60,
     fr.60.1,
     fr.60.1.mis.10.b,fr.60.1.mis.10.n,
     fr.60.1.mis.10.nb,fr.60.1.mis.20.b,
     fr.60.1.mis.20.n,fr.60.1.mis.20.nb,
     fr.60.1.mis.30.b,fr.60.1.mis.30.n,
     fr.60.1.mis.30.nb,fr.60.2,
     fr.60.2.sim.mis.10.b,
     fr.60.2,
     fr.60.2.sim.mis.10.n,fr.60.2.sim.mis.10.nb,
     fr.60.2.sim.mis.20.b,fr.60.2.sim.mis.20.n,
     fr.60.2.sim.mis.20.nb,fr.60.2.sim.mis.30.b,
     fr.60.2.sim.mis.30.n,fr.60.2.sim.mis.30.nb,
     toba.60,
     toba.60.1.mis.10.b,toba.60.1.mis.10.n,
     toba.60.1.mis.10.nb,toba.60.1.mis.20.b,
     toba.60.1.mis.20.n,toba.60.1.mis.20.nb,
     toba.60.1.mis.30.b,toba.60.1.mis.30.n,
     toba.60.1.mis.30.nb,toba.60.2.sim.mis.10.b,
     toba.60.2.sim.mis.10.n,toba.60.2.sim.mis.10.nb,
     toba.60.2.sim.mis.20.b,toba.60.2.sim.mis.20.n,
     toba.60.2.sim.mis.20.nb,toba.60.2.sim.mis.30.b,
     toba.60.2.sim.mis.30.n,toba.60.2.sim.mis.30.nb,
     sex.F.60,
     file = "./data/simulated/Data60_2waves.RData")
