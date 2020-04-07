library(RSiena) # or RSienaTest
source("./simulation/SimulateNetworksBehavior_v3.R")
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
friendship.3 <- friendship.3[not_missing, not_missing]

alcohol <- alcohol[not_missing, ]
tobacco <- tobacco[not_missing, ]
age <- age[not_missing]
sex.F <- sex.F[not_missing] - 1 # to 0s and 1s 

############################# Smaller dataset ##################################

fr.30.1 <- friendship.1[1:30, 1:30]
fr.30.2 <- friendship.2[1:30, 1:30]
fr.30.3 <- friendship.3[1:30, 1:30]

alco.30 <- alcohol[1:30, ]
toba.30 <- tobacco[1:30, ]
age.30 <- age[1:30]
sex.F.30 <- sex.F[1:30]

############################## Bigger dataset ##################################

fr.60.1 <- friendship.1[70:129, 70:129]
fr.60.2 <- friendship.2[70:129, 70:129]
fr.60.3 <- friendship.3[70:129, 70:129]

alco.60 <- alcohol[70:129, ]
toba.60 <- tobacco[70:129, ]
age.60 <- age[70:129]
sex.F.60 <- sex.F[70:129]

################################################################################
#######                                                                   ######
#######             Smaller dataset model estimation                      ######
#######                                                                   ######
################################################################################

friendship <- sienaDependent(array(c(fr.30.1, fr.30.2, fr.30.3),
                                     dim = c(30, 30, 3)))

drinkingbeh <- sienaDependent(alco.30, type = "behavior" )
gender <- coCovar(sex.F.30)
smoke1 <- coCovar(toba.30[, 1])

myCoEvolutionData <- sienaDataCreate(friendship, gender, smoke1, drinkingbeh)
myCoEvolutionEff <- getEffects(myCoEvolutionData)

print01Report(myCoEvolutionData, modelname = "complete_data_model_30")

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, transTrip, cycle3)


myCoEvolutionEff <- includeEffects(myCoEvolutionEff, simX,
                                    interaction1 = "smoke1" )
myCoEvolutionEff <- includeEffects(myCoEvolutionEff, sameX,
                                   interaction1 = "gender" )

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, altX, egoX, egoXaltX,
                                   interaction1 = "drinkingbeh" )
# If we want to parse out whether there is a selection or influence (or both)
# effect for drinking behaviour,
# we need to also include sender, receiver and homophily effects
# of drinking for friendship formation:
myCoEvolutionEff <- includeEffects(myCoEvolutionEff,
                                    name = "drinkingbeh", avAlt, avSim,
                                    interaction1 = "friendship")

myCoEvolutionEff <- includeEffects(myCoEvolutionEff,
                                    name = "drinkingbeh",
                                    quad, avAlt, avSim,
                                    interaction1 = c("", "friendship"))

# Check what effects you have decided to include:

myCoEvolutionEff

# Now we have to define the algorithm settings.
# The defaults are adequate. You only have to specify the filename
# that will receive the results in text format.

myCoEvAlgorithm <- sienaAlgorithmCreate(projname = "model_30", seed = 500)

# Finally, estimate the model; the whole command is put in parentheses
# to have the results printed directly to the screen.

(ans30 <- siena07ToConvergence(alg = myCoEvAlgorithm, dat = myCoEvolutionData,
                             eff = myCoEvolutionEff, threshold = 0.20))

################################################################################
#######                                                                   ######
#######             Bigger dataset model estimation                       ######
#######                                                                   ######
################################################################################

friendship <- sienaDependent(array(c(fr.60.1, fr.60.2, fr.60.3),
                                   dim = c(60, 60, 3)))

drinkingbeh <- sienaDependent(alco.60, type = "behavior")
smoke1 <- coCovar(toba.60[, 1])
gender <- coCovar(sex.F.60)

myCoEvolutionData <- sienaDataCreate(friendship, gender, smoke1, drinkingbeh)
myCoEvolutionEff <- getEffects(myCoEvolutionData)

print01Report(myCoEvolutionData, modelname = "complete_data_model_60")

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, transTrip, cycle3)


myCoEvolutionEff <- includeEffects(myCoEvolutionEff, simX,
                                   interaction1 = "smoke1" )
myCoEvolutionEff <- includeEffects(myCoEvolutionEff, sameX,
                                   interaction1 = "gender")

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, altX, egoX, egoXaltX,
                                   interaction1 = "drinkingbeh")
# If we want to parse out whether there is a selection or influence (or both)
# effect for drinking behaviour,
# we need to also include sender, receiver and homophily effects
# of drinking for friendship formation:

myCoEvolutionEff <- includeEffects(myCoEvolutionEff,
                                   name = "drinkingbeh", avAlt, avSim,
                                   interaction1 = "friendship")

myCoEvolutionEff <- includeEffects(myCoEvolutionEff,
                                   name = "drinkingbeh",
                                   quad, avAlt, avSim,
                                   interaction1 = c("", "friendship"))

# Check what effects you have decided to include:

myCoEvolutionEff


# Now we have to define the algorithm settings.
# The defaults are adequate. You only have to specify the filename
# that will receive the results in text format.

myCoEvAlgorithm <- sienaAlgorithmCreate(projname = "model_60", seed = 300)

(ans60 <- siena07ToConvergence(alg = myCoEvAlgorithm, dat = myCoEvolutionData,
                 eff = myCoEvolutionEff, threshold = 0.20))

################################################################################
#######                                                                   ######
#######             Smaller dataset simulation                            ######
#######                                                                   ######
################################################################################

n <- 30
M <- 3
c <- 5

rate <- ans30$theta[2]
dens <- ans30$theta[3]
rec <- ans30$theta[4]
tt <- ans30$theta[5]
c3 <- ans30$theta[6]

Vasame  <- ans30$theta[7]
Vbsim <- ans30$theta[8]

altX.b <- ans30$theta[9]
egoX.b <- ans30$theta[10]
egoXaltX.b <- ans30$theta[11]
rate.b <- ans30$theta[13]
lin.b <- ans30$theta[14]
qu.b <- ans30$theta[15]
avalt.b <- ans30$theta[17]
avsim.b <- ans30$theta[16]

SN_30 <- SimulateNetworksBehavior(net.w1 = fr.30.1, covara = age.30,
                                  covarb = toba.30[,1], b.w1 = alco.30[,1],
                                  n, M, c,
                                  rate, dens, rec, tt, c3,
                                  Vasame, Vbsim, altX.b, egoX.b, egoXaltX.b,
                                  rate.b, lin.b, qu.b, avalt.b, avsim.b)


################################################################################
#######                                                                   ######
#######             Bigger dataset simulation                             ######
#######                                                                   ######
################################################################################

n <- 60
M <- 3
c <- 5

rate <- ans60$theta[2]
dens <- ans60$theta[3]
rec <- ans60$theta[4]
tt <- ans60$theta[5]
c3 <- ans60$theta[6]

Vasame  <- ans60$theta[7]
Vbsim <- ans60$theta[8]

altX.b <- ans60$theta[9]
egoX.b <- ans60$theta[10]
egoXaltX.b <- ans60$theta[11]
rate.b <- ans60$theta[13]
lin.b <- ans60$theta[14]
qu.b <- ans60$theta[15]
avalt.b <- ans60$theta[17]
avsim.b <- ans60$theta[16]

SN_60 <- SimulateNetworksBehavior(net.w1 = fr.60.1, covara = age.60,
                                  covarb = toba.60[,1], b.w1 = alco.60[,1],
                                  n, M, c,
                                  rate, dens, rec, tt, c3,
                                  Vasame, Vbsim, altX.b, egoX.b, egoXaltX.b,
                                  rate.b, lin.b, qu.b, avalt.b, avsim.b)
