# install.packages('RSiena', repos='http://cran.us.r-project.org')
# conda install -c conda-forge r-mice

library(RSiena) # or RSienaTest
source('./simulation/siena07ToConvergence.R')
source('./smaller_experiment/simulateNetworkBehavior.R')

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
Nnodes = 31

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


net.w1 = fr.1
b1.w1 = alco[,1]
n = 60
c1 = 5
rate = ans60$theta[1] #2
dens = ans60$theta[2] 
rec = ans60$theta[3]
tt = ans60$theta[4]
tRt = ans60$theta[5]
c3 = ans60$theta[6]
inPopSq = ans60$theta[7]
outActSq = ans60$theta[8]
altX.b1 = ans60$theta[9] 
egoX.b1 = ans60$theta[10] #+ 1
egoXaltX.b1 = ans60$theta[11] #+ 1
rate.b1 = ans60$theta[12]
lin.b1 = ans60$theta[13]
qu.b1 = ans60$theta[14]
avalt.b1 = ans60$theta[15]

S = 100

fr.60.2.sim <- array(rep(0, n*n*S), c(n, n, S))
alco.60.2.sim <- array(rep(0, n*S), c(n, S))

for (i in 1:S) { 
  SN <- SimulateNetworksBehavior(net.w1, b1.w1, n, c1,
                                 rate, dens, rec, tt, tRt, c3, outActSq, inPopSq, 
                                 altX.b1, egoX.b1, egoXaltX.b1,
                                 rate.b1, lin.b1, qu.b1, avalt.b1)
  
  fr.60.2.sim[,,i] <- SN$networks
  alco.60.2.sim[,i] <- SN$behavior
  print(i)
  print(which(rowSums(fr.60.2.sim[,,i]) == 0 & colSums(fr.60.2.sim[,,i]) == 0)) 
}

save(fr.1, alco,
fr.60.2.sim, alco.60.2.sim,
ans60, file = "./data/simulated/smaller_exp.RData")