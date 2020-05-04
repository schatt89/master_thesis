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

S = 100
N = 60
M = 2
Nnodes = 31
c = 6

friendship <- sienaDependent(array(c(fr.1, fr.2),
                                   dim = c(N, N, M)))

drinkingbeh <- sienaDependent(alco[,1:M], type = "behavior" )

myCoEvolutionData <- sienaDataCreate(friendship,
                                     drinkingbeh)

myCoEvolutionEff <- getEffects(myCoEvolutionData)

print01Report(myCoEvolutionData, modelname = "complete_data_model_60")

myCoEvolutionEff <- includeEffects(myCoEvolutionEff, gwespFF,
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

#myCoEvolutionEff <- setEffect(myCoEvolutionEff, outTrunc, parameter = 7,
#                                            test = FALSE, fix = TRUE,          #initialValue = -1000)

# Check what effects you have decided to include:

myCoEvolutionEff

# Now we have to define the algorithm settings.
# The defaults are adequate. You only have to specify the filename
# that will receive the results in text format.

myCoEvAlgorithm <- sienaAlgorithmCreate(projname = "model", seed = 300,
MaxDegree = c(friendship = 6))

# Finally, estimate the model; the whole command is put in parentheses
# to have the results printed directly to the screen.

(ans60 <- siena07ToConvergence(alg = myCoEvAlgorithm, dat = myCoEvolutionData,
                               eff = myCoEvolutionEff, threshold = 0.25,
                               nodes = Nnodes))

################################################################################
#######                                                                   ######
#######                   Goodness of fit check                           ######
#######                                                                   ######
################################################################################

gof1.id <- sienaGOF(ans60, verbose=TRUE,
                    varName="friendship", IndegreeDistribution)
gof1.id
plot(gof1.id)

gof1.od <- sienaGOF(ans60, verbose=TRUE, varName="friendship",
                    OutdegreeDistribution)
gof1.od
plot(gof1.od)

TriadCensus <- function(i, data, sims, wave, groupName, varName, levls=1:16){
  unloadNamespace("igraph") # to avoid package clashes
  require(sna)
  require(network)
  x <- networkExtraction(i, data, sims, wave, groupName, varName)
  if (network.edgecount(x) <= 0){x <- symmetrize(x)}
  # because else triad.census(x) will lead to an error
  tc <- sna::triad.census(x)[1,levls]
  # triad names are transferred automatically
  tc
}

gof1.tc <- sienaGOF(ans60, verbose=TRUE,
                    varName="friendship", TriadCensus)

gof1.tc

plot(gof1.tc, scale=TRUE, center=TRUE)

################################################################################
#######                                                                   ######
#######                   Datasets simulation                             ######
#######                                                                   ######
################################################################################

getNet <- function(observedNet,edgeList) {
  # observedNet = observed network as adjacency with missing data
  # edgeList = edgelist that is returned by siena07(...)$sims
  observedNet[is.na(observedNet)] <- 0
  for (i in 1:nrow(edgeList)) {
    observedNet[edgeList[i,1],edgeList[i,2]] <- 1
  }
  return(observedNet)
}


sims = ans60$sims

fr.60.2.sim <- array(rep(NA, N*N*S), c(N,N,S))
alco.60.2.sim <- array(rep(NA, N*S), c(N, S))
edgelists <- list()
shape <- fr.1
shape[1:N,] <- 0

i <- 1
j <- 100
while (i <= 100) {
  edgelist = sims[[j]]$Data1$friendship$`1`
  avdegree = mean(table(edgelist[,1]))
  alcosim = sims[[j]]$Data1$drinkingbeh$`1`
  avalco = mean(alcosim)
  alcounique = length(unique(alcosim))
  print(j)
  if (avdegree >= 3 & avdegree < 3.2 &
      avalco >= 2.7 & avalco <= 3 & alcounique  == 5) {
    fr.60.2.sim[,,i] = getNet(shape, edgelist)
    alco.60.2.sim[,i] = alcosim
    i <- i + 1
    j <- j + 1
  } else {
    j <- j + 1
  }
}

save(fr.1, alco,
     fr.60.2.sim, alco.60.2.sim,
     ans60, file = "./data/simulated/smaller_exp.RData")
