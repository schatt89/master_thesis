library(RSiena) # or RSienaTest
source('./simulation/SimulateNetworksBehavior.R')
source('./simulation/siena07ToConvergence.R')

################################################################################
#######                                                                   ######
#######                      Data preparation                             ######
#######                                                                   ######
################################################################################

load('./data/Glasgow_data/Glasgow-friendship.RData')
load('./data/Glasgow_data/Glasgow-substances.RData')
load('./data/Glasgow_data/Glasgow-demographic.RData')

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
tobacco <- tobacco[not_missing, ]

sex.F <- sex.F[not_missing]

#################################################################################
##################   Select only boys ###########################################

friendship.1 <- friendship.1[sex.F == 0, sex.F == 0]
friendship.2 <- friendship.2[sex.F == 0, sex.F == 0]
friendship.3 <- friendship.3[sex.F == 0, sex.F == 0]

alcohol <- alcohol[sex.F == 0, ]
tobacco <- tobacco[sex.F == 0, ]
