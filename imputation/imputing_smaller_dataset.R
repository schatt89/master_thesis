library(RSiena) # or RSienaTest
library("mice")
source("./simulation/siena07ToConvergence.R")

getNet <- function(observedNet,edgeList) {
  # observedNet = observed network as adjacency with missing data
  # edgeList = edgelist that is returned by siena07(...)$sims
  observedNet[is.na(observedNet)] <- 0
  for (i in 1:nrow(edgeList)) {
    observedNet[edgeList[i,1],edgeList[i,2]] <- 1
  }
  return(observedNet)
}

load("./data/simulated/Data30.RData")


################################################################################
#########                                                             ##########
#########                   Imputing behavior with mice               ##########
#########                                                             ##########
################################################################################

miceImpAlco30.10.n <- array(rep(NA, 30*50*100), c(30, 50, 100))

for (i in 1:100) {
  indegree1 <- colSums(fr.30.1.mis.10.n, na.rm = TRUE)
  indegree2 <- colSums(fr.30.2.sim.mis.10.n[,,i], na.rm = TRUE)
  indegree3 <- colSums(fr.30.3.sim.mis.10.n[,,i], na.rm = TRUE)
  
  avgAltA1 <- rowSums(sweep(t(fr.30.1.mis.10.n),
                            MARGIN = 2, alco.30.1.mis.10.n,'*'), na.rm = TRUE) /
                            rowSums(t(fr.30.1.mis.10.n), na.rm = TRUE)
  
  avgAltA2 <- rowSums(sweep(t(fr.30.2.sim.mis.10.n[,,i]),
                    MARGIN = 2, alco.30.sim.mis.10.n[,1,i],'*'), na.rm = TRUE) /
                            rowSums(t(fr.30.2.sim.mis.10.n[,,i]), na.rm = TRUE)
  
  avgAltA3 <- rowSums(sweep(t(fr.30.3.sim.mis.10.n[,,i]),
                     MARGIN = 2, alco.30.sim.mis.10.n[,2,i],'*'),na.rm = TRUE) /
                            rowSums(t(fr.30.3.sim.mis.10.n[,,i]), na.rm = TRUE)
  
  avgAltA1[is.nan(avgAltA1)] <- NA
  avgAltA2[is.nan(avgAltA2)] <- NA
  avgAltA3[is.nan(avgAltA3)] <- NA
  
  miceData <- cbind(alco.30.1.mis.10.n, alco.30.sim.mis.10.n[,,i], indegree1,
                    indegree2, indegree3, avgAltA1, avgAltA2, avgAltA3)
  
  set.seed(11019)
  miceImp <- mice(miceData, m = 50, defaultMethod = "pmm", maxit = 20)
  for (d in 1:50) {
    miceImpAlco30.10.n[,d,i] <- complete(miceImp, d)$alco.30.1.mis.10.n
  }
  
  }
