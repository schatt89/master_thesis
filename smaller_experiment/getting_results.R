################################################################################
##############                                                     #############
##############                 general                             #############
##############                                                     #############
################################################################################

npar <- sum(effects.imputed$include)

rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

load("./data/results/result-20-mcar.RData")
load("./data/results/result-30-mcar.RData")
load("./data/results/result-20-b.RData")
load("./data/results/result-30-b.RData")
load("./data/results/result-20-n.RData")
load("./data/results/result-30-n.RData")
load("./data/results/result-20-nb.RData")
load("./data/results/result-30-nb.RData")

################################################################################
##############                                                     #############
##############                 dataset 20 b                        #############
##############                                                     #############
################################################################################

library(tidyr)
D = 50

N = length(saom.results.20.b.t1)

thetas.t1 <- c()
covthetas.t1 <- c()
thetas.t2 <- c()
covthetas.t2 <- c()
for (i in 1:length(saom.results.20.b.t1)) {
  all_thetas.t1 = saom.results.20.b.t1[[i]][[1]]
  all_covthetas.t1 = saom.results.20.b.t1[[i]][[2]]
  all_thetas.t2 = saom.results.20.b.t2[[i]][[1]]
  all_covthetas.t2 = saom.results.20.b.t2[[i]][[2]]
  for (j in 1:D) {
    thetas.t1 <- c(thetas.t1, list(all_thetas.t1[[j]]))
    covthetas.t1 <- c(covthetas.t1, list(all_covthetas.t1[[j]]))
    thetas.t2 <- c(thetas.t2, list(all_thetas.t2[[j]]))
    covthetas.t2 <- c(covthetas.t2, list(all_covthetas.t2[[j]]))
  }
}

SE.t1 = map(covthetas.t1, ~ sqrt(diag(.x)))
SE.t2 = map(covthetas.t2, ~ sqrt(diag(.x)))

bad.SE.t1 = rep(FALSE, length(SE.t1))
for(i in 1:length(SE.t1)) {
  bad.SE.t1[i] = sum(SE.t1[[i]] > 10) > 0
}

bad.SE.t2 = rep(FALSE, length(SE.t2))
for(i in 1:length(SE.t2)) {
  bad.SE.t2[i] = sum(SE.t2[[i]] > 10) > 0
}

thetas.t1 = thetas.t1[!bad.SE.t1]
thetas.t2 = thetas.t2[!bad.SE.t2]

covthetas.t1 = covthetas.t1[!bad.SE.t1]
covthetas.t2 = covthetas.t2[!bad.SE.t2]


################################################################################
##############                                                     #############
##############                 for theta 1, 20 b                   #############
##############                                                     #############
################################################################################

D = length(thetas.t1)

MIResults <- as.data.frame(matrix(NA,npar,(2 * D)))

for (i in 1:D) {
  names(MIResults)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(MIResults)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  MIResults[,i * 2 - 1] <- thetas.t1[[i]]  # estimates
  MIResults[,i * 2] <-  sqrt(diag(covthetas.t1[[i]]))  # standard errors
}

# Now we get the average covariance structure between the parameters
WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + covthetas.t1[[i]]
}

WDMIs <- (1/D) * WDMIs

# Using Rubin's Rules we combine the parameters and standard errors and
# complete the procedure
res20mis.b.t1 <- as.data.frame(matrix(NA,npar,2))
names(res20mis.b.t1) <- c("20misB.t1.Est.", "20misB.t1.SE")
rownames(res20mis.b.t1) <- effects.imputed$effectName[effects.imputed$include]

res20mis.b.t1[,1]<- rowMeans(MIResults[,seq(1,2*D,2)])
res20mis.b.t1[,2] <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                                  rowVar(MIResults[,seq(1,2*D,2)]))

################################################################################
##############                                                     #############
##############                 for theta 2, 20 b                   #############
##############                                                     #############
################################################################################

D = length(thetas.t2)

MIResults <- as.data.frame(matrix(NA,npar,(2 * D)))

for (i in 1:D) {
  names(MIResults)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(MIResults)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  MIResults[,i * 2 - 1] <- thetas.t2[[i]]  # estimates
  MIResults[,i * 2] <-  sqrt(diag(covthetas.t2[[i]]))  # standard errors
}

# Now we get the average covariance structure between the parameters
WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + covthetas.t2[[i]]
}

WDMIs <- (1/D) * WDMIs

# Using Rubin's Rules we combine the parameters and standard errors and
# complete the procedure
res20mis.b.t2 <- as.data.frame(matrix(NA,npar,2))
names(res20mis.b.t2) <- c("20misB.t2.Est.", "20misB.t2.SE")
rownames(res20mis.b.t2) <- effects.imputed$effectName[effects.imputed$include]

res20mis.b.t2[,1]<- rowMeans(MIResults[,seq(1,2*D,2)])
res20mis.b.t2[,2] <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                           rowVar(MIResults[,seq(1,2*D,2)]))

res20mis.b <- cbind(res20mis.b.t1, res20mis.b.t2)

################################################################################
##############                                                     #############
##############                 dataset 30 b                        #############
##############                                                     #############
################################################################################

D = 50

N = length(saom.results.30.b.t1)

thetas.t1 <- c()
covthetas.t1 <- c()
thetas.t2 <- c()
covthetas.t2 <- c()
for (i in 1:length(saom.results.30.b.t1)) {
  all_thetas.t1 = saom.results.30.b.t1[[i]][[1]]
  all_covthetas.t1 = saom.results.30.b.t1[[i]][[2]]
  all_thetas.t2 = saom.results.30.b.t2[[i]][[1]]
  all_covthetas.t2 = saom.results.30.b.t2[[i]][[2]]
  for (j in 1:D) {
    thetas.t1 <- c(thetas.t1, list(all_thetas.t1[[j]]))
    covthetas.t1 <- c(covthetas.t1, list(all_covthetas.t1[[j]]))
    thetas.t2 <- c(thetas.t2, list(all_thetas.t2[[j]]))
    covthetas.t2 <- c(covthetas.t2, list(all_covthetas.t2[[j]]))
  }
}

SE.t1 = map(covthetas.t1, ~ sqrt(diag(.x)))
SE.t2 = map(covthetas.t2, ~ sqrt(diag(.x)))

bad.SE.t1 = rep(FALSE, length(SE.t1))
for(i in 1:length(SE.t1)) {
  bad.SE.t1[i] = sum(SE.t1[[i]] > 10) > 0
}

bad.SE.t2 = rep(FALSE, length(SE.t2))
for(i in 1:length(SE.t2)) {
  bad.SE.t2[i] = sum(SE.t2[[i]] > 10) > 0
}

thetas.t1 = thetas.t1[!bad.SE.t1]
thetas.t2 = thetas.t2[!bad.SE.t2]

covthetas.t1 = covthetas.t1[!bad.SE.t1]
covthetas.t2 = covthetas.t2[!bad.SE.t2]


################################################################################
##############                                                     #############
##############                 for theta 1, 30 b                   #############
##############                                                     #############
################################################################################

D = length(thetas.t1)

MIResults <- as.data.frame(matrix(NA,npar,(2 * D)))

for (i in 1:D) {
  names(MIResults)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(MIResults)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  MIResults[,i * 2 - 1] <- thetas.t1[[i]]  # estimates
  MIResults[,i * 2] <-  sqrt(diag(covthetas.t1[[i]]))  # standard errors
}

# Now we get the average covariance structure between the parameters
WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + covthetas.t1[[i]]
}

WDMIs <- (1/D) * WDMIs

# Using Rubin's Rules we combine the parameters and standard errors and
# complete the procedure
res30mis.b.t1 <- as.data.frame(matrix(NA,npar,2))
names(res30mis.b.t1) <- c("30misB.t1.Est.", "30misB.t1.SE")
rownames(res30mis.b.t1) <- effects.imputed$effectName[effects.imputed$include]

res30mis.b.t1[,1]<- rowMeans(MIResults[,seq(1,2*D,2)])
res30mis.b.t1[,2] <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                            rowVar(MIResults[,seq(1,2*D,2)]))

################################################################################
##############                                                     #############
##############                 for theta 2, 30 b                   #############
##############                                                     #############
################################################################################

D = length(thetas.t2)

MIResults <- as.data.frame(matrix(NA,npar,(2 * D)))

for (i in 1:D) {
  names(MIResults)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(MIResults)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  MIResults[,i * 2 - 1] <- thetas.t2[[i]]  # estimates
  MIResults[,i * 2] <-  sqrt(diag(covthetas.t2[[i]]))  # standard errors
}

# Now we get the average covariance structure between the parameters
WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + covthetas.t2[[i]]
}

WDMIs <- (1/D) * WDMIs

# Using Rubin's Rules we combine the parameters and standard errors and
# complete the procedure
res30mis.b.t2 <- as.data.frame(matrix(NA,npar,2))
names(res30mis.b.t2) <- c("30misB.t2.Est.", "30misB.t2.SE")
rownames(res30mis.b.t2) <- effects.imputed$effectName[effects.imputed$include]

res30mis.b.t2[,1]<- rowMeans(MIResults[,seq(1,2*D,2)])
res30mis.b.t2[,2] <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                            rowVar(MIResults[,seq(1,2*D,2)]))

res30mis.b <- cbind(res30mis.b.t1, res30mis.b.t2)

################################################################################
##############                                                     #############
##############                 dataset 20 n                        #############
##############                                                     #############
################################################################################

D = 50

N = length(saom.results.20.n.t1)

thetas.t1 <- c()
covthetas.t1 <- c()
thetas.t2 <- c()
covthetas.t2 <- c()
for (i in 1:length(saom.results.20.n.t1)) {
  all_thetas.t1 = saom.results.20.n.t1[[i]][[1]]
  all_covthetas.t1 = saom.results.20.n.t1[[i]][[2]]
  all_thetas.t2 = saom.results.20.n.t2[[i]][[1]]
  all_covthetas.t2 = saom.results.20.n.t2[[i]][[2]]
  for (j in 1:D) {
    thetas.t1 <- c(thetas.t1, list(all_thetas.t1[[j]]))
    covthetas.t1 <- c(covthetas.t1, list(all_covthetas.t1[[j]]))
    thetas.t2 <- c(thetas.t2, list(all_thetas.t2[[j]]))
    covthetas.t2 <- c(covthetas.t2, list(all_covthetas.t2[[j]]))
  }
}

SE.t1 = map(covthetas.t1, ~ sqrt(diag(.x)))
SE.t2 = map(covthetas.t2, ~ sqrt(diag(.x)))

bad.SE.t1 = rep(FALSE, length(SE.t1))
for(i in 1:length(SE.t1)) {
  bad.SE.t1[i] = sum(SE.t1[[i]] > 10) > 0
}

bad.SE.t2 = rep(FALSE, length(SE.t2))
for(i in 1:length(SE.t2)) {
  bad.SE.t2[i] = sum(SE.t2[[i]] > 10) > 0
}

thetas.t1 = thetas.t1[!bad.SE.t1]
thetas.t2 = thetas.t2[!bad.SE.t2]

covthetas.t1 = covthetas.t1[!bad.SE.t1]
covthetas.t2 = covthetas.t2[!bad.SE.t2]


################################################################################
##############                                                     #############
##############                 for theta 1, 20 n                   #############
##############                                                     #############
################################################################################

D = length(thetas.t1)

MIResults <- as.data.frame(matrix(NA,npar,(2 * D)))

for (i in 1:D) {
  names(MIResults)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(MIResults)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  MIResults[,i * 2 - 1] <- thetas.t1[[i]]  # estimates
  MIResults[,i * 2] <-  sqrt(diag(covthetas.t1[[i]]))  # standard errors
}

# Now we get the average covariance structure between the parameters
WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + covthetas.t1[[i]]
}

WDMIs <- (1/D) * WDMIs

# Using Rubin's Rules we combine the parameters and standard errors and
# complete the procedure
res20mis.n.t1 <- as.data.frame(matrix(NA,npar,2))
names(res20mis.n.t1) <- c("20misN.t1.Est.", "20misN.t1.SE")
rownames(res20mis.n.t1) <- effects.imputed$effectName[effects.imputed$include]

res20mis.n.t1[,1]<- rowMeans(MIResults[,seq(1,2*D,2)])
res20mis.n.t1[,2] <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                            rowVar(MIResults[,seq(1,2*D,2)]))

################################################################################
##############                                                     #############
##############                 for theta 2, 20 n                   #############
##############                                                     #############
################################################################################

D = length(thetas.t2)

MIResults <- as.data.frame(matrix(NA,npar,(2 * D)))

for (i in 1:D) {
  names(MIResults)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(MIResults)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  MIResults[,i * 2 - 1] <- thetas.t2[[i]]  # estimates
  MIResults[,i * 2] <-  sqrt(diag(covthetas.t2[[i]]))  # standard errors
}

# Now we get the average covariance structure between the parameters
WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + covthetas.t2[[i]]
}

WDMIs <- (1/D) * WDMIs

# Using Rubin's Rules we combine the parameters and standard errors and
# complete the procedure
res20mis.n.t2 <- as.data.frame(matrix(NA,npar,2))
names(res20mis.n.t2) <- c("20misN.t2.Est.", "20misN.t2.SE")
rownames(res20mis.n.t2) <- effects.imputed$effectName[effects.imputed$include]

res20mis.n.t2[,1]<- rowMeans(MIResults[,seq(1,2*D,2)])
res20mis.n.t2[,2] <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                            rowVar(MIResults[,seq(1,2*D,2)]))

res20mis.n <- cbind(res20mis.n.t1, res20mis.n.t2)



################################################################################
##############                                                     #############
##############                 dataset 30 n                        #############
##############                                                     #############
################################################################################

D = 50

N = length(saom.results.30.n.t1)

thetas.t1 <- c()
covthetas.t1 <- c()
thetas.t2 <- c()
covthetas.t2 <- c()
for (i in 1:length(saom.results.30.n.t1)) {
  all_thetas.t1 = saom.results.30.n.t1[[i]][[1]]
  all_covthetas.t1 = saom.results.30.n.t1[[i]][[2]]
  all_thetas.t2 = saom.results.30.n.t2[[i]][[1]]
  all_covthetas.t2 = saom.results.30.n.t2[[i]][[2]]
  for (j in 1:D) {
    thetas.t1 <- c(thetas.t1, list(all_thetas.t1[[j]]))
    covthetas.t1 <- c(covthetas.t1, list(all_covthetas.t1[[j]]))
    thetas.t2 <- c(thetas.t2, list(all_thetas.t2[[j]]))
    covthetas.t2 <- c(covthetas.t2, list(all_covthetas.t2[[j]]))
  }
}

SE.t1 = map(covthetas.t1, ~ sqrt(diag(.x)))
SE.t2 = map(covthetas.t2, ~ sqrt(diag(.x)))

bad.SE.t1 = rep(FALSE, length(SE.t1))
for(i in 1:length(SE.t1)) {
  bad.SE.t1[i] = sum(SE.t1[[i]] > 10) > 0
}

bad.SE.t2 = rep(FALSE, length(SE.t2))
for(i in 1:length(SE.t2)) {
  bad.SE.t2[i] = sum(SE.t2[[i]] > 10) > 0
}

thetas.t1 = thetas.t1[!bad.SE.t1]
thetas.t2 = thetas.t2[!bad.SE.t2]

covthetas.t1 = covthetas.t1[!bad.SE.t1]
covthetas.t2 = covthetas.t2[!bad.SE.t2]


################################################################################
##############                                                     #############
##############                 for theta 1, 30 n                   #############
##############                                                     #############
################################################################################

D = length(thetas.t1)

MIResults <- as.data.frame(matrix(NA,npar,(2 * D)))

for (i in 1:D) {
  names(MIResults)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(MIResults)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  MIResults[,i * 2 - 1] <- thetas.t1[[i]]  # estimates
  MIResults[,i * 2] <-  sqrt(diag(covthetas.t1[[i]]))  # standard errors
}

# Now we get the average covariance structure between the parameters
WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + covthetas.t1[[i]]
}

WDMIs <- (1/D) * WDMIs

# Using Rubin's Rules we combine the parameters and standard errors and
# complete the procedure
res30mis.n.t1 <- as.data.frame(matrix(NA,npar,2))
names(res30mis.n.t1) <- c("30misN.t1.Est.", "30misN.t1.SE")
rownames(res30mis.n.t1) <- effects.imputed$effectName[effects.imputed$include]

res30mis.n.t1[,1]<- rowMeans(MIResults[,seq(1,2*D,2)])
res30mis.n.t1[,2] <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                            rowVar(MIResults[,seq(1,2*D,2)]))

################################################################################
##############                                                     #############
##############                 for theta 2, 30 n                   #############
##############                                                     #############
################################################################################

D = length(thetas.t2)

MIResults <- as.data.frame(matrix(NA,npar,(2 * D)))

for (i in 1:D) {
  names(MIResults)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(MIResults)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  MIResults[,i * 2 - 1] <- thetas.t2[[i]]  # estimates
  MIResults[,i * 2] <-  sqrt(diag(covthetas.t2[[i]]))  # standard errors
}

# Now we get the average covariance structure between the parameters
WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + covthetas.t2[[i]]
}

WDMIs <- (1/D) * WDMIs

# Using Rubin's Rules we combine the parameters and standard errors and
# complete the procedure
res30mis.n.t2 <- as.data.frame(matrix(NA,npar,2))
names(res30mis.n.t2) <- c("30misN.t2.Est.", "30misN.t2.SE")
rownames(res30mis.n.t2) <- effects.imputed$effectName[effects.imputed$include]

res30mis.n.t2[,1]<- rowMeans(MIResults[,seq(1,2*D,2)])
res30mis.n.t2[,2] <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                            rowVar(MIResults[,seq(1,2*D,2)]))

res30mis.n <- cbind(res30mis.n.t1, res30mis.n.t2)

################################################################################
##############                                                     #############
##############                 dataset 20 nb                       #############
##############                                                     #############
################################################################################

D = 50

N = length(saom.results.20.nb.t1)

thetas.t1 <- c()
covthetas.t1 <- c()
thetas.t2 <- c()
covthetas.t2 <- c()
for (i in 1:length(saom.results.20.nb.t1)) {
  all_thetas.t1 = saom.results.20.nb.t1[[i]][[1]]
  all_covthetas.t1 = saom.results.20.nb.t1[[i]][[2]]
  all_thetas.t2 = saom.results.20.nb.t2[[i]][[1]]
  all_covthetas.t2 = saom.results.20.nb.t2[[i]][[2]]
  for (j in 1:D) {
    thetas.t1 <- c(thetas.t1, list(all_thetas.t1[[j]]))
    covthetas.t1 <- c(covthetas.t1, list(all_covthetas.t1[[j]]))
    thetas.t2 <- c(thetas.t2, list(all_thetas.t2[[j]]))
    covthetas.t2 <- c(covthetas.t2, list(all_covthetas.t2[[j]]))
  }
}

SE.t1 = map(covthetas.t1, ~ sqrt(diag(.x)))
SE.t2 = map(covthetas.t2, ~ sqrt(diag(.x)))

bad.SE.t1 = rep(FALSE, length(SE.t1))
for(i in 1:length(SE.t1)) {
  bad.SE.t1[i] = sum(SE.t1[[i]] > 10) > 0
}

bad.SE.t2 = rep(FALSE, length(SE.t2))
for(i in 1:length(SE.t2)) {
  bad.SE.t2[i] = sum(SE.t2[[i]] > 10) > 0
}

thetas.t1 = thetas.t1[!bad.SE.t1]
thetas.t2 = thetas.t2[!bad.SE.t2]

covthetas.t1 = covthetas.t1[!bad.SE.t1]
covthetas.t2 = covthetas.t2[!bad.SE.t2]


################################################################################
##############                                                     #############
##############                 for theta 1, 20 nb                  #############
##############                                                     #############
################################################################################

D = length(thetas.t1)

MIResults <- as.data.frame(matrix(NA,npar,(2 * D)))

for (i in 1:D) {
  names(MIResults)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(MIResults)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  MIResults[,i * 2 - 1] <- thetas.t1[[i]]  # estimates
  MIResults[,i * 2] <-  sqrt(diag(covthetas.t1[[i]]))  # standard errors
}

# Now we get the average covariance structure between the parameters
WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + covthetas.t1[[i]]
}

WDMIs <- (1/D) * WDMIs

# Using Rubin's Rules we combine the parameters and standard errors and
# complete the procedure
res20mis.nb.t1 <- as.data.frame(matrix(NA,npar,2))
names(res20mis.nb.t1) <- c("20misNB.t1.Est.", "20misNB.t1.SE")
rownames(res20mis.nb.t1) <- effects.imputed$effectName[effects.imputed$include]

res20mis.nb.t1[,1]<- rowMeans(MIResults[,seq(1,2*D,2)])
res20mis.nb.t1[,2] <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                            rowVar(MIResults[,seq(1,2*D,2)]))

################################################################################
##############                                                     #############
##############                 for theta 2, 20 nb                  #############
##############                                                     #############
################################################################################

D = length(thetas.t2)

MIResults <- as.data.frame(matrix(NA,npar,(2 * D)))

for (i in 1:D) {
  names(MIResults)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(MIResults)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  MIResults[,i * 2 - 1] <- thetas.t2[[i]]  # estimates
  MIResults[,i * 2] <-  sqrt(diag(covthetas.t2[[i]]))  # standard errors
}

# Now we get the average covariance structure between the parameters
WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + covthetas.t2[[i]]
}

WDMIs <- (1/D) * WDMIs

# Using Rubin's Rules we combine the parameters and standard errors and
# complete the procedure
res20mis.nb.t2 <- as.data.frame(matrix(NA,npar,2))
names(res20mis.nb.t2) <- c("20misNB.t2.Est.", "20misNB.t2.SE")
rownames(res20mis.nb.t2) <- effects.imputed$effectName[effects.imputed$include]

res20mis.nb.t2[,1]<- rowMeans(MIResults[,seq(1,2*D,2)])
res20mis.nb.t2[,2] <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                            rowVar(MIResults[,seq(1,2*D,2)]))

res20mis.nb <- cbind(res20mis.nb.t1, res20mis.nb.t2)


################################################################################
##############                                                     #############
##############                 dataset 30 nb                       #############
##############                                                     #############
################################################################################

D = 50

N = length(saom.results.30.nb.t1)

thetas.t1 <- c()
covthetas.t1 <- c()
thetas.t2 <- c()
covthetas.t2 <- c()
for (i in 1:length(saom.results.30.nb.t1)) {
  all_thetas.t1 = saom.results.30.nb.t1[[i]][[1]]
  all_covthetas.t1 = saom.results.30.nb.t1[[i]][[2]]
  all_thetas.t2 = saom.results.30.nb.t2[[i]][[1]]
  all_covthetas.t2 = saom.results.30.nb.t2[[i]][[2]]
  for (j in 1:D) {
    thetas.t1 <- c(thetas.t1, list(all_thetas.t1[[j]]))
    covthetas.t1 <- c(covthetas.t1, list(all_covthetas.t1[[j]]))
    thetas.t2 <- c(thetas.t2, list(all_thetas.t2[[j]]))
    covthetas.t2 <- c(covthetas.t2, list(all_covthetas.t2[[j]]))
  }
}

SE.t1 = map(covthetas.t1, ~ sqrt(diag(.x)))
SE.t2 = map(covthetas.t2, ~ sqrt(diag(.x)))

bad.SE.t1 = rep(FALSE, length(SE.t1))
for(i in 1:length(SE.t1)) {
  bad.SE.t1[i] = sum(SE.t1[[i]] > 10) > 0
}

bad.SE.t2 = rep(FALSE, length(SE.t2))
for(i in 1:length(SE.t2)) {
  bad.SE.t2[i] = sum(SE.t2[[i]] > 10) > 0
}

thetas.t1 = thetas.t1[!bad.SE.t1]
thetas.t2 = thetas.t2[!bad.SE.t2]

covthetas.t1 = covthetas.t1[!bad.SE.t1]
covthetas.t2 = covthetas.t2[!bad.SE.t2]


################################################################################
##############                                                     #############
##############                 for theta 1, 30 nb                  #############
##############                                                     #############
################################################################################

D = length(thetas.t1)

MIResults <- as.data.frame(matrix(NA,npar,(2 * D)))

for (i in 1:D) {
  names(MIResults)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(MIResults)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  MIResults[,i * 2 - 1] <- thetas.t1[[i]]  # estimates
  MIResults[,i * 2] <-  sqrt(diag(covthetas.t1[[i]]))  # standard errors
}

# Now we get the average covariance structure between the parameters
WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + covthetas.t1[[i]]
}

WDMIs <- (1/D) * WDMIs

# Using Rubin's Rules we combine the parameters and standard errors and
# complete the procedure
res30mis.nb.t1 <- as.data.frame(matrix(NA,npar,2))
names(res30mis.nb.t1) <- c("30misNB.t1.Est.", "30misNB.t1.SE")
rownames(res30mis.nb.t1) <- effects.imputed$effectName[effects.imputed$include]

res30mis.nb.t1[,1]<- rowMeans(MIResults[,seq(1,2*D,2)])
res30mis.nb.t1[,2] <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                             rowVar(MIResults[,seq(1,2*D,2)]))

################################################################################
##############                                                     #############
##############                 for theta 2, 30 nb                  #############
##############                                                     #############
################################################################################

D = length(thetas.t2)

MIResults <- as.data.frame(matrix(NA,npar,(2 * D)))

for (i in 1:D) {
  names(MIResults)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(MIResults)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  MIResults[,i * 2 - 1] <- thetas.t2[[i]]  # estimates
  MIResults[,i * 2] <-  sqrt(diag(covthetas.t2[[i]]))  # standard errors
}

# Now we get the average covariance structure between the parameters
WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + covthetas.t2[[i]]
}

WDMIs <- (1/D) * WDMIs

# Using Rubin's Rules we combine the parameters and standard errors and
# complete the procedure
res30mis.nb.t2 <- as.data.frame(matrix(NA,npar,2))
names(res30mis.nb.t2) <- c("30misNB.t2.Est.", "30misNB.t2.SE")
rownames(res30mis.nb.t2) <- effects.imputed$effectName[effects.imputed$include]

res30mis.nb.t2[,1]<- rowMeans(MIResults[,seq(1,2*D,2)])
res30mis.nb.t2[,2] <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                             rowVar(MIResults[,seq(1,2*D,2)]))

res30mis.nb <- cbind(res30mis.nb.t1, res30mis.nb.t2)


################################################################################
##############                                                     #############
##############                 dataset 20 mcar                     #############
##############                                                     #############
################################################################################

D = 50
N = length(saom.results.20.mcar)

thetas <- c()
covthetas <- c()
for (i in 1:length(saom.results.20.mcar)) {
  all_thetas = saom.results.20.mcar[[i]][[1]]
  all_covthetas = saom.results.20.mcar[[i]][[2]]
  for (j in 1:D) {
    thetas <- c(thetas, list(all_thetas[[j]]))
    covthetas <- c(covthetas, list(all_covthetas[[j]]))
  }
}

SE = map(covthetas, ~ sqrt(diag(.x)))

bad.SE = rep(FALSE, length(SE))
for(i in 1:length(SE)) {
  bad.SE[i] = sum(SE[[i]] > 10) > 0
}

thetas = thetas[!bad.SE]
covthetas = covthetas[!bad.SE]


################################################################################
##############                                                     #############
##############                      20 mcar                        #############
##############                                                     #############
################################################################################

D = length(thetas)

MIResults <- as.data.frame(matrix(NA,npar,(2 * D)))

for (i in 1:D) {
  names(MIResults)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(MIResults)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  MIResults[,i * 2 - 1] <- thetas.t1[[i]]  # estimates
  MIResults[,i * 2] <-  sqrt(diag(covthetas.t1[[i]]))  # standard errors
}

# Now we get the average covariance structure between the parameters
WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + covthetas.t1[[i]]
}

WDMIs <- (1/D) * WDMIs

# Using Rubin's Rules we combine the parameters and standard errors and
# complete the procedure
res20mis.mcar <- as.data.frame(matrix(NA,npar,2))
names(res20mis.mcar) <- c("20mis.mcar.Est.", "20mis.mcar.SE")
rownames(res20mis.mcar) <- effects.imputed$effectName[effects.imputed$include]

res20mis.mcar[,1]<- rowMeans(MIResults[,seq(1,2*D,2)])
res20mis.mcar[,2] <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                             rowVar(MIResults[,seq(1,2*D,2)]))

################################################################################
##############                                                     #############
##############                 dataset 30 mcar                     #############
##############                                                     #############
################################################################################

D = 50
N = length(saom.results.30.mcar)

thetas <- c()
covthetas <- c()
for (i in 1:length(saom.results.30.mcar)) {
  all_thetas = saom.results.30.mcar[[i]][[1]]
  all_covthetas = saom.results.30.mcar[[i]][[2]]
  for (j in 1:D) {
    thetas <- c(thetas, list(all_thetas[[j]]))
    covthetas <- c(covthetas, list(all_covthetas[[j]]))
  }
}

SE = map(covthetas, ~ sqrt(diag(.x)))

bad.SE = rep(FALSE, length(SE))
for(i in 1:length(SE)) {
  bad.SE[i] = sum(SE[[i]] > 10) > 0
}

thetas = thetas[!bad.SE]
covthetas = covthetas[!bad.SE]


################################################################################
##############                                                     #############
##############                      30 mcar                        #############
##############                                                     #############
################################################################################

D = length(thetas)

MIResults <- as.data.frame(matrix(NA,npar,(2 * D)))

for (i in 1:D) {
  names(MIResults)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(MIResults)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  MIResults[,i * 2 - 1] <- thetas.t1[[i]]  # estimates
  MIResults[,i * 2] <-  sqrt(diag(covthetas.t1[[i]]))  # standard errors
}

# Now we get the average covariance structure between the parameters
WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + covthetas.t1[[i]]
}

WDMIs <- (1/D) * WDMIs

# Using Rubin's Rules we combine the parameters and standard errors and
# complete the procedure
res30mis.mcar <- as.data.frame(matrix(NA,npar,2))
names(res30mis.mcar) <- c("30mis.mcar.Est.", "30mis.mcar.SE")
rownames(res30mis.mcar) <- effects.imputed$effectName[effects.imputed$include]

res30mis.mcar[,1]<- rowMeans(MIResults[,seq(1,2*D,2)])
res30mis.mcar[,2] <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                            rowVar(MIResults[,seq(1,2*D,2)]))

res.mcar <- cbind(res20mis.mcar, res30mis.mcar)


################################################################################
##############                                                     #############
##############                 final results                       #############
##############                                                     #############
################################################################################


finalResults <- cbind(res.mcar, res20mis.n, res30mis.n, res20mis.b, res30mis.b,
                      res20mis.nb, res30mis.nb)


################################################################################
##############                                                     #############
##############                 plotting results                    #############
##############                                                     #############
################################################################################


finalResults$parameter <- c("Friend rate 1",
                            "Density",
                            "Reciprocity",
                            "gwespFF",
                            "Indeg. Pop. sqrt.",
                            "Outdeg. Act. sqrt.",
                            "Alter drinking",
                            "Ego drinking",
                            "Ego x alt drinking",
                            "Drinking rate 1",
                            "Drink. linear",
                            "Drink. quadratic",
                            "Avg. Alter Drink.")

ld <- as.data.frame(matrix(NA,13*14,2))

ld[,1] <- rep(finalResults$parameter,2)
ld[1:13,c(2,3)] <- finalResults[,c(1,2)]
ld[14:26,c(2,3)] <- finalResults[,c(3,4)]
ld[27:39,c(2,3)] <- finalResults[,c(5,6)]
ld[40:52,c(2,3)] <- finalResults[,c(7,8)]
ld[53:65,c(2,3)] <- finalResults[,c(9,10)]
ld[66:78,c(2,3)] <- finalResults[,c(11,12)]
ld[79:91,c(2,3)] <- finalResults[,c(13,14)]
ld[92:104,c(2,3)] <- finalResults[,c(15,16)]
ld[105:117,c(2,3)] <- finalResults[,c(17,18)]
ld[118:130,c(2,3)] <- finalResults[,c(19,20)]
ld[131:143,c(2,3)] <- finalResults[,c(21,22)]
ld[144:156,c(2,3)] <- finalResults[,c(23,24)]
ld[157:169,c(2,3)] <- finalResults[,c(25,26)]
ld[170:182,c(2,3)] <- finalResults[,c(27,28)]

ld[,4] <- c(rep("mis20mcar",13),
            rep("mis30mcar",13),
            rep("mis20n.t1",13),
            rep("mis20n.t2",13),
            rep("mis30n.t1",13),
            rep("mis30n.t2",13),
            rep("mis20b.t1",13),
            rep("mis20b.t2",13),
            rep("mis30b.t1",13),
            rep("mis30b.t2",13),
            rep("mis20nb.t1",13),
            rep("mis20nb.t2",13),
            rep("mis30nb.t1",13),
            rep("mis30nb.t2",13)
)

names(ld) <- c("parameter","Theta","SE","model")

ld$model = factor(ld$model, levels = c(unique(ld[,4]))[14:1])

ld$parameter = factor(ld$parameter,
                      levels = c("Friend rate 1",
                                 "Density",
                                 "Reciprocity",
                                 "gwespFF",
                                 "Indeg. Pop. sqrt.",
                                 "Outdeg. Act. sqrt.",
                                 "Alter drinking",
                                 "Ego drinking",
                                 "Ego x alt drinking",
                                 "Drinking rate 1",
                                 "Drink. linear",
                                 "Drink. quadratic",
                                 "Avg. Alter Drink.")[13:1])

g1 <- c("Friend rate 1",
        "Density",
        "Reciprocity",
        "gwespFF",
        "Indeg. Pop. sqrt.",
        "Outdeg. Act. sqrt.")


g2 <- c("Alter drinking",
       "Ego drinking",
       "Ego x alt drinking")

g3 <- c("Drinking rate 1",
        "Drink. linear",
        "Drink. quadratic",
        "Avg. Alter Drink.")

ld$g[ld$parameter %in% g1] <- "Structure"
ld$g[ld$parameter %in% g2] <- "Selection"
ld$g[ld$parameter %in% g3] <- "Behavior"

ld$g <- factor(ld$g, levels = c("Structure","Selection","Behavior"))

library(ggplot2)
p6 <- ggplot(ld, aes(x = parameter, y = Theta, group = model,
                     colour = model)) +
  facet_wrap(facets = ~g, scales = "free", nrow = 1) +
  geom_point( size = 2.2 , position = position_dodge(width = .75), shape=18) +
  geom_hline(yintercept = 0, size = .25, color = "black") +
  geom_errorbar(size = .7, aes(ymin = Theta - SE*2, ymax = Theta + SE*2),
                position = position_dodge(width = .75)) +
  theme_classic() + theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.text.y = element_text(angle = 0, hjust = 1),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.height = unit(3,"line"), legend.key.width = unit(3,"line"),
    text = element_text(size = 20),
    axis.title.x = element_blank()) + ylab("") + xlab("") + coord_flip()
p6
