library(tidyr)
D = 50

N = length(saom.results.20.b.t1)

load("./data/results/result-20-b.RData")

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
##############                general functions                    #############
##############                                                     #############
################################################################################

npar <- sum(effects.imputed$include)

rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

################################################################################
##############                                                     #############
##############                 for theta 1                         #############
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
##############                 for theta 2                         #############
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

################################################################################
##############                                                     #############
##############                 plotting results                    #############
##############                                                     #############
################################################################################

res20mis.b <- cbind(res20mis.b.t1, res20mis.b.t2)

res20mis.b$parameter <- c("Friend rate 1",
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

ld <- as.data.frame(matrix(NA,13*2,2))

ld[,1] <- rep(res20mis.b$parameter,2)
ld[1:13,c(2,3)] <- res20mis.b[,c(1,2)]
ld[14:26,c(2,3)] <- res20mis.b[,c(3,4)]

ld[,4] <- c(rep("mis20beh.t1",13),
            rep("mis20beh.t2",13)
)

names(ld) <- c("parameter","Theta","SE","model")

ld$model = factor(ld$model, levels = c("mis20beh.t1", "mis20beh.t2")[2:1])

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
p6 <- ggplot(ld, aes(x = parameter, y = Theta, group = model, colour = model,
                     shape = model)) +
  facet_wrap(facets = ~g, scales = "free", nrow = 1) +
  geom_point( size = 2.2 , position = position_dodge(width = .75)) +
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
