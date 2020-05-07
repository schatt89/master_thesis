load('./data/results/result-20-b.RData')

D = 50
# Now we have D RSiena results and all that is left is to combine them

# We need to extract all parameter and standard error estimates from the models
# and combine them using Rubin's Rules.

npar <- sum(effects.imputed$include)

rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

MIResults <- as.data.frame(matrix(NA,npar,(2 * D)))

for (i in 1:D) {
  names(MIResults)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(MIResults)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  MIResults[,i * 2 - 1] <- saom.results.20.b.t2[[1]][[1]][[i]]  # estimates
  MIResults[,i * 2] <-  sqrt(diag(saom.results.20.b.t2[[1]][[2]][[i]]))  # standard errors
}

# Now we get the average covariance structure between the parameters
WDMIs <- matrix(0,npar,npar)

for (i in 1:D) {
  WDMIs <- WDMIs + saom.results.20.b.t2[[1]][[2]][[i]]
}

WDMIs <- (1/D) * WDMIs

# Using Rubin's Rules we combine the parameters and standard errors and
# complete the procedure
finalResults <- as.data.frame(matrix(NA,npar,2))
names(finalResults) <- c("combinedEstimate", "combinedSE")
rownames(finalResults) <- effects.imputed$effectName[effects.imputed$include]

finalResults$combinedEstimate <- rowMeans(MIResults[,seq(1,2*D,2)])
finalResults$combinedSE <- sqrt(diag(WDMIs) + ((D + 1)/D) *
                                  rowVar(MIResults[,seq(1,2*D,2)]))

finalResults$completeTheta <- model.1$theta
finalResults$completeSE <- sqrt(diag(model.1$covtheta))

finalResults$defaultTheta <- model.2$theta
finalResults$defaultSE <- sqrt(diag(model.2$covtheta))


round(finalResults, 4)
# the names of user specified interactions have to be set manually

print(xtable(round(finalResults, 4), type = "latex"), file = "BdTdRound.tex")

save.image('mi.RData') # 100 MB
```



```{r}
###############################################################################
# proportion of variance due to imputation
BdTd <- as.data.frame(matrix(,npar,1))

names(BdTd) <- c("parameter")

BdTd$parameter <- finalResults$parameter

BdTd$MIsaom <- round(((51/50) * rowVar(MIResults[,seq(1,100,2)])) / (finalResults$combinedSE ** 2),2)

View(BdTd)

print(xtable(BdTd, type = "latex"), file = "BdTdRound.tex")

```


```{r}
################################################################################
#
# Plots
#
finalResults$parameter <- c("Friend rate 1",
                          "Density",
                          "Reciprocity",
                          "transTrip",
                          "transRecTrip",
                          "3-cycles",
                          "Indeg. Pop. sqrt.",
                          "Indeg. Act. sqrt.",
                          "Same sex",
                          "Alter drinking",
                          "Ego drinking",
                          "Ego x alt drinking",
                          "Alter smokimg",
                          "Ego smoking",
                          "Ego x alt smoking",
                          "Drinking rate 1",
                          "Drink. linear",
                          "Drink. quadratic",
                          "Avg. Alter Drink.",
                          "Drink: eff. from sex",
                          "Drink: eff. from smoke.",
                          "Smoking rate 1",
                          "Smoke. linear",
                          "Smoke. quadratic",
                          "Avg. Alter Smoke.",
                          "Smoke: eff. from sex",
                          "Smoke: eff. from drink.")



# lognformat - dirty

ld <- as.data.frame(matrix(NA,81,3))

ld[,1] <- rep(finalResults$parameter,3)
ld[1:27,c(2,3)] <- finalResults[,c(1,2)]
ld[28:54,c(2,3)] <- finalResults[,c(3,4)]
ld[55:81,c(2,3)] <- finalResults[,c(5,6)]


ld[,4] <- c(rep("MI-SAOM",27),
            rep("Complete Data",27),
            rep("Default SIENA",27)
            )

names(ld) <- c("parameter","Theta","SE","model")

ld$model = factor(ld$model, levels = c("MI-SAOM", "Default SIENA", "Complete Data" )[3:1]) # ,"Default SIENA"

ld$parameter = factor(ld$parameter,
                          levels = c("Friend rate 1",
                          "Density",
                          "Reciprocity",
                          "transTrip",
                          "transRecTrip",
                          "3-cycles",
                          "Indeg. Pop. sqrt.",
                          "Indeg. Act. sqrt.",
                          "Same sex",
                          "Alter drinking",
                          "Ego drinking",
                          "Ego x alt drinking",
                          "Alter smokimg",
                          "Ego smoking",
                          "Ego x alt smoking",
                          "Drinking rate 1",
                          "Drink. linear",
                          "Drink. quadratic",
                          "Avg. Alter Drink.",
                          "Drink: eff. from sex",
                          "Drink: eff. from smoke.",
                          "Smoking rate 1",
                          "Smoke. linear",
                          "Smoke. quadratic",
                          "Avg. Alter Smoke.",
                          "Smoke: eff. from sex",
                          "Smoke: eff. from drink.")[27:1])



g1 <- c("Friend rate 1",
                          "Density",
                          "Reciprocity",
                          "transTrip",
                          "transRecTrip",
                          "3-cycles",
                          "Indeg. Pop. sqrt.",
                          "Indeg. Act. sqrt.",
                          "Same sex")


g2 <- c("Alter drinking",
                          "Ego drinking",
                          "Ego x alt drinking",
                          "Alter smokimg",
                          "Ego smoking",
                          "Ego x alt smoking")
g3 <- c("Drinking rate 1",
                          "Drink. linear",
                          "Drink. quadratic",
                          "Avg. Alter Drink.",
                          "Drink: eff. from sex",
                          "Drink: eff. from smoke.",
                          "Smoking rate 1",
                          "Smoke. linear",
                          "Smoke. quadratic",
                          "Avg. Alter Smoke.",
                          "Smoke: eff. from sex",
                          "Smoke: eff. from drink.")

ld$g[ld$parameter %in% g1] <- "Structure"
ld$g[ld$parameter %in% g2] <- "Selection"
ld$g[ld$parameter %in% g3] <- "Behavior"

ld$g <- factor(ld$g, levels = c("Structure","Selection","Behavior"))

library(ggplot2)
p5 <- ggplot(ld, aes(x = parameter, y = Theta, group = model, colour = model,
                     shape = model)) +
         geom_point( size = 2.2 , position = position_dodge(width = .75)) +
         geom_hline(yintercept = 0, size = 1, color = "black") +
         geom_errorbar(size = .7, aes(ymin = Theta - SE, ymax = Theta + SE),
                        position = position_dodge(width = .75)) +
         theme_classic() +
         theme(axis.text.x = element_text(angle = 0, hjust = 1),
                axis.text.y = element_text(angle = 0, hjust = 1),
                legend.title = element_blank(),
                legend.position = "bottom",
                legend.direction = "vertical",
                legend.key.height = unit(1,"line"),
                legend.key.width = unit(1,"line"),
                text = element_text(size = 20),
                axis.title.x = element_blank()) +
         ylab("") + xlab("") + coord_flip()
p5


p6 <- ggplot(ld, aes(x = parameter, y = Theta, group = model, colour = model,
                     shape = model)) + scale_color_manual(values = c(
                    "grey75","grey60","grey45","grey30","grey15","grey0" )) +
  facet_wrap(facets = ~g, scales = "free", nrow = 3) +
  geom_point( size = 2.2 , position = position_dodge(width = .75)) +
  geom_hline(yintercept = 0, size = .25, color = "black") +
  geom_errorbar(size = .7, aes(ymin = Theta - SE, ymax = Theta + SE),
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


## 2 standard errors
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



a <- saom.results.20.b.t2[[2]]
same <- rep(NA, D)

for (i in 1:D) {
     if (i > 1) {
        if (sum(a[[1]][[i]] ==  a[[1]][[i-1]]) > 0) {
            same[i-1] = TRUE
            same[i] = TRUE
        } else {
            same[i] = FALSE
        }
     }
 }

print(which(same))
not_same = which(!same)

b <- list()
for (i in 1:length(not_same)) {
    b[[i]] = list(a[[1]][[not_same[i]]], a[[2]][[not_same[i]]])
}



npar <- sum(effects.imputed$include)

rowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

MIResults <- as.data.frame(matrix(NA,npar,(2 * D)))

for (i in 1:29) {
  names(MIResults)[i * 2 - 1] <- paste("imp" , "mean", sep = as.character(i))
  names(MIResults)[i * 2] <- paste("imp" , "se", sep = as.character(i))
  MIResults[,i * 2 - 1] <- b[[i]][[1]]  # estimates
  MIResults[,i * 2] <-  sqrt(diag(b[[i]][[2]]))  # standard errors
}

# Now we get the average covariance structure between the parameters
WDMIs <- matrix(0,npar,npar)

for (i in 1:29) {
  WDMIs <- WDMIs + b[[i]][[2]]
}

WDMIs <- (1/29) * WDMIs

# Using Rubin's Rules we combine the parameters and standard errors and
# complete the procedure
finalResults <- as.data.frame(matrix(NA,npar,2))
names(finalResults) <- c("combinedEstimate", "combinedSE")
rownames(finalResults) <- effects.imputed$effectName[effects.imputed$include]

finalResults$combinedEstimate <- rowMeans(MIResults[,seq(1,2*29,2)])
finalResults$combinedSE <- sqrt(diag(WDMIs) + ((29 + 1)/29) *
                                  rowVar(MIResults[,seq(1,2*29,2)]))



clean_results <- function(soam.result.list) {

}

library(sm)


diag(b[[i]][[2]])