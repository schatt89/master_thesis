SimulateNetworksBehavior <- function(net.w1, b1.w1, n, c1,
                                     rate, dens, rec, tt, tRt, c3, 
                                     outActSq, inPopSq,
                                     altX.b1, egoX.b1, egoXaltX.b1,
                                     rate.b1, lin.b1, qu.b1, avalt.b1){
  b1.w1[is.na(b1.w1)] <- 1
  
  # Create initial 2-wave data to get a suitable data structure.
  # arbitrarily, this initial network has an expected average degree of 3
  X0 <- matrix(rbinom(n * n, 1, 3 / (n - 1)), n, n)
  diag(X0) <- 0
  X1 <- X0
  # but X0 and X1 should not be identical for use in sienaDependent
  X0[1, 2] <- 0
  X0[2, 1] <- 1
  X1[1, 2] <- 1
  X1[2, 1] <- 0
  XX <- array(NA, c(n, n, 2))
  XX[, , 1] <- X0
  XX[, , 2] <- X1
  # Create behavior variable; initial distribution uniform on {1, ..., c}.
  ZZ1 <- pmin(matrix(trunc(c1 * runif(n * 2)) + 1, n, 2), c1)
  # With this data structure, we now can create the data.
  X  <- sienaDependent(XX, allowOnly = FALSE)
  Z1  <- sienaDependent(ZZ1, type = "behavior", allowOnly = FALSE)
  InitData <- sienaDataCreate(X, Z1)
  InitEff0 <- getEffects(InitData)
  # sink to avoid printing to the screen
  sink("eff.txt")
  # Specify the parameters.
  # The rate parameters are first multiplied by 10,
  # which will be used only to get from the totally random network XX[,,1] = X0
  # to the network that will be the simulated first wave.
  InitEff0 <- setEffect(InitEff0, Rate, type = "rate", initialValue = 10 * rate)
  InitEff0 <- setEffect(InitEff0, density, initialValue = dens)
  InitEff0 <- setEffect(InitEff0, recip, initialValue = rec)
  InitEff0 <- setEffect(InitEff0, transTrip, initialValue = tt)
  InitEff0 <- setEffect(InitEff0, cycle3, initialValue = c3)
  InitEff0 <- setEffect(InitEff0, outActSqrt, initialValue = outActSq)  
  InitEff0 <- setEffect(InitEff0, inPopSqrt, initialValue = inPopSq)
  
  # altX, egoX, egoXaltX
  
  InitEff0 <- setEffect(InitEff0, altX, initialValue = altX.b1,
                        interaction1 = "Z1")
  InitEff0 <- setEffect(InitEff0, egoX, initialValue = egoX.b1,
                        interaction1 = "Z1")
  InitEff0 <- setEffect(InitEff0, egoXaltX, initialValue = egoXaltX.b1,
                        interaction1 = "Z1")
  
  InitEff0 <- setEffect(InitEff0, name = "Z1", Rate, type = "rate",
                        initialValue = 10 * rate.b1)
  InitEff0 <- setEffect(InitEff0, name = "Z1", linear, initialValue = lin.b1)
  InitEff0 <- setEffect(InitEff0, name = "Z1", quad, initialValue = qu.b1)
  InitEff0 <- setEffect(InitEff0, name = "Z1", avAlt, interaction1 = "X",
                        initialValue = avalt.b1)
  
  
  ## The parameter given for n3 should be larger than sum(InitEff0$include)
  nthree <- sum(InitEff0$include)	+ 5
  InitAlg <- sienaAlgorithmCreate(projname = "Init", useStdInits = FALSE,
                                  cond = FALSE, nsub = 0, n3 = nthree,
                                  simOnly = TRUE)
  # Simulate the first wave.
  # InitSim   <- siena07(InitAlg, data = InitData, eff = InitEff0,
  #                      returnDeps = TRUE, batch = TRUE, silent = TRUE)
  # Now prepare for simulating waves 2 to M.
  # Create empty result network and behavior matrices
  Xs <- array(NA, dim = c(n, n, M))
  Z1s <- array(NA, dim = c(n, M))
  # The rate parameter values from the function call are reinstated in InitEff.
  InitEff <- InitEff0
  InitEff <- setEffect(InitEff, Rate, type = "rate", initialValue = rate)
  InitEff <- setEffect(InitEff, name = "Z1", Rate, type = "rate",
                       initialValue = rate.b1)
  
  sink()
  for (m in 1:M) {
    # Note that we start this loop with a previously simulated network.
    # Transform the previously simulated network
    # from edge list into adjacency matrix
    if (m == 1) {
      XXsim <- net.w1
      Z1sim <- b1.w1
      # Put simulated network and behavior into the result matrix.
      Xs[, , 1] <- XXsim
      Z1s[, 1] <- Z1sim
      
      # Put simulated network in desired places for the next simulation
      XX[, , 2] <- XX[, , 1] # used only to get the data structure
      XX[, , 1] <- XXsim
      ZZ1[, 2] <- ZZ1[, 1]
      ZZ1[, 1] <- Z1sim
    } else {
      XXsim <- matrix(0, n, n)
      nsim  <- InitAlg$n3
      XXsim[InitSim$sims[[nsim]][[1]]$X[[1]][, 1:2]]  <-
        InitSim$sims[[nsim]][[1]]$X[[1]][, 3]
      Z1sim <- InitSim$sims[[nsim]][[1]]$Z1[[1]]
      # Put simulated network and behavior into the result matrix.
      Xs[, , m] <- XXsim
      Z1s[, m] <- Z1sim
      # Put simulated network in desired places for the next simulation
      XX[, , 2] <- XX[, , 1] # used only to get the data structure
      XX[, , 1] <- XXsim
      ZZ1[, 2] <- ZZ1[, 1]
      ZZ1[, 1] <- Z1sim
    }
    
    if (m < M) {
      # The following is only to prevent the error that would occur
      # in the very unlikely event XX[,,1] == XX[,,2] or ZZ[,1] == ZZ[,2].
      if (identical(XX[, , 1], XX[, , 2])) {
        XX[1, 2, 1] <- 1 - XX[1, 2, 2]
      }
      if (identical(ZZ1[, 1], ZZ1[, 2])) {
        ZZ1[1, 1] <- ifelse((ZZ1[1, 1] == 1), 2, 1)
      }
      # Specify the two-wave network data set starting with XX[,,1].
      X <- sienaDependent(XX, allowOnly = FALSE)
      Z1 <- sienaDependent(ZZ1, type = "behavior", allowOnly = FALSE)
      # Simulate wave m+1 starting at XX[,,1] which is the previous XXsim
      InitData  <- sienaDataCreate(X, Z1)
      InitSim <- siena07(InitAlg, data = InitData, eff = InitEff,
                         returnDeps = TRUE, batch = TRUE, silent = TRUE)
    }
  }
  # Present the average degrees to facilitate tuning the outdegree parameter
  # to achieve a desired average value for the average degrees.
  closeAllConnections()
  cat("Average degrees ", round(colSums(Xs, dims = 2) / n, digits = 2), "\n")
  cat("Average behavior ", round(colSums(Z1s) / n, digits = 2), "\n")
  # Result: simulated data set; covara and covarb are vectors of length n;
  # networks is an array of dimension nxnxM;
  # behaviors is a matrix of dimension nxM
  list(networks = Xs[,,2:M], behavior = Z1s[,2])
}