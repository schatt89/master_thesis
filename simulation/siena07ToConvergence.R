siena07ToConvergence <- function(alg, dat, eff, ans0=NULL, threshold, nodes = 3,
                                 cluster = TRUE, n2startPrev = 1000, ...) {
  # parameters are:
  # alg, dat, eff: Arguments for siena07: algorithm, data, effects object.
  # ans0: previous answer, if available; used as prevAns in siena07.
  # threshold: largest satisfactory value
  #            for overall maximum convergence ratio (indicating convergence).
  # nodes: number of processes for parallel processing.
  numr <- 0
  if (is.null(ans0)) {
    ans <- siena07(alg, data = dat, effects = eff, prevAns = ans0,nbrNodes = nodes,
                   returnDeps = TRUE, useCluster = cluster, ...) # the first run
  } else {
    alg$nsub <- 1
    alg$n2start <- n2startPrev
    ans <- siena07(alg, data = dat, effects = eff, prevAns = ans0,nbrNodes = nodes,
                   returnDeps = TRUE, useCluster = cluster, ...)
  }
  repeat {
    #save(ans, file = paste("ans",numr,".RData",sep = "")) # to be safe
    numr <- numr + 1           # count number of repeated runs
    tm <- ans$tconv.max      # convergence indicator
    cat(numr,"tconv  max:", round(tm,3),"\n")       # report how far we are
    if (tm < threshold) {break}   # success
    if (tm > 3) {stop()}     # divergence without much hope
    # of returning to good parameter values
    if (numr > 100) {stop()}  # now it has lasted too long
    alg$nsub <- 1
    alg$n2start <- 1000 + numr * 1000
    alg$n3 <- 2000 + numr * 1000
    ans <- siena07(alg, data = dat,effects = eff,prevAns = ans,nbrNodes = nodes,
                   returnDeps = TRUE, useCluster = cluster, ...)
  }
  if (tm > threshold) {
    stop("Warning: convergence inadequate.\n")
  }
  ans
}