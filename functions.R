# om.R - DESC
# /om.R

# Copyright European Union, 2019
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the GPL 3.0


# lhPar {{{

lhPar <- function(...,
    m=list(model="gislason", params=c(m1=0.55, m2=-1.61, m3=1.44))) {

  # DEFAULT defined parameter values + m
  args <- list(a=0.0003, b=3, ato95=1, sel2=1, sel3=5000,
    s=0.9, v=1000, asym=1, m=m)
  
  # PARSE ...
  input <- list(...)
  
  # FIND FLPar(s) and data.frame(s) in input
  flp <- unlist(lapply(input, function(x) is(x, "FLPar")))
  dtf <- unlist(lapply(input, function(x) is(x, "data.frame")))

  # CONVERT all to single list of vectors
  input <- c(input[!flp & !dtf],
      Reduce("c", lapply(input[flp], as, "list")),
      unlist(lapply(input[dtf], unlist)))

  # DROP any NAs
  input <- input[!unlist(lapply(input, function(x) any(is.na(x))))]

  # MERGE input and default args
  args[names(input)] <- input

  # PARSE m params
  margs <- args$m[unlist(lapply(args$m, is, "numeric"))][[1]]

  # ENSURE m params are named
  if(is.null(names(margs)))
    names(margs) <- paste0("m", seq(length(margs)))

  # STORE m model name
  mmodel <- args$m[unlist(lapply(args$m, is, "character"))][[1]]

  # ADD args and m params to form params
  params <- c(args[names(args) != "m"], margs)

  # EXPAND to max iters
  its <- max(unlist(lapply(params, length)))

  # CREATE output FLPar
  params <- do.call("FLPar", lapply(params, rep, length.out=its))
  
  # DERIVED parameters
  #
  # Gislason, H., J.G. Pope, J.C. Rice, and N. Daan. 2008. Coexistence in
  # North Sea fish communities: Implications for growth and natural mortality.
  # ICES J. Mar. Sci. 65 (4): 514–30.

  # k
  if(!"k" %in% dimnames(params)$params)
    params <- rbind(params, FLPar(k=3.15 * params$linf ^ (-0.64)))

  # t0
  if(!"t0" %in% dimnames(params)$params)
    params <- rbind(params, FLPar(t0=-exp(-0.3922 - 0.2752 *
      log(params$linf) %-% (1.038 * log(params$k)))))

  # l50 - a50
  if(!"l50" %in% dimnames(params)$params) {
    if("a50" %in% dimnames(params)$params) {
      params <- rbind(params,
        FLPar(l50=vonB(age=c(params$a50), params[c("k", "t0", "linf"),])))
    } else {
      params <- rbind(params, FLPar(l50=0.72 * params$linf ^ 0.93))
    }
  }
  if(!"a50" %in% dimnames(params)$params) {
    params <- rbind(params, FLPar(a50=log(1-(params$l50 %/%
      params$linf)) %/% (-params$k) %+% params$t0))
  }

  # sel1
  if(!"sel1" %in% dimnames(params)$params)
    params <- rbind(params, FLPar(sel1=params$a50 + params$ato95))
  
  # bg
  if(!"bg" %in% dimnames(params)$params)
    params <- rbind(params, FLPar(bg=params$b))

  # SORT params
  order <- c("linf", "l50", "a50", "ato95", "k", "t0", "a", "b",
    "m1", "m2", "m3", "sel1", "sel2", "sel3", "asym", "s", "v", "bg")

  params <- params[order,]

  # KEEP mmodel as attribute
  attr(params, "mmodel") <- mmodel

  return(params)
} # }}}

# initiate {{{

initiate <- function(biom, f0=0, biol, catch.sel,
  minfbar=dims(biol)$min, maxfbar=dims(biol)$max) {

  # EXTRACT to vectors
  wt <- c(wt(biol)[, 1])
  m <- c(m(biol)[, 1])
  n <- c(n(biol)[, 1])
  sel <- c(catch.sel[, 1])

  # DEBUG GET fbar selectivity
  idx <- do.call(seq,
    as.list(match(ac(c(minfbar, maxfbar)), dimnames(n(biol))$age)))
  psel <- c(catch.sel[idx, 1])

  # DEBUG
  if(length(sel) != length(psel)) {
  
    # SOLVE for difference between F0 and Fbar
    foo <- function(f) {
      abs(mean(sel * f0) - mean(psel * f))
    }

    ff0 <- optim(f0, foo, method="Brent", lower=0, upper=f0 * 2)$par
  }

  # SOLVE R0 for biom, WT, M + F
  foo <- function(R0) {

    n[1] <- R0

    for(a in seq(2, length(n)))
      n[a] <- n[a-1] * exp(-(m[a-1] + f0 * sel[a-1]))

    return(abs(biom - sum(wt * n)))
  }

  res <- optim(2, foo, method="Brent", lower=1, upper=1e36)

  # INITIAL population
  n(biol)[1, 1] <- res$par
  for(a in seq(2, length(n)))
    n(biol)[a, 1] <- n(biol)[a-1, 1] * exp(-(m(biol)[a-1, 1] +
      f0 * catch.sel[a-1, 1]))

  # RETURN FLBiol

  return(biol)
}

# }}}

# simulator {{{

simulator <- function(biol, fisheries, v, d=0.5, control, ...) {

  # INITIATE N0
  biol <- initiate(biom=v, f0=0, biol, catch.sel(fisheries[[1]][[1]])[,1])

  # FWD to depletion
  eq <- fwd(biol, fisheries,
    control=fwdControl(year=2:30, quant="ssb_spawn",
    value=seq(c(ssb(biol)[,1]), c(ssb(biol)[,1]) * d, length=29), biol=1))

  n(biol)[,1] <- n(eq$biol)[,30]

  # FWD w/control
  res <- fwd(biol, fisheries, control=control, ...)

  return(res[c("biol", "fisheries")])
}

# }}}

# abcMcMC(biol, fisheries, catch, survey, priors, iters=nrow(priors)) {{{
# TODO ADD depletion, beta prior
abcMcMC <- function(biol, fisheries, catch, survey, surveySigma, iter=500,
  burnin=round(iter * 0.10), priors=list(v=c(log(400), 0.3), d=c(2, 2)),
  vars=c(v=100, d=1)) {

  # CONTROL from FLQuants
  control <- as(FLQuants(catch=catch), 'fwdControl')

  # DIMS
  nyears <- dims(biol)$year
  
  # PRIORS and PI(theta)
  v <- rlnorm(iter, priors$v[1], priors$v[2])
  pv <- dlnorm(v, priors$v[1], priors$v[2])

  d <- rbeta(iter, priors$d[1], priors$d[2])
  pd <- dbeta(d, priors$d[1], priors$d[2])
  
  ptheta <- pv * pd

  # CALCULATE mvtnorm as matrix normal
  nyx <- dim(index(survey))[2] - 1
  ky <- diag(1, nrow=nyx, ncol=nyx)
  ka <- surveySigma
  
  # Kronecker matrix product for matrix normal distribution
  sigmax <- ka %x% ky

  # CHAIN: theta (v, f0), ptheta, pdeviation, accept
  chain <- data.frame(v=v, d=d, ptheta=ptheta, pdeviation=as.numeric(NA),
    accept=NA)

  runs <- vector('list', length=iter)

  for(i in seq(iter)) {

    # RUN model
    run <- simulator(biol=biol, fisheries=fisheries, v=v[i], d=d[i],
      control=control)# , deviances=rlnorm(1, rec(biol)/rec(biol), 0.25))

    # EXTRACT objects
    flb <- run$biol
    flfs <- run$fisheries

    # KEEP all simulator$biols
    runs[[i]] <- flb

    # RESULT object
    if(i == 1) {
      omb <- propagate(flb, iter, fill.iter=FALSE)
      omfs <- lapply(flfs, propagate, iter, fill.iter=FALSE)
      # ACCEPT 1st run
      chain[1, "accept"] <- TRUE
    }

    # GET survey timing from spawning
    survt <- mean(range(survey, c("startf", "endf")))
    
    # GENERATE survey
    survn <- n(flb) * exp(-m(flb) * survt - harvest(flb, flfs) * survt) *
      sel.pattern(survey)

    # TEST
    # survn <- rlnorm(1, log(index(survey)), 0.02)

    # COMPUTE deviation log(D - X') (roc(survey) - roc(N*selex))
    deviation <- log(index(survey)[, -c(1)] / index(survey)[, -c(nyears)]) -
      log(survn[, -c(1)] / survn[, -c(nyears)])
    devx <- c(t(as.matrix(deviation@.Data[,, 1, 1, 1, 1])))

    # pi_E(D - X')
    chain[i, "pdeviation"] <- dmvnorm(x=devx, sigma=sigmax)

    # ACCEPTANCE probability
    if(i > 1) {
      
      # q(theta', theta_t)
      qthetanew <- mapply(function(x, y, z) dnorm(x, y, z),
        chain[i, c("v", "d")], chain[i - 1, c("v", "d")],
        vars)
    
      # q(theta_t, theta')
      qthetaold <- mapply(function(x, y, z) dnorm(x, y, z),
        chain[i - 1, c("v", "d")], chain[i, c("v", "d")],
        vars)

      # BUG: DROP runs with large deviations?
      if(chain[i, "pdeviation"] == 0)
        pprop <- 0
      else
        pprop <- min(0,
          sum(log(c(unlist(chain[i, c("pdeviation", "ptheta")]), qthetanew))) -
          sum(log(c(unlist(chain[i - 1, c("pdeviation", "ptheta")]),
            qthetaold))))

      # REJECT / ACCEPT
      if(pprop < log(runif(1, 0, 1))) {
        iter(omb, i) <- flb
        iter(omfs[[1]], i) <- flfs[[1]]

        chain[i, c("v", "d")] <- c(v[i], d[i])
        chain[i, "accept"] <- TRUE
        
      } else {
        iter(omb, i) <- iter(omb, i-1)
        iter(omfs[[1]], i) <- iter(omfs[[1]], i-1)

        chain[i, c("v", "d")] <- c(v[i-1], d[i-1])
        chain[i, "accept"] <- FALSE
      }
    }

  # PROGRESS
  cat("\r", i, "of", iter) 
  flush.console()
  }
  cat("\n")

  # TODO ADD burnin
  res <- list(biol=omb, fisheries=omfs, chain=data.table(chain),
    runs=runs)
  
  return(res)
}

# }}}
