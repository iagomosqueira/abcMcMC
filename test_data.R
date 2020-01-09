# test.R - CREATE test data set for abcMcMC

# Copyright Iago MOSQUEIRA (WMR), 2019
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(FLife)
library(FLasher)
library(ggplotFL)
library(mvtnorm)
library(data.table)

source("functions.R")

# --- CREATE test OM 

nyears <- 30

# GENERATE linked LH params from Linf, v=400
lhpar <- lhPar(FLPar(linf=100, v=400, a50=1))

# DERIVE equilibrium population
equil <- lhEql(lhpar, spwn=0.5, 
  range = c(min = 0, max = 8, minfbar = 0, maxfbar = 1, plusgroup = 8))

# COERCE into stock object
fpp <- as(equil, "FLStock")

# SUBSET for single initial F level F=0.0077
stk <- fpp[, 2]
dimnames(stk) <- list(year=1, age=1:9)
range(stk, c('minfbar', 'maxfbar')) <- c(2, 8)
om <- fwdWindow(stk, equil, end=nyears)

# DEBUG OM values: f0 = 0.00766, K = 400
ompar <- FLPar(f0=0.00766, K=400)

# SET SRR
bhm <- as(equil, "predictModel")

# R1: PROJECT to 3*FMSY, then 1.3*FMSY
trajectory <- c(
  # FROM low to FMSY * 1.4
  seq(c(fbar(stk)), c(fmsy(equil)) * 3, length=9),
  # FROM FMSY * 1.3 to FMSY
  seq(c(fmsy(equil)) * 3, c(fmsy(equil)), length=20))
control <- fwdControl(year=seq(2, nyears), quant="f", value=trajectory)
om <- fwd(om, sr=bhm, control=control)

plot(om)

# CONVERT om to FLB + FLFs
biol <- as(om, "FLBiol")
n(biol)[,-1] <- NA
rec(biol) <- bhm
fisheries <- FLFisheries(F=as(om, "FLFishery"))

# DEBUG ADD to coertion
names(fisheries[[1]]) <- "TES"
name(biol) <- "TES"

# SAVE test objects
save(biol, fisheries, om, file="test.RData", compress="xz")
