# test.R - DESC
# test.R

# Copyright Iago MOSQUEIRA (WMR), 2019
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(FLasher)
library(ggplotFL)
library(mvtnorm)
library(data.table)

source('functions.R')

load("test.RData")


# --- TEST initiate

# CALL initiate w/test biomass and F

tei <- initiate(c(quantSums(n(biol)[,1] * wt(biol)[,1])),
  f0=0, biol, catch.sel=catch.sel(fisheries[[1]][[1]])[,1])

# COMPARE B0 at F0=0
quantSums(n(tei)[,1] * wt(tei)[,1])
quantSums(n(biol)[,1] * wt(biol)[,1])
quantSums(stock.n(om)[,1] * stock.wt(om)[,1])

# COMPARE stock.n

n(tei)[,1] / stock.n(om)[,1]
wt(tei)[,1] / stock.wt(om)[,1]


# --- TEST simulator

# FOLLOW catch series

control <- as(FLQuants(catch=catch(om)[,-1]), "fwdControl")

sim <- simulator(biol, fisheries,
  v=c(quantSums(n(biol)[,1] * wt(biol)[,1])), d=0.5, control)


sims <- lapply(seq(0.1, 0.9, by=0.1), function(x)
  simulator(biol, fisheries, v=220, d=x, control)$biol)

lapply(sims, function(x) 220/ssb(x)[,1])

plot(FLQuants(lapply(sims, tsb))) + facet_grid(qname~.)


tes <- fwd(om, sr=rec(biol, FALSE), control=control)

stock.n(tes)[,1]
n(sim$biol)[,1]

stock.n(tes)[,2]
n(sim$biol)[,2]

plot(FLBiols(OM=biol, SIM=sim$biol))

# FOLLOW SSB (~I) series
control <- as(FLQuants(ssb_spawn=ssb(biol)[,-1]), "fwdControl")

sim <- simulator(biol, fisheries, v=ompar$K, d=0.5, control)

plot(FLBiols(OM=biol, SIM=sim$biol))

# FOLLOW F
control <- as(FLQuants(f=fbar(om)[,-1]), "fwdControl")
  control@target[,c('minAge')] <- 2
  control@target[,c('maxAge')] <- 8
sim <- simulator(biol, fisheries, v=ompar$K, f0=ompar$f0, control)

plot(FLBiols(OM=biol, SIM=sim$biol))


# --- TEST FLBiol

bi <- propagate(biol, 2, fill.iter=FALSE)

all.equal(iter(bi,1), iter(bi,2))

fec(iter(bi,1))
fec(iter(bi,2))

ssb(iter(bi,1))
ssb(iter(bi,2))

iter(bi,2) <- iter(bi,1)

iter(mat(bi, FALSE), 2) <- iter(mat(bi, FALSE), 1)



# References
