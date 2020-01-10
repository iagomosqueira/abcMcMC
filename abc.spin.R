# test.R - DESC
# /test.R

# Copyright Iago MOSQUEIRA, 2019
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the GPL 3.0

#' # Running abcMcMC on test dataset

library(FLasher)
library(ggplotFL)
library(mvtnorm)
library(data.table)

source("functions.R")

load("test.RData")

# EXTRACT catch series
catch <- catch(om)

# SURVEY selectivity
survsel <- FLQuant(c(0.2, 0.80, 1, 1, 0.98, 0.90, 0.80, 0.76, 0.72),
  dimnames=dimnames(m(om)))

# GENERATE survey
survey <- FLIndex(index=stock.n(om) * survsel * exp(-z(om) * 0.5),
  sel.pattern=survsel)
range(survey, c("startf", "endf")) <- c(0.4,0.6)

# SURVEY variance at age
surveysigma <- diag(c(0.2, 0.2, 0.2, 0.2, 0.25, 0.3, 0.35, 0.40, 45), ncol=9)

# COMPARE survey and OM
plot(FLQuants(OM=stock.n(om), SURVEY=index(survey)))

# TEST abcMcMC
system.time(
  t1 <- abcMcMC(biol, fisheries, catch[,-1], survey, surveysigma, iter=1000,
    vars=c(100, 2), verbose=TRUE)
)

# CATCH is matched perfectly, if possible
plot(FLQuants(C=catch(t1[[2]][[1]]), NC=catch)) + ylim(c(0,NA))

# COMPUTE acceptance rate
sum(t1$chain$accept, na.rm=TRUE) / dim(t1$chain)[1]

# PLOT trajectories SSB
ggplot(ssb(t1$biol), aes(x=year, y=data, group=iter)) + geom_line()

# PLOT corr(v, ssb[end]) & corr(d,ssb[end])
chain <- t1$chain[accept==TRUE,]
chain[, ssb:=c(ssb(t1$biol)[,'30'])]
ggplot(chain, aes(x=v, y=ssb)) + geom_point()
ggplot(chain, aes(x=d, y=ssb)) + geom_point()
ggplot(chain, aes(x=v*d, y=ssb)) + geom_point()

# PLOT v prior vs. posterior
chain[, ini:=d*v]
pps <- melt(chain[, .(ini, ssb)])
ggplot(pps, aes(value, group=variable, fill=variable)) + geom_histogram() + facet_wrap(~variable)

# COMPARE SSB OM ~ T1
plot(FLQuants(OM=ssb(om), T1=ssb(t1$biol)))

# PLOT parameters
ggplot(melt(t1$chain[,-5]), aes(x=value)) +
  geom_histogram() + xlab("") + ylab("") +
  facet_wrap(~variable, scales='free')

# SAVE
save(t1, om, file="t1.RData", compress="xz")
