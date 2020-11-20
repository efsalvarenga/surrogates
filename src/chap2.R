# Four motivating datasets

# setup
library(dplyr)
library(akima)
library(ggplot2)

# Rocket booster dynamics
# The first, oldest version of the data was derived from a less reliable code implementing the solver.
lgbb1 <- read.table("./data/lgbb/lgbb_original.txt", header=TRUE)
names(lgbb1)

# interpolates for plotting
g <- interp(x = lgbb1$mach,
            y = lgbb1$alpha,
            z = lgbb1$lift,
            duplicate = "mean")

# visual of interpolated data
image(g, col=heat.colors(128), xlab="mach", ylab="alpha")
points(lgbb1$mach, lgbb1$alpha, cex=0.25, pch=18)

# very few unique 'beta' variables
apply(lgbb1, 2, function(x) { length(unique(x)) })

# subset where 'alpha' == 1, and beta is either equal or non-equal to zero
a1b0 <- which(lgbb1$alpha == 1 & lgbb1$beta == 0)
a1bn0 <- which(lgbb1$alpha == 1 & lgbb1$beta != 0)
a1b0 <- a1b0[order(lgbb1$mach[a1b0])]
a1bn0 <- a1bn0[order(lgbb1$mach[a1bn0])]

# clearly there are some issues with lift outputs when beta != 0
# Also note the relatively low resolution, with each “curve” being traced
# out by just a handful of values – fewer than fifty in both cases. 
plot(lgbb1$mach[a1b0], lgbb1$lift[a1b0], type="l", xlab="mach", 
     ylab="lift", ylim=range(lgbb1$lift[c(a1b0, a1bn0)]), lwd=2)
lines(lgbb1$mach[a1bn0], lgbb1$lift[a1bn0], col=2, lty=2, lwd=2)
text(4, 0.3, paste("length(a1b0) =", length(a1b0)))
text(4, 0.25, paste("length(a1bn0) =", length(a1bn0)))
legend("topright", c("beta = 0", "beta != 0"), col=1:2, lty=1:2, lwd=2)

# Improved simulations were paired with model-based sequential design (§6.2)
# under a treed Gaussian process (TGP, §9.2.2) in order to obtain a more
# adaptive, automatic input “grid”
lgbb2 <- read.table("./data/lgbb/lgbb_as.txt", header=TRUE)

# compare frist and second study sampling
par(mfrow=c(1,2))
plot(lgbb1$mach, lgbb1$alpha, xlab="mach", ylab="alpha", pch=18, cex=0.5)
plot(lgbb2$mach, lgbb2$alpha, xlab="mach", ylab="alpha", pch=18, cex=0.5)

# more resolution on 'match'
apply(lgbb2[, -1], 2, function(x) { length(unique(x)) })

# (...) slice over mach and alpha when beta=1. The data behind this visual
# comes from an .RData file containing evaluations of the TGP surrogate,
# trained to 780 lift evaluations, on a dense predictive grid in the input space.
load("./data/lgbb/lgbb_fill.RData")

# much more resolution
apply(lgbb.fill, 2, function(x) { length(unique(x)) })
par(mfrow=c(1,3))
plot(lgbb1$mach, lgbb1$alpha, xlab="mach", ylab="alpha", pch=18, cex=0.5)
plot(lgbb2$mach, lgbb2$alpha, xlab="mach", ylab="alpha", pch=18, cex=0.5)
plot(lgbb.fill$mach, lgbb.fill$alpha, xlab="mach", ylab="alpha", pch=18, cex=0.5)

# smoother x = 'mach' y = 'alpha' z = 'lift'
# Notice how this view reveals a nice smooth surface with simple dynamics in
# high-speed regions, and a more complex relationship near the sound barrier – in
# particular for low speeds (mach) and high angles of attack (alpha).
par(mfrow=c(1,1))
lgbb.b1 <- lgbb.fill[lgbb.fill$beta == 1, ]
g <- interp(lgbb.b1$mach, lgbb.b1$alpha, lgbb.b1$lift)
image(g, col=heat.colors(128), xlab="mach [beta=1]", ylab="alpha [beta=1]")

# This view conveys nuance along the sound barrier more clearly than the previous
# image plot did
plot(lgbb.b1$mach, lgbb.b1$lift, type="n", xlab="mach", 
     ylab="lift [beta=1]")
for(ub in unique(lgbb.b1$alpha)) {
  colour <- which(ub == unique(lgbb.b1$alpha))
  a <- which(lgbb.b1$alpha == ub)
  a <- a[order(lgbb.b1$mach[a])]
  lines(lgbb.b1$mach[a], lgbb.b1$lift[a], type="l", lwd=2, col = colour) 
}

# University of Michigan’s Center for Radiative Shock Hydrodynamics (CRASH)

# load field experiment data
crash <- read.csv("./data/crash/CRASHExpt_clean.csv")
crash$BeThickness <- 21 #adding missing variable
names(crash)
str(crash)

# Load computer experiments data
ce1 <- read.csv("./data/crash/RS12_SLwithUnnormalizedInputs.csv")
ce2 <- read.csv("./data/crash/RS13Minor_SLwithUnnormalizedInputs.csv")
ce2$ElectronFluxLimiter <- 0.06 # correcting wrongly recorded variable
str(ce1)
str(ce2)
summary(ce1)
summary(ce2)

# expanding ce2 to match ce1 variables
sfmin <- ce2$EffectiveLaserEnergy/5000
sflen <- 10
ce2.sf <- matrix(NA, nrow = sflen * nrow(ce2), ncol = ncol(ce2) + 2)
for(i in 1:sflen) {
  sfi <- sfmin + (1 - sfmin) * (i / sflen)
  ce2.sf[(i - 1) * nrow(ce2) + (1 : nrow(ce2)),] <- 
    cbind(as.matrix(ce2), sfi, ce2$EffectiveLaserEnergy/sfi)
}
ce2.sf <- as.data.frame(ce2.sf)
names(ce2.sf) <- c(names(ce2), "EnergyScaleFactor", "LaserEnergy")

# difference in laser energy and energy sf between ce's
plot(ce2.sf$LaserEnergy, ce2.sf$EnergyScaleFactor, ylim=c(0.4, 1.1), 
     xlab="Laser Energy", ylab="Energy Scale Factor")
points(ce1$LaserEnergy, ce1$EnergyScaleFactor, col=2, pch=19)
legend("bottomleft", c("CE2", "CE1"), col=1:2, pch=c(21,19))

# lm on field data
u <- apply(crash, 2, function(x) { length(unique(x)) })
u
fit <- lm(ShockLocation ~ ., data=crash[, u > 1])
summary(fit)

# Time mops up nearly all of the variability in these data 
fit.time <- lm(ShockLocation ~ Time, data=crash)
summary(fit.time)
plot(crash$Time, crash$ShockLocation, xlab="time", ylab="location")
abline(fit.time)

# lm on ce1 data
ce1 <- ce1[,-1] ## first col is FileNumber
u.ce1 <- apply(ce1, 2, function(x) { length(unique(x)) })
u.ce1
fit.ce1 <- lm(ShockLocation ~ ., data=ce1[,u.ce1 > 1])
summary(fit.ce1)