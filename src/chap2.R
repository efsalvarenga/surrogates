# Four motivating datasets

# setup
library(dplyr)
library(akima)
library(ggplot2)
library(laGP)
library(olsrr)
library(lhs)

# -------------------------------------------------------------------------------------------------
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

# Homework exercises

# Our rocket booster example in §2.1 emphasized the lift output. Repeat similar slice visuals,
# for example like Figure 2.6, for the other five outputs. In each case you’ll need to choose
# a value of the third input, beta, to hold fixed for the visualization.

# a. Begin with the choice of beta=1 following the lift example. Comment on any trends or
# variations across the five (or six including lift) outputs.

# b. Experiment with other beta choices. In particular what happens when beta=0 for the latter
# three outputs: side, yaw and roll? How about with larger beta settings for those outputs?
# Explain what you think might be going on.
par(mfrow=c(2,3))

varNameVec <- names(lgbb.b1)[4:9]

for (chosenBeta in unique(lgbb.fill$beta)){
  lgbb.bX <- lgbb.fill[lgbb.fill$beta == chosenBeta, ]
  for (varName in varNameVec) {
    plot(lgbb.bX$mach, lgbb.bX[[varName]], type="n", xlab="mach", 
         ylab=paste0(varName, " [beta=", chosenBeta, ']'))
    for(ub in unique(lgbb.bX$alpha)) {
      colour <- which(ub == unique(lgbb.bX$alpha))
      a <- which(lgbb.bX$alpha == ub)
      a <- a[order(lgbb.bX$mach[a])]
      lines(lgbb.bX$mach[a], lgbb.bX[[varName]][a], type="l", lwd=2, col = colour) 
    }
  }
}
par(mfrow=c(1, 1))

# beta=0 means there is no relative wind to the longitudinal axis of the aircraft
# for this reason side force, yawl and roll are decreased to close to 0


# -------------------------------------------------------------------------------------------------
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
fit.time <- lm(ShockLocation ~ Time, data=crash[, u > 1])
summary(fit.time)
plot(crash$Time, crash$ShockLocation, xlab="time", ylab="location")
abline(fit.time)

# lm on ce1 data
ce1 <- ce1[,-1] ## first col is FileNumber
u.ce1 <- apply(ce1, 2, function(x) { length(unique(x)) })
u.ce1
fit.ce1 <- lm(ShockLocation ~ ., data=ce1[,u.ce1 > 1])
summary(fit.ce1)

# Figure  offers a view into how shock location varies with time
# and (scaled) laser energy. The heat plot in the figure is
# examining a linear interpolation of raw CE1 data
x <- ce1$Time
y <- ce1$LaserEnergy * ce1$EnergyScaleFactor
g <- interp(x/max(x), y/max(y), ce1$ShockLocation, dupl="mean")
image(g, col=heat.colors(128), xlab="scaled time", ylab="scaled energy")

# Time and energy are useful for estimating location
fit.te <- lm(ShockLocation ~ Time + LaserEnergy + Time * LaserEnergy, data=crash[, u > 1])
summary(fit.te)

# Homework exercises

# Revisit the CRASH simulation linear model/curve fitting analysis nearby Figure 2.10 by
# expanding the data and the linear basis.

# a. Form a data.frame combing CE1 data (ce1) and scale factor expanded CE2 data (ce2.sf),
# and don’t forget to drop the FileNumber column.
ceAll <- bind_rows(ce1, ce2.sf) %>%
  select(-FileNumber, -EffectiveLaserEnergy)
summary(ceAll)

# b. Fit a linear model with ShockLocation as the response and the other columns as predictors.
# Which predictors load (i.e., have estimated slope coefficients which are statistically
# different from zero)?
fit_CeAll01 <- lm(ShockLocation ~ ., data = ceAll)
summary(fit_CeAll01)

# c. Consider interactions among the predictors. Which load? Are there any collinearity
# concerns? Fix those if so. You might try stepwise regression.
ols_vif_tol(fit_CeAll01)

fit_CeAll02 <- lm(ShockLocation ~ (.)^2, data = ceAll)
summary(fit_CeAll02)
ols_vif_tol(fit_CeAll02)


fit_CeAll02_step1 <- stats::step(fit_CeAll02, scope = formula(fit_CeAll02), direction = 'backward',
                         k = log(nrow(ceAll)), trace = 1)
summary(fit_CeAll02_step1)
ols_vif_tol(fit_CeAll02_step1)

fit_CeAll02_step2 <- stats::step(fit_CeAll02, scope = formula(fit_CeAll02), direction = 'backward')
summary(fit_CeAll02_step2)
ols_vif_tol(fit_CeAll02_step2)

fit_CeAll02_step3 <- stats::step(fit_CeAll02, scope = formula(fit_CeAll02))
summary(fit_CeAll02_step3)
ols_vif_tol(fit_CeAll02_step3)

# c. Consider quadratic feature expansion (i.e., augment columns with squared terms) with
# and without interactions. Again watch out for collinearity; try step and comment on what loads.

# d. Contemplate higher-order polynomial terms as features. Does this represent a sensible,
# parsimonious approach to nonlinear surrogate modeling?


# -------------------------------------------------------------------------------------------------
# Predicting satellite drag
# tpm repo available at https://bitbucket.org/gramacylab/tpm/src/master/

# generating duplicate sample dataset
n <- 8
X <- data.frame(Umag=runif(n, 5500, 9500), Ts=runif(n, 100, 500), 
                Ta=runif(n, 200, 2000), theta=runif(n, -pi, pi), 
                phi=runif(n, -pi/2, pi/2), alphan=runif(n), sigmat=runif(n))
X <- rbind(X,X)

# mesh location for GRACE
mesh <- "../gramacylab-tpm-5f42493808f6/tpm/Mesh_Files/GRACE_A0_B0_ascii_redone.stl"

# sets up the atmospheric chemical composition as a unit-vector isolating pure helium (He)
moles <- c(0,0,0,0,1,0) 

# source the tpm function (I manually changed line 126 of the script)
# compiling it is required before continuing
source("../gramacylab-tpm-5f42493808f6/tpm/R/tpm.R")

system.time(y <- tpm(X, moles=moles, stl=mesh, verb=0))

mean((y[1:n] - y[(n + 1):length(y)])^2)

# LANL’s GRACE runs for pure He
train <- read.csv("../gramacylab-tpm-5f42493808f6/data/GRACE/CD_GRACE_1000_He.csv")
test <- read.csv("../gramacylab-tpm-5f42493808f6/data/GRACE/CD_GRACE_100_He.csv")
r <- apply(rbind(train, test)[,1:7], 2, range)
r

# Before fitting models, it helps to first convert to coded inputs.
X <- train[,1:7]
XX <- test[,1:7]
for(j in 1:ncol(X)) {
  X[,j] <- X[,j] - r[1,j]
  XX[,j] <- XX[,j] - r[1,j]
  X[,j] <- X[,j]/(r[2,j] - r[1,j])
  XX[,j] <- XX[,j]/(r[2,j] - r[1,j])
}

fit.gp <- newGPsep(X, train[,8], 2, 1e-6, dK=TRUE)
mle <- mleGPsep(fit.gp)
p <- predGPsep(fit.gp, XX, lite=TRUE)
rmspe <- sqrt(mean((100*(p$mean - test[,8])/test[,8])^2))
rmspe

# Homework exercises

# a. Get familiar with using tpm to calculate satellite drag coefficients. Begin by compiling
# the C code to produce a shared object for use with the R wrapper in tpm-git/tpm/R/tpm.R.
# Run R CMD SHLIB -o tpm.so *.c in the tpm-git/tpm/src directory. Then double-check you have a
# working tpm library by trying some of the examples in §2.3.2.

train[1,]
exa01 <- tpm(train[1, 1:7], moles = moles, stl = mesh, verb = 1)
exa01

train[12,]
exa02 <- tpm(train[12, 1:7], moles = moles, stl = mesh, verb = 1)
exa02

# b. Generate a random one-hundred run, 7d testing design uniformly in the ranges provided by
# Table 2.5 with restricted angle inputs as described in Table 2.6. Evaluate these inputs for
# GRACE under the mole fractions provided in §2.3.2 for a particular LEO position.

# create unit cube lhs
dimSize <- 7
exb01 <- data.frame(rbind(rep(0, dimSize),
                          randomLHS(98, 7),
                          rep(1, dimSize)))
names(exb01) <- names(train)[1:7]

# put cube lhs in right ranges
exb01 <- sapply(names(exb01), function(col){
  exb01[[col]] * (r[2, col] - r[1, col]) + r[1, col]
})
apply(exb01[,1:7], 2, range)
r

# run TPM
if (!file.exists('./data/exb01Y.csv')) {
  exb01Y <- tpm(exb01, moles = moles, stl = mesh, verb = 1)
  write.csv(exb01Y, './data/exb01Y.csv', row.names = FALSE)
} else {
  read.csv('./data/exb01Y.csv')
}

# c. Train six GP surrogates on the pure species data provided in the directory tpm-git/data/GRACE.
# A few helpful notes: pure He was done already (see §2.3.2), so there are only five left to do;
# be careful to use the same coding scheme for all six sets of inputs; you may wish to double-check
# your pure species predictors on the pure species testing sets provided.

# get train and test datasets for pure species
pureSpeciesFileNames <- list.files("../gramacylab-tpm-5f42493808f6/data/GRACE", ".csv")
pureSpeciesFileNames <- grep('CD', pureSpeciesFileNames, value = TRUE)
pureSpeciesFileNamesTrain <- grep('1000_', pureSpeciesFileNames, value = TRUE)
pureSpeciesFileNamesTest <- grep('100_', pureSpeciesFileNames, value = TRUE)

pureSpeciesTrainList <- sapply(pureSpeciesFileNamesTrain, function(filename) {
  importedData <- read.csv(paste0("../gramacylab-tpm-5f42493808f6/data/GRACE/", filename)) %>%
    select(Umag, Ts, Ta, theta, phi, alphan, sigmat, Cd, Cd_old)
  return(importedData)
}, simplify = FALSE, USE.NAMES = TRUE)

pureSpeciesTrainListCoded <- sapply(pureSpeciesFileNamesTrain, function(filename) {
  browser()
  
  pureSpeciesTrainList[[filename]]
  
  importedData <- read.csv(paste0("../gramacylab-tpm-5f42493808f6/data/GRACE/", filename)) %>%
    select(Umag, Ts, Ta, theta, phi, alphan, sigmat, Cd, Cd_old)
  return(importedData)
}, simplify = FALSE, USE.NAMES = TRUE)



Collect predictions from all six GP surrogates at the testing locations from #b. Combine them into a single prediction for that LEO position and calculate RMSPE to the true simulation outputs from #b. Is the RMSPE close to 1%?
Repeat #b–d for HST with an identical setup except that you’ll need to augment your design with 100 random panel angles in  
{
  0
  ,
  10
  ,
  20
  ,
  …
  ,
  80
  ,
  90
}
and utilize the appropriate mesh files in your simulations.

# -------------------------------------------------------------------------------------------------
# Groundwater remediation
bvv <- read.csv("runlock/pato_results.csv")
cols <- 3:14
nc <- length(cols)
matplot(bvv[1:1000,1], bvv[1:1000,cols], xlim=c(1,1500), type="l", bty="n",
        xlab="number of evaluations", ylab="pump-and-treat cost to remediate",
        lty=1:nc, col=1:nc)
legend("topright", names(bvv)[cols], lty=1:nc, col=1:nc, bty="n")
