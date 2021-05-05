# Steepest ascent and ridge analysis

# setup
library(dplyr)
library(knitr)

# -------------------------------------------------------------------------------------------------
# Chemical process
# Measurements of the response variables were subsequently collected on a central composite design
# whose center point was replicated four times.
chem <- read.table("data/chemical.txt", header=TRUE)
uchem <- unique(chem[,1:2])
reps <- apply(uchem, 1, function(x)
  {
  sum(apply(chem[,1:2], 1, function(y)
    { all(y == x) }))
  })

plot(uchem, type="n")
text(uchem, labels=reps)

# After expanding out into squared and interaction features, the second-order model may be fit to
# these data (in coded units) with the following code.
X <- data.frame(x1=chem$x1,
                x2=chem$x2,
                x11=chem$x1^2,
                x22=chem$x2^2, 
                x12=chem$x1*chem$x2)
y <- chem$y
fit <- lm(y ~ ., data=X)
kable(summary(fit)$coefficients, caption="Chemical process data.")

# fit model and predict on large grid
r <- cbind(c(200, 250), c(15,25))
d <- (r[2,] - r[1,])/2
xi1 <- seq(min(chem$temp), max(chem$temp), length=100)
xi2 <- seq(min(chem$conc), max(chem$conc), length=100)
xi <- expand.grid(xi1, xi2)
x <- cbind((xi[,1] - r[2,1] + d[1])/d[1], (xi[,2] - r[2,2] + d[2])/d[2])
XX <- data.frame(x1=x[,1],
                 x2=x[,2], 
                 x11=x[,1]^2,
                 x22=x[,2]^2,
                 x12=x[,1]*x[,2])
p <- predict(fit, newdata=XX)

# plot predictions
xlab <- "Temperature (°C)"
ylab <- "Concentration (%)"
cols <- heat.colors(128)
image(xi1, xi2, matrix(p, nrow=length(xi1)), 
      col=cols, xlab=xlab, ylab=ylab)
contour(xi1, xi2, matrix(p, nrow=length(xi1)), add=TRUE)

# canonical analysis
b <- coef(fit)[2:3]
B <- matrix(NA, nrow=2, ncol=2)
diag(B) <- coef(fit)[4:5]
B[1,2] <- B[2,1] <- coef(fit)[6]/2
xs <- -(1/2)*solve(B, b)
xs
xis <- xs*d + (r[2,] - d)
xis
E <- eigen(B)
E
# Both eigenvalues of B are negative, so the stationary point is a maximum (as is obvious
# by inspecting the surface in Figure 3.9).

# The code below saves those values, re-ordering them by magnitude, extracts eigenvectors
# (in that order), and converts them to their natural scale for visualization.
lambda <- E$values
o <- order(abs(lambda), decreasing=TRUE)
lambda <- lambda[o]
V <- E$vectors[,o]
Vxi <- V
for(j in 1:ncol(Vxi)) Vxi[,j] <- Vxi[,j]*d*10

image(xi1, xi2, matrix(p, nrow=length(xi1)), col=cols, xlab=xlab, ylab=ylab)
contour(xi1, xi2, matrix(p, nrow=length(xi1)), add=TRUE)
lines(c(-Vxi[1,1], Vxi[1,1])+xis[1], c(-Vxi[2,1], Vxi[2,1])+xis[2], lty=2)
lines(c(-Vxi[1,2], Vxi[1,2])+xis[1], c(-Vxi[2,2], Vxi[2,2])+xis[2], lty=2)
points(xis[1], xis[2])
text(xis[1], xis[2], "xs", pos=4)

# -------------------------------------------------------------------------------------------------
# Ridge analysis

# get data and fit model
rr <- read.table("data/risingridge.txt", header=TRUE)
rr$A2 <- rr$A^2
rr$B2 <- rr$B^2
rr$AB <- rr$A*rr$B
fit <- lm(y ~ ., data=rr)
kable(summary(fit)$coefficients, 
      caption="Linear model summary for rising ridge data.")

# extract b and B, find xs
b <- coef(fit)[2:3]
B <- matrix(NA, nrow=2, ncol=2)
diag(B) <- coef(fit)[4:5]
B[1,2] <- B[2,1] <- coef(fit)[6]/2
xs <- -(1/2)*solve(B, b)
xs

# bounds of design
apply(rr[,1:2], 2, range)
# so xs is outside the design

# Eigenvalues analysis estimate xs as a maximum, but one is too close to zero
E <- eigen(B)
lambda <- E$values
o <- order(abs(lambda), decreasing=TRUE)
V <- E$vectors[,o]*20
lambda <- lambda[o]
lambda

# canonical representation
ys <- coef(fit)[1] + 0.5*t(xs) %*% b
ys

# visuals to support findings so far
x <- seq(-6, 6, length=100)
xx <- expand.grid(x, x)
XX <- data.frame(A=xx[,1], B=xx[,2], A2=xx[,1]^2, 
                 B2=xx[,2]^2, AB=xx[,1]*xx[,2])
p <- predict(fit, newdata=XX)
image(x, x, matrix(p, nrow=length(x)), col=cols, xlab="A", ylab="B")
contour(x, x, matrix(p, nrow=length(x)), add=TRUE)
lines(c(-V[1,1], V[1,1]) + xs[1], c(-V[2,1], V[2,1]) + xs[2], lty=2)
lines(c(-V[1,2], V[1,2]) + xs[1], c(-V[2,2], V[2,2]) + xs[2], lty=2)
polygon(c(1,1,-1,-1), c(1,-1,-1,1), lty=3)
text(0,-0.5, "design region", cex=0.5)
points(xs[1], xs[2])
text(xs[1], xs[2], "xs", pos=4)

# example of chemical process (saddle topology)
# load and massage data
saddle <- read.table("data/saddle.txt", header=TRUE)
saddle <- cbind(saddle[,-5]^2, 
                model.matrix(~ .^2 - 1, saddle[,-5]), y=saddle[,5])
names(saddle)[1:4] <- paste("x", 1:4, 1:4, sep="")
names(saddle)[9:14] <- sub(":x", "", names(saddle)[9:14])
saddle <- saddle[c(5:8,1:4,9:15)]

# few useful features, proceed with caution
fit <- lm(y ~ ., data=saddle)
summary(fit)

# extract coefficient main effects vector b
# and matrix  B of second-order terms
b <- coef(fit)[2:5]
B <- matrix(NA, nrow=4, ncol=4)
diag(B) <- coef(fit)[6:9]
i <- 10
for(j in 1:3) {
  for(k in (j+1):4) {
    B[j,k] <- B[k,j] <- coef(fit)[i]/2
    i <- i + 1
  }
}

# Using those coefficients, stationary point xs
# may be calculated as follows.
xs <- -(1/2)*solve(B, b)
xs
apply(saddle[,1:4], 2, range)

# calculate eigenvalues
E <- eigen(B)
lambda <- E$values
o <- order(abs(lambda), decreasing=TRUE)
lambda <- lambda[o]
lambda

# choose mu
mul <- max(lambda)
x <- solve(B - mul*diag(4), -b/2)
x

# function to assist the search of a better mu
f <- function(mu, R2=1) 
{
  x <- solve(B - mu*diag(4), -b/2)
  R2 - t(x) %*% x
}

# find good mu
mu <- uniroot(f, c(mul, 10*mul), R2=1.4^2)$root
mu

# define new x
x <- solve(B - mu*diag(4), -b/2)
x

drop(sqrt(t(x) %*% x))

mus <- rs <- seq(0.1, 2, length=20)
xp <- matrix(NA, nrow=length(rs), ncol=4)
colnames(xp) <- c("x1", "x2", "x3", "x4")
for(i in 1:length(rs)) {
  mus[i] <- uniroot(f, c(mul, 100*mul), R2=rs[i]^2)$root
  xp[i,] <- solve(B - mus[i]*diag(4), -b/2)
}
xp <- rbind(rep(0,4), xp)
rs <- c(0, rs)
mus <- c(Inf, mus)

Xp <- data.frame(xp)
Xp <- cbind(Xp^2, model.matrix(~ .^2 - 1, Xp))
names(Xp)[1:4] <- paste("x", 1:4, 1:4, sep="")
names(Xp)[9:14] <- sub(":x", "", names(Xp)[9:14])
Xp <- Xp[c(5:8,1:4,9:14)]

p <- predict(fit, newdata=Xp, se.fit=TRUE)

cbind(R=rs, mu=mus, data.frame(pred=p$fit, se=p$se.fit), round(xp,6))

plot(rs, p$fit, type="b", ylim=c(20,100), xlab="radius (R)", 
     ylab="y.hat(x) & 95% CIs")
lines(rs, p$fit + 2*p$se, col=2, lty=2)
lines(rs, p$fit - 2*p$se, col=2, lty=2)

# investigating sampling properties
# creating a 2d subset of the chemical saddle data
x12 <- seq(-2,2,length=100)
g <- expand.grid(x12, x12)
Xp <- data.frame(x1=g[,1], x2=g[,2], 
                 x11=g[,1]^2, x22=g[,2]^2, x12=g[,1]*g[,2])

## x3 = x4 = 0
Xp$x3 <- Xp$x33 <- Xp$x4 <- Xp$x44 <- Xp$x34 <- 0
Xp$x13 <- Xp$x14 <- Xp$x23 <- Xp$x24 <- Xp$x34 <- 0
p0 <- predict(fit, newdata=Xp, se.fit=TRUE)

## x3 = x4 = 1
Xp$x3 <- Xp$x33 <- Xp$x4 <- Xp$x44 <- Xp$x34 <- 1
Xp$x13 <- Xp$x14 <- Xp$x1; Xp$x23 <- Xp$x24 <- Xp$x2
p1 <- predict(fit, newdata=Xp, se.fit=TRUE)

# produce plot
par(mfrow=c(1,2))
bs <- seq(min(p0$se.fit), max(p0$se.fit), length=129)
image(x12, x12, matrix(p0$se.fit, nrow=length(x12)), col=cols, breaks=bs,
      xlab="x1", ylab="x2", main="se.fit, x3 = x4 = 0")
contour(x12, x12, matrix(p0$se.fit, nrow=length(x12)), add=TRUE)
points(saddle[,1:2])
points(saddle[apply(saddle[,3:4] == c(0,0), 1, all),1:2], pch=19)
image(x12, x12, matrix(p1$se.fit, nrow=length(x12)), col=cols, breaks=bs,
      xlab="x1", ylab="x2", main="se.fit, x3 = x4 = 1")
contour(x12, x12, matrix(p1$se.fit, nrow=length(x12)), add=TRUE)
points(saddle[,1:2])
points(saddle[apply(saddle[,3:4] == c(1,1), 1, all), 1:2], pch=19)

# confidence in the sationaty point
crdat <- read.table("data/confreg.txt", header=TRUE)

# code below combines several steps: expanding to second-order features,
# model fitting and extracting b and B, differentiating and solving to
# estimate stationary point xs
crdat$x11 <- crdat$x1^2
crdat$x22 <- crdat$x2^2
crdat$x12 <- crdat$x1 * crdat$x2
fit <- lm(y ~ ., data=crdat)
summary(fit)
b <- coef(fit)[2:3]
B <- matrix(NA, nrow=2, ncol=2)
diag(B) <- coef(fit)[4:5]
B[1,2] <- B[2,1] <- coef(fit)[6]/2
xs <- -(1/2)*solve(B, b)
xs
eigen(B)$values # stationary point is a max

# Next evaluate predictive equations on a grid in order to visualize the
# response surface, design and stationary point, which is well within
# the interior of the experimental region
par(mfrow=c(1,1))
xx <- seq(-1.6, 1.6, length=200)
g <- expand.grid(xx, xx)
XX <- data.frame(x1=g[,1], x2=g[,2], 
                 x11=g[,1]^2, x22=g[,2]^2, x12=g[,1]*g[,2])
p <- as.numeric(predict(fit, newdata=XX))
image(xx, xx, matrix(p, nrow=length(xx)), col=cols, xlab="x1", ylab="x2")
contour(xx, xx, matrix(p, nrow=length(xx)), add=TRUE)
points(crdat$x1, crdat$x2, pch=20)
points(xs[1], xs[2])
text(xs[1], xs[2], "xs", pos=4)

# Combining  d(t) = b + 2Bt
d <- rbind(c(b[1], 2*B[,1]), c(b[2], 2*B[,2]))
colnames(d) <- c("(Intercept)", "t1", "t2")
d

# variance of d(t)
s2 <- summary(fit)$sigma^2
s2

# on n - p = 7 DoF
X <- cbind(1, as.matrix(crdat[,-3]))
XtXi <- solve(t(X) %*% X)
XtXi

# function to determine the covariance matrix as a function of its elements,
# and coordinates t
Vard <- function(t1, t2, s2, XtXi) 
{
  v11 <- XtXi[2,2] + 4*t1^2*XtXi[4,4] + t2^2*XtXi[6,6]
  v22 <- XtXi[3,3] + t1^2*XtXi[6,6] + 4*t2^2*XtXi[5,5]
  v12 <- v21 <- 4*t1*t2*XtXi[4,5] + t1*t2*XtXi[6,6]
  v <- s2 * matrix(c(v11,v12,v21,v22), ncol=2, byrow=TRUE)
  return(v)
}

# quadratic form of Eq 3.7
CIqf <- function(t1, t2, s2, XtXi, b, B)
{
  dt <- b + 2*B %*% c(t1, t2)
  V <- Vard(t1, t2, s2, XtXi)
  Vi <- solve(V)
  t(dt) %*% Vi %*% dt
}

# evaluate Eq 3.8
quadform <- rep(NA, nrow(g))
for(i in 1:nrow(g)) quadform[i] <- CIqf(g[i,1], g[i,2], s2, XtXi, b, B)

ci90 <- quadform <= 2*qf(0.9, 2, nrow(X)-ncol(X))
ci95 <- quadform <= 2*qf(0.95, 2, nrow(X)-ncol(X))
image(xx, xx, matrix(ci90 + ci95, ncol=length(xx)), xlab="x1", ylab="x2", 
      col=c("white", "lightgray", "darkgray"))
text(c(-0.2,-1), c(-0.5,1), c("90%","95%"))

summary(fit)$r.squared

# artificially improving the dataset to see response of CI
crdat$y <- c(87.6, 86.5, 85.7, 86.9, 86.7, 86.8, 87.4, 86.7, 90.3, 
             91.0, 90.8) 
fit2 <- lm(y ~ ., data=crdat)
summary(fit2)$r.squared # much better fit

# updating calculations
b <- coef(fit)[2:3]
B <- matrix(NA, nrow=2, ncol=2)
diag(B) <- coef(fit)[4:5]
B[1,2] <- B[2,1] <- coef(fit)[6]/2
xs <- -(1/2)*solve(B, b)
s2 <- summary(fit2)$sigma^2
for(i in 1:nrow(g)) quadform[i] <- CIqf(g[i,1], g[i,2], s2, XtXi, b, B)

ci95 <- quadform <= 2*qf(0.95, 2, nrow(X)-ncol(X))
image(xx, xx, matrix(ci95, ncol=length(xx)), xlab="x1", ylab="x2",
      col=c("white", "lightgray"))

# IMPORTANT REMARK
# This state of affairs is merely a reflection of the reality that the response
# surface is locally flat. A positive spin may be that such situations imply a
# degree of flexibility in choosing the optimum. Communicating that effectively
# to stakeholders requires care, of course, but is definitely of interest.
