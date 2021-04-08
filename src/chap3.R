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
