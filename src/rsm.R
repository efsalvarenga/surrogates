# RSM (Chap 1)

# chemical yield model
yield <- function(xi1, xi2) 
{
  xi1 <- 3*xi1 - 15
  xi2 <- xi2/50 - 13
  xi1 <- cos(0.5)*xi1 - sin(0.5)*xi2
  xi2 <- sin(0.5)*xi1 + cos(0.5)*xi2
  y <- exp(-xi1^2/80 - 0.5*(xi2 + 0.03*xi1^2 - 40*0.03)^2)
  return(100*y)
}

# generate values and evaluates chemical yield model
xi1 <- seq(1, 8, length=100)
xi2 <- seq(100, 1000, length=100)
g <- expand.grid(xi1, xi2)
y <- yield(g[,1], g[,2])

# 3d view
persp(xi1, xi2, matrix(y, ncol=length(xi2)), theta=45, phi=45, 
      lwd=0.5, xlab="xi1 : time", ylab="xi2 : temperature", 
      zlab="yield", expand=0.4)

# contour plot
cols <- heat.colors(128)
image(xi1, xi2, matrix(y, ncol=length(xi2)), col=cols, 
      xlab="xi1 : time", ylab="xi2 : temperature")
contour(xi1, xi2, matrix(y, ncol=length(xi2)), nlevels=4, add=TRUE)

# 1st order model
first.order <- function(x1, x2)
{
  50 + 8*x1 + 3*x2  
}

# generates values and evaluates 1st order model
x1 <- x2 <- seq(-1, 1, length=100)
g <- expand.grid(x1, x2)
eta1 <- matrix(first.order(g[,1], g[,2]), ncol=length(x2))

# views
par(mfrow=c(1,2))
persp(x1, x2, eta1, theta=30, phi=30, zlab="eta", expand=0.75, lwd=0.25)
image(x1, x2, eta1, col=heat.colors(128))
contour(x1, x2, matrix(eta1, ncol=length(x2)), add=TRUE)

# 1st order with interaction model
first.order.i <- function(x1, x2)
{
  50 + 8*x1 + 3*x2 - 4*x1*x2
}
eta1i <- matrix(first.order.i(g[,1], g[,2]), ncol=length(x2))

# views
par(mfrow=c(1,2))
persp(x1, x2, eta1i, theta=30, phi=30, zlab="eta", expand=0.75, lwd=0.25)
image(x1, x2, eta1i, col=heat.colors(128))
contour(x1, x2, matrix(eta1i, ncol=length(x2)), add=TRUE)

# 2nd order model
simple.max <- function(x1, x2)
{
  50 + 8*x1 + 3*x2 - 7*x1^2 - 3*x2^2 - 4*x1*x2
}
eta2sm <- matrix(simple.max(g[,1], g[,2]), ncol=length(x2))
par(mfrow=c(1,2))
persp(x1, x2, eta2sm, theta=30, phi=30, zlab="eta", expand=0.75, lwd=0.25)
image(x1, x2, eta2sm, col=heat.colors(128))
contour(x1, x2, eta2sm, add=TRUE)

# stationary ridge
stat.ridge <- function(x1, x2)
{
  80 + 4*x1 + 8*x2 - 3*x1^2 - 12*x2^2 - 12*x1*x2
}
eta2sr <- matrix(stat.ridge(g[,1], g[,2]), ncol=length(x2))
par(mfrow=c(1,2))
persp(x1, x2, eta2sr, theta=30, phi=30, zlab="eta", expand=0.75, lwd=0.25)
image(x1, x2, eta2sr, col=heat.colors(128))
contour(x1, x2, eta2sr, add=TRUE)

# rising ridge
rise.ridge <- function(x1, x2) 
{
  80 - 4*x1 + 12*x2 - 3*x1^2 - 12*x2^2 - 12*x1*x2
}
eta2rr <- matrix(rise.ridge(g[,1], g[,2]), ncol=length(x2))
par(mfrow=c(1,2))
persp(x1, x2, eta2rr, theta=30, phi=30, zlab="eta", expand=0.75, lwd=0.25)
image(x1, x2, eta2rr, col=heat.colors(128))
contour(x1, x2, eta2rr, add=TRUE)

# saddle
saddle <- function(x1, x2) 
{
  80 + 4*x1 + 8*x2 - 2*x1 - 12*x2 - 12*x1*x2 
}
eta2s <- matrix(saddle(g[,1], g[,2]), ncol=length(x2))
par(mfrow=c(1,2))
persp(x1, x2, eta2s, theta=30, phi=30, zlab="eta", expand=0.75, lwd=0.25)
image(x1, x2, eta2s, col=heat.colors(128))
contour(x1, x2, eta2s, add=TRUE)

# Aircraft wing weight example
wingwt <- function(Sw=0.48, Wfw=0.4, A=0.38, L=0.5, q=0.62, l=0.344, 
                   Rtc=0.4, Nz=0.37, Wdg=0.38)
{
  ## put coded inputs back on natural scale
  Sw <- Sw*(200 - 150) + 150 
  Wfw <- Wfw*(300 - 220) + 220 
  A <- A*(10 - 6) + 6 
  L <- (L*(10 - (-10)) - 10) * pi/180
  q <- q*(45 - 16) + 16 
  l <- l*(1 - 0.5) + 0.5 
  Rtc <- Rtc*(0.18 - 0.08) + 0.08
  Nz <- Nz*(6 - 2.5) + 2.5
  Wdg <- Wdg*(2500 - 1700) + 1700
  
  ## calculation on natural scale
  W <- 0.036*Sw^0.758 * Wfw^0.0035 * (A/cos(L)^2)^0.6 * q^0.006 
  W <- W * l^0.04 * (100*Rtc/cos(L))^(-0.3) * (Nz*Wdg)^(0.49)
  return(W)
}

x <- seq(0, 1, length=100)
g <- expand.grid(x, x)

W.A.Nz <- wingwt(A=g[,1], Nz=g[,2])

cs <- heat.colors(128)
bs <- seq(min(W.A.Nz), max(W.A.Nz), length=129)

image(x, x, matrix(W.A.Nz, ncol=length(x)), col=cs, breaks=bs, 
      xlab="A", ylab="Nz")
contour(x, x, matrix(W.A.Nz, ncol=length(x)), add=TRUE)

W.l.Wfw <- wingwt(l=g[,1], Wfw=g[,2])

image(x, x, matrix(W.l.Wfw,ncol=length(x)), col=cs, breaks=bs, 
      xlab="l", ylab="Wfw")
contour(x,x, matrix(W.l.Wfw,ncol=length(x)), add=TRUE)

library(lhs)
n <- 1000
X <- data.frame(randomLHS(n, 9))
names(X) <- names(formals(wingwt))

plot(X[,1:2], pch=19, cex=0.5)
abline(h=c(0.6, 0.8), col=2, lwd=2)

inbox <- X[,1] > 0.6 &  X[,1] < 0.8 
sum(inbox)/nrow(X)

Y <- wingwt(X[,1], X[,2], X[,3], X[,4], X[,5], X[,6], X[,7], X[,8], X[,9])

fit.lm <- lm(log(Y) ~ .^2, data=data.frame(Y,X))
fit.lmstep <- step(fit.lm, scope=formula(fit.lm), direction="backward", 
                   k=log(length(Y)), trace=0)

library(laGP)
fit.gp <- newGPsep(X, Y, 2, 1e-6, dK=TRUE)
mle <- mleGPsep(fit.gp)

baseline <- matrix(rep(as.numeric(formals(wingwt)), nrow(g)), 
                   ncol=9, byrow=TRUE)
XX <- data.frame(baseline)
names(XX) <- names(X)
XX$A <- g[,1]
XX$Nz <- g[,2]

p <- predGPsep(fit.gp, XX, lite=TRUE)

image(x, x, matrix(p$mean, ncol=length(x)), col=cs, breaks=bs, 
      xlab="A", ylab="Nz")
contour(x, x, matrix(p$mean, ncol=length(x)), add=TRUE)

meq1 <- meq2 <- me <- matrix(NA, nrow=length(x), ncol=ncol(X))
for(i in 1:ncol(me)) {
  XX <- data.frame(baseline)[1:length(x),]
  XX[,i] <- x
  p <- predGPsep(fit.gp, XX, lite=TRUE)
  me[,i] <- p$mean
  meq1[,i] <- qt(0.05, p$df)*sqrt(p$s2) + p$mean
  meq2[,i] <- qt(0.95, p$df)*sqrt(p$s2) + p$mean
}

matplot(x, me, type="l", lwd=2, lty=1, col=1:9, xlab="coded input")
matlines(x, meq1, type="l", lwd=2, lty=2, col=1:9)
matlines(x, meq2, type="l", lwd=2, lty=2, col=1:9)
legend("topleft", names(X)[1:5], lty=1, col=1:5, horiz=TRUE, 
       bty="n", cex=0.5)
legend("bottomright", names(X)[6:9], lty=1, col=6:9, horiz=TRUE, 
       bty="n", cex=0.5)

# Exercise 1
libray(dplyr)

wireRaw <- read.csv('http://bobby.gramacy.com/surrogates/wire.csv')

wireUnitCoefs <- sapply(names(wireRaw), function(selCol){
  linCoef <- min(wireRaw[[selCol]])
  angCoef <- max(wireRaw[[selCol]] - linCoef)
  return(c(linCoef, angCoef))
}) %>% data.frame()
row.names(wireUnitCoefs) <- c('linCoef', 'angCoef')

wireUnit <- sapply(names(wireRaw), function(selCol){
  return((wireRaw[[selCol]] - wireUnitCoefs['linCoef', selCol]) / wireUnitCoefs['angCoef', selCol])
}) %>% data.frame()

wireZeroOne <- sapply(names(wireRaw), function(selCol){
  linMean <- mean(wireRaw[[selCol]])
  angCoef <- max(linMean - min(wireRaw[[selCol]]), max(wireRaw[[selCol]]) - linMean)
  return((wireRaw[[selCol]] - linMean) / angCoef)
}) %>% data.frame()

wireY <- 'pstren'
pstren <- wireRaw[[wireY]]
wireX <- names(wireRaw)[!(names(wireRaw) %in% wireY)]

wireFit.lm <- lm(pstren ~ ., data = data.frame(pstren, wireUnit[wireX]))
wireFit.qlm <- lm(pstren ~ . ^ 2, data = data.frame(pstren, wireUnit[wireX]))
wireFit.lmstep <- step(wireFit.lm,
                       scope = formula(wireFit.qlm),
                       direction = "both", 
                       k = log(length(pstren)))


newValues <- data.frame(t(c(6, 20, 30, 90, 2, 2)))
names(newValues) <- wireX
newValuesUnit <- sapply(names(newValues), function(selCol){
  return((newValues[[selCol]] - wireUnitCoefs['linCoef', selCol]) / wireUnitCoefs['angCoef', selCol])
}) %>% t() %>% data.frame()

wirePred <- predict(wireFit.lm, newdata = newValuesUnit, interval  = 'confidence')
