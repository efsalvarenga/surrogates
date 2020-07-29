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
