rm(list=ls())
source("R/load_all.R")
set.seed(1)

n <- 200

x1 <- cbind(
  rnorm(n, mean = -2, sd = 0.5),
  rnorm(n, mean =  0, sd = 50)
)

x2 <- cbind(
  rnorm(n, mean =  2, sd = 0.5),
  rnorm(n, mean =  0, sd = 0.5)
)

X <- rbind(x1, x2)
true_lab <- rep(1:2, each = n)

km1 <- kmeans(X, centers = 2, nstart = 20)


X_scaled <- scale(X)
km2 <- kmeans(X_scaled, centers = 2, nstart = 20)

par(mfrow = c(1, 2))

plot(
  X,
  col = km1$cluster,
  pch = 17+true_lab,
  main = "K-means unscaled",
  xlab = "",
  ylab = "",
  axes=FALSE, cex.main=1.5, cex.lab=1.3
)
points(km1$centers, pch = 9, cex = 2, lwd = 2, col=5)
box()
plot(
  X_scaled,
  col = km2$cluster,
  pch = 17+true_lab,
  main = "K-means scaled",
  xlab = "",
  ylab = "",
  axes=FALSE, cex.main=1.5, cex.lab=1.3
)
points(km2$centers, pch = 9, cex = 2, lwd = 2, col=5)
box()




