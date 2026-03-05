hello_kmeans <- function(seed = 1, n = 200) {
  set.seed(seed)
  X <- cbind(rnorm(n), rnorm(n))
  km <- kmeans(X, centers = 2, nstart = 10)
  list(X = X, km = km, seed = seed, n = n)
}