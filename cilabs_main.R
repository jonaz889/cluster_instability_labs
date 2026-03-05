# WELCOME TO JONAS CLUSTER INSTABILITY LAB - CILAB 
#This file containts initial K-means experiments and investigations
#

## Reset
rm(list=ls())
setwd("C:/Users/jonas/Desktop/Learning/DTU/Term6/bachelor/algorithms/data")


# Libraries
library(MASS)  # mvnorm
library(wesanderson) # Fun colors
library(paletteer) # MORE COLORS
library(aricode) # Adjusted rand_index
library(clue) # Label alignment
library(gtools) # For permutations (instead of clue package)
library(usethis) # For version control

## Helper functions
# Computes the RI of two clusterings
rand_index <- function(c1,c2)
{
  n = length(c1)
  
  a <- 0 # Same Same
  b <- 0 # Different Different
  c <- 0 # Same Different
  d <- 0 # Different Same
  for (i in 1:(n-1)) 
  {
    for (j in (i+1):n) 
    {
      
      if(c1[i]==c1[j])
      {
        if(c2[i]==c2[j]){a = a + 1} # Same Same
        else {c = c + 1} # Same Different
      }
      else 
      {
        if(c2[i]==c2[j]){d = d + 1} # Different Same
        else {b = b + 1} # Different Different
      }
    }
  }
  ri <- (a+b)/(a+b+c+d) # Rand Index
  return(ri)
}

## Data generations
generate_gaussian_cluster_data_2D <- function (n_list, mu_list, sigma_list)
{
  k = length(n_list)
  if ( k != length(mu_list) || k != length(sigma_list) )
  {
    stop("Given means and std. deviations have different lenghts!")
  }
  if ( k == 0)
  {
    stop("Means has length 0 - cannot cluster")
  }
  
  samples <- Map(function(n, mu, Sigma) mvrnorm(n = n, mu = mu, Sigma = Sigma),
                 n_list, mu_list, sigma_list) # Returns list of function calls
  
  
  # Add cluster labels (Different for each list element, remember each list 
  #element is a new cluster)
  X_mat = lapply(seq(k) ,function(x) { cbind(samples[[x]], rep(x, n_list[x]))})
  # "flatten" the list
  X_mat <-do.call(rbind, X_mat)
  
  # Dataframe  and return it
  X_df = data.frame(
    x1 = X_mat[,1],
    x2 = X_mat[,2],
    label = X_mat[,3]
  )
  return(X_df)
}

## WHY DID I DO THIS  ( DELETE ? ? ? ?)
# Calculates WCSS given a df w/ col1=x1, col2=x2, col3=label
wcss_given_labels <- function (X)
{
  cluster_centers <- aggregate(cbind(x1,x2) ~ label, data=X, mean)
  sum(apply(X, 1, function(x) {
    center <- cluster_centers[x[["label"]], 2:3]
    sum((x[1:2] - center)^2)
  }))
}

## Experiments
# Graphical setup
pchs <- c(19, 17, 15, 18, 8, 3, 4,1, 2, 0, 5, 6, 7, 9,10, 11, 12, 13, 14,16)
cols = paletteer_d("ggthemes::Tableau_20")

####################################################
## Stability under different starts#################
## #################################################
# Question: Does Hartigan-Wong's k-means heuristic converge to different 
# solutions w/ similar WCSS.
#
# Inspired by supervision 1: What if we have symetri e.g. 3 symmetric gaussian 
# clouds w/ 2 center k-means initialization
#
# How: (1) Make generate a data scenario 
#      (2) Cluster it N times w/ random starts
#      (3) Check WCSS vs some measure (MAYBE ARI??)
set.seed(532)  # reproducibility

## Helper functions¨

##### Function for plotting near minima solutions
# nearness - Takes how close to optimum in %
# wcss - list of wcss scores
# clusterings - actual clusterings
plot_near_minima_sols <- function(clusterings, wcss, nearness) 
{
  best_wcss_sol <- clusterings[[which.min(wcss)]]$cluster
  near_best_wcss_sol <- which(wcss <= min(wcss) * (1+nearness/100))
  N <- length(clusterings)
  ## So this loop works by finding one 
  #  candidate from each group of solutions
  cand <- c() # Just uses index
  samesies <- c()
  for (i in near_best_wcss_sol) 
  {
    in_cand = FALSE
    for(j in cand)
    {
      if (ARI(clusterings[[i]]$cluster, clusterings[[j]]$cluster)==1)
      {
        in_cand = TRUE
        break
      }
    }
    if(!in_cand){cand <- c(cand, i)}
  } # OPTIM. IDEA: USE UNIQUE INSTEAD OF ABOVE
  
  # Now count #
  # Loop over each clustering and see if they are the same and count them
  for(i in cand)
  {
    samesies_counter <- 0 # counter of identicals (NOTE: 0 because it also counts itself)
    for (j in 1:N) 
    {
      if (ARI(clusterings[[i]]$cluster, clusterings[[j]]$cluster) == 1) 
      {
        samesies_counter <- samesies_counter + 1 # If same level up counter
      }
    }
    samesies <- c(samesies, samesies_counter)
  }
  ## Sort candidates by WCSS
  ord <- order(wcss[cand])
  cand <- cand[ord]
  samesies <- samesies[ord]
  
  cat("Total w/in wcss ball:", length(near_best_wcss_sol), "\n")
  cat("#distinct solutions:", length(cand), "\n")
  
  
  #Plotting
  k <- length(unique(best_wcss_sol)) # Gets #centers for this run
  par(mfrow=c(ceiling(length(cand) / ceiling(sqrt(length(cand)))), ceiling(sqrt(length(cand)))),mar=c(2,2,2.5,1))
  
  # plot each candidate and highlight different points
  for (i in 1:length(cand)) {
    cl <- clusterings[[cand[i]]]$cluster
    plot(X$x1, X$x2,
         col  = cols[cl],
         pch  = pchs[X$label],
         cex  = ifelse((cl != best_wcss_sol), 2, 1),
         main = paste0("#=",samesies[i],"/",N,"//wcss=", round(wcss[cand[i]], 1),"//ari=",round(ARI_to_best[cand[i]],3),"//#d=", sum(cl != best_wcss_sol)),
         xlab = expression(x[1]), ylab = expression(x[2]))
  }
  return(list(cand=cand, counts=samesies, near=near_best_wcss_sol, best=which.min(wcss)))
}


# Experiment 1 - two close clusters

## Settings
# Data settings
dist <- 10
n_list  <- c(200, 200)
mu_list<- list(c(-dist/2, 0), c(dist/2, 0))
sigma_list <- list(diag(2)*1, diag(2)*1)

# Clustering settings
k <- 2
N <- 1000

## Data generation
X <- generate_gaussian_cluster_data_2D(n_list, mu_list, sigma_list)

## Run k-means many times w/ random start
clusterings <- list()
wcss <- numeric(N)

kmeans_1 <- kmeans(X[, 1:2], centers = k)

kmeans_runs <- list()
X_list <- list()
for(i in 10:1)
{
  dist <- i
  n_list  <- c(200, 200)
  mu_list<- list(c(-dist/2, 0), c(dist/2, 0))
  
  X_list[[10-i+1]] <- generate_gaussian_cluster_data_2D(n_list, mu_list, sigma_list)
  kmeans_runs[[10 - i +1]] <- kmeans(X_list[[10-i+1]][, 1:2], centers=k) 
}

par(mfrow=c(5, 1))
for(i in 6:10)
{
  plot(X_list[[i]][,1:2], col=kmeans_runs[[i]]$cluster, pch=X_list[[i]][,3])
? kmean}


plot(X[,1:2], col=kmeans_1$cluster, pch=X[,3])



for (i in 1:N) {

  clusterings[[i]] <- kmeans(X[, 1:2], centers = k, nstart = 1)
  wcss[i] <- clusterings[[i]]$tot.withinss
}

## Best WCSS solution
best_wcss_sol <- clusterings[[which.min(wcss)]]$cluster

ARI_to_best <- numeric(N)
for (i in 1:N) {
  ARI_to_best[i] <- ARI(best_wcss_sol, clusterings[[i]]$cluster)
}


## Diagnostics
# ARI vs. wcss plot
par(mfrow=c(1,1))
plot(wcss, ARI_to_best,
     xlab = "WCSS",
     ylab = "ARI vs best solution",
     pch  = 19)

abline(v = min(wcss), lty = 2)

## Best clustering
plot(X$x1, X$x2,
     col  = cols[best_wcss_sol],
     pch  = pchs[X$label],
     xlab = expression(x[1]),
     ylab = expression(x[2]),
     main = paste("Best solution - #best=", sum(wcss==min(wcss) & ARI_to_best==1), "/",N, "- WCSS =", round(wcss[which.min(wcss)], 2)))

res<-plot_near_minima_sols(clusterings = clusterings, wcss = wcss, 5)


# Circle small radius 


# Experiment 2 - Equilateral triangle ----

## Settings
# Data settings
side_length = 2.5
n_list     <- c(200, 200, 200)
mu_list    <- list(c(-side_length/2, 0), c(side_length/2, 0), c(0,side_length*sqrt(3)/2))
sigma_list <- list(diag(2)*1, diag(2)*1, diag(2)*1)

# Clustering settings
k <- 2
N <- 5000

## Data generation
X <- generate_gaussian_cluster_data_2D(n_list, mu_list, sigma_list)

## Run k-means many times w/ random start
clusterings <- list()
wcss <- numeric(N)

for (i in 1:N) {
  clusterings[[i]] <- kmeans(X[, 1:2], centers = k, nstart = 1)
  wcss[i] <- clusterings[[i]]$tot.withinss
}

## Best WCSS solution
best_wcss_sol <- clusterings[[which.min(wcss)]]$cluster

ARI_to_best <- numeric(N)
for (i in 1:N) {
  ARI_to_best[i] <- ARI(best_wcss_sol, clusterings[[i]]$cluster)
}

## Diagnostics
# ARI vs. wcss plot
plot(wcss, ARI_to_best,
     xlab = "WCSS",
     ylab = "ARI vs best solution",
     pch  = 19)

abline(v = min(wcss), lty = 2)

## Best clustering
plot(X$x1, X$x2,
     col  = cols[best_wcss_sol],
     pch  = pchs[X$label],
     xlab = expression(x[1]),
     ylab = expression(x[2]),
     main = paste("Best solution - #best=", sum(wcss==min(wcss) & ARI_to_best==1), "/",N, "- WCSS =", round(wcss[which.min(wcss)], 2)))

res<-plot_near_minima_sols(clusterings = clusterings, wcss = wcss, 5)



# Experiment 3 - Isosceles triangle
## Settings
# Data settings
base_length = 10
height = 0
n_list     <- c(200, 200, 200)
mu_list    <- list(c(-base_length/2, 0), c(base_length/2, 0), c(0,height))
sigma_list <- list(diag(2)*1, diag(2)*1, diag(2)*1)

# Clustering settings
k <- 2
N <- 5000

## Data generation
X <- generate_gaussian_cluster_data_2D(n_list, mu_list, sigma_list)

## Run k-means many times w/ random start
clusterings <- list()
wcss <- numeric(N)

for (i in 1:N) {
  clusterings[[i]] <- kmeans(X[, 1:2], centers = k, nstart = 1)
  wcss[i] <- clusterings[[i]]$tot.withinss
}

## Best WCSS solution
best_wcss_sol <- clusterings[[which.min(wcss)]]$cluster

ARI_to_best <- numeric(N)
for (i in 1:N) {
  ARI_to_best[i] <- ARI(best_wcss_sol, clusterings[[i]]$cluster)
}

## Diagnostics
# ARI vs. wcss plot
plot(wcss, ARI_to_best,
     xlab = "WCSS",
     ylab = "ARI vs best solution",
     pch  = 19)

abline(v = min(wcss), lty = 2)

## Best clustering
plot(X$x1, X$x2,
     col  = cols[best_wcss_sol],
     pch  = pchs[X$label],
     xlab = expression(x[1]),
     ylab = expression(x[2]),
     main = paste("Best solution - #best=", sum(wcss==min(wcss) & ARI_to_best==1), "/",N, "- WCSS =", round(wcss[which.min(wcss)], 2)))

res<-plot_near_minima_sols(clusterings = clusterings, wcss = wcss, 30)



# Rest ----
sum(RIs != 1)
plot(RIs)
# WCSS plot
print(clustering$tot.withinss)

## Cluster plot - Color by cluster, point by label
k = 2
cols = paletteer_d("ggthemes::Tableau_20")
#wes_palette(n=k, "Darjeeling1")# wes_palette(n=5, "Darjeeling1") # See color palette
par(mfrow=c(1,1))

plot(cbind(X$x1, X$x2),
       col = cols[clusterings[[60]]$cluster],
       pch = pchs[X$label],
       xlab = expression(x[1]),
       ylab = expression(x[2]),
       main = paste("k =", k))


# E1 <- ELBOW FUN
n_list = c(100, 50, 50)
mu_list = list(c(0,0), c(-5, 10), c(5, 10))
sigma_list = list(diag(2)*3,diag(2)*3,diag(2)*3)


cdat <- generate_gaussian_cluster_data_2D(n_list, mu_list, sigma_list)

## Clustering
X = generate_gaussian_cluster_data_2D(n_list, mu_list, sigma_list)

clusterings = list()
Ks = 1:20
for (k in Ks)
{
  clusterings[[k]] = kmeans(X[,1:2],centers = k) 
}
# WCSS plot
wcss_dat = unlist(lapply(X=clusterings, FUN=function(x){return(x$tot.withinss)}))

plot(Ks, wcss_dat, type="b", main="Elbow-plot", xlab="k", ylab="Total WCSS")

## Cluster plot - Color by cluster, point by label
k = 6
cols = paletteer_d("ggthemes::Tableau_20")
#wes_palette(n=k, "Darjeeling1")# wes_palette(n=5, "Darjeeling1") # See color palette
par(mfrow=c(2,3))
for (k in 1:6) {
  plot(cbind(X$x1, X$x2),
       col = cols[clusterings[[k]]$cluster],
       pch = pchs[X$label],
       xlab = expression(x[1]),
       ylab = expression(x[2]),
       main = paste("k =", k))
  print(k)
}

## E2 - Clustering Noise
n_list = c(1000)
mu_list = list(c(0,0))
sigma_list = list(diag(2)*3)

## Clustering
X = generate_gaussian_cluster_data_2D(n_list, mu_list, sigma_list)

clusterings = list()
Ks = 1:8
for (k in Ks)
{
  clusterings[[k]] = kmeans(X[,1:2],centers = k, nstart=10000) 
}
# WCSS plot
wcss_dat = unlist(lapply(X=clusterings, FUN=function(x){return(x$tot.withinss)}))

plot(Ks, wcss_dat, type="b", main="Elbow-plot", xlab="k", ylab="Total WCSS")

## Cluster plot - Color by cluster, point by label
k = 6
cols = paletteer_d("ggthemes::Tableau_20")
#wes_palette(n=k, "Darjeeling1")# wes_palette(n=5, "Darjeeling1") # See color palette
par(mfrow=c(2,3))
for (k in 1:6) {
  plot(cbind(X$x1, X$x2),
       col = cols[clusterings[[k]]$cluster],
       pch = pchs[X$label],
       xlab = expression(x[1]),
       ylab = expression(x[2]),
       main = paste("k =", k))
  print(k)
}



## E3 - Symmetric cluster fun
v=12
n_list = c(100, 100, 100)
mu_list = list(c(-10,0), c(10,0), c(0,20))
sigma_list = list(diag(2)*v,diag(2)*v, diag(2)*v)

## Clustering
X = generate_gaussian_cluster_data_2D(n_list, mu_list, sigma_list)

clusterings = list()
Ks = 1:6
for (k in Ks)
{
  clusterings[[k]] = kmeans(X[,1:2],centers = k, nstart=1000) 
}
# WCSS plot
wcss_dat = unlist(lapply(X=clusterings, FUN=function(x){return(x$tot.withinss)}))

plot(Ks, wcss_dat, type="b", main="Elbow-plot", xlab="k", ylab="Total WCSS")

## Cluster plot - Color by cluster, point by label
k = 6
cols = paletteer_d("ggthemes::Tableau_20")
#wes_palette(n=k, "Darjeeling1")# wes_palette(n=5, "Darjeeling1") # See color palette
par(mfrow=c(2,3))
for (k in 1:6) {
  plot(cbind(X$x1, X$x2),
       col = cols[clusterings[[k]]$cluster],
       pch = pchs[X$label],
       xlab = expression(x[1]),
       ylab = expression(x[2]),
       main = paste("k =", k), asp=1)
  print(k)
}



(1/2)/(1-1/2)



  #create a hypothetical clustering outcome with 2 distinct clusters
g1 <- sample(1:2, size=10, replace=TRUE)
g2 <- sample(1:2, size=10, replace=TRUE)
rand.index(g1, g2)






 
######### Humble rebeginnings ######### 
set.seed(1)
# Basic setup ----
n_per_cluster <- 1000
k <- 3
N <- n_per_cluster * k

x <- rnorm(N); y <- rnorm(N)
g <- rep(1:k, each = n_per_cluster)

z1 <- c(3, 1, 3)
z2 <- c(10, 0, -10)
z3 <- c(0, 10, 0)

X <- cbind(
  x * z1[g] + z2[g],
  y * z1[g] + z3[g]
)

set.seed(1)
km <- kmeans(X, centers = k, nstart = 20)

# ---- Align labels (k=3: brute-force all permutations) ----
perms <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
acc <- sapply(perms, function(p) mean(p[km$cluster] == g))
best_p <- perms[[which.max(acc)]]

km_aligned <- best_p[km$cluster]
mis <- km_aligned != g

plot(X,
     pch = g,
     col = km_aligned,
     asp = 1,
     xlab = "x1", ylab = "x2",
     cex = 1)

# Optional: put a black ring around misclassified points
points(X[mis, 1], X[mis, 2], pch = 1, cex = 2.5, lwd = 2, col = "black")


# Experiment 1 ----
## Q: What is the relationship between misclassification and distance between  
## two clusters

# Reproduceability
set.seed(1)

e1_func <- function(sigma)
{
  gc()
  n <- 1000
  k <- 2
  N <- n * k
  exper_amount <- 1000
  # Base noise
  x_r <- rnorm(N)
  y_r <- rnorm(N)
  
  # True labels
  g <- rep(1:k, each = n)
  
  dist <- seq(0, 3, length.out=exper_amount)
  # Cluster parameters
  sx <- lapply(1:exper_amount, function(x){c(sigma,sigma)}) # x scale
  sy <- lapply(1:exper_amount, function(x){c(sigma,sigma)})     # y scale
  mx <- lapply(1:exper_amount, function(x){c(dist[x]/2, -dist[x]/2)}) # x shift
  my <- lapply(1:exper_amount, function(x){c(0, 0)})  # y shift
  
  # Construct dataset
  X <- lapply(1:exper_amount, function(x){
      cbind(
      x_r * sx[[x]][g] + mx[[x]][g],
      y_r * sy[[x]][g] + my[[x]][g]
      )
    }
  )
  
  
  ## cluster that dataK-means
  
  kms <- lapply(1:exper_amount, function(x){kmeans(X[[x]], centers = k, nstart = 40)})
  
  # Avoid label confusionism
  
  # Confusion matrix
  C <- lapply(1:exper_amount, function(x){as.matrix(table(g, kms[[x]]$cluster))})   # confusion counts (nonnegative)
  
  M <- lapply(1:exper_amount, function(x){max(C[[x]])})                            # largest count
  cost <- lapply(1:exper_amount, function(x){M[[x]] - C[[x]]})                          # nonnegative "cost" matrix
  
  assignment <- lapply(1:exper_amount, function(x){solve_LSAP(cost[[x]])})         # minimizes cost
  
  km_aligned <- lapply(1:exper_amount, function(i) {
    map_pred_to_true <- rep(NA_integer_, k)
    map_pred_to_true[assignment[[i]]] <- 1:k
    map_pred_to_true[kms[[i]]$cluster]
  })
  
  mis <- lapply(1:exper_amount, function(i) km_aligned[[i]] != g)
  mis_num <- unlist(lapply(1:exper_amount, function(x){sum(mis[[x]])}))
  mis_rate <- mis_num / N
  
  return(list(X = X,dist = dist ,g = g,kms = kms,km_aligned = km_aligned,mis = mis,mis_num = mis_num,mis_rate = mis_rate))
}
# Make plots
plot_range <- c(1, 1)

par(mfrow = c(diff(plot_range) + 1, 1), mar = c(4,4,2,1))
## Fisrst dist index w/ mis_num > 1
max(which(mis_num >= 1))
##
tail(which(dist <= 0.1), 1)
# Note all smaller dist has atleast 1 misclassification
all(mis_num[1:4651] >= 1)
#[1] TRUE

for (i in plot_range[1]:plot_range[2]) {
  plot(X[[i]],
       col = cols[km_aligned[[i]]+1],
       pch = 16+g^2,
       asp = 1,
       xlab = expression(x[1]),
       ylab = expression(x[2]),
       cex = ifelse(mis[[i]], 1.5, 1),
       main = sprintf("Ex1 Clustering - dist = %.3f  mis_rate = %.3f%%", dist[i], mis_rate[i] * 100))
  
  points(X[[i]][mis[[i]], 1], X[[i]][mis[[i]], 2],
         pch = 1, cex = 1.5, lwd = 2)
}
plot(X[[1]])

## Diagnostics
cap <- dist>=0
mis_num <- unlist(lapply(1:exper_amount, function(x){sum(mis[[x]])}))
mis_rate <- mis_num / N
plot(dist[cap],mis_rate[cap])
plot(mis_rate[cap]/abs(dist[cap]), dist[cap])
plot(abs(dist[cap])/mis_rate[cap], dist[cap])

abline(h=n_per_cluster)


# --- plotting only (uses your existing dist, mis_num, N) ---

mis_rate <- mis_num / N   # fraction misclassified
par(mfrow = c(1,1))
par(mar = c(5,4,4,2) + 0.1)   # default margins
plot(dist, mis_rate,
     type = "l", pch = 19,
     xlab = "Distance",
     ylab = "Misclassification Rate",
     main = "Miscalssification rate vs. distance",
     ylim = c(0, max(mis_rate, na.rm=TRUE)))

abline(h = 0.01, lty = 2)  # guideline: 5% misclass

plot(dist, mis_rate,
     type="l",
     log="y",
     xlab="separation",
     ylab="misclassification rate (log scale)")


bayes <- pnorm(-dist/2)
lines(dist, bayes, col="red")
legend("topright", legend=c("kmeans", "Bayes"),
       col=c("black","red"), lty=1)

slope <- diff(mis_rate) / diff(dist)
plot(dist[-1], slope, type="l")


grid()
plot(dist, mis_num,
     type = "b", pch = 19,
     xlab = "separation (dist)",
     ylab = "# misclassified",
     ylim = c(0, max(mis_num, na.rm=TRUE)))

abline(h = n, lty = 2)     
grid()


# Experiment 2 - variance and distance
# c(0.1, 0.5, 0.75, 1, 1.5, 2, 4, 6, 8, 10)
sigmas = c(.1,.5,1,2,3)
runs = list()
for (i in 1:length(sigmas))
{
  runs[[i]] = e1_func(sigmas[i])
}

# plot first line
plot(runs[[1]]$dist, runs[[1]]$mis_rate,
     col = cols[1],
     type = "l",
     lwd = 2,
     xlab = expression(paste("Dist (", Delta, ")")),
     ylab = "Misclassification rate",
     xlim = c(min(runs[[1]]$dist), max(runs[[1]]$dist) * 1.05),
     main = "Distance vs. misclassification rate", cex.lab=1.3, cex.main=1.5)

for(i in 2:length(sigmas)) {
  lines(runs[[i]]$dist, runs[[i]]$mis_rate, col = cols[i], lwd=2)
}

legend("topright",
       legend = paste0("σ = ", sigmas),
       col = cols,
       lty = 1,
       lwd = 2,
       bty = "n", cex=1.3)   


plot(runs[[1]]$dist/sigmas[1], runs[[1]]$mis_rate, type="l")
for(i in 2:length(sigmas)) lines(runs[[i]]$dist/sigmas[i], runs[[i]]$mis_rate)

plot(NULL, xlim=c(0,6), ylim=c(0,0.5),xlab=expression(Delta/sigma), ylab="Misclassification rate", main=("Misclassification rate as a constant function"), cex.lab=1.3, cex.main=1.5)

for(i in 1:length(sigmas)) {
  lines(runs[[i]]$dist/sigmas[i],
        runs[[i]]$mis_rate,
        col=cols[i], lwd=2)
}

legend("topright",
       legend = paste0("σ = ", sigmas),
       col = cols,
       lty = 1,
       lwd = 2,
       bty = "n", cex = 1.3)

# Experiment 2.1 - Investigating drops in mis_rate
# ============================================================
# Full runnable code (ONE sigma):
# - Same base data across all distances (only mean separation changes)
# - For each distance: run kmeans R times (nstart=1)
# - Bin solutions by ARI >= threshold
# - Store per-run: labels, WCSS (tot.withinss), mis_rate (aligned to true labels)
# - For each bin: keep representative, best-WCSS representative, and vectors of WCSS/mis
# - Also compute "optimal" curve at each distance = best WCSS among runs (proxy for nstart=200)
# ============================================================

# Packages
if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust")
if (!requireNamespace("clue", quietly = TRUE)) install.packages("clue")
library(mclust)
library(clue)

# ------------------------------------------------------------
# 1) Fixed data path across distances (same noise, only shift changes)
# ------------------------------------------------------------
e1_data_fixed <- function(sigma, seed = 1,
                          n = 1000,
                          exper_amount = 100,
                          dist_min = 0,
                          dist_max = 3) {
  set.seed(seed)
  k <- 2
  N <- n * k
  
  # fixed noise draw (same points across all distances)
  x_r <- rnorm(N)
  y_r <- rnorm(N)
  
  # true labels
  g <- rep(1:k, each = n)
  
  dist <- seq(dist_min, dist_max, length.out = exper_amount)
  
  # only the shift changes with dist[t]
  X <- lapply(seq_len(exper_amount), function(t) {
    mx <- c(dist[t] / 2, -dist[t] / 2)  # +/- shift in x
    cbind(
      x_r * sigma + mx[g],
      y_r * sigma
    )
  })
  
  list(X = X, dist = dist, g = g, sigma = sigma, seed = seed)
}

# ------------------------------------------------------------
# 2) Misclassification with optimal label alignment (k=2)
# ------------------------------------------------------------
mis_rate_aligned_2 <- function(g, pred) {
  # confusion counts
  C <- as.matrix(table(g, pred))
  M <- max(C)
  cost <- M - C
  assignment <- solve_LSAP(cost)  # maps TRUE rows -> PRED cols
  
  # assignment is a permutation of predicted labels:
  # predicted label = assignment[i] corresponds to true label i
  map_pred_to_true <- rep(NA_integer_, 2)
  map_pred_to_true[assignment] <- 1:2
  
  pred_aligned <- map_pred_to_true[pred]
  mean(pred_aligned != g)
}

# ------------------------------------------------------------
# 3) Basin/bin experiment:
#    For each distance t:
#      - run kmeans R times (nstart=1)
#      - compute ARI to bin reps
#      - store WCSS + mis_rate per run
#      - each bin stores vectors (wcss, mis_rate) + best representative by WCSS
#    Also returns "optimal" curve across all runs = min WCSS per distance
# ------------------------------------------------------------
bin_kmeans_path_with_metrics <- function(dat,
                                         R = 100,
                                         ari_thresh = 0.99,
                                         iter.max = 50,
                                         algorithm = "Hartigan-Wong") {
  X_list <- dat$X
  dist <- dat$dist
  g <- dat$g
  exper_amount <- length(dist)
  
  bins_per_dist <- vector("list", exper_amount)
  
  # these summarize the "optimal" run (min WCSS among the R runs) at each distance
  best_wcss <- numeric(exper_amount)
  best_mis  <- numeric(exper_amount)
  best_bin  <- integer(exper_amount)
  
  for (t in seq_len(exper_amount)) {
    
    # each bin keeps:
    # - representative clustering (first seen)
    # - counts
    # - vectors of wcss and mis_rate from runs assigned to this bin
    # - best-wcss representative (for "optimal" plotting / inspecting)
    reps_clusters <- list()
    reps_X <- list()
    bin_counts <- integer(0)
    
    bin_wcss <- list()
    bin_mis  <- list()
    
    # store best solution within each bin
    bin_best_wcss <- numeric(0)
    bin_best_run_cluster <- list()  # cluster labels for best wcss within bin
    bin_best_run_mis <- numeric(0)
    
    # global best among all bins
    global_best_wcss <- Inf
    global_best_mis  <- NA_real_
    global_best_bin  <- NA_integer_
    
    for (r in seq_len(R)) {
      # different initialization each run, same data
      set.seed(1000 + 10000 * t + r)
      
      km <- kmeans(X_list[[t]],
                   centers = 2,
                   nstart = 1,
                   iter.max = iter.max,
                   algorithm = algorithm)
      
      cl <- km$cluster
      wcss <- km$tot.withinss
      mis  <- mis_rate_aligned_2(g, cl)
      
      # assign run to a bin by ARI threshold
      if (length(reps_clusters) == 0) {
        reps_clusters[[1]] <- cl
        reps_X[[1]] <- X_list[[t]]
        bin_counts[1] <- 1
        bin_wcss[[1]] <- wcss
        bin_mis[[1]]  <- mis
        
        bin_best_wcss[1] <- wcss
        bin_best_run_cluster[[1]] <- cl
        bin_best_run_mis[1] <- mis
        
        j <- 1
      } else {
        ari_vals <- vapply(
          reps_clusters,
          function(repcl) mclust::adjustedRandIndex(repcl, cl),
          numeric(1)
        )
        j <- which.max(ari_vals)
        
        if (ari_vals[j] >= ari_thresh) {
          bin_counts[j] <- bin_counts[j] + 1
          bin_wcss[[j]] <- c(bin_wcss[[j]], wcss)
          bin_mis[[j]]  <- c(bin_mis[[j]], mis)
          
          # update best-within-bin if improved
          if (wcss < bin_best_wcss[j]) {
            bin_best_wcss[j] <- wcss
            bin_best_run_cluster[[j]] <- cl
            bin_best_run_mis[j] <- mis
          }
        } else {
          reps_clusters[[length(reps_clusters) + 1]] <- cl
          reps_X[[length(reps_X) + 1]] <- X_list[[t]]
          bin_counts[length(bin_counts) + 1] <- 1
          
          newj <- length(bin_counts)
          bin_wcss[[newj]] <- wcss
          bin_mis[[newj]]  <- mis
          
          bin_best_wcss[newj] <- wcss
          bin_best_run_cluster[[newj]] <- cl
          bin_best_run_mis[newj] <- mis
          
          j <- newj
        }
      }
      
      # track global best among all runs at this distance
      if (wcss < global_best_wcss) {
        global_best_wcss <- wcss
        global_best_mis  <- mis
        global_best_bin  <- j
      }
    }
    
    best_wcss[t] <- global_best_wcss
    best_mis[t]  <- global_best_mis
    best_bin[t]  <- global_best_bin
    
    bins_per_dist[[t]] <- list(
      dist = dist[t],
      
      # bins
      n_bins = length(reps_clusters),
      counts = bin_counts,
      
      # representative (first seen) per bin (handy for plotting)
      reps_clusters = reps_clusters,
      reps_X = reps_X,
      
      # metrics per bin across assigned runs
      bin_wcss = bin_wcss,
      bin_mis  = bin_mis,
      
      # best-within-bin (proxy "what would be chosen if you searched hard but stayed in that basin")
      bin_best_wcss = bin_best_wcss,
      bin_best_cluster = bin_best_run_cluster,
      bin_best_mis = bin_best_run_mis
    )
  }
  
  list(
    sigma = dat$sigma,
    seed = dat$seed,
    g = dat$g,
    dist = dist,
    R = R,
    ari_thresh = ari_thresh,
    bins_per_dist = bins_per_dist,
    
    # global-best per distance (proxy for nstart very large)
    best_wcss = best_wcss,
    best_mis = best_mis,
    best_bin = best_bin
  )
}

# ------------------------------------------------------------
# 4) Plot helpers
# ------------------------------------------------------------
plot_mis_best <- function(res) {
  plot(res$dist, res$best_mis,
       type = "l", lwd = 2,
       xlab = expression(Delta),
       ylab = "Misclassification rate",
       main = paste("Mis-rate of best-WCSS solution vs distance (sigma =", res$sigma, ")"))
}

plot_bins_and_domfrac <- function(res) {
  n_bins <- sapply(res$bins_per_dist, function(z) z$n_bins)
  dom_frac <- sapply(res$bins_per_dist, function(z) max(z$counts) / sum(z$counts))
  
  par(mfrow = c(1,2))
  plot(res$dist, n_bins, type="l", lwd=2,
       xlab=expression(Delta), ylab="# bins",
       main="Number of ARI bins vs distance")
  plot(res$dist, dom_frac, type="l", lwd=2, ylim=c(0,1),
       xlab=expression(Delta), ylab="Dominant bin fraction",
       main="Dominant bin frequency vs distance")
  par(mfrow = c(1,1))
}

# Plot representative configuration for a chosen distance+bin
plot_rep_bin <- function(res, t, bin = NULL, use_best = FALSE) {
  info <- res$bins_per_dist[[t]]
  if (is.null(bin)) bin <- which.max(info$counts)
  
  X <- info$reps_X[[bin]]
  if (!use_best) {
    cl <- info$reps_clusters[[bin]]
    ttl <- paste0("Rep (first) at Δ=", round(info$dist,3),
                  " | bin ", bin, " (count=", info$counts[bin], ")")
  } else {
    cl <- info$bin_best_cluster[[bin]]
    ttl <- paste0("Rep (best WCSS) at Δ=", round(info$dist,3),
                  " | bin ", bin,
                  " (count=", info$counts[bin],
                  ", best_mis=", round(info$bin_best_mis[bin],3),
                  ", best_wcss=", round(info$bin_best_wcss[bin],1), ")")
  }
  plot(X, col = cl, pch = 16, xlab="x", ylab="y", main = ttl)
}

# ------------------------------------------------------------
# RUN IT
# ------------------------------------------------------------
sigma <- 1
dat <- e1_data_fixed(sigma = sigma, seed = 12, n = 1000, exper_amount = 2000)

res <- bin_kmeans_path_with_metrics(dat, R = 400, ari_thresh = 0.99)

# 1) Misclassification rate for the "optimal" (min WCSS among R) solution at each distance
plot_mis_best(res)

# 2) Bin diagnostics (to relate to your drop)
plot_bins_and_domfrac(res)
find_biggest_drop <- function(dist, y) {
  dy <- y[-length(y)] - y[-1]
  tstar <- which.max(dy)
  list(tstar=tstar, drop=dy[tstar],
       Delta_before=dist[tstar], Delta_after=dist[tstar+1],
       y_before=y[tstar], y_after=y[tstar+1])
}

dropinfo <- find_biggest_drop(res$dist, res$best_mis)
dropinfo
tstar <- dropinfo$tstar

t_before <- tstar      # "before" distance index
t_after  <- tstar + 1  # "after" distance index

find_biggest_drop <- function(dist, y) {
  dy <- y[-length(y)] - y[-1]
  tstar <- which.max(dy)
  list(tstar=tstar, drop=dy[tstar],
       Delta_before=dist[tstar], Delta_after=dist[tstar+1],
       y_before=y[tstar], y_after=y[tstar+1])
}

dropinfo <- find_biggest_drop(res$dist, res$best_mis)
dropinfo
tstar <- dropinfo$tstar

t_before <- tstar      # "before" distance index
t_after  <- tstar + 1  # "after" distance index

plot_bin_best_wcss <- function(res, t, main_prefix="") {
  info <- res$bins_per_dist[[t]]
  barplot(info$bin_best_wcss,
          names.arg = seq_along(info$bin_best_wcss),
          xlab = "Bin id",
          ylab = "Best WCSS in bin",
          main = paste0(main_prefix, " Δ=", round(info$dist,3),
                        " | global best WCSS=", round(res$best_wcss[t],1)))
}

par(mfrow=c(1,2))
plot_bin_best_wcss(res, t_before, "BEFORE:")
plot_bin_best_wcss(res, t_after,  "AFTER :")
par(mfrow=c(1,1))

plot_best_2x3 <- function(res, t_center) {
  idx <- (t_center-2):(t_center+3)          # 3 before (t-2,t-1,t) and 3 after (t+1,t+2,t+3)
  idx <- idx[idx >= 1 & idx <= length(res$dist)]
  
  par(mfrow=c(2,3), mar=c(3,3,3,1))
  for (t in idx) {
    info <- res$bins_per_dist[[t]]
    b <- res$best_bin[t]
    X <- info$reps_X[[1]]                   # same X for all bins at fixed t
    cl <- info$bin_best_cluster[[b]]
    
    plot(X, col=cl, pch=16, cex=0.6, xlab="x", ylab="y",
         main=paste0("Δ=", round(info$dist,3),
                     "\nbest_bin=", b,
                     " best_mis=", round(res$best_mis[t],3)))
  }
  par(mfrow=c(1,1))
}

plot_best_2x3(res, tstar)

plot_all_bins <- function(res, t, use_best=TRUE, ncol=3) {
  info <- res$bins_per_dist[[t]]
  K <- info$n_bins
  nrow <- ceiling(K / ncol)
  
  par(mfrow=c(nrow, ncol), mar=c(3,3,3,1))
  
  for (j in 1:K) {
    X <- info$reps_X[[1]]
    cl <- if (use_best) info$bin_best_cluster[[j]] else info$reps_clusters[[j]]
    
    ttl <- paste0("Δ=", round(info$dist,3),
                  " | bin ", j,
                  "\ncount=", info$counts[j],
                  " best_mis=", round(info$bin_best_mis[j],3),
                  " best_wcss=", round(info$bin_best_wcss[j],1))
    
    plot(X, col=cl, pch=16, cex=0.55, xlab="x", ylab="y", main=ttl)
  }
  
  # fill leftover panels
  for (j in (K+1):(nrow*ncol)) {
    if (j > nrow*ncol) break
    plot.new()
  }
  
  par(mfrow=c(1,1))
}
png("bins_before.png", width=1600, height=1600)
plot_all_bins(res, t_before, use_best=TRUE)
dev.off()

png("bins_after.png", width=1600, height=1600)
plot_all_bins(res, t_after, use_best=TRUE)
dev.off()# ------------------------------------------------------------
# Inspect bins around a suspected drop
# ------------------------------------------------------------
# Example: pick a distance index near where your mis curve drops
t <- 40
res$dist[t]
res$bins_per_dist[[t]]$n_bins
res$bins_per_dist[[t]]$counts
res$best_bin[t]          # which bin produced global best WCSS at this distance?
res$best_mis[t]          # mis-rate of global best
res$best_wcss[t]         # wcss of global best

# Plot dominant bin representative (first seen)
plot_rep_bin(res, t)

# Plot best-WCSS representative in the bin that actually wins (proxy "nstart huge")
plot_rep_bin(res, t, bin = res$best_bin[t], use_best = TRUE)


## STORY

# ============================================================
# Add/replace these plotting helpers
# - plot_bin_counts(): barplot of bin frequencies (counts), with best-bin and best-WCSS in title
# - plot_all_bins_safe(): plots ALL reps at a distance with count + best_mis + best_wcss
#   and ALSO marks which bin is chosen globally (best WCSS) in the title
# ============================================================

plot_bin_counts <- function(res, t, main_prefix="") {
  info <- res$bins_per_dist[[t]]
  bb <- res$best_bin[t]
  
  barplot(info$counts,
          names.arg = seq_along(info$counts),
          xlab = "Bin id",
          ylab = "Frequency (count)",
          main = paste0(main_prefix,
                        " Δ=", round(info$dist,3),
                        " | n_bins=", info$n_bins,
                        " | best_bin=", bb,
                        " (best_wcss=", round(info$bin_best_wcss[bb], 1),
                        ", best_mis=", round(info$bin_best_mis[bb], 3), ")"))
}

plot_all_bins_safe <- function(res, t, use_best=TRUE, ncol=NULL) {
  info <- res$bins_per_dist[[t]]
  K <- info$n_bins
  bb <- res$best_bin[t]
  
  if (is.null(ncol)) {
    ncol <- if (K <= 3) K else if (K <= 6) 3 else 4
  }
  nrow <- ceiling(K / ncol)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  par(mfrow = c(nrow, ncol),
      mar = c(1.4, 1.4, 2.4, 0.6),
      oma = c(0, 0, 0, 0))
  
  X <- info$reps_X[[1]]  # same X for all bins at fixed t
  
  for (j in 1:K) {
    cl <- if (use_best) info$bin_best_cluster[[j]] else info$reps_clusters[[j]]
    
    chosen_tag <- if (j == bb) "  <-- chosen (min WCSS)" else ""
    ttl <- paste0("Δ=", round(info$dist, 3),
                  " bin ", j, chosen_tag,
                  "\ncount=", info$counts[j],
                  " mis=", round(info$bin_best_mis[j], 3),
                  " wcss=", round(info$bin_best_wcss[j], 1))
    
    plot(X, col=cl, pch=16, cex=0.5, axes=FALSE, xlab="", ylab="", main=ttl)
    box()
  }
}


# Find drop index (you already have this)
tstar <- which.max(res$best_mis[-length(res$best_mis)] - res$best_mis[-1])
t_before <- tstar
t_after  <- tstar + 1

# 1) Counts before/after + best-bin info in title
png("story_counts_before_after.png", width=1800, height=900)
par(mfrow=c(1,2), mar=c(4,4,4,1))
plot_bin_counts(res, t_before, "BEFORE:")
plot_bin_counts(res, t_after,  "AFTER :")
par(mfrow=c(1,1))
dev.off()

# 2) All representatives BEFORE (with WCSS + which one is chosen)
png("bins_before.png", width=1800, height=1200)
plot_all_bins_safe(res, t_before, use_best=TRUE)
dev.off()

# 3) All representatives AFTER
png("bins_after.png", width=1800, height=1200)
plot_all_bins_safe(res, t_after, use_best=TRUE)
dev.off()

# Experiment 3 - equilateral triangles
e3_func <- function(sigma){
  gc()
  n <- 300
  k_true <- 3
  k <- 2
  N <- n * k_true
  exper_amount <- 200
  # Base noise
  x_r <- rnorm(N)
  y_r <- rnorm(N)
  
  # True labels
  g <- rep(1:k_true, each = n)
  
  #dist <- seq(0, 2000, length.out=exper_amount)
  side_length <- 20
  pull <- seq(0,5, length.out = exper_amount)
  # Cluster parameters
  sx <- lapply(1:exper_amount, function(x){c(sigma,sigma, sigma)}) # x scale
  sy <- lapply(1:exper_amount, function(x){c(sigma,sigma, sigma)})     # y scale
  mx <- lapply(1:exper_amount, function(x){c(side_length/2, -side_length/2, 0)}) # x shift
  my <- lapply(1:exper_amount, function(x){c(-pull[x], 0, (side_length * sqrt(3))/2)})  # y shift
  
  # Construct dataset
  X <- lapply(1:exper_amount, function(x){
    cbind(
      x_r * sx[[x]][g] + mx[[x]][g],
      y_r * sy[[x]][g] + my[[x]][g]
    )
  }
  )
  ## cluster that dataK-means
  kms <- lapply(1:exper_amount, function(x){kmeans(X[[x]], centers = k, nstart = 20)})
 
  ari <- unlist(lapply(1:exper_amount, function(t)
    mclust::adjustedRandIndex(g, kms[[t]]$cluster)
  ))
  
  return(list(
    X = X, side_length = side_length, pull = pull, g = g,
    kms = kms,
    clusters = lapply(1:exper_amount, function(t) kms[[t]]$cluster),
    ari = ari
  ))
}

# Count cluster bins
bin_clusterings_by_ari <- function(sigma, R = 50, efunc) {
  runs <- vector("list", R)
  for (r in 1:R) runs[[r]] <- efunc(sigma)
  
  pull <- runs[[1]]$pull
  exper_amount <- length(pull)
  
  bins_per_pull <- vector("list", exper_amount)
  
  for (t in 1:exper_amount) {
    
    reps_clusters <- list()   # representative cluster label vectors
    reps_X <- list()          # corresponding data matrix
    bin_counts <- integer(0)
    
    for (r in 1:R) {
      
      cl <- runs[[r]]$clusters[[t]]   # predicted labels
      X_t <- runs[[r]]$X[[t]]         # dataset
      
      if (length(reps_clusters) == 0) {
        reps_clusters[[1]] <- cl
        reps_X[[1]] <- X_t
        bin_counts[1] <- 1
      } else {
        
        ari_to_reps <- sapply(
          reps_clusters,
          function(repcl) mclust::adjustedRandIndex(repcl, cl)
        )
        
        j <- which(ari_to_reps == 1)
        
        if (length(j) > 0) {
          bin_counts[j[1]] <- bin_counts[j[1]] + 1
        } else {
          reps_clusters[[length(reps_clusters) + 1]] <- cl
          reps_X[[length(reps_X) + 1]] <- X_t
          bin_counts[length(bin_counts) + 1] <- 1
        }
      }
    }
    
    bins_per_pull[[t]] <- list(
      pull = pull[t],
      n_bins = length(reps_clusters),
      counts = bin_counts,
      reps_clusters = reps_clusters,
      reps_X = reps_X
    )
  }
  
  list(
    pull = pull,
    bins_per_pull = bins_per_pull
  )
}

plot_reps_for_pull <- function(func_dat, t, cols, pch = 16) {
  bp <- func_dat$bins_per_pull[[t]]
  reps_X <- bp$reps_X
  reps_cl <- bp$reps_clusters
  counts <- bp$counts
  pull_val <- bp$pull
  
  # grid layout (roughly square)
  ncol <- ceiling(sqrt(bp$n_bins))
  nrow <- ceiling(bp$n_bins / ncol)
  
  op <- par(mfrow = c(ceiling(bp$n_bins / ceiling(sqrt(bp$n_bins))), ceiling(sqrt(bp$n_bins))), mar = c(3.5, 3.5, 3, 1))
  on.exit(par(op), add = TRUE)
  
  for (i in 1:bp$n_bins) {
    X <- reps_X[[i]]
    cl <- reps_cl[[i]]
    plot(X,
         col = cols[cl], pch = pch, asp = 1,
         xlab = expression(x[1]),ylab = expression(x[2]),
         cex.main = 1.5, cex.lab = 1.3,
         main = sprintf("pull = %.3f   freq = %d", pull_val, counts[b]))
  }
  
}


out <- bin_clusterings_by_ari(sigma = 1, R = 400)
out$summary_df

plot_reps_for_pull(out, t = 130, cols=c(cols[1], cols[3]))


bin_count_list <- sapply(1:length(out$pull), function(p){out$bins_per_pull[[p]]$n_bins})
plot( out$pull, bin_count_list, pch = 17, col=cols[1], 
      main="Cluster bins vs. pull length",
      cex.lab=1.3, cex.main = 1.5,
      xlab="Pull", ylab = "# Different clusterings")
pull_trans_1 <- (out$pull[which(bin_count_list < 3)[1]]+out$pull[which(bin_count_list >= 3)[length(which(bin_count_list >= 3))]])/2
pull_trans_2 <- (out$pull[which(bin_count_list < 2)[1]]+out$pull[which(bin_count_list >= 2)[length(which(bin_count_list >= 2))]])/2

abline(v=pull_trans_1, lty = 3, col=cols[3], lwd=2)
abline(v=pull_trans_2, lty = 3, col=cols[5], lwd=2)
legend("topright",
       legend = c(paste0("pull = ", sprintf("%.4f", pull_trans_1)),paste0("pull = ", sprintf("%.4f", pull_trans_2))),
       col = c(cols[3], cols[4]),
       lty = 3,
       lwd = 2, cex=1)   


# TODO add some sort of identifier for clusterings such that I can get 
# probability of being in each cluster at a given pull length


# Experiment 4 - Square
e4_func <- function(sigma){
  gc()
  n <- 300
  k_true <- 4
  k <- 2
  N <- n * k_true
  exper_amount <- 200
  
  # Base noise (fixed across pull values)
  x_r <- rnorm(N)
  y_r <- rnorm(N)
  
  # True labels
  g <- rep(1:k_true, each = n)
  
  side_length <- 20
  pull <- seq(0, 5, length.out = exper_amount)
  

  sx <- lapply(1:exper_amount, function(x) rep(sigma, k_true))
  sy <- lapply(1:exper_amount, function(x) rep(sigma, k_true))
  mx <- lapply(1:exper_amount, function(x) c(-side_length/2,  side_length/2,  side_length/2, -side_length/2))
  my <- lapply(1:exper_amount, function(x) c(-side_length/2, -side_length/2 - pull[x],  side_length/2,  side_length/2))
  
  # Construct datasets
  X <- lapply(1:exper_amount, function(x){
    cbind(
      x_r * sx[[x]][g] + mx[[x]][g],
      y_r * sy[[x]][g] + my[[x]][g]
    )
  })
  
  # kmeans
  kms <- lapply(1:exper_amount, function(x) kmeans(X[[x]], centers = k, nstart = 20))
  
  ari <- unlist(lapply(1:exper_amount, function(t)
    mclust::adjustedRandIndex(g, kms[[t]]$cluster)
  ))
  
  list(
    X = X,
    side_length = side_length,
    pull = pull,
    g = g,
    kms = kms,
    clusters = lapply(1:exper_amount, function(t) kms[[t]]$cluster),
    ari = ari
  )
}


out4 <- bin_clusterings_by_ari(sigma = 1, R = 400, efunc= e4_func)
plot_reps_for_pull(out4, t = 19, cols = cols)



bin_count_list <- sapply(1:length(out4$pull), function(p){out4$bins_per_pull[[p]]$n_bins})
plot( out4$pull, bin_count_list, pch = 17, col=cols[1], 
      main="Cluster bins vs. pull length",
      cex.lab=1.3, cex.main = 1.5,
      xlab="Pull", ylab = "# Different clusterings")
pull_trans_1 <- (out4$pull[which(bin_count_list < 2)[1]]+out4$pull[which(bin_count_list >= 2)[length(which(bin_count_list >= 2))]])/2

abline(v=pull_trans_1, lty = 3, col=cols[3], lwd=2)
legend("topright",
       legend = c(paste0("pull = ", sprintf("%.4f", pull_trans_1))),
       col = c(cols[4]),
       lty = 3,
       lwd = 2, cex=1)   



# Experiment 5 - Data pertubation two clusters
stability_under_added_gauss_noise <- function(X0, k = 2, sig = seq(0, 1, length.out = 50),
                                     R = 200, nstart = 20) {

  cl0 <- kmeans(X0, centers = k, nstart = nstart)$cluster
  
  stab <- numeric(length(sig))  
  
  for (i in 1:length(sig)) {
    s <- sig[i]
    ari_vals <- numeric(R)
    
    for (r in 1:R) {
      Xp <- X0 + matrix(rnorm(length(X0), sd = s), ncol = ncol(X0))
      clp <- kmeans(Xp, centers = k, nstart = nstart)$cluster
      ari_vals[r] <- mclust::adjustedRandIndex(cl0, clp)
    }
    
    stab[i] <- mean(ari_vals)
  }
  
  list(sig = sig, stability = stab)
}

n <- 300
k_true <- 1
k <- 2
N = n * k_true
# Base noise
x_r <- rnorm(N)
y_r <- rnorm(N)


dist <- 5
# Cluster parameters
sigma = 1
sx <- c(sigma,sigma) # x scale
sy <- c(sigma,sigma)    # y scale
mx <- c(dist/2, -dist/2)# x shift
my <- c(0,0)  # y shift

# Construct dataset
X <- cbind(
    x_r * sx + mx,
    y_r * sy + my
  )

st <- stability_under_added_gauss_noise(X, k = 2, sig = seq(0, 10, length.out = 400), R = 200)

plot(dist(st$sig), st$stability, type="l", lwd=2,
     xlab=expression(sigma[added]), ylab="Mean ARI vs Noise",
     main="Stability plot")



# Ring clustering
sigma <- 0.2
n <- 3000

# mean of each cluster + noise
r1 <- 2 + rnorm(n, 0, sigma)
r2 <- 5 + rnorm(n, 0, sigma)

a1 <- runif(n,0, 2*pi)
x_inner <- r1 *cos(a1)
y_inner <- r1 * sin(a1)

a2 <- runif(n,0,2*pi)
x_outer <- r2 *cos(a2)
y_outer <- r2 *sin(a2)

X <- rbind(
  cbind(x_inner+10, y_inner+10),
  cbind(x_outer+10, y_outer+10)
)




X_alt <- sqrt((X[,1]-mean(X[,1]))^2 + (X[,2]-mean(X[,2]))^2)
#X_norm <- X / sqrt(rowSums(X^2))
#X_norm <- X - colMeans(X)
mu <- colMeans(X)
X_norm <- X - matrix(mu, nrow(X), 2, byrow = TRUE)



km_rings <- kmeans(X, centers=2)
km_rings_alt <- kmeans(X_alt, centers=2)
km_rings_norm <- kmeans(X_norm, centers=2)
g <- c(rep(1, n), rep(2, n))


par(mfrow=c(2,2))
plot(X,
     col = cols[km_rings$cluster+1],
     pch = g+16,
     asp = 1,
     main = "Rings", axes=FALSE)
box()

plot(X_alt,
     col = cols[km_rings_alt$cluster+1],
     pch = g+16,
     main = "Rings alt", axes=FALSE)
box()

plot(X,
     col = cols[km_rings_alt$cluster+1],
     pch = g+16,
     main = "Rings alt", asp=1, axes=FALSE)
box()



# Half circle fun

sigma <- 0.1


dist <- 8
a1 <- -runif(n,0, pi)
x_inner <- r1 *cos(a1)+ r2 + rnorm(n, 0, sigma)
y_inner <- r1 * sin(a1) + rnorm(n, 0, sigma) + r2- dist

a2 <- runif(n,0,pi)
x_outer <- r2 *cos(a2) + rnorm(n, 0, sigma)
y_outer <- r2 *sin(a2) +  rnorm(n, 0, sigma)

X <- rbind(
  cbind(x_inner, y_inner),
  cbind(x_outer, y_outer)
)
X_alt <- sqrt(X[,1]^2 + X[,2]^2)

g <- c(rep(1, n), rep(2, n))




km_rings <- kmeans(X, centers=3)
km_rings_alt <- kmeans(X_alt, centers=2)

plot(X,
     col = cols[km_rings_alt$cluster+1],
     pch = g+16,
     main = "Rings alt")
plot(X,
     col = cols[km_rings$cluster+1],
     pch = g+16,
     asp = 1,
     main = "Half rings")






### Recreation
n <- 3000
ks <- 2:7
x_r <- rnorm(n)
y_r <- rnorm(n)
mx <- 0
my <- 0
sx <- 1
sy <- 1

X <- cbind(x_r*sx + mx, y_r*sy + my)


kms <- lapply(ks, function(k){ kmeans(X, center=k, nstart = 50) })

par(mfrow=c(3,2))
for(k in ks)
  {
    plot(X, pch=19,col=cols[kms[[k+1 - min(ks)]]$cluster], main=paste("k = ",k))
  }

kmeans(X, center=k, nstart=50)
plot(X, col=cols[km$cluster])




## Silluette score and CH

set.seed(1)
X <- scale(iris[, 1:4])   # IMPORTANT: k-means + these indices assume sensible scaling

k <- 3
km <- kmeans(X, centers = k, nstart = 50)

km$tot.withinss   # W (within)
km$betweenss      # B (between)
library(cluster)

# For k-means: silhouette is typically computed using Euclidean distances
d <- dist(X, method = "euclidean")

sil <- silhouette(km$cluster, d)

# Average silhouette (overall score)
mean(sil[, "sil_width"])
plot(sil, main = paste("Silhouette plot, k =", k), border = NA)
abline(v = mean(sil[, "sil_width"]), lty = 2)
aggregate(sil[, "sil_width"], by = list(cluster = sil[, "cluster"]), FUN = mean)
install.packages("fpc")   # once
library(fpc)

ch <- calinhara(X, km$cluster)
ch
n <- nrow(X)
k <- length(unique(km$cluster))
W <- km$tot.withinss
B <- km$betweenss

CH <- (B / (k - 1)) / (W / (n - k))
CH


library(cluster)
library(fpc)

set.seed(1)
ks <- 2:149

sil_scores <- numeric(length(ks))
ch_scores  <- numeric(length(ks))
wcss <- numeric(length(ks))

for (i in seq_along(ks)) {
  k <- ks[i]
  km <- kmeans(X, centers = k, nstart = 200)
  
  # Elbow method
  wcss[i] <- km$tot.withinss
  
  # Silhouette
  sil_scores[i] <- mean(silhouette(km$cluster, dist(X))[, "sil_width"])
  
  # CH
  ch_scores[i] <- calinhara(X, km$cluster)
}

sil_scores
ch_scores
par(mfrow=c(1,3))
plot(wcss, type="b")
plot(sil_scores,type="b")
plot(ch_scores, type="b")

plot(-diff(wcss)/max(diff(wcss)))
abline(h=0.1)

t <- seq(0, 30, 0.001)
plot(cos(10*t), , ylim=c(-1,1), asp=1)




# Experiment 6 - Rim data replication 
library(gtools)

set.seed(510) # seed 510 gives desired clustering
# Generate data
n <- 10000
X <- cbind(rnorm(n), rnorm(n))

# kmeans
k <- 5
km_k  <- kmeans(X, centers = k,   nstart = 50)
km_k1 <- kmeans(X, centers = k+1, nstart = 50)

# Plot
par(mfrow=c(1,2))
plot(X, col=cols[km_k$cluster], axes=FALSE, main="k=5", pch=16, cex=0.4)
box()
plot(X, col=cols[km_k1$cluster], axes=FALSE, main="k=6", pch=16, cex=0.4)
box()

 # Make a matrix of all paired center distances 
dist_euclid <- function(a, b) sqrt(sum((a-b)^2))

dist_mat <- matrix(NA_real_, nrow = k, ncol = k+1)

for (i in 1:k) { # Loop over old centers
  for (j in 1:(k+1)) { # loop over new centers
    dist_mat[i, j] <- dist_euclid(km_k$centers[i,], km_k1$centers[j,]) # Entrance i,j is the dist between old center and new center 
  }
}

# WE will align the clusters and identify the "new cluster
# Use gtools permutations to get all possible permutations of pairings
# Note we leave one out which is the new cluster so if perms= 2 3 4 1 5
# then it means that:
# old 1 -> new 2, old 2 -> new 3, old 3 -> new 4, old 4 -> new 1, old 5 -> new 5 
perms <- permutations(k, k)

best_sum_distances <- numeric(k + 1) # Vector to store the best sum of distances between cluster centers
best_perm_for_i <- matrix(NA_integer_, nrow = k + 1, ncol = k)  # store best permutation for each i

for (i in 1:(k + 1)) {
  
  # build the k x k distance mat w/o col i
  dist_mat_wo_i <- dist_mat[, -i]
  
  # Find best permutation leaving out i
  best_sum <- Inf
  best_perm <- rep(NA_integer_, k)
  
  # loop over each permuatation  
  for (perm in 1:nrow(perms)) {
    total <- 0
    # loop over each element in old clustering
    for (j in 1:k) {
      # We add the distance from old center to the new center as per
      # current permutation (perm)
      total <- total + dist_mat_wo_i[j, perms[perm, j]]
    }
    
    # Update best to reflect 
    if (total < best_sum) {
      best_sum <- total
      best_perm <- perms[perm, ]
    }
  }
  
  best_sum_distances[i] <- best_sum
  best_perm_for_i[i, ] <- best_perm
}

best_removed_i <- which.min(best_sum_distances)

pair_old_to_new <- rep(NA_integer_, k)
perm <- best_perm_for_i[best_removed_i, ]# 

for (i in 1:k) { # pair the old clusters to the 5 new w/o i. 
  # note that the permutation is an "indexer", so if perm[1]=3 then
  # it means that the old 1 -> new 3, but new 3 might corrospond to something
  # different as we removed i. e.g. if i=2 then 3 corresponds to 4
  pair_old_to_new[i] <- (1:(k+1))[-best_removed_i][perm[i]]
}

for (r in 1:k) {cat("old cluster", r, "->new cluster", pair_old_to_new[r], "\n")}

# Plotting
new_cols<- rep(cols[15], k + 1)

# Permutate cols such that it aligns with pair_old_to_new
# s.t. if we have pairing old 1 -> new 3 then new_cols[3]=cols[1]
for (i in 1:k) {
  new_cols[pair_old_to_new[i]] <- cols[i]
}

par(mfrow = c(1, 2), mar = c(2, 2, 3, 1))
plot(X, col = cols[km_k$cluster], pch = 16, cex = 0.5,
     axes = FALSE, main = paste0("k = ", k, " (old)"))
box()
plot(X, col = new_cols[km_k1$cluster], pch = 16, cex = 0.5,
     axes = FALSE, main = paste0("k = ", k + 1, " (new)"))
box()

# Relabeling instead....
new_to_old_map <- numeric(1+k) + k+1

for (i in 1:k) {
  new_to_old_map[pair_old_to_new[i]] <- i
}
km_k1_aligned <- new_to_old_map[km_k1$cluster]

# So we only count rim data if points switches between assigned clusters which is 
# NOT the "new cluster"
is_rim <- (km_k1_aligned != km_k$cluster) & (km_k1_aligned != k + 1)

cat("rim rate:", sum(is_rim)/length(is_rim), "\n")

# Plot w/ new cluster highlighted
par(mfrow = c(1, 1))
plot(X, col = ifelse(km_k1_aligned == 6, 1, 2), pch = 16, cex = 0.4,
     main = sprintf("k = 6", k+1), axes=FALSE)
box()

# rim points highlighted
par(mfrow = c(1, 2))
plot(X, col = cols[km_k$cluster], pch = 16, cex = ifelse(is_rim, 1, 0.35),
     main = "k=5",xlab="",ylab="", axes=FALSE, asp=1, cex.main=1.5)
box()
plot(X, col = cols[km_k1_aligned] , pch = 16, cex = ifelse(is_rim, 1, 0.35),
     main = "k=6",xlab="",ylab="",axes=FALSE, asp=1,cex.main=1.5)
box()



# Experiment 7 - 3D rim data

library(gtools)

set.seed(510) # seed 510 gives desired clustering
par(mfrow=c(3,3))

# Generate data
n <- 10000
d <- 16
X <- matrix(rnorm(n * d), nrow = n, ncol = d)

# kmeans
ks <- 1:10

kms <- list()
for(k in ks)
{
  kms[[k]] <- kmeans(X, centers=k, nstart=50)
}
# Make a matrix of all paired center distances 
dist_euclid <- function(a, b) sqrt(sum((a-b)^2))

dist_mats <- list()
km_k1_aligned <- list()
is_rim <- list()
for(k in 1:(max(ks)-1))
{
  dist_mats[[k]] <- matrix(NA_real_, nrow = k, ncol = k+1)

  for (i in 1:k) { # Loop over old centers
    for (j in 1:(k+1)) { # loop over new centers
      dist_mats[[k]][i, j] <- dist_euclid(kms[[k]]$centers[i,], kms[[k+1]]$centers[j,]) # Entrance i,j is the dist between old center and new center 
    }
  }

  # WE will align the clusters and identify the "new cluster
  # Use gtools permutations to get all possible permutations of pairings
  # Note we leave one out which is the new cluster so if perms= 2 3 4 1 5
  # then it means that:
  # old 1 -> new 2, old 2 -> new 3, old 3 -> new 4, old 4 -> new 1, old 5 -> new 5 
  perms <- permutations(k, k)
  
  best_sum_distances <- numeric(k + 1) # Vector to store the best sum of distances between cluster centers
  best_perm_for_i <- matrix(NA_integer_, nrow = k + 1, ncol = k)  # store best permutation for each i

  for (i in 1:(k + 1)) {
    
    # build the k x k distance mat w/o col i
    dist_mat_wo_i <- as.matrix(dist_mats[[k]][, -i])
    
    # Find best permutation leaving out i
    best_sum <- Inf
    best_perm <- rep(NA_integer_, k)
    
    # loop over each permuatation  
    for (perm in 1:nrow(perms)) {
      total <- 0
      # loop over each element in old clustering
      for (j in 1:k) {
        # We add the distance from old center to the new center as per
        # current permutation (perm)
        total <- total + dist_mat_wo_i[j, perms[perm, j]]
      }
      
      # Update best to reflect 
      if (total < best_sum) {
        best_sum <- total
        best_perm <- perms[perm, ]
      }
    }
    
    best_sum_distances[i] <- best_sum
    best_perm_for_i[i, ] <- best_perm
  }
  
  best_removed_i <- which.min(best_sum_distances)
  
  pair_old_to_new <- numeric(k)
  perm <- best_perm_for_i[best_removed_i, ]
  
  for (i in 1:k) { # pair the old clusters to the 5 new w/o i. 
    # note that the permutation is an "indexer", so if perm[1]=3 then
    # it means that the old 1 -> new 3, but new 3 might corrospond to something
    # different as we removed i. e.g. if i=2 then 3 corresponds to 4
    pair_old_to_new[i] <- (1:(k+1))[-best_removed_i][perm[i]]
  }

  # Relabeling
  new_to_old_map <- numeric(1+k) + k+1
  
  for (i in 1:k) {
    new_to_old_map[pair_old_to_new[i]] <- i
  }
  km_k1_aligned[[k]] <- new_to_old_map[kms[[k+1]]$cluster]

  # So we only count rim data if points switches between assigned clusters which is 
  # NOT the "new cluster"
  is_rim[[k]] <- (km_k1_aligned[[k]] != kms[[k]]$cluster) & (km_k1_aligned[[k]] != k + 1)
}

rim_rates <- numeric(max(ks)-1)
for(k in 1:(max(ks)-1))
{
  rim_rates[k] <- sum(is_rim[[k]])/length(is_rim[[k]])
  cat("k = ", k, "k + 1 = ",k+1, "rim rate = ", rim_rates[k], "\n")
}

plot(rim_rates, type="b", main=paste("Rim data frequency by # cluster, d=",d), ylim=c(0,1))

