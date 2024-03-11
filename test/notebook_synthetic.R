source("utils.R")

############################################################
## In this notebook, we illustrate the use of the msPCA package on synthetic data.
## We generate a spiked-covariance matrix with 2 4-sparse spikes (or PCs).
##
## We generate datasets of increasing size consistent with this model, compute
## the empirical covariance matrix, and apply our msPCA algorithm on the data.
## We measure performance in terms of support recovery of the true PCs, fraction
## of variance explained by the PCs, and orthogonality violation.
##
## We compare with the spca method from the elasticnet package.
############################################################

library(mvtnorm)
library(msPCA)
library(elasticnet)

p = 50 #Dimension
r = 2 #Number of sparse PCs
k = 20 #Sparsity of each PC
b = 2 #Signal strength

#Construct random -1,0,+1 vectors with disjoint supports (use 1/sqrt(k) to ensure norm 1)
xtrue <- matrix(0,p,r)
xtrue[1:k,1] <- sign(runif(k) - .5)/sqrt(k)
xtrue[(k+1):(2*k),2] <- sign(runif(k) - .5)/sqrt(k)

shufflecoords <- sample(1:p) #Shuffle coordinates
xtrue <- xtrue[shufflecoords,]

#Compute support recover accuracy metric
accuracy <- function(xtrue, U){
  a1 <- sum( (abs(U) > 0)*(abs(xtrue) > 0) )
  a2 <- sum( (abs(U[,c(2,1)]) > 0)*(abs(xtrue) > 0) )
  max(a1,a2) / sum(abs(xtrue) > 0)
}

#True covariance matrix
Strue <- b*outer(xtrue[,1],xtrue[,1]) + b*outer(xtrue[,2],xtrue[,2]) +diag(p)


## We comare the algorithms on varying
resdf <- data.frame(n=numeric(), package=character(),
                   sparsity=numeric(), accuracy=numeric(),
                   ortho_viol=numeric(), varexplained=numeric())
icol <- 1
for (iter in 1:5){ #Number of replications

  Xfull <- mvtnorm::rmvnorm(2000, mean = numeric(p), sigma = Strue)

  for (n in c(25,50,75,100,150,200,250,500,1000,1500)){ #Size of the dataset used to recover the PCs
  #for (n in seq(25, 250, by = 25)){ #Size of the dataset used to recover the PCs
    print(n)

    X <- Xfull[1:n,]
    S <- cov(X)

    #msPCA algorithm
    mspca_results <- msPCA::iterativeDeflationHeuristic(S, 2, c(k,k), verbose=FALSE)
    resdf[icol,] <- list(n,"msPCA",
                   sum(abs(mspca_results$x_best) > 0), #Sparsity level
                   accuracy(xtrue, mspca_results$x_best), #Support recovery accuracy
                   orthogonality_violation(mspca_results$x_best), #Orthogonality violation
                   fraction_variance_explained(Strue,mspca_results$x_best)) #Fraction of variance explained

    icol <- icol +1

    #elasticnet algorithm
    enet_results <- elasticnet::spca(S, 2, sparse="varnum", para = c(k,k), type="Gram")
    resdf[icol,] <- list(n,"elasticnet",
                      sum(abs(enet_results$loadings) > 0),
                      accuracy(xtrue, enet_results$loadings),
                      orthogonality_violation(enet_results$loadings),
                      fraction_variance_explained(Strue,enet_results$loadings))
    icol <- icol +1

  }
 }

library(readr)
write_csv(resdf, "msPCA_synthetic_results.csv") #For saving results

resdf <- read_csv("msPCA_synthetic_results.csv")

library(dplyr)
library(ggplot2)
library(RColorBrewer)
std.error <- function(x){sd(x)/sqrt(length(x))}

sumdf <- resdf %>%
  group_by(n,package) %>%
  summarise(accuracy=mean(accuracy),
            ortho_viol_se = std.error(ortho_viol),
            ortho_viol_mean =mean(ortho_viol),
            varexplained_se = std.error(varexplained),
            varexplained_mean =mean(varexplained),
            nreplicates=n()
            ) %>%
  ungroup() %>%
  mutate(ortho_viol_min=ortho_viol_mean-2*ortho_viol_se,
         ortho_viol_max=ortho_viol_mean+2*ortho_viol_se,
         varexplained_min=varexplained_mean-2*varexplained_se,
         varexplained_max=varexplained_mean+2*varexplained_se
         )

sumdf %>%
  #filter(n <= 500) %>%
  ggplot() + aes(x=n, y=ortho_viol_mean, group=package, color=package,shape=package) +
  geom_line(size=0.7) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=ortho_viol_min,ymax=ortho_viol_max), width=0.1) +
  theme_minimal() +
  scale_color_brewer(palette="Set1") +
  labs(x="Sample size n", y="Orthogonality violation", group="Package", shape="Package", color="Package") +
  theme(legend.position="bottom",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor = element_line(colour = "grey90"))

sumdf %>%
  #filter(n <= 500) %>%
  ggplot() + aes(x=n, y=varexplained_mean, group=package, color=package, shape=package) +
  geom_line(size=0.7) + geom_point(size=3) +
  geom_errorbar(aes(ymin=varexplained_min,ymax=varexplained_max), width=0.1) +
  theme_minimal() +
  geom_hline(aes(yintercept=6/54, linetype="Information theoretic upper bound"), color="black") +
  scale_color_brewer(palette="Set1") +
  labs(x="Sample size n", y="Out-of-sample fraction of variance explained",
       linetype="", group="Package", shape="Package", color="Package") +
  theme(legend.position="bottom",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor = element_line(colour = "grey90"))

