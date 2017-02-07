
#generate sequence of seeds of length equal to the number of parallel tasks.
#make sure that ntasks and the sbatch script agree!
ntasks <- 100L
common_start <- 6415L
library(rngtools)
rng <-rngtools::RNGseq(ntasks,common_start)
# set RNG seed specific to the 100th run
rngtools::setRNG(rng[[ntasks]])



#################################
# Distinguising Simulation Paramters of mod1_hub_p500_*
masterDegree <- 14L ## minimum out-degree of master predictors 
indepX <- 10L ## number of independent X variables. 
n <- 150L ## sample size
#sig.fit is currently  estimated and requires iterative estimates. Known to be 1's in this case.
iter <- 3   ## number of outer-iterations (for sigma^{ii} update) of the space step
# Logical status of assigning (x,y) nodes randomly or not. 
randxy <- FALSE 
#Specify what algorithms to run.
runAlgo <- list(scggm=FALSE,spmap=TRUE,spmapG0=FALSE,spmapEquals=FALSE,space=FALSE)
#Specify if cross-validation is being run
runCV <- TRUE

#Additional parameters 
tol <- 1e-6 ##convergence tolerance
cd_iter <- 1e7 ## max coordinate descent iterations 
#offsetRemMap <- 70 ##additional penalty on remMap (L1) to reduce false discovery rate.
sridge <- 0 ## ridge penalty for elastic net of space
kcv <- 10  ## number of cross-validation folds (if applicable)
#ncores <- 28 ## Number of cores for parallel functionality.

#################################
#generate data: ParCor -> Sig.Inv -> Sig

##get the template partial correlation matrix from space package example
library(space) 
data(spaceSimu)
ParCor.true=spaceSimu$ParCor.true

#get adjacency matrix
true.adj <- abs(ParCor.true) > tol 
diag(true.adj) <- FALSE
# 
# #get hub variables with degree > 11. Assign to X set.
# deg <- rowSums(true.adj)
# hist(deg)
# Xhubs <- which(deg >= 14)
# Xsubhubs <- sample(x = Xhubs, size = 5)
# 
# linkOut <- function(src) { 
#   #Assign Y set to all variables that have edge to X (i.e. targets of X)
#   targets <- sapply(src, function(i) which(true.adj[i,]))
#   targets <- sort(Reduce(union, targets))
#   #remove src in targets
#   targets <- setdiff(targets, src)
# }
# 
# Y1 <- linkOut(Xsubhubs)
# Y2 <- linkOut(Y1)
# subvertices <- sort(Reduce(base::union, list(Xhubs, Y1, Y2)))
# 
# #avg degree of each variable type
# degy <- rowSums(true.adj)[union(Y1,Y2)]
# mean(degy)
# degx <- rowSums(true.adj)[Xsubhubs]
# mean(degx)


#get the disconnected components
library(igraph) 
dc <-  function(x) { 
  sg <- igraph::graph.adjacency(adjmatrix = x)
  sc <- igraph::components(sg)
  scgl <- groups(sc)
  scgl
}
true.dc <- dc(true.adj)
#get size of disconnected component
true.dc.len <- sapply(true.dc, length)
#use 2 of the largest dc
dc.2.idx <- sample(which(true.dc.len > 80),2)
subvertices <- unlist(true.dc[dc.2.idx])

trueParCorXY <- ParCor.true[subvertices,subvertices]
dim(trueParCorXY)
true.adj2 <- abs(trueParCorXY) > tol 
diag(true.adj2) <- FALSE
deg2 <- rowSums(true.adj2)
mean(deg2)

## generate true covariance matrix, and the X, Y indices 
#evaluate model parameters
source("~/repos/spacemap/sim/utils/model1/simulateMod1.R")
model.c <- simSpaceMapModel1(ParCor.true = trueParCorXY, mdegree = masterDegree, 
                             randxy = randxy,
                             pnone = indepX, tol = tol)

# Data generation with the correlation matrix instead of the covariance: 
# The purpose is to make the diagonal elements of the precision matrix to be different from 
# 1 so that the estimates, sig.fit, are not set to the truth by default. 

model.c$Sigma <- cov2cor(model.c$Sigma)

# diag(solve(model.c$Sigma))
# summary(diag(solve(cov2cor(model.c$Sigma))))
# identical(abs(solve(model.c$Sigma)) > tol, abs(solve(cov2cor(model.c$Sigma))) > tol)

## generate data
data.c=simSpaceMapData1(n=n, model=model.c)

#standardize the data
data.c$XY <- scale(data.c$XY)
data.c$X <- scale(data.c$X)
data.c$Y <- scale(data.c$Y)

#extract indices from data generation
Yindex <- data.c$Yindex
Xindex <- data.c$Xindex     
#useful for repeated evaluation of performance against truth
trueParCor <- list(xy = model.c$trueXY, yy = model.c$trueYY)

#The particular model 1 parameter settings generated a true graph comprised of: 
#  ```{r model_dims} 
length(Yindex) ## number of Y nodes
length(Xindex) ## number of X nodes
length(Xindex) - indepX ## number of master predictors.
#```

# 
# tmp <- new.env()
# load(file = "~/../Downloads/cvScore_cov2cor_iter3_sparse_rep_mod1_Hub_p500_11_20_250_dataset_12_results_spmap2.RData", envir = tmp)
# all.equal(tmp$data.c$X[1,], data.c$X[1,])
# tmp$cvSpmap$tuneGrid

sim1 <- data.c
sim1$trueParCor <- trueParCor
library(devtools)
devtools::use_data(sim1)
