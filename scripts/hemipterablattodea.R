# This script is used to estimate rates of chromosome
# number evolution in Hemiptera and Blattodea. It is amended from Ruckman 2020.

# load packages
library(ape) # basic phylo tools
library(chromePlus) # chromosome models
library(diversitree) # basic likelihood functions
library(doSNOW)
library(phytools)
library(maps)


iter <- 1000
# for troubleshooting we might run on just a couple of trees
# usually ntree should be set equal to the number of trees being
# analyzed.
ntree <- 100

# load custom functions
source("functions_hemipterablattodea.R")

# read in the data and the tree
trees <- read.nexus("../data/trees/misof.backbone.nex")
dat <- read.csv("../data/chrom.data/data.invert.csv", as.is = T)[, -c(6:9)]
#Change Isoptera to Blattodea
dat$Order <- sub("Isoptera", "Blattodea", dat$Order)



# lets find our orders first
# foo <- getData(trees, dat)
# tree <- foo[[1]]
# genera <- row.names(foo[[2]])
# orders <- dat$Order[dat$Genus %in% genera]
# orders <- as.data.frame(table(orders))



#####Prune trees and data
orders <- c("Blattodea", "Hemiptera")

Blattodea_trees_pruned <- getDataOrder(trees, dat, order=orders[1])[[1]]
Hemiptera_trees_pruned <- getDataOrder(trees, dat, order=orders[2])[[1]]


Blattodea_dat <- dat[dat$Order == "Blattodea",]
Hemiptera_dat <- dat[dat$Order == "Hemiptera",]


#####Blattodea

# from prelim analyses we have w and w2
#w <- c(0.192, 0.146)
#w2 <- c(0.4323, 1.1948, 0.0778)

###
pruned.trees <- list()
listdat <- list()
for (i in 1:length(Blattodea_trees_pruned)){
  print(i)
  overlap <- Blattodea_trees_pruned[[i]]$tip.label %in% dat$Genus
  keep <- Blattodea_trees_pruned[[i]]$tip.label[overlap]
  pruned.trees[[i]] <- keep.tip(Blattodea_trees_pruned[[i]], keep)
  hapcount <- vector(length=length(pruned.trees[[i]]$tip.label))
  names(hapcount) <- pruned.trees[[i]]$tip.label
  for(j in 1:length(pruned.trees[[i]]$tip.label)){
    hit <- which(dat$Genus == names(hapcount)[j])
    if(length(hit)>1){
      hit <- sample(hit, 1)
    }
    hapcount[j] <- dat$haploid.num[hit]
  }
  hapcount <- data.frame(names(hapcount), hapcount)
  listdat[[i]] <- hapcount
}


# prep some data use datatoMatrix
chrom.matlist <- list()
for (k in 1:length(listdat)){
  chrom.matlist[[k]] <- datatoMatrix(x = listdat[[k]],
                                     range=c(7, max(listdat[[k]]$hapcount+1)),
                                     hyper = F)
}

nClust <- 20
cl <- makeCluster(nClust, outfile = "")
registerDoSNOW(cl)

x <- foreach (i = 1:100,
              .verbose = T,
              .packages = c("phytools", "ape", "maps", "chromePlus", "diversitree")) %dopar% {
                lik <- make.mkn(tree = pruned.trees[[i]], 
                                states = chrom.matlist[[i]],
                                k = ncol(chrom.matlist[[i]]),
                                strict = F,
                                control = list(method="ode"))
                conlik1 <- constrainMkn(chrom.matlist[[i]], lik, F, F,
                                        constrain=list(drop.poly=T, drop.demi=T ))
                results <- mcmc(lik = conlik1, x.init = runif(2), nsteps = 1000, w = w,
                                prior = make.prior.exponential(2))
              }
stopCluster(cl)
saveRDS(x, "blattodea-mod1.RData")

###polyploidy
nClust <- 20
cl <- makeCluster(nClust, outfile = "")
registerDoSNOW(cl)

x <- foreach (i = 1:100,
              .verbose = T,
              .packages = c("phytools", "ape", "maps", "chromePlus", "diversitree")) %dopar% {
                lik <- make.mkn(tree = pruned.trees[[i]], 
                                states = chrom.matlist[[i]],
                                k = ncol(chrom.matlist[[i]]),
                                strict = F,
                                control = list(method="ode"))
                conlik1 <- constrainMkn(chrom.matlist[[i]], lik, F, F,
                                        constrain=list(drop.poly=F, drop.demi=T ))
                results <- mcmc(lik = conlik1, x.init = runif(3), nsteps = 1000, w = w2,
                                prior = make.prior.exponential(2))
              }
stopCluster(cl)
saveRDS(x, "blattodea-mod2.RData")


##########Hemiptera
pruned.trees <- list()
listdat <- list()
for (i in 1:length(Hemiptera_trees_pruned)){
  print(i)
  overlap <- Hemiptera_trees_pruned[[i]]$tip.label %in% dat$Genus
  keep <- Hemiptera_trees_pruned[[i]]$tip.label[overlap]
  pruned.trees[[i]] <- keep.tip(Hemiptera_trees_pruned[[i]], keep)
  hapcount <- vector(length=length(pruned.trees[[i]]$tip.label))
  names(hapcount) <- pruned.trees[[i]]$tip.label
  for(j in 1:length(pruned.trees[[i]]$tip.label)){
    hit <- which(dat$Genus == names(hapcount)[j])
    if(length(hit)>1){
      hit <- sample(hit, 1)
    }
    hapcount[j] <- dat$haploid.num[hit]
  }
  hapcount <- data.frame(names(hapcount), hapcount)
  listdat[[i]] <- hapcount
}

# prep some data use datatoMatrix
chrom.matlist <- list()
for (k in 1:length(listdat)){
  chrom.matlist[[k]] <- datatoMatrix(x = listdat[[k]],
                                     range=c(1, max(listdat[[k]]$hapcount+1)),
                                     hyper = F)
}
#these w's were from previous analyses
#w <- c(0.0967, 0.099)
#w2 <- c(0.6602, 0.7299, 0.0959)

#######
nClust <- 20
cl <- makeCluster(nClust, outfile = "")
registerDoSNOW(cl)

x <- foreach (i = 1:100,
              .verbose = T,
              .packages = c("phytools", "ape", "maps", "chromePlus", "diversitree")) %dopar% {
                lik <- make.mkn(tree = pruned.trees[[i]], 
                                states = chrom.matlist[[i]],
                                k = ncol(chrom.matlist[[i]]),
                                strict = F,
                                control = list(method="ode"))
                conlik1 <- constrainMkn(chrom.matlist[[i]], lik, F, F,
                                        constrain=list(drop.poly=T, drop.demi=T ))
                results <- mcmc(lik = conlik1, x.init = runif(2), nsteps = 1000, w = w,
                                prior = make.prior.exponential(2))
              }
stopCluster(cl)
saveRDS(x, "hemiptera-mod1.RData")


###polyploidy
nClust <- 20
cl <- makeCluster(nClust, outfile = "")
registerDoSNOW(cl)

x <- foreach (i = 1:100,
              .verbose = T,
              .packages = c("phytools", "ape", "maps", "chromePlus", "diversitree")) %dopar% {
                lik <- make.mkn(tree = pruned.trees[[i]], 
                                states = chrom.matlist[[i]],
                                k = ncol(chrom.matlist[[i]]),
                                strict = F,
                                control = list(method="ode"))
                conlik1 <- constrainMkn(chrom.matlist[[i]], lik, F, F,
                                        constrain=list(drop.poly=F, drop.demi=T ))
                results <- mcmc(lik = conlik1, x.init = runif(3), nsteps = 1000, w = w2,
                                prior = make.prior.exponential(2))
              }
stopCluster(cl)
saveRDS(x, "hemiptera-mod2.RData")
