#orthoptera
library(ape)
library(chromePlus)
library(diversitree)
library(phytools)


dat <- read.csv("../data/chrom.data/orthoptera.csv")[,c(4,10)]

trees <- read.nexus("../data/trees/sylvester-2020-orthoptera.nex")

#remove nas from dat
#dat <- na.omit(dat)




pruned.trees <- list()
listdat <- list()
for (i in 1:length(trees)){
  print(i)
  overlap <- trees[[i]]$tip.label %in% dat$species
  keep <- trees[[i]]$tip.label[overlap]
  pruned.trees[[i]] <- keep.tip(trees[[i]], keep)
  hapcount <- vector(length=length(pruned.trees[[i]]$tip.label))
  names(hapcount) <- pruned.trees[[i]]$tip.label
  for(j in 1:length(pruned.trees[[i]]$tip.label)){
    hit <- which(dat$species == names(hapcount)[j])
    if(length(hit)>1){
      hit <- sample(hit, 1)
    }
    hapcount[j] <- dat$hap[hit]
  }
  hapcount <- data.frame(names(hapcount), hapcount)
  listdat[[i]] <- hapcount
}

# prep some data use datatoMatrix
chrom.matlist <- list()
for (k in 1:length(listdat)){
  chrom.matlist[[k]] <- datatoMatrix(x = listdat[[k]],
                                     range=c(3, max(listdat[[k]]$hapcount+1)), hyper = F)
}

conlik <- list()
for (j in 1:length(pruned.trees)){
  lik <- make.mkn(tree = pruned.trees[[j]], 
                  states = chrom.matlist[[j]],
                  k = ncol(chrom.matlist[[j]]),
                  strict = F,
                  control = list(method="ode"))
  conlik[[j]] <- constrainMkn(chrom.matlist[[j]], lik, F, F,
                               constrain=list(drop.poly=T, drop.demi=T ))
}


# fit with diversitree

fit <- mcmc(lik = conlik[[i]], x.init = runif(2), nsteps = 100, w = 1,
            prior = make.prior.exponential(2))
w <- diff(sapply(fit[,2:3], quantile, c(.05,.95)))

fitlist <- list()

for (i in 1:100){
  fitlist[[i]] <- mcmc(lik = conlik[[i]], x.init = runif(2), nsteps = 1000,
                       w = w, prior = make.prior.exponential(2))
  
}

saveRDS(fitlist, "orthoptera-mod1_fitlist.RData")

######Polyploidy
conlik2 <- list()
for (j in 1:length(pruned.trees)){
  lik <- make.mkn(tree = pruned.trees[[j]], 
                  states = chrom.matlist[[j]],
                  k = ncol(chrom.matlist[[j]]),
                  strict = F,
                  control = list(method="ode"))
  conlik2[[j]] <- constrainMkn(chrom.matlist[[j]], lik, F, F,
                              constrain=list(drop.poly=F, drop.demi=T ))
}


# fit with diversitree

fit2 <- mcmc(lik = conlik2[[i]], x.init = runif(3), nsteps = 100, w = 1,
            prior = make.prior.exponential(2))
w2 <- diff(sapply(fit2[,2:4], quantile, c(.05,.95)))

fitlist2 <- list()

for (i in 1:100){
  fitlist2[[i]] <- mcmc(lik = conlik2[[i]], x.init = runif(3), nsteps = 1000,
                       w = w2, prior = make.prior.exponential(2))
  
}

saveRDS(fitlist2, "orthoptera-mod2_fitlist.RData")




