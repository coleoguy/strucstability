#coleoptera
library(ape)
library(chromePlus)
library(diversitree)
library(stringr)
library(phytools)
library(evobiR)

dat <- read.csv("../data/chrom.data/coleoptera.csv")
trees <- read.nexus("../data/trees/blackmon-2014-coleoptera.nexus")

#need to remove first capital letters from trees - all characters until first
#underscore (and including that first underscore)
named.trees <- list()
for(i in 1:100){
  curtree <- trees[[i]]
  for(j in 1:length(trees[[1]]$tip.label)){
    curtip <- trees[[i]]$tip.label[j]
    tipvec <- strsplit(curtip, "_", fixed=F)[[1]]
    curtree$tip.label[j] <- paste(tipvec[-1], collapse="_")
  }
  named.trees[[i]] <- curtree
}
class(named.trees) <- class(trees)

#sub spaces for underscores in the dat column
dat$species <- sub(" ", "_", x=dat$species)

#delete spaces at beginning and end of species name
dat$species <- str_trim(dat$species, "both")

#get number of matches
sum(dat$species %in% named.trees[[1]]$tip.label)
rm(curtree,trees,curtip,i,j,tipvec)




pruned.trees <- list()
listdat <- list()
for (i in 1:length(named.trees)){
  print(i)
  keep <- named.trees[[i]]$tip.label[named.trees[[i]]$tip.label %in% dat$species]
  pruned.trees[[i]] <- keep.tip(named.trees[[i]], keep)
  hapcount <- vector(length=length(pruned.trees[[i]]$tip.label))
  names(hapcount) <- pruned.trees[[i]]$tip.label
  for(j in 1:length(pruned.trees[[i]]$tip.label)){
    hit <- which(dat$species == names(hapcount)[j])
    if(length(hit)>1){
      hit <- sample(hit, 1)
    }
    hapcount[j] <- dat$haploid.homo[hit]
  }
  hapcount <- data.frame(names(hapcount), hapcount)
  listdat[[i]] <- hapcount
}


# prep some data use datatoMatrix
chrom.matlist <- list()
for (k in 1:length(listdat)){
  chrom.matlist[[k]] <- datatoMatrix(x = listdat[[k]],
                                     range=c(1, max(listdat[[k]]$hapcount+1)), hyper = F)
  print(ncol(chrom.matlist[[k]]))
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

# fit <- diversitree::mcmc(lik = conlik[[i]], x.init = runif(2), nsteps = 100, w = 1,
#             prior = make.prior.exponential(2))
# w <- diff(sapply(fit[50:100,2:3], quantile, c(.05,.95)))
w <- c(0.02582563, 0.01621471)


library(doMC)
registerDoMC(8)

x <- foreach(i=1:length(conlik)) %dopar% {
  mcmc(lik = conlik[[i]], x.init = runif(2), nsteps = 1000,
       w = w, prior = make.prior.exponential(2))
}


saveRDS(x, "coleopteraresults_mod1.RData")

####polyploidy model
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

fit2 <- diversitree::mcmc(lik = conlik2[[i]], x.init = runif(3, 0,.1), nsteps = 100, w = .1,
            prior = make.prior.exponential(2))
w2 <- diff(sapply(fit2[50:100,2:4], quantile, c(.05,.95)))
w2 <- c(0.01513523, 0.01193033, 0.001211298)

x2 <- foreach(i=1:length(conlik2), .verbose = T) %dopar% {
  print(paste("working on tree:", i))
  foo <- diversitree::mcmc(lik = conlik2[[i]], x.init = runif(3,0,.05), nsteps = 1000,
       w = w2, prior = make.prior.exponential(2))
  write.csv(foo, file=paste("coleoptera-mod2-", i, ".csv",sep=""), row.names = F)
}
saveRDS(x2, "coleopteraresults_mod2.RData")
