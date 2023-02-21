#lepidoptera

library(ape)
library(devtools)
#install_github('coleoguy/chromePlus')
library(chromePlus)
library(diversitree)
library(phytools)
library(maps)
library(doSNOW)

dat <- read.csv("../data/chrom.data/lepidoptera.csv")[,c(5,16)]
#calculate haploid number
dat$haploid <- dat$male2N/2
#remove NA rows
dat <- na.omit(dat)
#remove subspecies information
dat$binomial <- sub("^(\\S*\\s+\\S+).*", "\\1", dat$binomial)
#remove perfectly duplicatedrows
dat <- dat[!duplicated(dat),]



tree <- read.tree("../data/trees/wiemers-2020-lepidoptera.trees")[sample(1:1000, 100)]

#sub spaces for underscores in the dat column
dat$binomial <- sub(" ", "_", x=dat$binomial)

# get number of matches for grins
sum(tree[[1]]$tip.label %in% dat$binomial)

##there are duplicated species. randomly sample duplicated values
#shuffle data
#dat1 <- dat[sample(1:nrow(dat)),]
#dat1 <- dat1[!duplicated(dat1$binomial),]
#dat2 <- dat[sample(1:nrow(dat)),]
#dat2 <- dat2[!duplicated(dat2$binomial),]



###
pruned.trees <- list()
listdat <- list()
for (i in 1:length(tree)){
  print(i)
  overlap <- tree[[i]]$tip.label %in% dat$binomial
  keep <- tree[[i]]$tip.label[overlap]
  pruned.trees[[i]] <- keep.tip(tree[[i]], keep)
  hapcount <- vector(length=length(pruned.trees[[i]]$tip.label))
  names(hapcount) <- pruned.trees[[i]]$tip.label
  for(j in 1:length(pruned.trees[[i]]$tip.label)){
    hit <- which(dat$binomial == names(hapcount)[j])
    if(length(hit)>1){
      hit <- sample(hit, 1)
    }
    hapcount[j] <- dat$haploid[hit]
  }
  hapcount <- data.frame(names(hapcount), hapcount)
  listdat[[i]] <- hapcount
}


# prep some data use datatoMatrix
chrom.matlist <- list()
for (k in 1:length(listdat)){
  chrom.matlist[[k]] <- datatoMatrix(x = listdat[[k]],
                                     range=c(7, max(listdat[[k]]$hapcount+1)), hyper = F)
}


lik <- make.mkn(tree = pruned.trees[[1]], 
                states = chrom.matlist[[1]],
                k = ncol(chrom.matlist[[1]]),
                strict = F,
                control = list(method="ode"))
conlik1 <- constrainMkn(chrom.matlist[[1]], lik, F, F,
                             constrain=list(drop.poly=F, drop.demi=T ))
argnames(conlik1)
#fit <- mcmc(lik = conlik1, x.init = runif(3), nsteps = 100, w = 1,
#            prior = make.prior.exponential(2))
#w <- diff(sapply(fit[,2:4], quantile, c(.05,.95)))
# based on original run
 w <- c(0.32381,0.20206,0.00558)
 
for (j in 1:100){
   lik <- make.mkn(tree = pruned.trees[[j]], 
                         states = chrom.matlist[[j]],
                         k = ncol(chrom.matlist[[j]]),
                         strict = F,
                         control = list(method="ode"))
   conlik1 <- constrainMkn(chrom.matlist[[j]], lik, F, F,
                                constrain=list(drop.poly=F, drop.demi=T ))
   argnames(conlik1)
   conlik1(c(.1,.1,.1))
   fit <- mcmc(lik = conlik1, x.init = runif(3), nsteps = 1000, w = w,
               prior = make.prior.exponential(2))
   write.csv(fit, file=paste("lepmod1-",j,".csv",sep=""))
}

 
 
##it failed on tree3 - try starting at tree 4
nClust <- 25
cl <- makeCluster(nClust, outfile = "")
registerDoSNOW(cl)

x <- foreach (i = 4:100,
              .verbose = T,
              .packages = c("phytools", "ape", "maps", "chromePlus", "diversitree")) %dopar% {
                lik <- make.mkn(tree = pruned.trees[[j]], 
                                states = chrom.matlist[[j]],
                                k = ncol(chrom.matlist[[j]]),
                                strict = F,
                                control = list(method="ode"))
                conlik1 <- constrainMkn(chrom.matlist[[j]], lik, F, F,
                                        constrain=list(drop.poly=F, drop.demi=T ))
                fit <- mcmc(lik = conlik1, x.init = runif(3), nsteps = 1000, w = w,
                            prior = make.prior.exponential(2))
                write.csv(fit, file=paste("lepmod2-", i, sep=""))
                
} 
stopCluster(cl)
saveRDS(x, "lepidoptera-mod2.RDS")


#write.csv(fit, file=paste("lepmod1-",j,".csv",sep=""))

# do a short initial run to geta  sense of appropriate values for w
#fit <- mcmc(lik = conlik1[[i]], x.init = runif(2), nsteps = 100, w = 1,
#            prior = make.prior.exponential(2))
#w <- diff(sapply(fit[,2:3], quantile, c(.05,.95)))
# based on original run
#w <- c(6.99, 6.68)

#saveRDS(conlik1, "conlik1")
#saveRDS(conlik1, "conlik1_real.RData")
#saveRDS(w, "w")
#saveRDS(w, "w_real.RData")



#conlik1 <- readRDS("conlik1_real.RData")
#w <- readRDS("w_real.RData")
#library(doMC)
#registerDoMC(6)


# define number of clusters for parallel computing
# number of CPUs
# make connections for each CPU
# setting outfile = "" would print the output on the console
# however in Windows OS you will not see an output unless you
# run the R script through command prompt,powershell, or
# windows subsystem for linux


nClust <- 60
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
saveRDS(x, "lepidopteraresults.RData")
