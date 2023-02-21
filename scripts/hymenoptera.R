#hymenoptera
library(ape)
library(chromePlus)
library(diversitree)
library(doMC)
dat <- read.csv("../data/chrom.data/hymenoptera.csv")[,c(3:6,16)]
dat$species <- paste(dat$Genus, dat$species)
tree <- read.tree("../data/trees/ross-2015-hymenoptera.tre")

# get number of matches for grins
sum(sub("_", " ", x=tree$tip.label) %in% dat$species)

overlap <- sub("_", " ", x=tree$tip.label) %in% dat$species
keep <- tree$tip.label[overlap]
pruned.tree <- keep.tip(tree, keep)
dat <- dat[dat$species %in% sub("_", " ", x=tree$tip.label), 4:5]

##scale tree depth
#treedepth <- max(branching.times(pruned.tree))
#pruned.tree$edge.length <- pruned.tree$edge.length / max(branching.times(pruned.tree))

# prep some data use datatoMatrix

chrom.mat <- datatoMatrix(x = dat, range=c(1, max(dat$Chromosome.number..male..2N)+1), hyper = F)
rownames(chrom.mat) <- sub(" ", "_", x=rownames(chrom.mat))
# make an mkn model in diversitree
lik <- make.mkn(tree = pruned.tree, 
                states = chrom.mat,
                k = 58,
                strict = F,
                control = list(method="ode"))

# constrain it with chromeplus
argnames(lik)
conlik <- constrainMkn(chrom.mat, lik, F, F, constrain=list(drop.poly=T, drop.demi=T ))
argnames(conlik)
# fit with diversitree
fit <- mcmc(lik = conlik, x.init = runif(2), nsteps = 100, w = 1,
            prior = make.prior.exponential(2))
w <- diff(sapply(fit[,2:3], quantile, c(.05,.95)))

plot(fit$p)


fit <- mcmc(lik = conlik, x.init = runif(2), nsteps = 1000, w = w,
            prior = make.prior.exponential(2))
plot(fit$p)

#rescale to millions of years
#fit[,2:3] <- fit[,2:3]/treedepth


write.csv(fit,"../results/hymenoptera-mod1.csv", row.names = F)


####polyploidy model
conlik2 <- constrainMkn(chrom.mat, lik, F, F, constrain=list(drop.poly=F, drop.demi=T ))
argnames(conlik2)
# fit with diversitree
fit2 <- mcmc(lik = conlik2, x.init = runif(3), nsteps = 100, w = 1,
            prior = make.prior.exponential(2))
w <- diff(sapply(fit2[,2:4], quantile, c(.05,.95)))

plot(fit2$p)


#registerDoMC(6)
# x <- foreach(i = 1:6) %dopar% {
#   fit2 <- mcmc(lik = conlik2, x.init = runif(3), nsteps = 200, w = w,
#                prior = make.prior.exponential(2))
#   
# }
fit2 <- mcmc(lik = conlik2, x.init = runif(3), nsteps = 1000, w = w,
             prior = make.prior.exponential(2))


write.csv(fit2,"../results/hymenoptera-mod2.csv", row.names = F)


