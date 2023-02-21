#odonata 
library(ape)
library(chromePlus)
library(diversitree)
library(phytools)

dat <- read.csv("../data/chrom.data/odonata.csv")
colnames(dat)[7] <- "species"
dat$haploid <- dat$dipliod.autosome.count/2
tree <- read.tree("../data/trees/letsch-2016-odonata.tre")

# get number of matches for grins
sum(sub("_", " ", x=tree$tip.label) %in% dat$species)

overlap <- sub("_", " ", x=tree$tip.label) %in% dat$species
keep <- tree$tip.label[overlap]
pruned.tree <- keep.tip(tree, keep)
dat <- dat[dat$species %in% sub("_", " ", x=tree$tip.label), ]
dat <- data.frame(dat$species, dat$haploid)
colnames(dat) <- c("species", "haploid")
##there are duplicated species. randomly sample duplicated values
#shuffle data
dat1 <- dat[sample(1:nrow(dat)),]
dat1 <- dat1[!duplicated(dat1$species),]

# this data set was run and provides similar results with no qualitative changes
dat2 <- dat[sample(1:nrow(dat)),]
dat2 <- dat2[!duplicated(dat2$species),]

# prep some data use datatoMatrix

chrom.mat1 <- datatoMatrix(x = dat1, range=c(3, max(dat1$haploid)+1), hyper = F)
rownames(chrom.mat1) <- sub(" ", "_", x=rownames(chrom.mat1))


pruned.tree <- force.ultrametric(pruned.tree)


pruned.tree <- multi2di(pruned.tree)
# make an mkn model in diversitree
lik1 <- make.mkn(tree = pruned.tree, 
                states = chrom.mat1,
                k = ncol(chrom.mat1),
                strict = F,
                control = list(method="ode"))

# constrain it with chromeplus
argnames(lik1)
conlik1 <- constrainMkn(chrom.mat1, lik1, F, F, constrain=list(drop.poly=T, drop.demi=T ))



argnames(conlik1)
# fit with diversitree
fit_sample1 <- mcmc(lik = conlik1, x.init = runif(2), nsteps = 100, w = 1,
            prior = make.prior.exponential(2))


w1 <- diff(sapply(fit_sample1[,2:3], quantile, c(.05,.95)))

plot(fit_sample1$p)

fit_sample1 <- mcmc(lik = conlik1, x.init = runif(2), nsteps = 1000, w = w1,
            prior = make.prior.exponential(2))
plot(fit_sample1$p)

write.csv(fit_sample1,"../results/odonata-mod1_sample1.csv")




#now do it again for second sample

chrom.mat2 <- datatoMatrix(x = dat2, range=c(1, max(dat2$haploid)+1), hyper = F)
rownames(chrom.mat2) <- sub(" ", "_", x=rownames(chrom.mat2))


pruned.tree <- force.ultrametric(pruned.tree)


pruned.tree <- multi2di(pruned.tree)
# make an mkn model in diversitree
lik2 <- make.mkn(tree = pruned.tree, 
                 states = chrom.mat2,
                 k = ncol(chrom.mat2),
                 strict = F,
                 control = list(method="ode"))

# constrain it with chromeplus
argnames(lik2)
conlik2 <- constrainMkn(chrom.mat2, lik2, F, F, constrain=list(drop.poly=T, drop.demi=T ))



argnames(conlik2)
# fit with diversitree
fit_sample2 <- mcmc(lik = conlik2, x.init = runif(2), nsteps = 100, w = 1,
                    prior = make.prior.exponential(2))


w2 <- diff(sapply(fit_sample2[,2:3], quantile, c(.05,.95)))

plot(fit_sample2$p)

fit_sample2 <- mcmc(lik = conlik2, x.init = runif(2), nsteps = 1000, w = w2,
                    prior = make.prior.exponential(2))
plot(fit_sample2$p)

write.csv(fit_sample2,"../results/odonata-mod1_sample2.csv")


###samples are the same - going with sample 1

#####polyploidy
conlik3 <- constrainMkn(chrom.mat1, lik1, F, F, constrain=list(drop.poly=F, drop.demi=T ))
argnames(conlik3)
# fit with diversitree
fit3 <- mcmc(lik = conlik3, x.init = runif(3), nsteps = 100, w = 1,
             prior = make.prior.exponential(2))
w3 <- diff(sapply(fit3[,2:4], quantile, c(.05,.95)))

plot(fit3$p)

fit3 <- mcmc(lik = conlik3, x.init = runif(3), nsteps = 1000, w = w3,
               prior = make.prior.exponential(2))

write.csv(fit3,"../results/odonata-mod2.csv", row.names = F)







