#diptera
library(ape)
library(chromePlus)
library(diversitree)
library(phytools)

dat <- read.csv("../data/chrom.data/diptera-2022-12-08.csv")

tree <- read.nexus("../data/trees/flytree.chrono.tre")


unique(dat$Family[dat$Family %in% tree$tip.label])

matches <- c()
for(i in 1:nrow(dat)){
  hit <- dat$Family[i] %in% tree$tip.label
  if(hit){
    matches <- c(matches, i)
  }
}

pruned.data <- dat[matches,]

pruned.tree <- keep.tip(tree, pruned.data$Family)

plotTree(pruned.tree)

smpvals <- c()
for(i in 1:length(pruned.tree$tip.label)){
  curfam <- pruned.tree$tip.label[i]
  possvals <- pruned.data$Haploid.Number[pruned.data$Family == curfam]
  if(length(possvals)>1){
    smpvals[i] <- sample(possvals, 1)
  }else{
    smpvals[i] <- possvals
  }
}

names(smpvals) <- pruned.tree$tip.label
print(smpvals)


chrom.mat <- datatoMatrix(x = data.frame(pruned.tree$tip.label, smpvals),
                          range = c(2, 14),
                          hyper = F)

pruned.tree <- force.ultrametric(pruned.tree)
pruned.tree <- multi2di(pruned.tree)

lik <- make.mkn(tree = pruned.tree, 
                states = chrom.mat,
                k = ncol(chrom.mat),
                strict = F)

con.lik <- constrainMkn(data = chrom.mat,
                        lik = lik,
                        hyper = F,
                        polyploidy = F,
                        constrain = list(drop.poly = T, drop.demi = T))

argnames(con.lik)

fit <- mcmc(lik = con.lik,
            x.init = runif(2),
            nsteps = 100,
            w = 1,
            prior = make.prior.exponential(2))
w <- diff(sapply(fit[,2:3], quantile, c(.05,.95)))

fit <- mcmc(lik = con.lik,
            x.init = runif(2),
            nsteps = 1000,
            w = w,
            prior = make.prior.exponential(2))


write.csv(fit,"../results/diptera-mod1.csv", row.names = F)

###polyploidy model
conlik2 <- constrainMkn(chrom.mat, lik, F, F, constrain=list(drop.poly=F, drop.demi=T ))
argnames(conlik2)
# fit with diversitree
fit2 <- mcmc(lik = conlik2, x.init = runif(3), nsteps = 100, w = 1,
             prior = make.prior.exponential(2))
w <- diff(sapply(fit2[,2:4], quantile, c(.05,.95)))
fit2 <- mcmc(lik = conlik2, x.init = runif(3), nsteps = 1000, w = w,
             prior = make.prior.exponential(2))
write.csv(fit2,"../results/diptera-mod2.csv", row.names = F)

