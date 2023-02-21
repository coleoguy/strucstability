library(chromePlus)
library(diversitree)

###### Lepidoptera #######
# Read and prep chromosome data
dat <- read.csv("../data/chrom.data/lepidoptera.csv")[,c(5,16)]

#calculate haploid number
dat$haploid <- dat$male2N/2

#remove NA rows
dat <- na.omit(dat)

#remove subspecies information
dat$binomial <- sub("^(\\S*\\s+\\S+).*", "\\1", dat$binomial)

#remove perfectly duplicatedrows
dat <- dat[!duplicated(dat),]

#sub spaces for underscores in the dat column
dat$binomial <- sub(" ", "_", x=dat$binomial)

# read in one randomly drawn tree
tree <- read.tree("../data/trees/wiemers-2020-lepidoptera.trees")[sample(1:1000, 1)][[1]]
pruned.tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% dat$binomial])
hapcount <- vector(length=length(pruned.tree$tip.label))
names(hapcount) <- pruned.tree$tip.label
for(j in 1:length(pruned.tree$tip.label)){
  hit <- which(dat$binomial == names(hapcount)[j])
  if(length(hit)>1){
    hit <- sample(hit, 1)
  }
  hapcount[j] <- dat$haploid[hit]
}
hapcount <- data.frame(names(hapcount), hapcount)

# remove polymattus (highest variance genus)
hapcount2 <- hapcount[-71:-91,]

# reprune tree
pruned.tree2 <- keep.tip(tree, tree$tip.label[tree$tip.label %in% row.names(hapcount2)])

# prep some data use datatoMatrix
chrom.mat2 <- datatoMatrix(x = hapcount2,
                          range=c(min(hapcount2$hapcount-1), 
                                  max(hapcount2$hapcount+1)), 
                          hyper = F)

lik2 <- make.mkn(tree = pruned.tree2, 
                states = chrom.mat2,
                k = ncol(chrom.mat2),
                strict = F,
                control = list(method="ode"))
conlik2 <- constrainMkn(chrom.mat2, lik2, F, F,
                       constrain=list(drop.poly=T, drop.demi=T ))
result2 <- mcmc(lik = conlik2, x.init = runif(2), nsteps = 500, w = 1,
               prior = make.prior.exponential(2))

##### Done with Leps #####


#####  Beetles  ########

dat <- read.csv("../../data/chrom.data/coleoptera.csv")
tree <- read.nexus("../../data/trees/blackmon-2014-coleoptera.nexus")[[1]]

#need to remove first capital letters from trees - all characters until first
#underscore (and including that first underscore)
for(j in 1:length(tree$tip.label)){
  curtip <- tree$tip.label[j]
  tipvec <- strsplit(curtip, "_", fixed=F)[[1]]
  tree$tip.label[j] <- paste(tipvec[-1], collapse="_")
}

#sub spaces for underscores in the dat column
dat$species <- sub(" ", "_", x=dat$species)

#delete spaces at beginning and end of species name
dat$species <- str_trim(dat$species, "both")

#get number of matches
sum(dat$species %in% tree$tip.label)

keep <- tree$tip.label[tree$tip.label %in% dat$species]
pruned.tree <- keep.tip(tree, keep)
hapcount <- vector(length=length(pruned.tree$tip.label))
names(hapcount) <- pruned.tree$tip.label
for(j in 1:length(pruned.tree$tip.label)){
  hit <- which(dat$species == names(hapcount)[j])
  if(length(hit)>1){
    hit <- sample(hit, 1)
  }
  hapcount[j] <- dat$haploid.homo[hit]
}
hapcount <- data.frame(names(hapcount), hapcount)

# record genera names and remove Calathus high variance clade
# with more than 10 species to keep it similar ot Leps
genera <- c()
for(i in 1:nrow(hapcount)){
  genera[i] <- strsplit(row.names(hapcount)[i], "_")[[1]][1]
}
hapcount2 <- hapcount[genera != "Calathus", ]
keep <- tree$tip.label[tree$tip.label %in% hapcount2$names.hapcount.]
pruned.tree <- keep.tip(tree, keep)

chrom.mat <- datatoMatrix(x = hapcount2,
                          range=c(1, max(hapcount2$hapcount+1)), hyper = F)
lik <- make.mkn(tree = pruned.tree, 
                states = chrom.mat,
                k = ncol(chrom.mat),
                strict = F,
                control = list(method="ode"))
conlik <- constrainMkn(chrom.mat, lik, F, F,
                       constrain=list(drop.poly=T, drop.demi=T ))
result <- diversitree::mcmc(lik = conlik, x.init = runif(2), nsteps = 500,
                            w = 1, prior = make.prior.exponential(2))


##### Done with Coleoptera #######


