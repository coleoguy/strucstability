library(ape)
library(phytools)
tree <- read.tree(text="(((((Zoraptera,Dermaptera),(((((Phasmatodea,Embiidina),(Mantophasmatodea,Grylloblata)),Orthoptera),Plecoptera),(Mantodea,(Blattodea,Isoptera)))),((Thysanoptera,Hemiptera),(Psocoptera,(Hymenoptera,(((Strepsiptera,Coleoptera),(Raphidioptera,(Megaloptera,Neuroptera))),((Lepidoptera,Trichoptera),(Diptera,(Mecoptera,Siphonaptera)))))))),(Odonata,Ephemeroptera)),Zygentoma);")
tree$edge.lengths <- rep(1,54)
tree <- force.ultrametric(tree, method="extend")
plot(tree)
dat <- read.csv("data.invert.csv")
unique(dat$Order)[!(unique(dat$Order)) %in% tree$tip.label]


unique(tree$tip.label)[!tree$tip.label%in% unique(dat$Order) ]
range(dat$haploid.homo)
dat2 <- matrix(NA, length(unique(tree$tip.label)), 223)
row.names(dat2) <- tree$tip.label
for(i in 1:nrow(dat2)){
  vals <- dat$haploid.homo[dat$Order == row.names(dat2)[i]] 
  for(j in 1:ncol(dat2)){
    dat2[i, j] <- sum(vals == j)
  }
}
dat3 <- dat2[,1:100]
dat3[21,100] <-21


new_palette <- c("white", viridis(n = 78, begin = 0, end = 1, option="B"))
phylo.heatmap(tree,log1p(dat3), fsize = c(0.6, 0.6, 1), colors = new_palette, grid = TRUE,
              split = c(0.2, 0.8), legend = TRUE)


