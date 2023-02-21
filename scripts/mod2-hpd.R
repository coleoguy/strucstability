library(ggplot2)
library(coda)

######Mod 2
dip <- read.csv("../results/posterior-samples/dipteramod2.csv")
coleo <- read.csv("../results/posterior-samples/coleoptera-mod2.csv")
blatt <- read.csv("../results/posterior-samples/blattodea-mod2.csv")
hemi <- read.csv("../results/posterior-samples/hemiptera-mod2.csv")
hymen <- read.csv("../results/posterior-samples/hymenopteramod2.csv")
lep <- read.csv("../results/posterior-samples/lepidoptera-mod2.csv")
odo <- read.csv("../results/posterior-samples/odonatamod2.csv")
ortho <- read.csv("../results/posterior-samples/orthoptera-mod2.csv")

HPDinterval(as.mcmc(blatt$desc1))
HPDinterval(as.mcmc(blatt$asc1))
HPDinterval(as.mcmc(blatt$pol1))

HPDinterval(as.mcmc(coleo$desc1))
HPDinterval(as.mcmc(coleo$asc1))
HPDinterval(as.mcmc(coleo$pol1))

HPDinterval(as.mcmc(dip$desc1))
HPDinterval(as.mcmc(dip$asc1))
HPDinterval(as.mcmc(dip$pol1))

HPDinterval(as.mcmc(hemi$desc1))
HPDinterval(as.mcmc(hemi$asc1))
HPDinterval(as.mcmc(hemi$pol1))

HPDinterval(as.mcmc(hymen$desc1))
HPDinterval(as.mcmc(hymen$asc1))
HPDinterval(as.mcmc(hymen$pol1))

HPDinterval(as.mcmc(lep$desc1))
HPDinterval(as.mcmc(lep$asc1))
HPDinterval(as.mcmc(lep$pol1))

HPDinterval(as.mcmc(odo$desc1))
HPDinterval(as.mcmc(odo$asc1))
HPDinterval(as.mcmc(odo$pol1))

HPDinterval(as.mcmc(ortho$desc1))
HPDinterval(as.mcmc(ortho$asc1))
HPDinterval(as.mcmc(ortho$pol1))
