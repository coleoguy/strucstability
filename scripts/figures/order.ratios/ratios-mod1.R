library(ggplot2)
library(coda)
library(viridis)
library(ape)

tree <- read.tree (text="((((((Diptera,Lepidoptera),Coleoptera),Hymenoptera),Hemiptera),(Blattodea,Orthoptera)),Odonata);")

ggtheme <- theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.background = element_blank(),
                              panel.border=element_blank(),
                              axis.line = element_line(colour="grey30"),
                              axis.title = element_text(colour="grey20"),
                              axis.text = element_text(colour="grey30"),
                              legend.title = element_text(colour="grey20"),
                              legend.text = element_text(colour="grey30"))

dip <- read.csv("../../../results/posterior-samples/dipteramod1.csv")
coleo <- read.csv("../../../results/posterior-samples/coleoptera-mod1.csv")
blatt <- read.csv("../../../results/posterior-samples/blattodea-mod1.csv")
hemi <- read.csv("../../../results/posterior-samples/hemiptera-mod1.csv")
hymen <- read.csv("../../../results/posterior-samples/hymenopteramod1.csv")[,-1]
lep <- read.csv("../../../results/posterior-samples/lepidoptera-mod1.csv")
odo <- read.csv("../../../results/posterior-samples/odonatamod1.csv")[,-1]
ortho <- read.csv("../../../results/posterior-samples/orthoptera-mod1.csv")

dip$order <- "Diptera 50"
coleo$order <- "Coleoptera 663"
blatt$order <- "Blattodea 40"
hemi$order <- "Hemiptera 102"
hymen$order <- "Hymenoptera 301"
lep$order <- "Lepidoptera 238"
odo$order <- "Odonata 228"
ortho$order <- "Orthoptera 36"

dat <- rbind(dip, coleo, blatt, hemi, hymen, lep, odo, ortho)
dat$ratio <- dat$desc1 / (dat$asc1 + dat$desc1)
rm(list=ls()[-c(3,5)])

dat$order <- as.factor(dat$order)
dat$order <- factor(dat$order, levels=rev(c("Hemiptera 102", "Coleoptera 663","Blattodea 40",
                                        "Hymenoptera 301", 
                                        "Lepidoptera 238", "Diptera 50",
                                        "Orthoptera 36", "Odonata 228")))

ode.hpd <- ort.hpd <- bla.hpd <- hem.hpd <- hym.hpd <- 
  col.hpd <- lep.hpd <- dip.hpd <- as.data.frame(matrix(NA,2,2))
colnames(ode.hpd) <- colnames(ort.hpd) <- colnames(bla.hpd) <- 
  colnames(hem.hpd) <- colnames(hym.hpd) <- colnames(col.hpd) <- 
  colnames(lep.hpd) <- colnames(dip.hpd) <- c("X","Y")
ode.hpd$X <- as.numeric(HPDinterval(as.mcmc(dat$ratio[dat$order == "Odonata 228"])))
ort.hpd$X <- as.numeric(HPDinterval(as.mcmc(dat$ratio[dat$order == "Orthoptera 36"])))
bla.hpd$X <- as.numeric(HPDinterval(as.mcmc(dat$ratio[dat$order == "Blattodea 40"])))
hem.hpd$X <- as.numeric(HPDinterval(as.mcmc(dat$ratio[dat$order == "Hemiptera 102"])))
hym.hpd$X <- as.numeric(HPDinterval(as.mcmc(dat$ratio[dat$order == "Hymenoptera 301"])))
col.hpd$X <- as.numeric(HPDinterval(as.mcmc(dat$ratio[dat$order == "Coleoptera 663"])))
lep.hpd$X <- as.numeric(HPDinterval(as.mcmc(dat$ratio[dat$order == "Lepidoptera 238"])))
dip.hpd$X <- as.numeric(HPDinterval(as.mcmc(dat$ratio[dat$order == "Diptera 50"])))

ode.hpd$Y <- c(1,1)
ort.hpd$Y <- c(2,2)
bla.hpd$Y <- c(6,6)
hem.hpd$Y <- c(8,8)
hym.hpd$Y <- c(5,5)
col.hpd$Y <- c(7,7)
lep.hpd$Y <- c(4,4)
dip.hpd$Y <- c(3,3)

means <- as.data.frame(aggregate(dat$ratio, list(dat$order), FUN=mean))
means$y <- 1:8
#Plotting
ggplot(dat, aes(y=order, x = ratio)) +
  geom_vline(xintercept = .5, linetype = "dashed", color = "gray")+
  geom_jitter(size=.4, height=.2,color=viridis(10,alpha=.4,option="B")[6])+
  xlab("Fusion / (Fusion + Fission)")+ 
  geom_line(data=ode.hpd, aes(x=X, y=Y), linewidth=1.2, 
            lineend="round") +
  geom_line(data=ort.hpd, aes(x=X, y=Y), linewidth=1.2, 
            lineend="round") +
  geom_line(data=bla.hpd, aes(x=X, y=Y), linewidth=1.2, 
            lineend="round") +
  geom_line(data=hem.hpd, aes(x=X, y=Y), linewidth=1.2, 
            lineend="round") +
  geom_line(data=hym.hpd, aes(x=X, y=Y), linewidth=1.2, 
            lineend="round") +
  geom_line(data=col.hpd, aes(x=X, y=Y), linewidth=1.2, 
            lineend="round") +
  geom_line(data=lep.hpd, aes(x=X, y=Y), linewidth=1.2, 
            lineend="round") +
  geom_line(data=dip.hpd, aes(x=X, y=Y), linewidth=1.2, 
            lineend="round") +
  geom_point(data=means, aes(x=x, y=y), size=3) +
  ggtheme
# export 6"x6"
