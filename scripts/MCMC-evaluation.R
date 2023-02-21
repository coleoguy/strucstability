#####  Diptera #### 
dat1 <- read.csv("../results/raw-results/diptera-mod1.csv")
dat2 <- read.csv("../results/raw-results/diptera-mod2.csv")
plot(dat1$p, pch=16, cex=.2)
plot(dat2$p, pch=16, cex=.2)
write.csv(dat1[701:1000,], file="../results/posterior-samples/dipteramod1.csv", row.names = F)
write.csv(dat2[701:1000,], file="../results/posterior-samples/dipteramod2.csv", row.names = F)
rm(list=ls())
#####  Diptera #### 



#####  Hymenoptera #### 
dat1 <- read.csv("../results/raw-results/hymenoptera-mod1.csv")
dat2 <- read.csv("../results/raw-results/hymenoptera-mod2.csv")
plot(dat1$p, pch=16, cex=.2)
plot(dat2$p, pch=16, cex=.2)
abline(v=700)
write.csv(dat1[701:1000,], file="../results/posterior-samples/hymenopteramod1.csv", row.names = F)
write.csv(dat2[701:1000,], file="../results/posterior-samples/hymenopteramod2.csv", row.names = F)
rm(list=ls())
#####  Hymenoptera #### 





#####  Hemiptera #### 
dat1 <- readRDS("../results/raw-results/hemiptera-mod1.RData")
dat2 <- readRDS("../results/raw-results/hemiptera-mod2.RData")

hem1 <- dat1[[1]][sample(701:1000, 3), ]
hem2 <- dat2[[1]][sample(701:1000, 3), ]
for(i in 2:100){
  hem1 <- rbind(hem1, dat1[[i]][sample(701:1000, 3), ])
  hem2 <- rbind(hem2, dat2[[i]][sample(701:1000, 3), ])
}
write.csv(hem1,  "../results/posterior-samples/hemiptera-mod1.csv", row.names = F)
write.csv(hem2,  "../results/posterior-samples/hemiptera-mod2.csv", row.names = F)


rm(list=ls())
#####  Hemiptera #### 




#####  Odonata #### 
dat1 <- read.csv("../results/raw-results/odonata-mod1-sample1.csv")
dat2 <- read.csv("../results/raw-results/odonata-mod2.csv")
plot(dat1$p, pch=16, cex=.2)
abline(v=700)
plot(dat2$p, pch=16, cex=.2)
write.csv(dat1[701:1000,], file="../results/posterior-samples/odonatamod1.csv", row.names = F)
write.csv(dat2[701:1000,], file="../results/posterior-samples/odonatamod2.csv", row.names = F)
rm(list=ls())
#####  Odonata #### 





#####  Orthoptera #### 
dat1 <- readRDS("../results/raw-results/orthoptera-mod1.RData")
dat2 <- readRDS("../results/raw-results/orthoptera-mod2.RData")
for(i in 1:100){
  plot(dat1[[i]]$p, pch=16, type="l")
}
for(i in 1:100){
  plot(dat2[[i]]$p, pch=16, type="l")
}
orth1 <- dat1[[1]][sample(901:1000, 3), ]
for(i in 2:100){
  orth1 <- rbind(orth1, dat1[[i]][sample(701:1000, 3), ])
}
orth2 <- dat2[[1]][sample(901:1000, 3), ]
for(i in 2:100){
  orth2 <- rbind(orth2, dat2[[i]][sample(701:1000, 3), ])
}
write.csv(orth1, "../results/posterior-samples/orthoptera-mod1.csv", row.names = F)
write.csv(orth2, "../results/posterior-samples/orthoptera-mod2.csv", row.names = F)
rm(list=ls())
#####  Orthoptera #### 



#### Blattodea ####
dat1 <- readRDS("../results/raw-results/blattodea-mod1.RData")
blat1 <- dat1[[1]][sample(901:1000, 3), ]
for(i in 2:100){
  blat1 <- rbind(blat1, dat1[[i]][sample(901:1000, 3), ])
}
blat2files <- list.files("../results/raw-results/blattodea-mod2/")
blat2 <- read.csv(paste("../results/raw-results/blattodea-mod2/", blat2files[1], sep=""))
blatmod2 <- blat2[sample(701:1000, 3), ]
for(i in 2:100){
  blat2 <- read.csv(paste("../results/raw-results/blattodea-mod2/", blat2files[i], sep=""))
  blatmod2 <- rbind(blatmod2, blat2[sample(701:1000, 3), ])
}
write.csv(blat1,  "../results/posterior-samples/blattodea-mod1.csv", row.names = F)
write.csv(blatmod2[,-1],  "../results/posterior-samples/blattodea-mod2.csv", row.names = F)
rm(list=ls())
#### Blattodea ####



#####  Coleoptera #### 
dat1 <- readRDS("../results/raw-results/coleoptera-mod1.RData")
col1 <- dat1[[1]][sample(901:1000, 3), ]
for(i in 2:100){
  col1 <- rbind(col1, dat1[[i]][sample(901:1000, 3), ])
}
col2files <- list.files("../results/raw-results/coleoptera-mod2/")
col2 <- read.csv(paste("../results/raw-results/coleoptera-mod2/", col2files[1], sep=""))
colmod2 <- col2[sample(701:1000, 3), ]
for(i in 2:100){
  col2 <- read.csv(paste("../results/raw-results/coleoptera-mod2/", col2files[i], sep=""))
  colmod2 <- rbind(colmod2, col2[sample(701:1000, 3), ])
}
write.csv(col1,  "../results/posterior-samples/coleoptera-mod1.csv", row.names = F)
write.csv(colmod2,  "../results/posterior-samples/coleoptera-mod2.csv", row.names = F)
rm(list=ls())
#####  Coleoptera #### 

#####  Lepidoptera #### 
dat1 <- readRDS("../results/raw-results/lepidoptera-mod1.RData")
lep1 <- dat1[[1]][sample(701:1000, 3), ]
for(i in 2:100){
  lep1 <- rbind(lep1, dat1[[i]][sample(701:1000, 3), ])
}
lep2files <- list.files("../results/raw-results/lepidoptera-mod2/")
lep2 <- read.csv(paste("../results/raw-results/lepidoptera-mod2/", lep2files[1], sep=""))
lepmod2 <- lep2[sample(701:1000, 3), ]
for(i in 2:100){
  lep2 <- read.csv(paste("../results/raw-results/lepidoptera-mod2/", lep2files[i], sep=""))
  lepmod2 <- rbind(lepmod2, lep2[sample(701:1000, 3), ])
}
write.csv(lep1,  "../results/posterior-samples/lepidoptera-mod1.csv", row.names = F)
write.csv(lepmod2[,-1],  "../results/posterior-samples/lepidoptera-mod2.csv", row.names = F)
rm(list=ls())
#####  Lepidoptera #### 




