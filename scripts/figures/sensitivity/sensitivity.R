library(viridis)

coleo.full <- read.csv("../../../results/posterior-samples/coleoptera-mod1.csv")
# drop 200 for burnin
coleo.red <- read.csv("../../../results/sensitivity/coleoptera.sensitivity.csv")[201:500,]
lep.full <- read.csv("../../../results/posterior-samples/lepidoptera-mod1.csv")
# drop 200 for burnin
lep.red <- read.csv("../../../results/sensitivity/lepidoptera.sensitivity.csv")[201:500,]

cols <- viridis(10, option="B", alpha=.5)[c(3,8)]
plot(density(coleo.red$asc1), xlim=c(0,.1), xlab="rate (MY)", main="")
polygon(density(coleo.red$asc1), col=cols[1])
lines(density(coleo.full$asc1))
polygon(density(coleo.full$asc1), col=cols[2])
points(0,78, pch=16, col=cols[1], cex=2)
text(0,77.6, "Reduced data", pos=4, cex=.7)
points(0,74, pch=16, col=cols[2], cex=2)
text(0,73.6, "Full data", pos=4, cex=.7)

plot(density(lep.full$asc1), xlim=c(0,21), xlab="rate (MY)", main="", ylim=c(0,1.5))
polygon(density(lep.full$asc1), col=cols[2])
lines(density(lep.red$asc1))
polygon(density(lep.red$asc1), col=cols[1])
points(0,1.5, pch=16, col=cols[1], cex=2)
text(0,1.5, "Reduced data", pos=4, cex=.7)
points(0,1.4, pch=16, col=cols[2], cex=2)
text(0,1.4, "Full data", pos=4, cex=.7)
