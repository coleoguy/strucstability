# Michelle Jonika
# December 20th
# This code makes a plot of the order data for each model - with and without 
# polploidy

library(ggplot2)
library(ggdist)
library(viridis)
ggtheme <- theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.background = element_blank(),
                              panel.border=element_blank(),
                              axis.line = element_line(colour="grey30"),
                              axis.title = element_text(colour="grey20"),
                              axis.text = element_text(colour="grey30"),
                              legend.title = element_text(colour="grey20"),
                              legend.text = element_text(colour="grey30"))
#read in the data needed
wpoly <- read.csv("../../../results/posterior-samples/final.orderdata.wpoly.csv")
wopoly <- read.csv("../../../results/posterior-samples/final.orderdata.wopoly.csv")

# reorder so taxa are in right order
wpoly$order <- factor(wpoly$order, levels=rev(c("odonata", "orthoptera", "coleoptera", 
                                                "hymenoptera", 
                                                "blattodea", 
                                                "hemiptera", "lepidoptera", "diptera")))
wopoly$order <- factor(wopoly$order, levels=rev(c("odonata", "hemiptera", "coleoptera", "orthoptera", 
                                                  "blattodea", 
                                                  "hymenoptera", "diptera", 
                                                  "lepidoptera")))

# we used this code to just get a legend we dont actually use
# or care about the data plotted here
cols <- viridis(10,alpha=.1,option="B")[c(4,6)]
ggplot(wopoly, aes(x=log(rate), y=order, color=factor(type))) +
  geom_jitter(position= position_jitterdodge(jitter.height = .5,
                                             dodge.width = 0.75,
                                             seed = NA)) +
  scale_color_manual(values = c("fission" = cols[1], "fusion" = cols[2])) +
  stat_summary(aes(y = order, x = log(rate), fill=factor(type)),
               fun.data = "mean_hdci", 
               fun.args = list(mult=1),
               size = .4,
               position = position_jitterdodge(0),
               inherit.aes = FALSE)+
  ggtheme


cols <- viridis(10,alpha=.1,option="B")[c(4,6,8)]
ggplot(wpoly, aes(x=log(rate), y=order, color=factor(type))) +
  geom_jitter(position= position_jitterdodge(jitter.height = .4,
                                             dodge.width = 0.65,
                                             seed = NA)) +
  scale_color_manual(values = c("fission" = cols[1], "fusion" = cols[2], "polyploidy" = cols[3])) +
  stat_summary(aes(y = order, x = log(rate), fill=factor(type)),
               fun.data = "mean_hdci", 
               fun.args = list(mult=1),
               size = .4,
               position = position_jitterdodge(0),
               inherit.aes = FALSE)+
  ggtheme


# exported both at 6x6 and then in illustrtor 
# change cirlce to match correct color
