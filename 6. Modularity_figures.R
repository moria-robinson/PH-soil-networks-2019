################################################
# plots from subsampling / resampling networks #
################################################
# Prior script (modularity_empirical.rarefied.R) looped through a resampling pipeline to 1) calculate modularity for many rarefied networks of different sizes, and 2) calculate the maximum possible modularity for each network size (to divide by for Qnorm). 
# Here, we build plots of that subsampling, & ask if/when the distributions are different.

#############
# Packages #
#############
library(ggplot2)

########
# Data #
########
PH <- read.csv("~/Github/data/2014-2015_network.csv") # just for soil/year info 

# Description of methods for each df:
subsamp.q <- read.csv("~/Github/data/Q_from.subsampling.csv") # Subsampled networks & calculated Q. We subsampled networks as follows: for each network, we created intervals of network size ranging from 25 interactions to the total number of interactions in the network. At each interval, we randomly subsampled that number of interactions from each network and calculated modularity of the resulting network (=mod.max). Modularity was calculated as described in main text methods (function “DIRTLPAwb+”, Beckett [2016]). We iterated this for 100 random networks at each size interval.

subsamp.q_max <- read.csv("~/Github/data/NormalizedQ_from.subsampling.csv") # For this, we did a similar procedure as above, but for each size interval (of 5) we made only 10 random subwebs (="rep.q"). For each of those subwebs, we then "modularizated" the subweb - that is, aggregated the total # of interactions per herbivore species in the web, and then randomly assigned each one (and all of its interactions) to one of the 4 plant genera. This keeps the total number of interactions/herbivore species the same, but makes all of the species specialists (removes all edges between modules). Then modularity was calculated for this modularized version. We made 15 "modularized" versions of each random subweb, and took the maximum modularity value across those 15 webs (= max.q) to represent the maxumum modularity possible for that random subweb. Therefore, the theoretical maximum modularity for each size interval is the maximum of the 10 max.q values, per size.

emp.q <- read.csv("~/Github/data/emp_q&q_max.csv") # These are the empirical Q values for each network, and the normalized Q values. In this case, q.max was gotten from 100 different "modularized" versions of the web (where "modularization" again means randomly assigning each species & all its interactions to a single host plant)

##########################################
# Clean up data, clarify variable names #
##########################################
# soil/year info to merge
PH$abu <- 1
PH$soil_year <- paste(PH$soil, PH$year, sep="_"); unique(PH$soil_year)
info <- aggregate(abu ~ soil_year + soil + year, data=PH, length); info$abu <- NULL
info
info$web <- info$soil_year

# Q - subsampling
head(subsamp.q)
subsamp.q$X <- NULL
names(subsamp.q) <- c("mod.q", "rep","size","web") # change mod.max to mod.q, so as not to confuse with the theoretical maximum mod in the other df
subsamp.q$soil_year <- subsamp.q$web
subsamp.q <- merge(subsamp.q, info, by="soil_year"); head(subsamp.q)
str(subsamp.q)
subsamp.q$year <- as.factor(subsamp.q$year)

# Single empirical Qs (n=4)
emp.q$year <- as.factor(emp.q$year)

# First, pull out single highest max.q per size & network:
head(subsamp.q_max)
subsamp.q_max$cat <- paste(subsamp.q_max$web, subsamp.q_max$size)
# function to subset the DF to the maximum q.max per network & size category
max <- do.call(rbind,by(subsamp.q_max,subsamp.q_max$cat,function(dat)dat[order(dat$max.q,decreasing=TRUE)[1:1],])); nrow(max) # length 423, which is one for each size category, for each network
head(max)
max2 <- subset(max, select=c("soil_year", "max.q", "size"))

# For each size, merge in the theoretical maximum Q
subsamp.q_norm <- merge(subsamp.q, max2, by=c("soil_year", "size")); nrow(subsamp.q); nrow(subsamp.q_norm); head(subsamp.q_norm)
# Calculate normalized Q
subsamp.q_norm$norm.q <- subsamp.q_norm$mod.q/subsamp.q_norm$max.q

head(subsamp.q_norm); str(subsamp.q_norm)
head(emp.q); str(emp.q)
emp.q$norm.q <- emp.q$q.norm # slightly different name; rename norm.q in the empirical/point est. df:

###############################################
# PLOT 1 : raw modularity Q across subsamples # THIS FIGURE NOT IN MS
###############################################
head(subsamp.q); str(subsamp.q)
head(emp.q); str(emp.q)

ggplot(subsamp.q, aes(x=size, y=mod.q)) +
  geom_point(size=.2, alpha=0.03, aes(colour=soil)) +
  geom_smooth(method="loess", size=.5, aes(colour=soil), se=TRUE, level=0.95) +
  geom_point(data=emp.q, size=2, aes(colour=soil)) +
  facet_grid(~year) +
  scale_colour_manual(values=c("s"="#E69F00", "ns"="#0072B2")) +
  scale_fill_manual(values=c("2014"="#E69F00", "2015"="#0072B2")) +
  theme_bw()+
  theme(panel.border = element_rect(colour="black"), panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), strip.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Network Size (# interactions)") +
  ylab("Modularity (Q)") +
  theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=12), strip.text.x=element_text(size=15, face="bold"), axis.title.y=element_text(size=15,face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x=element_text(size=15,face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  # theme(legend.title = element_text(size=12, face="bold")) +
  theme(legend.position="none") 
# Best export dimensions - Export PDF, 4 x 5 

##########################################################################################
# PLOT: Normalized modularity Q across subsamples + test for difference in distributions # THIS FIGURE IN MS
##########################################################################################
ggplot(subsamp.q_norm, aes(x=size, y=norm.q)) +
  geom_point(size=.2, alpha=0.03, aes(colour=soil)) +
  geom_smooth(method="loess", size=.5, aes(colour=soil), se=TRUE, level=0.95) +
  geom_point(data=emp.q, size=2, aes(colour=soil)) +
  geom_segment(data = subset(subsamp.q_norm, year == 2014), aes(x=140, y=.5, xend=140, yend=1), linetype=3, size=.5, alpha=0.5) +
  geom_segment(data = subset(subsamp.q_norm, year == 2015), aes(x=140, y=.5, xend=140, yend=1), linetype=3, size=.5, alpha=0.5) +
  facet_grid(~year) +
  scale_colour_manual(values=c("s"="#E69F00", "ns"="#0072B2")) +
  scale_fill_manual(values=c("2014"="#E69F00", "2015"="#0072B2")) +
  theme_bw()+
  theme(panel.border = element_rect(colour="black"), panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), strip.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Network Size (# interactions)") +
  ylab("Normalized modularity\n(Q/Qmax)") +
  theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=12), strip.text.x=element_text(size=15, face="bold"), axis.title.y=element_text(size=15,face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x=element_text(size=15,face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  # theme(legend.title = element_text(size=12, face="bold")) +
  theme(legend.position="none") +
scale_y_continuous(limits=c(0.5,1.0), breaks = c(.5, .6, .7, .8, .9, 1.0)) +
  scale_x_continuous(limits=c(0,1250), breaks = c(0, 250, 500, 750, 1000, 1250))
# Best export dimensions - Export PDF, 4 x 5 

# Test for differences @ thresold of 1/2 the smallest network size; = 1/2 of 286 = 143; so @ 140. 
webs_2014 <- subset(subsamp.q_norm, year == 2014 & size == 140); nrow(webs_2014)
t.test(norm.q ~ soil, data=webs_2014, paired=FALSE)
# t = -20.21, df=176.14, p < 0.001
webs_2015 <- subset(subsamp.q_norm, year == 2015 & size == 140); nrow(webs_2015)
t.test(norm.q ~ soil, data=webs_2015, paired=FALSE)
# t = -4.92, df=177.53, p < 0.001
