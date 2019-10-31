# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Shifts in generalist interaction strength + plant resistance index #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# 1) Calculate a Resistance Index for each plant species in each soil type
# 2) Calculate shifts in interaction strength by generalists in each soil type
# 3) Compare

#############
# Packages #
#############
library(ggplot2)

########
# Data #
########
all.traits <- read.csv("~/Github/data/plant.traits_RI.csv")
PH <- read.csv("~/Github/data/2014-2015_network.csv")

########################################
# 1) Calculate Fine's resistance index #
########################################
# The steps (from Fine 2006) would be...
# 1) Average each trait value across soils
# 2) For 'good' traits (ie nitrogen, water content), take the inverse of the species averages, because a larger amount = lower resistance
# 3) Z-transform to give the traits among the 4 pairs of plants a mean of 0 and SD of 1
# 4) Sum all standardized resistance variables to create a single value, a "resistance index", for each plant

# Fine (2006)'s analysis was...
# 5) For each genus, subtract the RI of the less-harsh-soil plant from the RI of the harsh-soil plant (RI[s] - RI[ns])

# 6) Fine used these "difference scores" (RI[s] - RI[ns]) to ask if S are more resistant than NS. I will compare the rank order of those difference scores (4 values) with the rank order of shifts in generalist interactions between the same serp / non-serp plants.

# 1) Average of each trait, per plant and soil, is in this DF:
all.traits

# 2) Inverse the "good" traits - N, LWC
all.traits$nitrogen.inverse <- all.traits$nitrogen*(-1)
all.traits$lwc.inverse <- all.traits$lwc*(-1)

# 3) Z-transform each trait
all.traits$nitrogen.ztrans <- scale(all.traits$nitrogen.inverse)
all.traits$carbon.ztrans <- scale(all.traits$carbon)
all.traits$lwc.ztrans <- scale(all.traits$lwc.inverse)
all.traits$nickel.ztrans <- scale(all.traits$nickel)
all.traits$tough_young.ztrans <- scale(all.traits$tough_young)
all.traits$tough_old.ztrans <- scale(all.traits$tough_old)

# 4) Sum them for each plant to get a "resistance index"
all.traits$RI <- all.traits$nitrogen.ztrans + all.traits$carbon.ztrans + all.traits$lwc.ztrans + all.traits$nickel.ztrans + all.traits$tough_young.ztrans + all.traits$tough_old.ztrans

RI_long <- subset(all.traits, select=c("genus", "soil", "RI"))
RI_wide <- spread(RI_long, soil, RI)
names(RI_wide) <- c("Genus", "RI_NS", "RI_S")

RI_long$abbrev <- ""
RI_long$abbrev[which(RI_long$genus=="arctostaphylos" & RI_long$soil== "s")] <- "ARCT"
RI_long$abbrev[which(RI_long$genus=="quercus" & RI_long$soil== "s")] <- "QUER"
RI_long$abbrev[which(RI_long$genus=="ceanothus" & RI_long$soil== "s")] <- "CEAN"
RI_long$abbrev[which(RI_long$genus=="adenostoma" & RI_long$soil== "s")] <- "ADEN"

# Make a simple figure
ggplot(data=RI_long, aes(x=soil, y=RI)) +
  geom_line(aes(group=genus), size=.5, colour="grey50", lty=2) +
  geom_point(aes(group=genus), size=5) +
  geom_text(aes(label = abbrev), hjust = 0, nudge_x = 0.07) +
  theme_bw()+
  theme(axis.line = element_line(color = 'black')) +
  theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  xlab(NULL) +
  ylab("Resistance Index (RI)") +
  scale_x_discrete(labels = c("Non-serpentine \n Soil","Serpentine \n Soil")) +
  scale_y_continuous(breaks = c(-8, -4, 0, 4, 8), limits=c(-8, 8)) +
  theme(axis.title.x=element_text(size=15,face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y=element_text(size=15,face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0))) + theme(axis.text.x = element_text(size = 12), axis.text.y=element_text(size=12)) + theme(legend.title = element_text(size=12, face="bold")) 

# Calculate "difference scores" (RI[s] - RI[ns])
RI_wide$difference.score <- RI_wide$RI_S - RI_wide$RI_NS
RI_wide

############################################################
# 2) Shift in generalist interaction strength across soils #
############################################################
# Subset to known generalists
PH1 <- subset(PH, hostbreadth == "gen")
# only use species with known diet breadth
PH2 <- subset(PH1, hostbreadth != "unk"); PH2$cat.genus.species <- factor(PH2$cat.genus.species) # re-factor to remove hidden levels
unique(PH2$cat.genus.species)
# Remove morphotypes that are not known to species (e.g. for which we could not confirm diet breadth in the literature)
remove <- c("Erannis.sp", "Noctuid.19")
PH3 <- PH2[ ! PH2$cat.genus.species %in% remove, ]
# Check if it worked...
PH3$cat.genus.species <- factor(PH3$cat.genus.species) # re-factor to remove hidden levels
unique(PH3$cat.genus.species)

# Calculate OVERALL interaction strength (generalists only) per host plant genus
PH3$abu <- 1
per.genus <- aggregate(abu ~ soil + year + genus, data=PH3, length)
# Average across years:
per.genus2 <- summarySE(per.genus, measurevar="abu", groupvars=c("soil", "genus"))

# Plot
ggplot(data=per.genus2, aes(x=soil, y=abu)) +
  geom_errorbar(aes(group=genus, ymin=abu-se, ymax=abu+se), width=.3, size=.5, position=position_dodge(width=0.4)) +
  geom_line(aes(group=genus), size=.5, colour="grey50", lty=2, position=position_dodge(width=0.4)) +
  geom_point(aes(group=genus), size=5, position=position_dodge(width=0.4)) + # add colour=genus to visualize which is which
  theme_bw()+
  theme(axis.line = element_line(color = 'black')) +
  theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  xlab(NULL) +
  ylab("Total Interaction Strength \n (fundamental generalists)") +
  scale_x_discrete(labels = c("Non-serpentine \n Soil","Serpentine \n Soil")) +
  theme(axis.title.x=element_text(size=15,face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y=element_text(size=15,face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0))) + theme(axis.text.x = element_text(size = 12), axis.text.y=element_text(size=12)) + theme(legend.title = element_text(size=12, face="bold")) 

# Interaction strength = abu column
per.genus2
