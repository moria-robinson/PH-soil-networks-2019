##########################################################################
# Shift in generalist & specialist species' abundance between soil types #
##########################################################################
# FIGURE
# 1) Average abundance of generalists and specialist species in non-serpentine and serpentine soils

# ANALYSIS
# 2) Model quantifying effect of soil on abundance

#############
# Packages #
#############

library(Rmisc)
library(ggplot2)
library(lme4)
library(bbmle)
library(MASS)

########
# Data #
########
PH <- read.csv("~/Github/data/2014-2015_network.csv")

#########################################################################
# Calculate abundance of each generalist species, in each soil and year #
# Remove species with unkown or unconfirmed diet breadth, & add 0's 
#########################################################################
# Get abundances for each caterpillar species, for each network (e.g. each soil X year)
PH$abu <- 1
PH2 <- aggregate(abu ~ cat.genus.species + hostbreadth + soil + year, data=PH, length)
# only use species with known diet breadth
PH3 <- subset(PH2, hostbreadth != "unk") 
PH3$hostbreadth <- factor(PH3$hostbreadth)
# Remove morphotypes that are not known to species (e.g. for which we could not confirm diet breadth in the literature)
remove <- c("Geometrid.09", "Geometrid.18", "Geometrid.19", "Noctuid.11", "Noctuid.19", "Geometrid.39", "Geometrid.40", "Geometrid.47", "Noctuid.17", "Erannis.sp")
PH4 <- PH3[ ! PH3$cat.genus.species %in% remove, ]

# Quick count of species
spec <- subset(PH4, hostbreadth=="spec"); factor(spec$cat.genus.species) # 31 specialists
gen <- subset(PH4, hostbreadth=="gen"); factor(gen$cat.genus.species) # 22 generalists

# Add 0s for species that are not found in a given web
species.levels <- as.character(factor(unique(PH4$cat.genus.species)))
year.levels <- c(2014, 2014, 2015, 2015)
soil.levels <- c("s", "ns", "s", "ns")
cat.genus.species <- sort(rep(species.levels, 4))
year <- as.integer(rep(year.levels, 53))
soil <- rep(soil.levels, 53)
merger <- data.frame(cat.genus.species, year, soil)

PH4_temp <- subset(PH4, select=c("cat.genus.species", "soil", "year", "abu"))
PH5 <- merge(PH4_temp, merger, by=c("cat.genus.species", "soil", "year"), all=TRUE); PH5$abu[is.na(PH5$abu)] <- 0

species_hostbreadth <- subset(PH4, select=c("cat.genus.species", "hostbreadth"))
temp <- species_hostbreadth[!duplicated(species_hostbreadth), ]

# OK - use these two slices. PH6 = zeros added. PH6_nozero = original data (non-zero data only)
PH6 <- merge(PH5, temp, by="cat.genus.species"); head(PH6)
PH6_nozero <- subset(PH6, abu != 0)

###########
# FIGURES #
###########

########################################################################################
# 1) Average abundance of generalists and specialists in non-serpentine and serpentine #
########################################################################################
# With zeros
summ <- summarySE(PH6, measurevar="abu", groupvars=c("hostbreadth", "soil", "year"))
summ$Year <- as.factor(summ$year)
summ$Hostbreadth <- "Generalist"
summ$Hostbreadth[which(summ$hostbreadth=="spec")] <- "Specialist"
summ$Hostbreadth <- as.factor(summ$Hostbreadth)

# Plot:
summ$linegroup <- paste(summ$Hostbreadth, summ$year, sep="_")

pos <- position_dodge(width=0.7) # set dodge amount

ggplot(summ, aes(x=soil, y=abu)) + # NOTE THAT YEAR ICONS in LEGEND ARE WRONG - OPEN ICONS = 2015; FIX MANUALLY in ILLUSTRATOR
  geom_errorbar(aes(group=linegroup, ymin=abu-se, ymax=abu+se), width=.3, size=.5, position=pos) +
  geom_line(aes(group=linegroup), size=.5, colour="grey50", lty=2, position=pos) +
  geom_point(aes(group=linegroup, shape=Hostbreadth, fill=Year), position=pos, colour="black", size=3.5, stroke=1.2) + 
  scale_fill_manual(name="Year", values=c("2015"= "White", "2014"="Black")) +
  scale_shape_manual(name="Hostbreadth", values=c("Generalist"=21, "Specialist"=24)) +
  theme_bw()+
  theme(axis.line = element_line(color = 'black')) +
  theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  xlab("Soil") +
  ylab("Average Abundance \n (per species)") +
  scale_x_discrete(labels = c("Non-\nserpentine","Serpentine")) +
  theme(axis.title.x=element_text(size=15,face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y=element_text(size=15,face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12), axis.text.y=element_text(size=12)) +
  theme(legend.title = element_text(size=12, face="bold")) 

############
# ANALYSES # 
############

##################################
# 2) Effect of soil on abundance #
##################################
str(PH6)
PH6$year <- as.factor(PH6$year)
hist(PH6$abu)
var(PH6$abu)/mean(PH6$abu) # Very overdispersed; use negative binomial model

m1 <- glmer.nb(abu ~ soil*hostbreadth + year + (1|cat.genus.species), data=PH6)
summary(m1)
coef_m1 <- data.frame(summary(m1)$coef); str(coef_m1)
coef_m1$beta_trans <- (exp(coef_m1$Estimate)-1)*100 
anova(m1) # Get F statistics

# Calculate residual DF. 
nrow(PH6) # 212 observations (each row is a species). "Total variance = N-1. [here n=212 - 1 = 211]. The model degrees of freedom corresponds to the number of coefficients estimated minus 1.  Including the intercept, there are 5 coefficients, so the model has 5-1=4 degrees of freedom. The Residual degrees of freedom is the DF total minus the DF model, 211 – 4 =207."

# test significance of species random effect. 
m1.without_species <- glm.nb(abu ~ soil*hostbreadth + year, data=PH6) # can't run glmer.nb without a random effect in the formula...
AICtab(m1, m1.without_species) # With species is much better (delta AIC = 60.5), but I would prefer not to compare models from different families.

# Second approach:
# If I drop species from the model, the resulting structure has no random effects; but, I could add a 'dummy' random effect to both models.
levels <- as.factor(c(1:20)) # 20 'species'
dummy <- sample(levels, size=nrow(PH6), replace=TRUE) # make this into a fake random effect

PH6$dummy <- dummy 
m1.with_species <- glmer.nb(abu ~ soil*hostbreadth + year + (1|cat.genus.species) + (1|dummy), data=PH6)
m1.without_species <- glmer.nb(abu ~ soil*hostbreadth + year + (1|dummy), data=PH6)
# Both are overfit with the dummy in addition to the species (boundary (singular) fit: see ?isSingular). Checked documentation; 
AICtab(m1.with_species, m1.without_species) # Actually, same difference in AIC - that is encouraging.
anova(m1.with_species, m1.without_species) # SO - adding species matters. Species is significant as a random effect (df = 1, x2 = 62.5, p < 0.001)

# Check residuals of model
plot(residuals(m1)) # look good
