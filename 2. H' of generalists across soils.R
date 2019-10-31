#############################################
# Shift in generalist H' between soil types #
#############################################
# FIGURES
# 1) Average H' of generalists in non-serpentine and serpentine 
# 2) Individual points/lines per generalist species (different shapes for years, colors for sp).
# 3) Correlation between H of each species/soil & their abundance in that soil

# ANALYSES
# 4) Model quantifying effect of soil on H, abundance as a covariate

# SIZE CONTROL
# 5) As an additional test to control for effects of network size on H, we rarify the networks to the same size & repeat Figure 1 + analysis.

#############
# Packages #
#############
library(bipartite)
library(RColorBrewer)
library(labdsv)
library(bbmle)
library(ggplot2)
library(Rmisc)
library(car)
library(lme4)
#library(lmerTest)
library(glmmADMB) # Run code chunk below to install, as this package will soon be deprecated

# To donwload glmmADMB (soon to be replaced by glmmTMB)
install.packages("R2admb")
install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos",getOption("repos")),type="source")

########
# Data #
########
PH <- read.csv("~/Github/data/2014-2015_network.csv")

##################################################################
# Calculate H for each generalist species, in each soil and year #
##################################################################
# a) Make a variable to loop over, since I want to do this separately for each soil, and also separately within each year (soil x year partitions)
PH$soil_year <- paste(PH$soil, PH$year, sep="_"); unique(PH$soil_year)

# b) Get abundances for each caterpillar species, for each plant genus X soil X year, for interaction matrices.
PH$abu <- 1
PH2 <- aggregate(abu ~ cat.genus.species + genus + soil_year, data=PH, length)

# c) Make a DF to merge with later; abundance of each species in each soil and year (to include as a covariate, and to make Figure 3)
species.abu <- aggregate(abu ~ cat.genus.species + soil_year + hostbreadth, data=PH, length)

# OK! Now, write a loop. This is what it does:
# 1) Subsets out each soil X year partition of the data
# 2) Computes the H for each species in that network

partner.div <- list() # empty list to fill with the four dataframes (each soil X year partition)

web <- unique(PH2$soil_year)
nweb <- length(unique(PH2$soil_year))

for(i in 1:nweb) {
  
  temp1 <- subset(PH2, soil_year == web[i])
  temp2 <- subset(temp1, select=c("genus", "cat.genus.species", "abu"))
  mat <- matrify(temp2)
  
  H <- specieslevel(mat, index="partner diversity", level="higher")
  
  H$species <- rownames(H)
  H$soil_year <- unique(temp1$soil_year)
  rownames(H) <- NULL

  partner.div[[i]] <- H # Save each output dataframe (species, soil, year, and H)
}

partner.div1 <- do.call("rbind", partner.div)

# Merge back in per-species info (abundance, hostbreadth)
species.abu$species <- species.abu$cat.genus.species
species.abu$cat.genus.species <- NULL

partner.div2 <- merge(partner.div1, species.abu, by=c("species", "soil_year"))

# Split soil_year in to two columns
partner.div2$soil <- factor(gsub("_.+$", "", partner.div2$soil_year)) # also makes it a factor
partner.div2$year <- factor(gsub("^.+_", "", partner.div2$soil_year)) # also makes it a factor
str(partner.div2)

###########
# FIGURES #
###########
# Subset to just generalists
gen.H <- subset(partner.div2, hostbreadth == "gen"); length(unique(gen.H$species)) 

# Check that all have species-level ID
unique(gen.H$species) # might need to remove Noct 19, if it isn't excluded in the next step

# To avoid spurious designations of specialization due to rarity, subset to just species with > 3 individuals per soil X year
gen.H.min3 <- gen.H[which(gen.H$abu > 3),]

# Check again that all have species-level ID
unique(gen.H.min3$species)

# Get sample sizes per soil & year
nrow(subset(gen.H.min3, soil == "s" & year == 2014)) # n = 4
nrow(subset(gen.H.min3, soil == "s" & year == 2015)) # n = 8

nrow(subset(gen.H.min3, soil == "ns" & year == 2014)) # n = 8
nrow(subset(gen.H.min3, soil == "ns" & year == 2015)) # n = 13

#################################################################
# 1) Average H' of generalists in non-serpentine and serpentine #
#################################################################
summ <- summarySE(gen.H.min3, measurevar="partner.diversity", groupvars=c("soil", "year"))
summ$Year <- summ$year # for legend title purposes

# Plot:
pos <- position_dodge(width=0.7) # set dodge amount

ggplot(summ, aes(x=soil, y=partner.diversity)) +
  geom_errorbar(aes(ymin=partner.diversity-se, ymax=partner.diversity+se, linetype=NULL, group=as.factor(Year)), width=.3, size=.7, position=pos) +
  geom_point(aes(group=Year, shape=Year), position=pos, size=3.5, fill="white", colour="black", stroke=1) + 
  scale_shape_manual(name="Year", values=c("2015"=21, "2014"=19)) +
  theme_bw()+
  theme(axis.line = element_line(color = 'black')) +
  theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  xlab("Soil") +
  ylab("Generalist partner diversity (H)") +
  scale_x_discrete(labels = c("Non-\nserpentine","Serpentine")) +
  theme(axis.title.x=element_text(size=15,face="bold", vjust=-.5), axis.title.y=element_text(size=15,face="bold", vjust=1.5)) +
  theme(axis.text.x = element_text(size = 12), axis.text.y=element_text(size=12)) +
  scale_y_continuous(limits=c(0,1)) +
  theme(legend.title = element_text(size=12, face="bold")) 

##########################################
# 2) Individual points/lines per species #
##########################################
pos <- position_dodge(width=0.2) # set dodge amount

gen.H.min3$species2 <- gsub("\\.", " ", gen.H.min3$species) # Make a new species name with a space rather than a period
gen.H.min3$linegroup <- paste(gen.H.min3$species2, gen.H.min3$year, sep="_")

# Generate a unique color for each species
colourCount = length(unique(gen.H.min3$species2)) # 13 unique species
all.colors <- colorRampPalette(brewer.pal(8,"Dark2"))(colourCount)

# Name the colors by the Species names
names(all.colors) <- levels(gen.H.min3$species2)

ggplot(gen.H.min3, aes(x=soil, y=partner.diversity)) +
  theme_bw()+
  geom_line(aes(group=linegroup), size=1, colour="grey50", position=pos) +
  geom_point(aes(group=linegroup, colour=species2, shape=year), size=6, position=pos) +
  scale_colour_manual(name = "Species", values = all.colors) +
  scale_shape_manual(name="Year", values=c("2015"=17,"2014"=16)) +
  theme_classic() +
  theme(axis.line.x=element_line(),axis.line.y=element_line()) +
  theme(legend.title = element_text(colour="grey23", size=10, face="bold")) +
  theme(legend.key = element_blank()) +
  xlab("Soil") +
  ylab("Partner diversity (H)") +
  scale_x_discrete(labels = c("Non-\nserpentine","Serpentine")) +
  theme(axis.title.x=element_text(size=19,face="bold", vjust=-.5), axis.title.y=element_text(size=19,face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 16), axis.text.y=element_text(size=16)) +
  theme(legend.title = element_text(size=16, face="bold")) 

##################################
# 3) Regression of H & abundance #
##################################
ggplot(gen.H.min3, aes(x=abu, y=partner.diversity)) +
  geom_point(aes(colour=species), alpha=0.8, size=3) + 
  scale_colour_manual(name = "Species", values = all.colors) +
  geom_smooth(method="lm") +
  theme_bw()+
  theme(axis.line = element_line(color = 'black')) +
  theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  xlab("Abundance") +
  ylab("Generalist partner diversity (H)") +
  theme(axis.title.x=element_text(size=15,face="bold", vjust=-.5), axis.title.y=element_text(size=15,face="bold", vjust=1.5)) +
  theme(axis.text.x = element_text(size = 12), axis.text.y=element_text(size=12)) +
  theme(legend.title = element_text(size=12, face="bold")) 

###############
# 4) ANALYSES #
###############
hist(gen.H.min3$partner.diversity) # not bad!
str(gen.H.min3)
gen.H.min3$species <- factor(gen.H.min3$species) # Make species a factor for the model
gen.H.min3$abu <- as.numeric(as.character(gen.H.min3$abu)) # Make abundance numeric

m1 <- glmmadmb(partner.diversity ~ soil + abu + year + (1|species), data=gen.H.min3, family="gaussian")
summary(m1) 
Anova(m1, type=2, test="F") 

# What's the % difference in partner diversity, caused by soil type? After accounting for effects of covariates, the group mean H for non-serpentine is 0.6675 (intercept), and for serpentine 0.6675-0.2190 = 0.4485. So, .2190/.6675=0.328 - 33%. Serpentine reduces partner diversity by 33%, on average. 

# What's the effect of herbivore abundance on partner diversity? After accounting for effects of covariates, the beta shows that Shannon increases by 0.0120 for every one-unit increase in herbivore abundance. Generalist herbivores were, on average, how different in abundance?
abu.diff <- summarySE(gen.H.min3, measurevar="abu", groupvar="soil") # OK - 13.38 in non-serp, 9.91 in serp. A difference of 3.47 individuals. So - 3.47*.0120 = 0.042. 

# test significance of species random effect
levels <- as.factor(c(1:20)) # 20 'species'
dummy <- sample(levels, size=nrow(gen.H.min3), replace=TRUE) # make this into a fake random effect
gen.H.min3$dummy <- as.factor(dummy)

m2 <- glmmadmb(partner.diversity ~ soil + abu + year + (1|species) + (1|dummy), data=gen.H.min3, family="gaussian")
m2.1 <- glmmadmb(partner.diversity ~ soil + abu + year + (1|dummy), data=gen.H.min3, family="gaussian")

AICtab(m2, m2.1)
anova(m2, m2.1) # No difference; no effect of species.

########################################################################
# 5) Size control. For each web, subsample the same # of interactions: #
########################################################################
# Let's subsample each web for 200 interactions.
web <- unique(PH$soil_year)
nweb <- length(unique(PH$soil_year))
all.H.allwebs <- list()

set.seed(1)

for(i in 1:nweb) {
  #  temp1 <- subset(PH, soil_year == "s_2015"); nrow(temp1)
  temp1 <- subset(PH, soil_year == web[i])
  all.H <- list() # empty list to fill
  
  for(j in 1:5000) { # repeat the follow subsampling & H calculation 500x...
    subweb <- temp1[sample(nrow(temp1), 200),]; nrow(subweb) # take random subset of x interactions
    # gather a list of species to remove due to low abundances by chance
    abu.test <- aggregate(abu ~ cat.genus.species, data=subweb, length)
    keep <- subset(abu.test, abu > 3) # can play with this threshold; maybe increase it?
    # This will only keep species with total abundance > 3
    subweb2 <- subweb[subweb$cat.genus.species %in% keep$cat.genus.species,]; nrow(subweb2)
    subweb3 <- aggregate(abu ~ cat.genus.species + genus, data=subweb2, length)
    # count up abundances in the subweb to put into the output
    abu.sp <- aggregate(abu ~ cat.genus.species, data=subweb3, length); nrow(abu.sp)
    abu.sp$species <- abu.sp$cat.genus.species; abu.sp$cat.genus.species <- NULL
    subweb4 <- subset(subweb3, select=c("genus", "cat.genus.species", "abu"))
    subweb5 <- matrify(subweb4)
    H <- specieslevel(subweb5, index="partner diversity", level="higher")
    H$species <- rownames(H)
    H$soil_year <- unique(subweb$soil_year)
    H$subweb <- j
    rownames(H) <- NULL
    H2 <- merge(H, abu.sp, by="species")
    all.H[[j]] <- H2
  }
  all.H.allwebs1 <- data.frame(do.call(rbind, all.H))
  all.H.allwebs[[i]] <- all.H.allwebs1
  # Save each output dataframe (species, soil, year, and H)
}

out <- data.frame(do.call(rbind,all.H.allwebs)); nrow(out)
out_avg <- summarySE(measurevar="partner.diversity", groupvars=c("species", "soil_year"), data=out); nrow(out_avg)
abu.avg <- summarySE(measurevar="abu", groupvars=c("species", "soil_year"), data=out); nrow(abu.avg)
abu.avg2 <- subset(abu.avg, select=c("species", "soil_year", "abu"))

H_abu <- merge(out_avg, abu.avg2, by=c("species", "soil_year")); nrow(H_abu); head(H_abu)

# Merge in per-species info (hostbreadth)
# hostbreadth info to merge in later
info <- aggregate(abu ~ cat.genus.species + hostbreadth, data=PH, length); info$abu <- NULL; info
info$species <- info$cat.genus.species; info$cat.genus.species <- NULL
out_avg2 <- merge(H_abu, info, by="species"); nrow(out_avg); nrow(out_avg2)
# Split soil_year in to two columns
out_avg2$soil <- factor(gsub("_.+$", "", out_avg2$soil_year)) # also makes it a factor
out_avg2$year <- factor(gsub("^.+_", "", out_avg2$soil_year)) # also makes it a factor
str(out_avg2)
# Subset to just generalists
gen.H <- subset(out_avg2, hostbreadth == "gen"); nrow(gen.H) 

summ <- summarySE(gen.H, measurevar="partner.diversity", groupvars=c("soil", "year"))
summ$Year <- as.factor(summ$year)

# Plot:
pos <- position_dodge(width=0.5) # set dodge amount
ggplot(summ, aes(x=soil, y=partner.diversity, group=Year)) +
  geom_errorbar(aes(ymin=partner.diversity-se, ymax=partner.diversity+se, linetype=NULL, group=as.factor(Year)), width=.3, size=.5, position=pos) +
  geom_line(position=pos) +
  geom_point(aes(group=Year, shape=Year), position=pos, size=3.5, stroke=1.2) + 
  scale_shape_manual(name="Year", values=c("2015"=19, "2014"=17)) +
  theme_bw()+
  theme(axis.line = element_line(color = 'black')) +
  theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  xlab("Soil") +
  ylab("Partner diversity (H)") +
  scale_x_discrete(labels = c("Non-\nserpentine","Serpentine")) +
  theme(axis.title.x=element_text(size=15,face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.title.y=element_text(size=15,face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12), axis.text.y=element_text(size=12)) +
  scale_y_continuous(limits=c(0,1)) +
  theme(legend.position="none") +
  theme(legend.title = element_text(size=12, face="bold"))

############
# Analysis #
############
hist(gen.H$partner.diversity) 
str(gen.H)
gen.H$species <- as.factor(gen.H$species)

m1 <- lmer(partner.diversity ~ soil + year + abu + (1|species), data=gen.H)
summary(m1); hist(resid(m1)); plot(resid(m1))
Anova(m1, type=2, test="F") # Even after subsampling networks to the same size, soil and abundance are still associated with for H' of generalists

