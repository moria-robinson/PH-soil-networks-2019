############################################
# Plant-herbivore modularity analyses - II #
############################################
# Address two issues. 

# 1) Sampling completeness. Build Interaction Accumulation Curves (IAC) for each network, in each year. Jordano (2016) has some good R package suggestions for this. 

#############
# Packages #
#############

library(data.table)
library(flextable)
library(vegan)

########
# Data #
########
PH <- read.csv("~/Github/data/2014-2015_network.csv")

#############################
# Prepare data for analyses #
#############################

# a) Make a variable to loop over, since I want to do this separately for each soil, and also separately within each year (soil x year partitions)
PH$soil_year <- paste(PH$soil, PH$year, sep="_"); unique(PH$soil_year)

# b) Get abundances for each caterpillar species, for each plant genus X soil X year. I will make matrices from these abundances.
PH$abu <- 1
PH2 <- aggregate(abu ~ cat.genus.species + genus + soil_year, data=PH, length)

#######################################
# Sampling completeness statistics.
# Modelled after Table 4 in Jordano (2016). I want one column per network (e.g. each soil & year).
#######################################
PH2$interaction <- paste(PH2$cat.genus.species, PH2$genus, sep="_")

web <- unique(PH2$soil_year)
nweb <- length(unique(PH2$soil_year))

#-------------------------------------#
# A = observed # of herbivore species #
#-------------------------------------#
web.id <- NA # empty vector to fill with the soil_year partition
A <- NA # empty vector to fill with counts of species
web <- unique(PH2$soil_year)
nweb <- length(unique(PH2$soil_year))

for(i in 1:nweb) {
  temp1 <- subset(PH2, soil_year == web[i])
  A[i] <- nrow(data.frame(unique(temp1$cat.genus.species)))
  web.id[i] <- web[i]
}

A <- data.frame(web.id, A) # count of unique herbivore species per web

#---------------------------------#
# P = observed # of plant species #
#---------------------------------#
web.id <- NA # empty vector to fill with the soil_year partition
P <- NA # empty vector to fill with counts of species
web <- unique(PH2$soil_year)
nweb <- length(unique(PH2$soil_year))

for(i in 1:nweb) {
  temp1 <- subset(PH2, soil_year == web[i])
  P[i] <- nrow(data.frame(unique(temp1$genus)))
  web.id[i] <- web[i]
}

P <- data.frame(web.id, P) # count of unique plant species per web

#----------------------------------#
# I(max) = Maximum potential links #
#----------------------------------#
# This is simply multiplication of A and P.
API.max <- merge(A, P, by="web.id"); API.max
API.max$I.max <- API.max$A*API.max$P; API.max

#--------------------------------------#
# N = Number of records / observations #
#--------------------------------------#
web.id <- NA # empty vector to fill with the soil_year partition
N <- NA # empty vector to fill with records (use the PH dataframe for this; each row is a record)
web <- unique(PH$soil_year)
nweb <- length(unique(PH$soil_year))

for(i in 1:nweb) {
  temp1 <- subset(PH, soil_year == web[i])
  N[i] <- nrow(temp1)
  web.id[i] <- web[i]
}

N <- data.frame(web.id, N) # count of # records per web

#-------------------------------------------------------#
# I = Observed number of unique pairwise  interactions  #
#-------------------------------------------------------#
web.id <- NA # empty vector to fill with the soil_year partition
I <- NA # empty vector to fill with records (use the PH dataframe for this; each row is a record)
web <- unique(PH2$soil_year)
nweb <- length(unique(PH2$soil_year))

for(i in 1:nweb) {
  temp1 <- subset(PH2, soil_year == web[i])
  I[i] <- nrow(data.frame(unique(temp1$interaction)))
  web.id[i] <- web[i]
}

I <- data.frame(web.id, I) # count of # links per web

#-------------------------------------------------------------------#
# **** Merge all these descriptors of the raw interaction data ***  #
#-------------------------------------------------------------------#
API.max.N <- merge(API.max, N, by="web.id")
API.max.N.I <-  merge(API.max.N, I, by="web.id")

raw.descript <- API.max.N.I
raw.descript

# transpose
raw.descript2 <- transpose(raw.descript)
# get row and colnames in order
colnames(raw.descript2) <- rownames(raw.descript)
raw.descript2$index <- colnames(raw.descript)

table <- autofit(theme_vanilla(regulartable(raw.descript2))); table

#------------------------------------------------------------------------------#
# Chao1 - asymptotic estimator for number of unique pairwise interactions (I)
#------------------------------------------------------------------------------#
# For each network, create a matrix with all possible interactions as the first column, and then the counts of those interactions in the abu column. I'll have this separately for each network (n=4), and calculate the estimated number of pairwise interactions for each network.
head(PH2)
PH2$interaction <- paste(PH2$cat.genus.species, PH2$genus, sep="_")
temp1 <- dcast(PH2, soil_year ~ interaction, value.var="abu"); head(temp1)
temp1[is.na(temp1)] <- 0
matrix <- temp1[, 2:ncol(temp1)] 
specpool(matrix) # across all soils & years
estimateR(matrix) # by soil/year network

