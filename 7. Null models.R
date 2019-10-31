#######################################
# Plant-herbivore modularity analyses #
#######################################
# Question: What is rthe relative importance of species' abundances and interaction structure in shaping network structure. 
# 1) Calculate expectations for modularity under 3 different expectations:
#       Null Model I) - Patefield = "r2dtable" (marginal totals as observed) 
#       Null Model II) - Shuffle = "shuffle.web" (relocates entries; connectance as observed)
#       Null Model III) - Swap = "swap.web" (marginal totals AND connectance as observed) 
# 2) Make a figure
# 3) Calculate z-scores

#############
# Packages #
#############
library(bipartite)
library(matrixStats)
library(reshape2)
library(Rmisc)
library(ggplot2)

########
# Data #
########

PH <- read.csv("~/Github/data/2014-2015_network.csv") # just for soil/year info 
PH <- read.csv("~/Google Drive/Manuscripts/Chapter 2 - Network structure and specialization/Data/2014-2015_network.csv")

#############################
# Prepare data for analyses #
#############################
# a) Make a variable to loop over, since I want to do this separately for each soil, and also separately within each year (soil x year partitions)
PH$soil_year <- paste(PH$soil, PH$year, sep="_"); unique(PH$soil_year)

# b) Get abundances for each caterpillar species, for each plant genus X soil X year. I will make matrices from these abundances.
PH$abu <- 1
PH2 <- aggregate(abu ~ cat.genus.species + genus + soil_year, data=PH, length)

######################################################################################
# 1) Calculate modularity (EMPIRICAL) for the serpentine and non-serpentine networks #              
######################################################################################
# OK! Now, write a loop. This is what it does:
# 1) Subsets out each soil X year partition of the data
# 2) Computes the modularity for that partition many times (NREPS)
# 3) Takes the maximum value from those 15 runs (this is the best modularity estimate)
# 4) Reports this single maximum value from those 15 independent runs 

web.id <- NA # empty vector to fill with the soil_year partition
modularity.est <- NA # empty vector to fill with modularity estimates
mod.reps <- list() # empty list to fill with the various estimates of modularity from the independent runs (n=NREPS)

web <- unique(PH2$soil_year)
nweb <- length(unique(PH2$soil_year))

NREPS <- 15 # modelled after Vizentin-Bugoni; they did 10 independent runs

for(i in 1:nweb) {
  
  temp1 <- subset(PH2, soil_year == web[i])
  temp2 <- subset(temp1, select=c("genus", "cat.genus.species", "abu"))
  mat <- matrify(temp2)
  
  mod.rep <- replicate(NREPS, {   # Better to do more. 100 (Gomez 2014, New Phytologist) and 10 (Vizentin-Bugoni 2016, Maruyama paper too) - I'll do 150
    tryCatch({
      mod <- computeModules(mat, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE, 
                            steps = 10000000, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)
      return(mod@likelihood)},
      error=function(e){})
  })
  
  mod.reps[[i]] <- mod.rep # shouldn't really need to refer to this, except to confirm that the 15 independent runs weren't arriving at variable modularities (the high number of swaps within the computeModules algorithm I think helps it find the same maximum every time, at least from what I have seen)
  modularity.est[i] <- max(unlist(mod.rep)) # take maximum likelihood computed; this is the best estimate of modularity
  web.id[i] <- temp1[3,3]
  
}

results.empirical.maxmod <- data.frame(web.id, modularity.est) # The single highest value from each of the 15 independent computeModules runs per web
results.empirical.reps <- data.frame(do.call(rbind, mod.reps)) # The 15 independent computeModules runs, per web (4 X 15 dataframe)
results.empirical.reps$web.id <- results.empirical.maxmod$web.id

###############################################################################################################
# 2)  Compare modularity of each of these four networks to the modularity we would expect from null structure. 
###############################################################################################################
# This loop...
# 1) Uses that empirical web to generate some number of null networks (n=NNULLS). For each of those null networks, it calculates the modularity score. 
# 2) Use different null models to ask different questions about mechanism. See Dormann (2009) - "Indices, graphs, and null models" (p12) for an awesome description of what these nulls are telling us, biologically, about the different roles of species abundances & interactions!

# Null Model I) - Patefield = "r2dtable" (marginal totals as observed) 
# Null Model II) - Shuffle = "shuffle.web" (relocates entries; connectance as observed)
# Null Model III) - Swap = "swap.web" (marginal totals AND connectance as observed) 

# 3) I can make a figure with the modularity distributions generated from various null models, and also plot my empirical modularities (point estimates).

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Null Model I - PATEFIELD (marginal totals as observed) #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
web.id <- NA # empty vector to fill with the soil_year partition
null.modularity <- list() # create an empty list to fill with the null matrices' modularity calculations (# NNULLS per soil_year nurtwurk)
max.mod.null <- NA # create an empty vector to fill with the maximum modularity acheived across all the null models (not really useful, since I will be using the distribution of null values, but just interesting to know)

web <- unique(PH2$soil_year)
nweb <- length(unique(PH2$soil_year))

NNULLS <- 1000 # Vizentin-Bugoni did 1000 null networks to look at other network network metrics, but only 100 for modularity, because too computationally intensive - no longer, thanks to Beckett! I can run 1000!

for(i in 1:nweb) {
  
  temp1 <- subset(PH2, soil_year==web[i])
  temp2 <- subset(temp1, select=c("genus", "cat.genus.species", "abu"))
  nurtwurk <- matrify(temp2)
  
  message(" *** generating null networks ***")
  nulls <- nullmodel(nurtwurk, N=NNULLS, method="r2dtable") # creates many null matrices from the given soil x year partition. 
  
  message(" *** starting computeModules() for null networks ***")
  
  likelihood <- unlist(lapply(nulls, function(null) {
    cmres <- tryCatch({computeModules(null, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE, 
                                      steps = 10000000, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)}, # For each of the null matrices (n = NNULLS), compute modularity
                      error=function(e){})
    if (is.null(cmres)) {
      warning("computeModules failed")
      return(NA)
    }
    return(cmres@likelihood) # Likelihood = a vector of null modularity values for the null matrices (should be the length of the NNULLS - might be a few less if some of the nulls dont converge) 
  }))
  
  message("   computeModules() for nulls finished")
  
  web.id[i] <- temp1[3,3]
  max.mod.null[i] <- max(likelihood)
  null.modularity[[i]] <- likelihood # vector of all the modularity outputs for the nulls (length NNULLS)
}

# Organize the null model runs into a DF (4 X 1000)
info.patefield <- data.frame(web.id); info.patefield$merger <- rownames(info.patefield)
null.mods.patefield <- data.frame(do.call(rbind, null.modularity)); null.mods.patefield$merger <- rownames(null.mods.patefield)
null.patefield <- merge(info.patefield, null.mods.patefield, by="merger") # OK - so this DF has 4 rows, one for each web, and many columns (n = NNULLS) - one for the output of each of the null webs' modularities.

# To understand what's up under the hood, let's look at the dimensionality and entries of a few of the nulls, from the last run:
nurtwurk # dimensions are 4 X 36
null.10 <- data.frame(nulls[10])
null.500 <- data.frame(nulls[500])
null.999 <- data.frame(nulls[999])

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Null Model II - SHUFFLE (connectance only as observed, marginal totals vary) # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
web.id <- NA # empty vector to fill with the soil_year partition
null.modularity <- list() # create an empty list to fill with the null matrices' modularity calculations (# NNULLS per soil_year nurtwurk)
max.mod.null <- NA # create an empty vector to fill with the maximum modularity acheived across all the null models (not really useful, since I will be using the distribution of null values, but just interesting to know)

web <- unique(PH2$soil_year)
nweb <- length(unique(PH2$soil_year))

NNULLS <- 1000 # Vizentin-Bugoni did 1000 null networks to look at other network network metrics, but only 100 for modularity, because too computationally intensive - no longer, thanks to Beckett! I can run 1000!

for(i in 1:nweb) {
  
  temp1 <- subset(PH2, soil_year==web[i])
  temp2 <- subset(temp1, select=c("genus", "cat.genus.species", "abu"))
  nurtwurk <- matrify(temp2)
  
  message(" *** generating null networks ***")
  nulls <- nullmodel(nurtwurk, N=NNULLS, method="shuffle.web") # creates many null matrices from the given soil x year partition. 
  
  message(" *** starting computeModules() for null networks ***")
  
  likelihood <- unlist(lapply(nulls, function(null) {
    cmres <- tryCatch({computeModules(null, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE, 
                                      steps = 10000000, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)}, # For each of the null matrices (n = NNULLS), compute modularity
                      error=function(e){})
    if (is.null(cmres)) {
      warning("computeModules failed")
      return(NA)
    }
    return(cmres@likelihood) # Likelihood = a vector of null modularity values for the null matrices (should be the length of the NNULLS - might be a few less if some of the nulls dont converge) 
  }))
  
  message("   computeModules() for nulls finished")
  
  web.id[i] <- temp1[3,3]
  max.mod.null[i] <- max(likelihood)
  null.modularity[[i]] <- likelihood # vector of all the modularity outputs for the nulls (length NNULLS)
}

# Organize the null model runs into a DF (4 X 1000)
info.shuffle <- data.frame(web.id); info.shuffle$merger <- rownames(info.shuffle)
null.mods.shuffle <- data.frame(do.call(rbind, null.modularity)); null.mods.shuffle$merger <- rownames(null.mods.shuffle)
null.shuffle <- merge(info.shuffle, null.mods.shuffle, by="merger") # OK - so this DF has 4 rows, one for each web, and many columns (n = NNULLS) - one for the output of each of the null webs' modularities.

# To understand what's up under the hood, let's look at the dimensionality and entries of a few of the nulls, from the last run:
nurtwurk # dimensions are 4 X 36
null.10 <- data.frame(nulls[10])
null.500 <- data.frame(nulls[500])
null.999 <- data.frame(nulls[999])

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Null Model III - SWAP (connectance AND marginal totals as observed) #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
web.id <- NA # empty vector to fill with the soil_year partition
null.modularity <- list() # create an empty list to fill with the null matrices' modularity calculations (# NNULLS per soil_year nurtwurk)
max.mod.null <- NA # create an empty vector to fill with the maximum modularity acheived across all the null models (not really useful, since I will be using the distribution of null values, but just interesting to know)

web <- unique(PH2$soil_year)
nweb <- length(unique(PH2$soil_year))

NNULLS <- 1000 # Vizentin-Bugoni did 1000 null networks to look at other network network metrics, but only 100 for modularity, because too computationally intensive - no longer, thanks to Beckett! I can run 1000!

for(i in 1:nweb) {
  
  temp1 <- subset(PH2, soil_year==web[i])
  temp2 <- subset(temp1, select=c("genus", "cat.genus.species", "abu"))
  nurtwurk <- matrify(temp2)
  
  message(" *** generating null networks ***")
  nulls <- nullmodel(nurtwurk, N=NNULLS, method="swap.web") # creates many null matrices from the given soil x year partition. 
  
  message(" *** starting computeModules() for null networks ***")
  
  likelihood <- unlist(lapply(nulls, function(null) {
    cmres <- tryCatch({computeModules(null, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE, 
                                      steps = 10000000, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)}, # For each of the null matrices (n = NNULLS), compute modularity
                      error=function(e){})
    if (is.null(cmres)) {
      warning("computeModules failed")
      return(NA)
    }
    return(cmres@likelihood) # Likelihood = a vector of null modularity values for the null matrices (should be the length of the NNULLS - might be a few less if some of the nulls dont converge) 
  }))
  
  message("   computeModules() for nulls finished")
  
  web.id[i] <- temp1[3,3]
  max.mod.null[i] <- max(likelihood)
  null.modularity[[i]] <- likelihood # vector of all the modularity outputs for the nulls (length NNULLS)
}

# Organize the null model runs into a DF (4 X 1000)
info.swap <- data.frame(web.id); info.swap$merger <- rownames(info.swap)
null.mods.swap <- data.frame(do.call(rbind, null.modularity)); null.mods.swap$merger <- rownames(null.mods.swap)
null.swap <- merge(info.swap, null.mods.swap, by="merger") # OK - so this DF has 4 rows, one for each web, and many columns (n = NNULLS) - one for the output of each of the null webs' modularities.

# To understand what's up under the hood, let's look at the dimensionality and entries of a few of the nulls, from the last run:
nurtwurk # dimensions are 4 X 36
null.10 <- data.frame(nulls[10])
null.500 <- data.frame(nulls[500])
null.999 <- data.frame(nulls[999])

######################
# 3)  Make a figure! #
######################
# To visualize: plot mean + 95% CI of values generated from nulls (see Vizentin-Bugoni 2016, p264), + my point estimates
# In text, to discuss statistical significance, pair with Z-score comparison (next section)

# Reformat data for plotting
# Null Model I
patefield.plot1 <- null.patefield
patefield.plot1$soil <- gsub("_.+$", "", patefield.plot1$web.id)
patefield.plot1$year <- gsub("^.+_", "", patefield.plot1$web.id)
patefield.plot1$null <- "Patefield"
patefield.plot2 <- patefield.plot1[ , -which(names(patefield.plot1) %in% c("merger","web.id"))] # remove merger and web.id - they'll mess up the melt
molten.patefield <- melt(patefield.plot2, variable.name = "run", value.names = "value", id.vars = c("soil", "year", "null")); head(molten.patefield)
molten.patefield$run <- NULL

# Null Model II
shuffle.plot1 <- null.shuffle
shuffle.plot1$soil <- gsub("_.+$", "", shuffle.plot1$web.id)
shuffle.plot1$year <- gsub("^.+_", "", shuffle.plot1$web.id)
shuffle.plot1$null <- "Shuffle"
shuffle.plot2 <- shuffle.plot1[ , -which(names(shuffle.plot1) %in% c("merger","web.id"))] # remove merger and web.id - they'll mess up the melt
molten.shuffle <- melt(shuffle.plot2, variable.name = "run", value.names = "value", id.vars = c("soil", "year", "null")); head(molten.patefield)
molten.shuffle$run <- NULL

# Null Model III
swap.plot1 <- null.swap
swap.plot1$soil <- gsub("_.+$", "", swap.plot1$web.id)
swap.plot1$year <- gsub("^.+_", "", swap.plot1$web.id)
swap.plot1$null <- "Swap"
swap.plot2 <- swap.plot1[ , -which(names(swap.plot1) %in% c("merger","web.id"))] # remove merger and web.id - they'll mess up the melt
molten.swap <- melt(swap.plot2, variable.name = "run", value.names = "value", id.vars = c("soil", "year", "null")); head(molten.swap)
molten.swap$run <- NULL

# Combine all the null model runs
all.nulls <- rbind(molten.patefield, molten.shuffle, molten.swap)
all.nulls$soil2 <- "Serpentine"
all.nulls$soil2[which(all.nulls$soil == "ns")] <- "Non-serpentine"

# Summarize the null models...
summ <- summarySE(data=all.nulls, measurevar = "value", groupvars = c("soil2", "year", "null")) 
summ$sd2 <- summ$sd*2

# Modify the dataframe with the empirical values to overplot...
mod.emp <- results.empirical.maxmod
mod.emp$soil <- gsub("_.+$", "", mod.emp$web.id)
mod.emp$year <- gsub("^.+_", "", mod.emp$web.id)
mod.emp$value <- modularity.est
mod.emp$soil2 <- "Serpentine"
mod.emp$soil2[which(mod.emp$soil == "ns")] <- "Non-serpentine"

# Make an alternative legend for paper clarity...
summ$null2 <- "Null Model I (Patefield)"
summ$null2[which(summ$null == "Shuffle")] <- "Null Model II (Shuffle)"
summ$null2[which(summ$null == "Swap")] <- "Null Model III (Swap)"

# Plot
ggplot(data=summ, aes(x=soil2, y=value)) +
  geom_errorbar(aes(ymin=value-sd2, ymax=value+sd2, color=null2), size=4, alpha=0.5, width=0) +
  scale_color_manual("Null Model", values=c("#0072B2", "#009E73", "#E69F00")) +
  geom_point(data = mod.emp, aes(x=soil2, y=value), size=4, color="grey20") +
# geom_line(data = mod.emp, aes(x=soil2, y=value, group=year), size=0.5, lty=3) +
  theme_bw()+
  facet_grid(~year) +
  theme(panel.border = element_rect(colour="black"), panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), strip.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Soil") +
  ylab("Modularity (Q)") +
  scale_x_discrete(labels = c("Non-\nserpentine","Serpentine")) +
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=10), strip.text.x=element_text(size=15, face="bold"), axis.title.y=element_text(size=15,face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.title.x=element_text(size=15,face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.title = element_text(size=12, face="bold")) 

########################
# 4)  Z-score approach #
########################
# See Dormann and Strauss (2013) PREPRINT (p12).

# "Furthermore, the absolute value of Q (like all network indices: Dormann et al., 2009) is dependent on network size (i.e. the number of species) as well as the number of links and the total number of interactions observed (see also Thébault, 2013). We would thus recommend a null model comparison (e.g. Vázquez & Aizen, 2003; Blüthgen et al., 2008; Dormann et al., 2009) to correct the observed value of Q by null model expectation (e.g. by standardising them to z-scores. The resultant value gives the number of standard deviations higher than expected the empirical web modularity is from the random networks. Since z-scores are assumed to be normally distributed, values above 2 are considered significant. In R, this could be achieved by the following code. 

# Are my empirical modularity values significantly different from the null values generated by....

# NOTE - as these null models rely on random resampling, the average modularities and Z-scores below will be slightly different each time this code is run.

##############
# Patefield? # YES - EXTREMELY DIFFERENT 
# Interpretation: modularity not recapitulated from species abundance distributions alone
##############
null.patefield$mean.null.mod <- rowMeans(as.matrix(subset(null.patefield, select=X1:X1000)))
null.patefield$sd.null.mod <- rowSds(as.matrix(subset(null.patefield, select=X1:X1000)))
patefield2 <- merge(null.patefield, results.empirical.maxmod, by="web.id")
patefield2$z <- (patefield2$modularity.est - patefield2$mean.null.mod)/(patefield2$sd.null.mod)
patefield <- subset(patefield2, select=c("web.id", "modularity.est", "mean.null.mod", "sd.null.mod", "z"))
patefield$null <- "Patefield (r2d)"

############
# SHUFFLE? # NO - MODULARITY RECAPITULATED
# Interpretation: moving cell entries around results in similar modularity (similarly sparse matrix) 
############
null.shuffle$mean.null.mod <- rowMeans(as.matrix(subset(null.shuffle, select=X1:X1000)))
null.shuffle$sd.null.mod <- rowSds(as.matrix(subset(null.shuffle, select=X1:X1000)))
shuffle2 <- merge(null.shuffle, results.empirical.maxmod, by="web.id")
shuffle2$z <- (shuffle2$modularity.est - shuffle2$mean.null.mod)/(shuffle2$sd.null.mod)
shuffle <- subset(shuffle2, select=c("web.id", "modularity.est", "mean.null.mod", "sd.null.mod", "z"))
shuffle$null <- "Shuffle"

#########
# Swap? # YES - DIFFERENT
# Interpretation: including abundances reduces the modularity; reinforces that abundances are not driving this pattern 
#########
null.swap$mean.null.mod <- rowMeans(as.matrix(subset(null.swap, select=X1:X1000)))
null.swap$sd.null.mod <- rowSds(as.matrix(subset(null.swap, select=X1:X1000)))
swap2 <- merge(null.swap, results.empirical.maxmod, by="web.id")
swap2$z <- (swap2$modularity.est - swap2$mean.null.mod)/(swap2$sd.null.mod)
swap <- subset(swap2, select=c("web.id", "modularity.est", "mean.null.mod", "sd.null.mod", "z"))
swap$null <- "Swap"