#######################################
# Plant-herbivore modularity analyses #
#######################################
# Goal: build curves showing relationship between network size (# of species, # of links) and modularity. Show that differences in modularity are maintained across a range of rarefaction targets

# 1) Do this for regular Q (output of DIRTLPAWb+)
# 2) Also calculate normalized Q (e.g. divide Q by the maximum achievable modularity given the network size (e.g. without any generalized interactions - per correspondence with Stephen Beckett)

#############
# Packages #
#############
library(bipartite)
library(labdsv)
library(ggplot2)

########
# Data #
########
PH <- read.csv("~/Github/data/2014-2015_network.csv")

#############################
# Prepare data for analyses #
#############################
# a) Make a variable to loop over, since I want to do this separately for each network in each year (soil x year partition)
PH$soil_year <- paste(PH$soil, PH$year, sep="_"); unique(PH$soil_year)

# b) Get abundances for each caterpillar species, for each plant genus X soil X year. I will make matrices from these abundances.
PH$abu <- 1
PH2 <- aggregate(abu ~ cat.genus.species + genus + soil_year, data=PH, length)

#### Will use this 'info' later , to merge in soil & year info #####
info <- aggregate(abu ~ soil_year + soil + year, data=PH, length); info$abu <- NULL

############################################################
# 1) Rarefy all networks & calculate modularity Q for each #           
############################################################
# This is what the loop does:
# 1) Counts the number of interactions in a particular soil_web; makes an interval of network sizes from 25 to that number, going by 5's. 
# 2) For each of those size intervals, assembles a random network with that # of interactions & calculates modularity. Does this multiple times, so that multiple random networks of the same size are assembled and modularity calculated.
# 3) Final desired output: for each soil_web network, plot a curve showing how modularity changes across different subsampled network sizes. Error/variance around the curve shows the variation in modularities achieved from different networks @ that same size, due to random assembly of interactions.
web <- unique(PH2$soil_year)
nweb <- length(unique(PH2$soil_year))
all.subwebs <- list()
# ns_2014: 429 interactions; 81 DFs, ends @ 425
# ns_2015: 1122 interactions; 220 DFs, ends @ 1120
# s_2014: 286 interactions; 53 DFs, ends @ 285
# s_2015: 367 interactions; 69 DFs, ends @ 365

for(i in 1:nweb) {
#  temp1 <- subset(PH, soil_year == "s_2015"); nrow(temp1)
  temp1 <- subset(PH, soil_year == web[i])
  count.ints <- nrow(temp1)
  sizes <- seq(from=25, to=count.ints, by=5) # REDUCE to 5
  nsizes <- length(sizes); nsizes 
  all.mod <- list() # clear all.mod list for each network 
  
    for(j in 1:nsizes) { # at each matrix size (# interactions).... FROM HERE
        mod.max <- NA; rep <- NA
      for(x in 1:100) { # repeat the follow subsampling & modularity calculation 100x...
#         subweb <- temp1[sample(nrow(temp1), 25),]
          subweb <- temp1[sample(nrow(temp1), sizes[j]),] # take random subset of j interactions 
          subweb$abu <- 1
          subweb2 <- aggregate(abu ~ cat.genus.species + genus, data=subweb, length)
          subweb3 <- matrify(subweb2)
      #   calculate modularity for this subweb
            mod.rep <- replicate(15, {   # Nreps: Better to do more, to find the highest modularity. 100 (Gomez 2014, New Phytologist) and 10 (Vizentin-Bugoni 2016, Maruyama paper too) - I'll do 15 for now.
            tryCatch({
            mod <- computeModules(subweb3, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE,
            steps = 10000000, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)
            return(mod@likelihood)},
            error=function(e){})
            })
          mod.max[x] <- max(mod.rep)
          rep[x] <- x
        }
      all <- data.frame(cbind(mod.max, rep))
      all$size <- sizes[j]
      all.mod[[j]] <- all # for each web, should be a list of j dfs; one df per size category
      print(all)
      print(temp1[1,22])
    }  
  subweb.together <- data.frame(do.call(rbind, all.mod)) # put all the j dfs together
  subweb.together$web <- temp1[1,22] # add a column to name the web
  all.subwebs[[i]] <- subweb.together
  print(subweb.together)
  flush.console()
}

str(all.subwebs) # should be a list of 4 dfs, one per soil/year network, each of different dimensions 
all.webs <- data.frame(do.call(rbind, all.subwebs)); nrow(all.webs) # for 100 random subwebs per slice, should be 42300 rows 
# write.csv(all.webs, "~/Github/data/Q_from.subsampling.csv")

################################################################
# 1.1) Calculate EMPIRICAL modularity Q for each network (n=4) #           
################################################################
# First, I might want to add the empirical modularities to the plot...
# Calculate single empirical modularity for each web
#------------------#
# Serpentine, 2014 #
#------------------#
s_2014 <- subset(PH, soil_year == "s_2014")
s_2014$abu <- 1
s_2014.2 <- aggregate(abu ~ cat.genus.species + genus, data=s_2014, length)
s_2014.3 <- matrify(s_2014.2)
mod.rep <- replicate(15, {   # Nreps: Better to do more, to find the highest modularity. 100 (Gomez 2014, New Phytologist) and 10 (Vizentin-Bugoni 2016, Maruyama paper too) - I'll do 15 for now.
        tryCatch({
          mod <- computeModules(s_2014.3, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE,
          steps = 10000000, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)
          return(mod@likelihood)},
          error=function(e){})})
mod.q1 <- data.frame(max(mod.rep))
mod.q1$size <- nrow(s_2014)
mod.q1$network <- "s_2014" 
names(mod.q1) <- c("mod.q", "size", "soil_year"); mod.q1
#------------------#
# Serpentine, 2015 #
#------------------#
s_2015 <- subset(PH, soil_year == "s_2015")
s_2015$abu <- 1
s_2015.2 <- aggregate(abu ~ cat.genus.species + genus, data=s_2015, length)
s_2015.3 <- matrify(s_2015.2)
mod.rep <- replicate(15, {   # Nreps: Better to do more, to find the highest modularity. 100 (Gomez 2014, New Phytologist) and 10 (Vizentin-Bugoni 2016, Maruyama paper too) - I'll do 15 for now.
  tryCatch({
    mod <- computeModules(s_2015.3, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE,
                          steps = 10000000, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)
    return(mod@likelihood)},
    error=function(e){})})
mod.q2 <- data.frame(max(mod.rep))
mod.q2$size <- nrow(s_2015)
mod.q2$network <- "s_2015" 
names(mod.q2) <- c("mod.q", "size", "soil_year"); mod.q2

#----------------------#
# Non-serpentine, 2014 #
#----------------------#
ns_2014 <- subset(PH, soil_year == "ns_2014")
ns_2014$abu <- 1
ns_2014.2 <- aggregate(abu ~ cat.genus.species + genus, data=ns_2014, length)
ns_2014.3 <- matrify(ns_2014.2)
mod.rep <- replicate(15, {   # Nreps: Better to do more, to find the highest modularity. 100 (Gomez 2014, New Phytologist) and 10 (Vizentin-Bugoni 2016, Maruyama paper too) - I'll do 15 for now.
  tryCatch({
    mod <- computeModules(ns_2014.3, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE,
                          steps = 10000000, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)
    return(mod@likelihood)},
    error=function(e){})})
mod.q3 <- data.frame(max(mod.rep))
mod.q3$size <- nrow(ns_2014)
mod.q3$network <- "ns_2014" 
names(mod.q3) <- c("mod.q", "size", "soil_year"); mod.q3

#----------------------#
# Non-serpentine, 2015 #
#----------------------#
ns_2015 <- subset(PH, soil_year == "ns_2015")
ns_2015$abu <- 1
ns_2015.2 <- aggregate(abu ~ cat.genus.species + genus, data=ns_2015, length)
ns_2015.3 <- matrify(ns_2015.2)
mod.rep <- replicate(15, {   # Nreps: Better to do more, to find the highest modularity. 100 (Gomez 2014, New Phytologist) and 10 (Vizentin-Bugoni 2016, Maruyama paper too) - I'll do 15 for now.
  tryCatch({
    mod <- computeModules(ns_2015.3, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE,
                          steps = 10000000, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)
    return(mod@likelihood)},
    error=function(e){})})
mod.q4 <- data.frame(max(mod.rep))
mod.q4$size <- nrow(ns_2015)
mod.q4$network <- "ns_2015" 
names(mod.q4) <- c("mod.q", "size", "soil_year"); mod.q4

all_single_mods <- rbind(mod.q1, mod.q2, mod.q3, mod.q4)
all_single_mods2 <- merge(all_single_mods, info, by="soil_year")

############################################################################
# 2.1) SUBSAMPLE all networks & calculate NORMALIZED modularity Q for each #           
############################################################################
# This is what the loop does:
# 1) Counts the number of interactions in a particular soil_web; makes an interval of network sizes from 25 to that number, going by 5's. 
# 2) For each of those size intervals, assembles a random network with that # of interactions & then:
#     a. Calculates modularity Q (this is the same as above)
#     b. Jiggles the network to remove any generalized interactions. Calculates modularity. Does this jiggling / modularity calculation a few times. Save the maximum modularity of those (=Qmax). 
#     c. Divide Q by Qmax. This is the normalized modularity for that random network.
# 3) Make curves, like above, but with the normalized modularity scores across sizes.

web <- unique(PH2$soil_year)
nweb <- length(unique(PH2$soil_year))

all.subwebs <- list()
# ns_2014: 429 interactions; 81 DFs, ends @ 425
# ns_2015: 1122 interactions; 220 DFs, ends @ 1120
# s_2014: 286 interactions; 53 DFs, ends @ 285
# s_2015: 367 interactions; 69 DFs, ends @ 365

for(i in 1:nweb) {
#  temp1 <- subset(PH, soil_year == "s_2015"); nrow(temp1) # 367
  temp1 <- subset(PH, soil_year == web[i])
  count.ints <- nrow(temp1)
  sizes <- seq(from=25, to=count.ints, by=5) # Intervals of 5
  nsizes <- length(sizes); nsizes 
  all.mod <- list() # clear all.mod list for each network 
  
  for(j in 1:nsizes) { # at each matrix size (# interactions)....
    mod.q <- NA; rep.q <- NA; max.q <- NA
#-#     
    for(x in 1:10) { 
#     subweb <- temp1[sample(nrow(temp1), 100),] # MAKE ONE SUBWEB
      subweb <- temp1[sample(nrow(temp1), sizes[j]),] # take random subset of j interactions 
      subweb$abu <- 1
      subweb2 <- aggregate(abu ~ cat.genus.species + genus, data=subweb, length)
      subweb3 <- matrify(subweb2)
      # calculate modularity Q for this subweb
      mod.rep.q <- replicate(15, {tryCatch({
          mod <- computeModules(subweb3, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE, steps = 10000000, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)
          return(mod@likelihood)},
          error=function(e){})
      })
      
      mod.q[x] <- max(mod.rep.q) # max modularity value for that subweb
      rep.q[x] <- x # the rep (1-5, etc.)
      
      # calculate maximum possible modularity for this subweb
      mod.max.q.max <- NA; rep.q.max <- NA    
      for (y in 1:15) {
            all.abu <- aggregate(abu ~ cat.genus.species, data=subweb, length) # get their abundances, aggregating across all host plants
            all.abu$ran.genus <- sample(1:4, size = nrow(all.abu), replace = TRUE)
            all.abu$ran.genus2 <- NA
            all.abu$ran.genus2[which(all.abu$ran.genus == 1)] <- "Ceanothus"
            all.abu$ran.genus2[which(all.abu$ran.genus == 2)] <- "Adenostoma"
            all.abu$ran.genus2[which(all.abu$ran.genus == 3)] <- "Arctostaphylos"
            all.abu$ran.genus2[which(all.abu$ran.genus == 4)] <- "Quercus"
            all.abu$ran.genus2 <- as.factor(all.abu$ran.genus2)
            subweb_modular1 <- subset(all.abu, select=c("cat.genus.species", "ran.genus2", "abu"))
            subweb_modular2 <- matrify(subweb_modular1)
            # Calculate modularity for this modularized subweb
            mod.rep.max <- replicate(15, {tryCatch({
            mod <- computeModules(subweb_modular2, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE, steps = 10000000, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)
          return(mod@likelihood)},
          error=function(e){})
      })
      mod.max.q.max[y] <- max(mod.rep.max) # single value that is the maximum value for that subset's dimensions
      rep.q.max[y] <- y # the rep (1-5 etc)
    }
    all.qmax <- data.frame(cbind(mod.max.q.max, rep.q.max)) # should be a df 2 X the # of modularized subwebs
    max.q[x] <- max(all.qmax$mod.max.q.max) # for each of the subwebs, this is the max modularity value across those 20 modularized versions of the subweb
    }
#-# 
    q.norm <- data.frame(cbind(mod.q, rep.q, max.q)) # for that size, 5 random webs & the corresponding maximum modularity value possible for that random subsample
    q.norm$size <- sizes[j]
    all.mod[[j]] <- q.norm # for each web, should be a list of j dfs; one df per size category
    print(q.norm)
    print(temp1[1,22])
  } #-# 
  subweb.together <- data.frame(do.call(rbind, all.mod)) # put all the j dfs together
  subweb.together$web <- temp1[1,22] # add a column to name the web
  all.subwebs[[i]] <- subweb.together
  print(subweb.together)
  flush.console()
}

str(all.subwebs) # should be a list of 4 dfs, one per soil/year network, each of different dimensions 
all.norm <- data.frame(do.call(rbind, all.subwebs)); nrow(all.norm) 
all.norm$soil_year <- all.norm$web
all.norm3 <- merge(all.norm, info, by="soil_year")
# Calculate NORMALIZED modularity for each web
all.norm3$q.norm <- all.norm3$mod.q/all.norm3$max.q
#-#-#-#-#-#-#-
# Write .csv #
#-#-#-#-#-#-#-
head(all.norm3)
write.csv(all.norm3, "~/Github/data/NormalizedQ_from.subsampling.csv")

####################################################################
# 2.2) Calculate MAXIMUM  modularity Q for EMPIRICAL network (n=4) #           
####################################################################
#------------------#
# serpentine, 2014 #
#------------------#
s_2014 <- subset(PH, soil_year == "s_2014")
# calculate maximum possible modularity for this subweb
q.max <- NA
for (y in 1:100) {
  all.abu <- aggregate(abu ~ cat.genus.species, data=s_2014, length) # get their abundances, aggregating across all host plants
  all.abu$ran.genus <- sample(1:4, size = nrow(all.abu), replace = TRUE)
  all.abu$ran.genus2 <- NA
  all.abu$ran.genus2[which(all.abu$ran.genus == 1)] <- "Ceanothus"
  all.abu$ran.genus2[which(all.abu$ran.genus == 2)] <- "Adenostoma"
  all.abu$ran.genus2[which(all.abu$ran.genus == 3)] <- "Arctostaphylos"
  all.abu$ran.genus2[which(all.abu$ran.genus == 4)] <- "Quercus"
  all.abu$ran.genus2 <- as.factor(all.abu$ran.genus2)
  subweb_modular1 <- subset(all.abu, select=c("cat.genus.species", "ran.genus2", "abu"))
  subweb_modular2 <- matrify(subweb_modular1)
  # Calculate modularity for this modularized subweb
  mod.rep.max <- replicate(15, {tryCatch({
    mod <- computeModules(subweb_modular2, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE, steps = 10000000, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)
    return(mod@likelihood)},
    error=function(e){})
  })
  q.max[y] <- max(mod.rep.max)}
max(q.max) # 0.7465402 after 50 modular webs; 0.7474449 after 100. So, 100 seems fine.
q.max1 <- data.frame(max(q.max))
q.max1$size <- nrow(s_2014)
q.max1$network <- "s_2014" 
names(q.max1) <- c("q.max", "size", "soil_year"); q.max1

#------------------#
# serpentine, 2015 #
#------------------#
s_2015 <- subset(PH, soil_year == "s_2015")
# calculate maximum possible modularity for this subweb
q.max <- NA
for (y in 1:100) {
  all.abu <- aggregate(abu ~ cat.genus.species, data=s_2015, length) # get their abundances, aggregating across all host plants
  all.abu$ran.genus <- sample(1:4, size = nrow(all.abu), replace = TRUE)
  all.abu$ran.genus2 <- NA
  all.abu$ran.genus2[which(all.abu$ran.genus == 1)] <- "Ceanothus"
  all.abu$ran.genus2[which(all.abu$ran.genus == 2)] <- "Adenostoma"
  all.abu$ran.genus2[which(all.abu$ran.genus == 3)] <- "Arctostaphylos"
  all.abu$ran.genus2[which(all.abu$ran.genus == 4)] <- "Quercus"
  all.abu$ran.genus2 <- as.factor(all.abu$ran.genus2)
  subweb_modular1 <- subset(all.abu, select=c("cat.genus.species", "ran.genus2", "abu"))
  subweb_modular2 <- matrify(subweb_modular1)
  # Calculate modularity for this modularized subweb
  mod.rep.max <- replicate(15, {tryCatch({
    mod <- computeModules(subweb_modular2, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE, steps = 10000000, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)
    return(mod@likelihood)},
    error=function(e){})
  })
  q.max[y] <- max(mod.rep.max)}
q.max2 <- data.frame(max(q.max))
q.max2$size <- nrow(s_2015)
q.max2$network <- "s_2015" 
names(q.max2) <- c("q.max", "size", "soil_year"); q.max2

#----------------------#
# Non-serpentine, 2014 #
#----------------------#
ns_2014 <- subset(PH, soil_year == "ns_2014")
# calculate maximum possible modularity for this subweb
q.max <- NA
for (y in 1:100) {
  all.abu <- aggregate(abu ~ cat.genus.species, data=ns_2014, length) # get their abundances, aggregating across all host plants
  all.abu$ran.genus <- sample(1:4, size = nrow(all.abu), replace = TRUE)
  all.abu$ran.genus2 <- NA
  all.abu$ran.genus2[which(all.abu$ran.genus == 1)] <- "Ceanothus"
  all.abu$ran.genus2[which(all.abu$ran.genus == 2)] <- "Adenostoma"
  all.abu$ran.genus2[which(all.abu$ran.genus == 3)] <- "Arctostaphylos"
  all.abu$ran.genus2[which(all.abu$ran.genus == 4)] <- "Quercus"
  all.abu$ran.genus2 <- as.factor(all.abu$ran.genus2)
  subweb_modular1 <- subset(all.abu, select=c("cat.genus.species", "ran.genus2", "abu"))
  subweb_modular2 <- matrify(subweb_modular1)
  # Calculate modularity for this modularized subweb
  mod.rep.max <- replicate(15, {tryCatch({
    mod <- computeModules(subweb_modular2, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE, steps = 10000000, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)
    return(mod@likelihood)},
    error=function(e){})
  })
  q.max[y] <- max(mod.rep.max)}
max(q.max) 
q.max3 <- data.frame(max(q.max))
q.max3$size <- nrow(ns_2014)
q.max3$network <- "ns_2014" 
names(q.max3) <- c("q.max", "size", "soil_year"); q.max3

#----------------------#
# Non-serpentine, 2015 #
#----------------------#
ns_2015 <- subset(PH, soil_year == "ns_2015")
# calculate maximum possible modularity for this subweb
q.max <- NA
for (y in 1:100) {
  all.abu <- aggregate(abu ~ cat.genus.species, data=ns_2015, length) # get their abundances, aggregating across all host plants
  all.abu$ran.genus <- sample(1:4, size = nrow(all.abu), replace = TRUE)
  all.abu$ran.genus2 <- NA
  all.abu$ran.genus2[which(all.abu$ran.genus == 1)] <- "Ceanothus"
  all.abu$ran.genus2[which(all.abu$ran.genus == 2)] <- "Adenostoma"
  all.abu$ran.genus2[which(all.abu$ran.genus == 3)] <- "Arctostaphylos"
  all.abu$ran.genus2[which(all.abu$ran.genus == 4)] <- "Quercus"
  all.abu$ran.genus2 <- as.factor(all.abu$ran.genus2)
  subweb_modular1 <- subset(all.abu, select=c("cat.genus.species", "ran.genus2", "abu"))
  subweb_modular2 <- matrify(subweb_modular1)
  # Calculate modularity for this modularized subweb
  mod.rep.max <- replicate(15, {tryCatch({
    mod <- computeModules(subweb_modular2, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE, steps = 10000000, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)
    return(mod@likelihood)},
    error=function(e){})
  })
  q.max[y] <- max(mod.rep.max)}
q.max4 <- data.frame(max(q.max))
q.max4$size <- nrow(ns_2015)
q.max4$network <- "ns_2015" 
names(q.max4) <- c("q.max", "size", "soil_year"); q.max4

all_single_max.mods <- rbind(q.max1, q.max2, q.max3, q.max4)
all_mods <- merge(all_single_mods, all_single_max.mods, by=c("soil_year", "size"))
all_mods
all_mods$q.norm <- all_mods$mod.q/all_mods$q.max
all_mods2 <- merge(all_mods, info, by="soil_year")

#-#-#-#-#-#-#-
# Write .csv #
#-#-#-#-#-#-#-
write.csv(all_mods2, "~/Github/data/emp_q&q_max.csv")
