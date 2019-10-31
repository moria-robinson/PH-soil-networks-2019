#######################################
# Plant-herbivore networks - visualize #
#######################################

# 1) Visualize Plant-Herbivore networks

#############
# Packages #
#############
library(bipartite)
library(labdsv)

########
# Data #
########
PH <- read.csv("~/Github/data/2014-2015_network.csv")

# For abundance of host plants
plants <- read.csv("~/Github/data/plants.csv")
soil <- read.csv("~/Github/data/soil.csv") 

#############################
# Prepare data for analyses #
#############################
# First, notice that PH has some presumed diet breadths of un-identified species. These aren't errors, they are intentional; they are from the 'species key' .csv, and presumed from large numbers of individuals collected from the same plant genyues). However, for these figures, we want to match our species-level analyses of diet breadth, which are only with species confirmed to species ID. SO:

# Manually change diet breadth to "unk" for any unconfirmed IDs:
PH$hostbreadth[which(PH$cat.genus.species=="Geometrid.09")] <- "unk"
PH$hostbreadth[which(PH$cat.genus.species=="Geometrid.16")] <- "unk"
PH$hostbreadth[which(PH$cat.genus.species=="Geometrid.18")] <- "unk"
PH$hostbreadth[which(PH$cat.genus.species=="Geometrid.19")] <- "unk"
PH$hostbreadth[which(PH$cat.genus.species=="Geometrid.22")] <- "unk"
PH$hostbreadth[which(PH$cat.genus.species=="Geometrid.27")] <- "unk"
PH$hostbreadth[which(PH$cat.genus.species=="Geometrid.39")] <- "unk"
PH$hostbreadth[which(PH$cat.genus.species=="Geometrid.40")] <- "unk"
PH$hostbreadth[which(PH$cat.genus.species=="Geometrid.47")] <- "unk"
PH$hostbreadth[which(PH$cat.genus.species=="Noctuid.11")] <- "unk"
PH$hostbreadth[which(PH$cat.genus.species=="Noctuid.17")] <- "unk"
PH$hostbreadth[which(PH$cat.genus.species=="Noctuid.19")] <- "unk"

# Get independent abundance estimates for the lower trophic level (host plants)
temp1 <- merge(plants, soil, by="transect")
temp2 <- subset(temp1, species != "wislizenii")
temp3 <- subset(temp2, select=c("plant.id", "genus", "species", "soil")); temp3$abundance <- 1

s <- subset(temp3, soil == "s")
ns <- subset(temp3, soil == "ns")

temp.s <- aggregate(abundance ~ genus + soil, data=s, length)
temp.ns <- aggregate(abundance ~ genus + soil, data=ns, length)

serpentine <- temp.s$abundance # counts of each genus
names(serpentine) <- temp.s$genus
nonserpentine <- temp.ns$abundance # counts of each genus
names(nonserpentine) <- temp.ns$genus

# Make matrices from serpentine and non-serpentine slices of the data to build networks
PH$abu <- 1

serp.net3 <- subset(PH, soil == "s")
serp.net2 <- aggregate(abu ~ cat.genus.species + genus, data=serp.net3, length)
serp.net1 <- subset(serp.net2, select=c("genus", "cat.genus.species", "abu"))
serp.net <- matrify(serp.net1)

nonserp.net3 <- subset(PH, soil == "ns")
nonserp.net2 <- aggregate(abu ~ cat.genus.species + genus, data=nonserp.net3, length)
nonserp.net1 <- subset(nonserp.net2, select=c("genus", "cat.genus.species", "abu"))
nonserp.net <- matrify(nonserp.net1)

# Hostbreadth slice to use later
hostbreadth <- aggregate(abu ~ cat.genus.species + hostbreadth, data=PH, length); hostbreadth$abu <- NULL
hostbreadth$Genus.species <- hostbreadth$cat.genus.species
hostbreadth$host.breadth <- "Generalist"
hostbreadth$host.breadth[which(hostbreadth$hostbreadth=="spec")] <- "Specialist"
hostbreadth$host.breadth[which(hostbreadth$hostbreadth=="unk")] <- "unk"
hostbreadth2 <- subset(hostbreadth, select=c("Genus.species", "host.breadth"))

########################
# Whole-soil Networks! #
########################
plotweb(serp.net, 
        abuns.type="independent",
        low.abun = serpentine,   
        method="normal",
        arrow="down.center",
        col.interaction="gray84",
        bor.col.interaction="black",
        low.abun.col="gray16",
        bor.low.abun.col="gray16",
        high.abun.col="gray16",
        bor.high.abun.col="gray16")

# order of serp caterpillar species:
s.species <- data.frame(colnames(serp.net))
s.species$Genus.species <- s.species$colnames.serp.net.
s.species$number_left_to_right <- rownames(s.species); nrow(s.species)
s.species2 <- merge(hostbreadth2, s.species, by="Genus.species"); nrow(s.species2)
s.species3 <- subset(s.species2, select=c("Genus.species", "host.breadth", "number_left_to_right"))

plotweb(nonserp.net, 
        abuns.type="independent",
        low.abun = nonserpentine,  
        method="normal",
        arrow="down.center",
        col.interaction="gray84",
        bor.col.interaction="black",
        low.abun.col="gray16",
        bor.low.abun.col="gray16",
        high.abun.col="gray16",
        bor.high.abun.col="gray16")

# order of nonserp caterpillar species:
ns.species <- data.frame(colnames(nonserp.net))
ns.species$Genus.species <- ns.species$colnames.nonserp.net.
ns.species$number_left_to_right <- rownames(ns.species); nrow(ns.species)
ns.species2 <- merge(hostbreadth2, ns.species, by="Genus.species"); nrow(ns.species2)
ns.species3 <- subset(ns.species2, select=c("Genus.species", "host.breadth", "number_left_to_right"))

#############################################################################################
# What about 'decomposing' each soil network into its specialist and generalist components? #
#############################################################################################

  # # # # # #  
# Serpentine #
 # # # # # # 
# Specialists and unknowns
serp.spec3 <- subset(serp.net3, hostbreadth == "spec" | hostbreadth == "unk")
serp.spec2 <- aggregate(abu ~ cat.genus.species + genus, data=serp.spec3, length)
serp.spec1 <- subset(serp.spec2, select=c("genus", "cat.genus.species", "abu"))
serp.spec <- matrify(serp.spec1)
# Confirmed 'fundamental generalists'
serp.gen3 <- subset(serp.net3, hostbreadth == "gen")
serp.gen2 <- aggregate(abu ~ cat.genus.species + genus, data=serp.gen3, length)
serp.gen1 <- subset(serp.gen2, select=c("genus", "cat.genus.species", "abu"))
serp.gen <- matrify(serp.gen1)

  # # # # # # # # 
# Non-Serpentine #
  # # # # # # # #
# Specialists and unknowns
nonserp.spec3 <- subset(nonserp.net3, hostbreadth == "spec" | hostbreadth == "unk")
nonserp.spec2 <- aggregate(abu ~ cat.genus.species + genus, data=nonserp.spec3, length)
nonserp.spec1 <- subset(nonserp.spec2, select=c("genus", "cat.genus.species", "abu"))
nonserp.spec <- matrify(nonserp.spec1)
# Confirmed 'fundamental generalists'
nonserp.gen3 <- subset(nonserp.net3, hostbreadth == "gen")
nonserp.gen2 <- aggregate(abu ~ cat.genus.species + genus, data=nonserp.gen3, length)
nonserp.gen1 <- subset(nonserp.gen2, select=c("genus", "cat.genus.species", "abu"))
nonserp.gen <- matrify(nonserp.gen1)

#------------------------------------------------------------------#
# Non - Serpentine: specialist web, generalist web, and entire web #
#------------------------------------------------------------------#
# Fundamental specialist (and unknown sp) component
plotweb(nonserp.spec, 
        abuns.type="independent",
        low.abun = nonserpentine,  
        method="normal",
        arrow="down.center",
        col.interaction="gray84",
        bor.col.interaction="black",
        low.abun.col="gray16",
        bor.low.abun.col="gray16",
        high.abun.col="gray16",
        bor.high.abun.col="gray16")

# Order of species across the sub-web
nonserp.spec.sp <- data.frame(colnames(nonserp.spec))
nonserp.spec.sp$Genus.species <- nonserp.spec.sp$colnames.nonserp.spec.
nonserp.spec.sp$number_left_to_right <- rownames(nonserp.spec.sp); nrow(nonserp.spec.sp)
nonserp.spec.sp2 <- merge(hostbreadth2, nonserp.spec.sp, by="Genus.species"); nrow(nonserp.spec.sp2)
nonserp.spec.sp3 <- subset(nonserp.spec.sp2, select=c("Genus.species", "host.breadth", "number_left_to_right"))

# Fundamental generalist component
plotweb(nonserp.gen, 
        abuns.type="independent",
        low.abun = nonserpentine,  
        method="normal",
        arrow="down.center",
        col.interaction="gray84",
        bor.col.interaction="black",
        low.abun.col="gray16",
        bor.low.abun.col="gray16",
        high.abun.col="gray16",
        bor.high.abun.col="gray16")

# Order of species across the sub-web
nonserp.gen.sp <- data.frame(colnames(nonserp.gen))
nonserp.gen.sp$Genus.species <- nonserp.gen.sp$colnames.nonserp.gen.
nonserp.gen.sp$number_left_to_right <- rownames(nonserp.gen.sp); nrow(nonserp.gen.sp)
nonserp.gen.sp2 <- merge(hostbreadth2, nonserp.gen.sp, by="Genus.species"); nrow(nonserp.gen.sp2)
nonserp.gen.sp3 <- subset(nonserp.gen.sp2, select=c("Genus.species", "host.breadth", "number_left_to_right"))

#------------------------------------------------------------#
# Serpentine: specialist web, generalist web, and entire web #
#------------------------------------------------------------#
# Fundamental specialist (and unknown sp) component
plotweb(serp.spec, 
        abuns.type="independent",
        low.abun = serpentine,  
        method="normal",
        arrow="down.center",
        col.interaction="gray84",
        bor.col.interaction="black",
        low.abun.col="gray16",
        bor.low.abun.col="gray16",
        high.abun.col="gray16",
        bor.high.abun.col="gray16")

# Order of species across the sub-web
serp.spec.sp <- data.frame(colnames(serp.spec))
serp.spec.sp$Genus.species <- serp.spec.sp$colnames.serp.spec.
serp.spec.sp$number_left_to_right <- rownames(serp.spec.sp); nrow(serp.spec.sp)
serp.spec.sp2 <- merge(hostbreadth2, serp.spec.sp, by="Genus.species"); nrow(serp.spec.sp2)
serp.spec.sp3 <- subset(serp.spec.sp2, select=c("Genus.species", "host.breadth", "number_left_to_right"))

# Fundamental generalist component
plotweb(serp.gen, 
        abuns.type="independent",
        low.abun = serpentine,  
        method="normal",
        arrow="down.center",
        col.interaction="gray84",
        bor.col.interaction="black",
        low.abun.col="gray16",
        bor.low.abun.col="gray16",
        high.abun.col="gray16",
        bor.high.abun.col="gray16")

# Order of species across the sub-web
serp.gen.sp <- data.frame(colnames(serp.gen))
serp.gen.sp$Genus.species <- serp.gen.sp$colnames.serp.gen.
serp.gen.sp$number_left_to_right <- rownames(serp.gen.sp); nrow(serp.gen.sp)
serp.gen.sp2 <- merge(hostbreadth2, serp.gen.sp, by="Genus.species"); nrow(serp.gen.sp2)
serp.gen.sp3 <- subset(serp.gen.sp2, select=c("Genus.species", "host.breadth", "number_left_to_right"))


