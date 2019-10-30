# PH-soil-networks-2019
Analyses of plant-herbivore ecological networks from serpentine and non-serpentine soils. Per-species abundance and H' analyses; empirical networks, subsampled networks, and null models

*Raw Data
2014-2015_network.csv
plants.csv
soil.csv

*Randomization outputs: subsampling code takes a long time to run. These are the outputs used in MS.  
emp_q&q_max.csv
- These are the empirical Q values for each network, and the normalized Q values. In this case, q.max was gotten from 100 different "modularized" versions of the web (where "modularization" again means randomly assigning each species & all its interactions to a single host plant

Q_from.subsampling.csv
- Subsampled networks & calculated Q. We subsampled networks as follows: for each network, we created intervals of network size ranging from 25 interactions to the total number of interactions in the network. At each interval, we randomly subsampled that number of interactions from each network and calculated modularity of the resulting network (=mod.max). Modularity was calculated as described in main text methods (function “DIRTLPAwb+”, Beckett [2016]). We iterated this for 100 random networks at each size interval.

NormalizedQ_from.subsampling.csv
- For this, we did a similar procedure as above, but for each size interval (of 5) we made only 10 random subwebs (="rep.q"). For each of those subwebs, we then "modularizated" the subweb - that is, aggregated the total # of interactions per herbivore species in the web, and then randomly assigned each one (and all of its interactions) to one of the 4 plant genera. This keeps the total number of interactions/herbivore species the same, but makes all of the species specialists (removes all edges between modules). Then modularity was calculated for this modularized version. We made 15 "modularized" versions of each random subweb, and took the maximum modularity value across those 15 webs (= max.q) to represent the maxumum modularity possible for that random subweb. Therefore, the theoretical maximum modularity for each size interval is the maximum of the 10 max.q values, per size.

*Script name followed by contents

1. Abundance of generalists & specialists across soils.R
      - Figure 3
      - Figure S5

2. H' of generalists across soils.R
      - Figure 4
      - Figure S3 & S4
      - Table S6 & S7
      
3. Sampling completeness.R
      - Table S4
      
4. Bipartite networks.R
      - Figure 5
      - Figure S5
