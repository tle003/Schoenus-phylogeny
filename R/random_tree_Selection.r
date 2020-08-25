## Title: Cyperaceae trees truncated
### Name: Tammy L. Elliott
### Date:August 4, 2020
#### R version 4.02

## Prepare data
rm(list = ls()) #Clear your work environment

#set working directory
setwd("/Users/tammy/PostDoc/Smuts/Schoeneae/Phylogeny/manuscript/RAxML-msc-sangor-prelim/concatenated/Bayesian-analysis/Bayesian-run-Cyperaceae-6-calib-final-12-AUG")


#Load relevant libraries
library(ape)
library(phytools)

#Read in tree
Cyperaceae.trees<-read.nexus("Cyperaceae-all-taxa-6calib-Aug12-comb.trees")



#randomly select 100 of the 18,018 trees 
Cyperaceae.trees.100<-sample(Cyperaceae.trees, 100, replace = FALSE, prob = NULL)


#Write the trees (100)
#write.tree(Cyperaceae.trees.100, file = "Cyperaceae.trees.100.trees", append = FALSE,
 #          digits = 10, tree.names = FALSE)


#make sure you can read them in
Cyperaceae.trees.100<-read.tree("Cyperaceae.trees.100.trees")

#Prune to Schoenus only
Schoenus.pruned.1<-lapply(Cyperaceae.trees.100,drop.tip,tip=c("Xyroschoenus_hornei","Gymnoschoenus_sphaerocephalus"))
"multiPhylo"->class(Schoenus.pruned.1)
Schoenus.pruned<-lapply(Schoenus.pruned.1, function(x) drop.tip(x,x$tip.label[!grepl("Schoenus", x$tip.label)]))
"multiPhylo"->class(Schoenus.pruned)

#Write the trees (100)
#write.tree(Schoenus.pruned, file = "Schoenus.trees.100.trees", append = FALSE,
 #          digits = 10, tree.names = FALSE)

#Schoeneae trees
Schoeneae_trees_prelim<-lapply(Cyperaceae.trees.100, drop.tip, tip=c("Schoenoplectus_tabernaemontani","Bulbostylis_barbata","Fimbristylis_complanata","Cyperus_papyrus","Calliscirpus_criniger","Calliscirpus_brachythrix",
    "Carpha_glomerata","Cladium_mariscus","Lagenocarpus_alboniger","Rhynchospora_rugosa","Mapania_cuspidata","Fuirena_hirsuta","Ficinia_indica","Carex_siderosticta","Dulichium_arundinaceum",
    "Eriophorum_vaginatum","Carex_rupestris","Blysmus_compressus","Scirpus_hattorianus","Rhynchospora_capitellata","Sumatroscirpus_rupestris", "Khaosokia_caricoides",
    "Scleria_brownii","Scleria_rugosa","Fimbristylis_ovata", "Scirpus_pendulus","Hypolytrum_nemorum","Eleocharis_acicularis", "Eleocharis_palustris","Cladium_chinense",
    "Diplacrum_caricinum","Hypolytrum_africanum","Exochogyne_amazonica","Calyptrocarya_glomerulata","Trianoptiles_capensis"))
"multiPhylo"->class(Schoeneae_trees_prelim)
#Save tree file
#write.tree(Schoeneae_trees_prelim, file = "Schoeneae.trees.100.trees", append = FALSE,digits = 10)