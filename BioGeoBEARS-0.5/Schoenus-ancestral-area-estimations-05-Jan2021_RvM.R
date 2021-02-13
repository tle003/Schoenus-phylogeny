## Title:  Schoeneae ancestoral areas - refined
### Name: Tammy Elliott
### Date: June 18, 2020


# Load the required libraries.
#install.packages("rexpokit")
#install.packages("cladoRcpp")
#install.packages("devtools")
#library(devtools)
#devtools::install_github(repo="nmatzke/BioGeoBEARS")
# Load the package (after installation, see above).
library(GenSA)    # GenSA is better than optimx (although somewhat slower)
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(rexpokit)
library(cladoRcpp)
library(devtools)
library(BioGeoBEARS)
library(ape)
library(FossilSim)
library(phangorn)
library(phytools)

#Get August bug fix
#devtools::install_github(repo="nmatzke/BioGeoBEARS", INSTALL_opts="--byte-compile")

#set working directory
setwd("BioGeoBEARS-0.5")  # [RvM]
getwd()




#Check timeperiods file
timeperiods<-read.table("timeperiods.txt")
dispersal.multipliers<-read.table("dispersal_multipliers-nums.txt",header = TRUE, fill = TRUE)
#Write .txt file
#write.table(Schoenus.cont, file = "Schoenus.cont.txt", quote=F)



#load geography data
#extdata loaad
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir

geo_file_name <- np("Schoeneae-DEC-9areas.txt")
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=geo_file_name)
tipranges

# Maximum range size observed:
max(rowSums(dfnums_to_numeric(tipranges@df)))
#Format Schoeneae tree
#Cyperaceae_tree<-read.nexus("/Users/tammy/PostDoc/Smuts/Schoeneae/Phylogeny/manuscript/RAxML-msc-sangor-prelim/Bayesian-results-Jan2021/Longer_runs/Cyperaceae-all-taxa-6-calib-st6-Dec26.tre")
#Schoeneae_tree_prelim<-drop.tip(Cyperaceae_tree, c("Schoenoplectus_tabernaemontani","Bulbostylis_barbata","Fimbristylis_complanata","Cyperus_papyrus","Calliscirpus_criniger","Calliscirpus_brachythrix",
#    "Carpha_glomerata","Cladium_mariscus","Lagenocarpus_alboniger","Rhynchospora_rugosa","Mapania_cuspidata","Fuirena_hirsuta","Ficinia_indica","Carex_siderosticta","Dulichium_arundinaceum",
#    "Eriophorum_vaginatum","Carex_rupestris","Blysmus_compressus","Scirpus_hattorianus","Rhynchospora_capitellata","Sumatroscirpus_rupestris", "Khaosokia_caricoides",
#    "Scleria_brownii","Scleria_rugosa","Fimbristylis_ovata", "Scirpus_pendulus","Hypolytrum_nemorum","Eleocharis_acicularis", "Eleocharis_palustris","Cladium_chinense",
#    "Diplacrum_caricinum","Hypolytrum_africanum","Exochogyne_amazonica","Calyptrocarya_glomerulata","Trianoptiles_capensis"))
#is.ultrametric(Schoeneae_tree_prelim)
#Make tree ultrametric
#Schoeneae_tree_ultrametric<-force.ultrametric(Schoeneae_tree_prelim)
#is.ultrametric(Schoeneae_tree_ultrametric)
#Save tree file
#write.tree(Schoeneae_tree_ultrametric, file = "/Users/tammy/PostDoc/Smuts/Schoeneae/Phylogeny/manuscript/RAxML-msc-sangor-prelim/Bayesian-results-Jan2021/Longer_runs/Schoeneae_tree_ultrametric-Jan2021.tre", append = FALSE,
 #         digits = 10)

#Read tree
Schoeneae_tree_ultrametric<-read.tree("Schoeneae_tree_ultrametric-Jan2021.tre")  # [RvM]

#Impose mimimum branch lengths
#Schoeneae_tree_min_branch<-impose_min_brlen(Schoeneae_tree_ultrametric, min_brlen = 0.01, leave_BL0_terminals = TRUE,
#    direct_ancestor_brlen = 1e-07, printlevel = 2)

#######################################################
# This is the example Newick file for Hawaiian Schoeneae
# (from Ree & Smith 2008)
# "trfn" = "tree file name"
trfn = "Schoeneae_tree_ultrametric-Jan2021.tre"  # [RvM]

# Look at the raw Newick file:
moref(trfn)

# Look at your phylogeny (plots to a PDF, which avoids issues with multiple graphics in same window):
pdffn = "tree.pdf"
pdf(file=pdffn, width=9, height=12)

tr = read.tree(trfn)
tr
plot(tr)
title("Schoeneae phylogeny")
nodelabels(cex=0.4)
axisPhylo() # plots timescale

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)


numstates_from_numareas(numareas=10, maxareas=2, include_null_range=TRUE)
####################################################
####################################################


######################
#DEC analysis

#Read tree;
tree_file_name <- np("Schoeneae_tree_ultrametric-Jan2021.tre")  # [RvM]
tree_file_name

tr <- read.tree(tree_file_name)
plot(tr, cex=0.3)
nodelabels(cex=0.4)

tree<-tr


#Get root age
tree.max(tree, root.edge = FALSE)

#Get tree table to get root ages
prt(tree)


####################################################
####################################################
# KEY HINT: The number of states (= number of different possible geographic ranges)
# depends on (a) the number of areas and (b) max_range_size.
# If you have more than about 500-600 states, the calculations will get REALLY slow,
# since the program has to exponentiate a matrix of e.g. 600x600.  Often the computer
# will just sit there and crunch, and never get through the calculation of the first
# likelihood.
#
# (this is also what is usually happening when LAGRANGE hangs: you have too many states!)
#
# To check the number of states for a given number of ranges, try:
numstates_from_numareas(numareas=4, maxareas=4, include_null_range=TRUE)
numstates_from_numareas(numareas=4, maxareas=4, include_null_range=FALSE)
numstates_from_numareas(numareas=4, maxareas=3, include_null_range=TRUE)
numstates_from_numareas(numareas=4, maxareas=2, include_null_range=TRUE)

# Large numbers of areas have problems:
numstates_from_numareas(numareas=9, maxareas=6, include_null_range=TRUE)

# ...unless you limit the max_range_size:


#Set-up BioGeoBEARS parametres
BioGeoBEARS_run_object <- define_BioGeoBEARS_run()

# Set tree and data set objects
BioGeoBEARS_run_object$trfn <- tree_file_name
BioGeoBEARS_run_object$geogfn <- geo_file_name

# Set parametres
BioGeoBEARS_run_object$max_range_size = 6
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE

# Set up a time-stratified analysis:
BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "dispersal_multipliers-nums.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.


# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 8
# (use more cores to speed it up; this requires
# library(parallel) and/or library(snow). The package "parallel"
# is now default on Macs in R 3.0+, but apparently still
# has to be typed on some Windows machines. Note: apparently
# parallel works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC model
# (nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# For a slow analysis, run once, then set runslow=FALSE to just
# load the saved result.
#runslow = TRUE
#resfn = "Schoeneae_DEC_M0_unconstrained_v1.Rdata"
#if (runslow)
    #{
    #res = bears_optim_run(BioGeoBEARS_run_object)
    #res

    #save(res, file=resfn)
    #resDEC = res
   # } else {
    # Loads to "res"
   # load(resfn)
  #  resDEC = res
 #   }

#Run analysis
#results_DEC <- bears_optim_run(BioGeoBEARS_run_object)  # [RvM]

#saveRDS(results_DEC, file = "results_DEC_constrained.rds")
#saveRDS(results_DEC, file = "results_DEC_constrained-reparamaterization.rds")
#saveRDS(results_DEC, file = "results_DEC_constrained-constrained05.rds")
#saveRDS(results_DEC, file = "results_DEC_constrained-constrained05-Jan2021.rds")
#######################################################
# PDF plots
#######################################################
#results_DEC<-readRDS(file = "results_DEC_constrained.rds")

results_DEC<-readRDS(file = "results_DEC_constrained-constrained05-Jan2021.rds")

#Plot ancestral states
#Plot the DEC and DEC+J models

tree_file_name <- np("Schoeneae_tree_ultrametric-Jan2021.tre")  # [RvM]
tree_file_name

tr <- read.tree(tree_file_name)

#extdata loaad
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir

geo_file_name <- np("Schoeneae-DEC-9areas.txt")
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=geo_file_name)
tipranges

pdffn = "Schoeneae-constrained.pdf"
pdf(pdffn, width=6, height=6)
analysis_titletxt = "DEC on Schoeneae"

scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
res2 = plot_BioGeoBEARS_results(results_DEC, analysis_titletxt, addl_params=list("j"), plotwhat="text",
label.offset=0.45, tipcex=0.4, statecex=0.4, splitcex=0.35, titlecex=0.8, plotsplits=TRUE,
cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

plot_BioGeoBEARS_results(results_DEC, analysis_titletxt, addl_params=list("j"), plotwhat="pie",
label.offset=0.45, tipcex=0.4, statecex=0.4, splitcex=0.35, titlecex=0.8, plotsplits=TRUE,
cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)


dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)



######################################################
# Get the inputs for Biogeographical Stochastic Mapping
# Note: this can be slow for large state spaces and trees, since
# the independent likelihoods for each branch are being pre-calculated
# E.g., for 10 areas, this requires calculation of a 1024x1024 matrix
# for each branch.  On a tree with ~800 tips and thus ~1600 branches, this was about 1.6 gigs
# for storage of "BSM_inputs_file.Rdata".
# Update: 2015-09-23 -- now, if you used multicore functionality for the ML analysis,
# the same settings will be used for get_inputs_for_stochastic_mapping().
#######################################################
#Read in data
results_DEC<-readRDS(file = "results_DEC_constrained-constrained05-Jan2021.rds")
results_DEC$inputs$wd <- getwd()  # [Rvm]
results_DEC$inputs$geogfn <- "Schoeneae-DEC-9areas.txt"  # [Rvm]
results_DEC$inputs$trfn <- "Schoeneae_tree_ultrametric.tre"  # [Rvm]
res<-results_DEC

tree_file_name <- np("Schoeneae_tree_ultrametric-Jan2021.tre")  # [RvM]
tree_file_name

tr <- read.tree(tree_file_name)

#extdata loaad
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir

geo_file_name <- np("Schoeneae-DEC-9areas.txt")
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=geo_file_name)
tipranges


# [Rvm]
#BSM_inputs_fn = "BSM_inputs_file.Rdata"
#runInputsSlow = TRUE
#if (runInputsSlow)
#    {
#    stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
#    save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
#    } else {
#    # Loads to "stochastic_mapping_inputs_list"
#    load(BSM_inputs_fn)
#    } # END if (runInputsSlow)
#
## Check inputs (doesn't work the same on unconstr)
#names(stochastic_mapping_inputs_list)
#stochastic_mapping_inputs_list$phy2
#stochastic_mapping_inputs_list$COO_weights_columnar
#stochastic_mapping_inputs_list$unconstr
#set.seed(seed=as.numeric(Sys.time()))
# [/Rvm]

runBSMslow = FALSE  # [RvM]
if (runBSMslow == TRUE)
    {
    # Saves to: RES_clado_events_tables.Rdata
    # Saves to: RES_ana_events_tables.Rdata
    BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=150, nummaps_goal=1000, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)

    RES_clado_events_tables = BSM_output$RES_clado_events_tables
    RES_ana_events_tables = BSM_output$RES_ana_events_tables
    } else {
    # Load previously saved...

    # Loads to: RES_clado_events_tables
    load(file="RES_clado_events_tables.Rdata")
    # Loads to: RES_ana_events_tables
    load(file="RES_ana_events_tables.Rdata")
    BSM_output = NULL
    BSM_output$RES_clado_events_tables = RES_clado_events_tables
    BSM_output$RES_ana_events_tables = RES_ana_events_tables
    } # END if (runBSMslow == TRUE)




# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables
head(clado_events_tables[[1]])
head(ana_events_tables[[1]])
length(clado_events_tables)
length(ana_events_tables)


include_null_range = TRUE
areanames = names(tipranges@df)
areas = areanames
max_range_size = 6

# Note: If you did something to change the states_list from the default given the number of areas, you would
# have to manually make that change here as well! (e.g., areas_allowed matrix, or manual reduction of the states_list)
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=6, include_null_range=include_null_range)

colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=5, plot_null_range=TRUE)

##Start here

############################################
# Setup for painting a single stochastic map
############################################
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = TRUE
clado_events_table = clado_events_tables[[1]]
ana_events_table = ana_events_tables[[1]]

# cols_to_get = names(clado_events_table[,-ncol(clado_events_table)])
# colnums = match(cols_to_get, names(ana_events_table))
# ana_events_table_cols_to_add = ana_events_table[,colnums]
# anagenetic_events_txt_below_node = rep("none", nrow(ana_events_table_cols_to_add))
# ana_events_table_cols_to_add = cbind(ana_events_table_cols_to_add, anagenetic_events_txt_below_node)
# rows_to_get_TF = ana_events_table_cols_to_add$node <= length(tr$tip.label)
# master_table_cladogenetic_events = rbind(ana_events_table_cols_to_add[rows_to_get_TF,], clado_events_table)

############################################
# Open a PDF
############################################
#pdffn = paste0("DEC", "Schoeneae_single_stochastic_map_n1.pdf")  # [RvM]
#pdf(file=pdffn, width=6, height=6)  # [RvM]

# Convert the BSM into a modified res object
master_table_cladogenetic_events = clado_events_tables[[1]]
resmod = stochastic_map_states_into_res(res=res, master_table_cladogenetic_events=master_table_cladogenetic_events, stratified=stratified)

#plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE)  # [RvM]

# Paint on the branch states
#paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, colors_list_for_states=colors_list_for_states, lwd=4, lty=par("lty"), root.edge=TRUE, stratified=stratified)  # [RvM]

#plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=TRUE)  # [RvM]

############################################
# Close PDF
############################################
#dev.off()  # [RvM]
#cmdstr = paste("open ", pdffn, sep="")  # [RvM]
#system(cmdstr)  # [RvM]

# [RvM] Save-out ancestral areas' relative probs at nodes ----------------------
# (Written by RvM)

library(tidyverse)
library(tidytree)
library(ggtree)

# NOTE: resmod = stochastic_map_states_into_res(res)
# where res = results_DEC
# where results_DEC = "results_DEC_constrained.rds"

nodes_state_probs <- as.data.frame(
  resmod$ML_marginal_prob_each_state_at_branch_top_AT_node
)

# Remove the tip-nodes from the matrix
nodes_state_probs <-  nodes_state_probs[-(1:Ntip(tr)), ]

# Check:
Nnode(tr) == nrow(nodes_state_probs)

nstates <- ncol(nodes_state_probs)

# Replace column names in state probability matrix with state names
for (i in 1:nstates) {
  colnames(nodes_state_probs)[[i]] <- paste(collapse = "_",
    areas[states_list_0based[[i]] + 1]  # because 9x areas indexed 0:8
  )
}
colnames(nodes_state_probs)[[1]] <- "na"

# Save state probability matrix
write.csv(
  nodes_state_probs,
  "ancestral_areas_relative_probs.csv"
)

# Tidy state probability matrix
nodes_state_probs_tidy <- nodes_state_probs %>%
  as_tibble(rownames = "node") %>%
  gather(state, prob, -node) %>%
  group_by(node) %>%
  # Only keep the most probable state
  arrange(node, desc(prob)) %>%
  slice(1) %>%
  mutate(node = as.numeric(node))

# Save tidied state probability matrix
write.csv(
  nodes_state_probs_tidy,
  "ancestral_areas_relative_probs_tidy.csv"
)

# Combine node ancestral area data with phylogeny proper
tr_w_ancestral_areas <- tr %>%
  as_tibble() %>%
  full_join(nodes_state_probs_tidy) %>%
  as.treedata()

# Save tree with data as BEAST-style NEXUS-file
treeio::write.beast(tr_w_ancestral_areas, "tr_w_ancestral_areas.tre")
