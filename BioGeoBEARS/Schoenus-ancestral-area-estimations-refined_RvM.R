# Tammy Elliott & Ruan van Mazijk, 2020

# [RvM] Setup ------------------------------------------------------------------

setwd("BioGeoBEARS")  # [Rvm]

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

# [RvM] BSM analysis -----------------------------------------------------------

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
results_DEC<-readRDS(file = "results_DEC_constrained-reparamaterization.rds")
results_DEC$inputs$wd <- getwd()  # [Rvm]
results_DEC$inputs$geogfn <- "Schoeneae-DEC-9areas.txt"  # [Rvm]
results_DEC$inputs$trfn <- "Schoeneae_tree_ultrametric.tre"  # [Rvm]
res<-results_DEC

tree_file_name <- np("Schoeneae_tree_ultrametric.tre")
tree_file_name

tr <- read.tree(tree_file_name)

#extdata loaad
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir

geo_file_name <- np("Schoeneae-DEC-9areas.txt")
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=geo_file_name)
tipranges



BSM_inputs_fn = "BSM_inputs_file.Rdata"
runInputsSlow = TRUE
if (runInputsSlow)
    {
    stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
    save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
    } else {
    # Loads to "stochastic_mapping_inputs_list"
    load(BSM_inputs_fn)
    } # END if (runInputsSlow)

# Check inputs (doesn't work the same on unconstr)
names(stochastic_mapping_inputs_list)
stochastic_mapping_inputs_list$phy2
stochastic_mapping_inputs_list$COO_weights_columnar
stochastic_mapping_inputs_list$unconstr
set.seed(seed=as.numeric(Sys.time()))

runBSMslow = TRUE
if (runBSMslow == TRUE)
    {
    # Saves to: RES_clado_events_tables.Rdata
    # Saves to: RES_ana_events_tables.Rdata
    BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=150, nummaps_goal=100, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)

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



############################################
# Setup for painting a single stochastic map
############################################
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = TRUE
clado_events_table = clado_events_tables[[1]]
ana_events_table = ana_events_tables[[1]]

# [RvM] BSM results ------------------------------------------------------------

# Convert the BSM into a modified res object
master_table_cladogenetic_events = clado_events_tables[[1]]
resmod = stochastic_map_states_into_res(res=res, master_table_cladogenetic_events=master_table_cladogenetic_events, stratified=stratified)

#plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE)

# Paint on the branch states
#paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, colors_list_for_states=colors_list_for_states, lwd=4, lty=par("lty"), root.edge=TRUE, stratified=stratified)

#plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=TRUE)

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
