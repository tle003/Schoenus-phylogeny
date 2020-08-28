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
results_DEC<-readRDS(file = "results_DEC_constrained.rds")
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

plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE)

# Paint on the branch states
paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, colors_list_for_states=colors_list_for_states, lwd=4, lty=par("lty"), root.edge=TRUE, stratified=stratified)

plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=TRUE)

# [RvM] Save-out ancestral areas' relative probs at nodes ----------------------
# (Written by RvM)

library(tidyverse)
library(tidytree)
library(ggtree)

# NOTE: resmod = stochastic_map_states_into_res(res)
# where res = results_DEC
# where results_DEC = "results_DEC_constrained.rds"

states_relative_probs_for_nodes <- as.data.frame(
  resmod$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS
)

# Check:
(Nnode(tr) + Ntip(tr)) == nrow(states_relative_probs_for_nodes)

nstates <- ncol(states_relative_probs_for_nodes)

# Replace column names in state probability matrix with state names
for (i in 1:nstates) {
  colnames(states_relative_probs_for_nodes)[[i]] <- paste(collapse = "_",
    areas[states_list_0based[[i]] + 1]
  )
}
colnames(states_relative_probs_for_nodes)[[1]] <- "na"

# Save state probability matrix
write.csv(
  states_relative_probs_for_nodes,
  "ancestral_areas_relative_probs.csv"
)

# Tidy state probability matrix
states_relative_probs_for_nodes_tidy <- states_relative_probs_for_nodes %>%
  as_tibble(rownames = "label") %>%
  gather(state, relative_prob, -label) %>%
  group_by(label) %>%
  # Only keep the most probable state
  arrange(desc(relative_prob)) %>%
  slice(1)

# Save tidied state probability matrix
write.csv(
  states_relative_probs_for_nodes_tidy,
  "ancestral_areas_relative_probs_tidy.csv"
)

# Combine node ancestral area data with phylogeny proper
tr_w_ancestral_areas <- tr %>%
  as_tibble() %>%
  full_join(states_relative_probs_for_nodes_tidy) %>%
  as.treedata()

# FIXME:
#treeio::write.beast(tr_w_ancestral_areas, "tr_w_ancestral_areas.tre")

# FIXME:
#plotTree(tr_w_ancestral_areas@phylo)
# (But this works??: `plot(tr_w_ancestral_areas@phylo)`)

# Check that as.treedata() didn't break the phylogeny?
tr2 <- as.phylo(tr_w_ancestral_areas@phylo)
tr2$node.label <- NULL
plotTree(tr2$edge.length)
write.tree(tr2, "foo.tre")
tr3 <- read.tree("foo.tre")
plotTree(tr3)
# Nope!

tr_w_ancestral_areas2 <- tr3 %>%
  force.ultrametric(method = "extend") %>%
  as_tibble() %>%
  full_join(states_relative_probs_for_nodes_tidy) %>%
  as.treedata()

# FIXME:
#tree_plot <- ggtree(tr_w_ancestral_areas2)
#ggsave("WIP.pdf", tree_plot, width = 10, height = 10)

# Try plotting with phytools::/ape::/base::?

my_palette <- scales::brewer_pal(palette = "Paired")(
  n = length(c(
    "Cape",
    "Africa",
    "Western Australia",
    "Australia",
    "New Zealand",
    "Neotropics",
    "Pacific",
    "Tropical Asia",
    "Holarctic"
  ))
)
# Darken purple
my_palette[[9]] <- "#AB71C7"

op <- par()
par(mar = c(0, 0, 0, 0))
pdf("WIP.pdf", width = 10, height = 10)
plot(ladderize(tr3), cex = 0.5)
nodelabels(
  pie = as.matrix(states_relative_probs_for_nodes),
  piecol = my_palette,  # FIXME: palette too small, no. states != no. regions
  cex = 0.5
)
dev.off()
par(op)
