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

results_DEC <- readRDS("BioGeoBEARS/results_DEC_constrained-constrained05.rds")
tr <- read.tree("BioGeoBEARS/Schoeneae_tree_ultrametric.tre")
tipranges <- getranges_from_LagrangePHYLIP("BioGeoBEARS/Schoeneae-DEC-9areas.txt")

nodes_state_probs <- as.data.frame(
  results_DEC$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS
)

# Remove the tip-nodes from the matrix
nodes_state_probs <- nodes_state_probs[-(1:Ntip(tr)), ]

# Check:
Nnode(tr) == nrow(nodes_state_probs)

include_null_range = TRUE
areanames = names(tipranges@df)
areas = areanames
max_range_size = 6
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=6, include_null_range=include_null_range)

nstates <- ncol(nodes_state_probs)

# Replace column names in state probability matrix with state names
for (i in 1:nstates) {
  colnames(nodes_state_probs)[[i]] <- paste(collapse = "_",
    areas[states_list_0based[[i]] + 1]  # because 9x areas indexed 0:8
  )
}
colnames(nodes_state_probs)[[1]] <- "na"

# Tidy state probability matrix
nodes_state_probs_tidy <- nodes_state_probs %>%
  as_tibble(rownames = "node") %>%
  gather(state, prob, -node) %>%
  group_by(node) %>%
  arrange(node, desc(prob)) %>%
  mutate(node = as.numeric(node))

relprobs_matrix <-
  results_DEC$ML_marginal_prob_each_state_at_branch_top_AT_node
ranges_list <- states_list_0based_to_ranges_txt_list(
  state_indices_0based = states_list_0based,
  areanames            = areas
)
statenames <- unlist(ranges_list)
MLstates <- get_ML_states_from_relprobs(
  relprobs_matrix,
  statenames,
  returnwhat = "states",
  if_ties    = "takefirst"
)
possible_ranges_list_txt <-areas_list_to_states_list_new(
  areas,
  maxareas           = max_range_size,
  split_ABC          = FALSE,
  include_null_range = results_DEC$inputs$include_null_range
)
colors_list_for_states <- mix_colors_for_states(
  get_colors_for_numareas(length(areas)),
  rcpp_areas_list_to_states_list(
    areas,
    maxareas = max_range_size,
    include_null_range = results_DEC$inputs$include_null_range
  ),
  plot_null_range = results_DEC$inputs$include_null_range
)
cols_byNode <- rangestxt_to_colors(
  possible_ranges_list_txt,
  colors_list_for_states,
  MLstates
)
colnums_to_keep_in_probs = NULL
probs = results_DEC$ML_marginal_prob_each_state_at_branch_top_AT_node
probs2 = probs
maxprob = rep(0, nrow(probs))
other = rep(0, nrow(probs))
num_to_keep = 1
cat(
  "\nSince simplify_piecharts==TRUE, reducing prob pie charts to (most probable, other)...\n"
)
nodes = (length(ladderize(tr, right = FALSE)$tip.label) + 1):(length(ladderize(tr, right = FALSE)$tip.label) +
  ladderize(tr, right = FALSE)$Nnode)
tips = 1:length(ladderize(tr, right = FALSE)$tip.label)
for (i in 1:nrow(probs)) {
  cat(i, " ", sep = "")
  tmprow = probs[i,]
  positions_highest_prob_to_lowest = rev(order(tmprow))
  positions_to_keep = positions_highest_prob_to_lowest[1:num_to_keep]
  colnums_to_keep_in_probs = c(colnums_to_keep_in_probs,
                               positions_to_keep)
  keepTF = rep(FALSE, length(tmprow))
  keepTF[positions_to_keep] = TRUE
  otherTF = keepTF == FALSE
  other[i] = sum(tmprow[otherTF])
  tmprow[otherTF] = 0
  probs2[i,] = tmprow
}
cat("\n")
colnums_to_keep_in_probs_in_order = sort(unique(colnums_to_keep_in_probs))
probs3 = cbind(probs2[, colnums_to_keep_in_probs_in_order], other)
probs3 = probs3[nodes,]
newcols = c(
  colors_list_for_states[colnums_to_keep_in_probs_in_order],
  "white"
)

pdf("Schoeneae-constrained_phytools.pdf", width = 10, height = 15)
plotTree(
  ladderize(tr, right = FALSE),
  fsize  = 0.6,
  ftype  = "i",
  offset = 2.25,
  mar    = c(4, 0.1, 0, 0.1)
)
axisPhylo()
title(xlab = "Ma")
nodelabels(
  text = MLstates[nodes],
  node = nodes,
  bg   = cols_byNode[nodes],
  cex  = 0.5
)
tiplabels(
  text   = MLstates[tips],
  tip    = tips,
  bg     = cols_byNode[tips],
  cex    = 0.5,
  offset = 2.25
)
dev.off()

pdf("Schoeneae-constrained.pdf", width = 10, height = 15)
plot_BioGeoBEARS_results(results_DEC,
  tr           = tr,
  tipranges    = tipranges,
  plotwhat     = "pie",
  addl_params  = list("j"),
  label.offset = 0.5,
  tipcex       = 0.5,
  statecex     = 0.5,
  splitcex     = 0.5,
  plotsplits   = TRUE,
  cornercoords_loc =
    np(system.file("extdata/a_scripts", package = "BioGeoBEARS"))
)
dev.off()
