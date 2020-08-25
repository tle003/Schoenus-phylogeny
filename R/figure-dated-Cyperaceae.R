# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(tidyverse)
library(phytools)
library(ggtree)   # Multi-phylo plots
                  # (Installed with BiocManager::install("ggtree"))
library(treeio)   # For ::read.beast()
library(jntools)  # For ::get_tips_in_ape_plot_order()
library(lemon)    # For ::coord_capped_cart()
library(strap)    # For ::geoscalePhylo()

# Import data ------------------------------------------------------------------

MCC_tree <- read.beast("data/phylogenies/Cyperaceae-all-taxa-6calib-max-clad-AUG12.tre")

subtribes <- read_csv("data/Schoeneae-subtribes.csv")

# Define helper functions ------------------------------------------------------

find_node <- function(tree, tip_pattern,
                            additional_taxa = NULL,
                            taxa_to_exclude = NULL) {
  pattern_matches <- tree$tip.label[str_detect(tree$tip.label, tip_pattern)]
  taxa <- c(pattern_matches, additional_taxa)
  taxa <- taxa[!(taxa %in% taxa_to_exclude)]
  getMRCA(tree, taxa)
}

vector2regex <- function(...) {
  x <- c(...)
  paste0("(",
    paste(x, collapse = "|"),
  ").+")
}

find_subtribe <- function(tree,
                          subtribe_name, subtribes_df,
                          additional_taxa_to_exclude = NULL) {
  genera <- subtribes_df %>%
    filter(subtribe == subtribe_name, is.na(species)) %>%
    pull(taxa)
  additional_taxa <- subtribes_df %>%
    filter(subtribe == subtribe_name, !is.na(species)) %>%
    pull(taxa)
  taxa_to_exclude <- subtribes_df %>%
    filter(subtribe != subtribe_name, genus %in% genera) %>%
    pull(taxa)
  find_node(tree,
    tip_pattern = vector2regex(genera),
    additional_taxa,
    c(taxa_to_exclude, additional_taxa_to_exclude)
  )
}

# Tidy data --------------------------------------------------------------------

MCC_tree@phylo <- force.ultrametric(MCC_tree@phylo, method = "extend")
MCC_tree@phylo$tip.label <- str_replace(MCC_tree@phylo$tip.label, "_", " ")

# Prune tree to Schoeneae-only
Schoeneae_node <-
  find_node(MCC_tree@phylo, "Schoenus", "Gymnoschoenus sphaerocephalus")
Schoenoid_taxa <- treeio::offspring(MCC_tree, Schoeneae_node)
non_Schoenoid_taxa <- MCC_tree@data$node %>%
  {.[!(. %in% Schoenoid_taxa)]} %>%
  as.numeric()
Schoeneae_tree <- drop.tip(MCC_tree, non_Schoenoid_taxa)

# Export tree with spaces in tip names for FigTree
Schoeneae_tree2 <- Schoeneae_tree
Schoeneae_tree2@phylo$tip.label <- paste0("'",
  Schoeneae_tree2@phylo$tip.label,
"'")
write.beast(Schoeneae_tree2, "data/phylogenies/Cyperaceae-all-taxa-6calib-max-clad-AUG12_Schoeneae.tre")

Schoeneae_tree3 <- Schoeneae_tree
Schoeneae_tree3@phylo$root.time <- max(nodeHeights(Schoeneae_tree3@phylo))

# Duplicate this strap::geoscalePhylo() but modify so that grey and white boxes
# follow user specified timescale:
my_geoscalePhylo <- function (tree, ages, direction = "rightwards", units = c("Period",
    "Epoch", "Age"), boxes = "Age", tick.scale = "myr", user.scale,
    cex.age = 0.3, cex.ts = 0.3, cex.tip = 0.3, width = 1, label.offset,
    ts.col = TRUE, vers = "ICS2013", x.lim, quat.rm = FALSE,
    erotate, arotate, urotate, ...)
{
    options <- as.list(match.call())
    if (any(names(options) == "type")) {
        if (all(options$type != c("phylogram", "cladogram", "p",
            "c"))) {
            return(cat("type must be either 'phylogram' or 'cladogram'."))
        }
    }
    if (all(direction != c("rightwards", "upwards"))) {
        return(cat("direction must be either 'rightwards' or 'upwards', here set to 'rightwards'."))
    }
    if (is.null(tree$root.time)) {
        return(cat("\n tree$root.time is missing, check tree is time scaled."))
    }
    else {
        root.age <- tree$root.time
    }
    # Remove the following in my version of thus function: ---------------------
    #if (boxes == "User" && any(units != "User")) {
    #    boxes <- "no"
    #}
    #if (all(boxes != units)) {
    #    boxes <- "no"
    #}
    if (tick.scale == "User" && all(units != "User")) {
        tick.scale <- "myr"
    }
    if (missing(ages) == FALSE) {
        ranges <- TRUE
    }
    else {
        ranges <- FALSE
    }
    if (missing(user.scale) & any(units == "User")) {
        units <- units[units != "User"]
        cat("\n user.scale not provided, 'Other' removed from units.")
    }
    if (missing(ages) == FALSE) {
        ages <- ages[tree$tip.label, ]
    }
    if (any(units == "User") & !missing(user.scale)) {
        Midpoint <- matrix(ncol = 1, nrow = length(user.scale[,
            1]))
        Midpoint[, 1] <- (user.scale[, "Start"] + user.scale[,
            "End"])/2
        user.scale <- cbind(user.scale, Midpoint)
    }
    if (all(units != "Age") && boxes == "Age") {
        boxes <- "no"
    }
    units <- paste(toupper(substring(units, 1, 1)), substring(units,
        2), sep = "")
    units[units == "Eonothem"] <- "Eon"
    units[units == "Erathem"] <- "Era"
    units[units == "Series"] <- "Epoch"
    units[units == "System"] <- "Period"
    units[units == "Stage"] <- "Age"
    units <- unique(units)
    boxes[boxes == "Eonothem"] <- "Eon"
    boxes[boxes == "Erathem"] <- "Era"
    boxes[boxes == "Series"] <- "Epoch"
    boxes[boxes == "System"] <- "Period"
    boxes[boxes == "Stage"] <- "Age"
    if (length(units) == 1) {
        ts.width = 0.15
    }
    else if (length(units) == 2) {
        ts.width = 0.2
    }
    else if (length(units) >= 3) {
        ts.width = 0.25
    }
    if (ranges == TRUE && missing(ages) == FALSE) {
        missing.tip.names <- setdiff(tree$tip.label, row.names(ages))
        if (length(missing.tip.names) > 0) {
            cat(paste("\n", missing.tip.names, "not present in ages file, ranges set to FALSE"))
            cat("\n ranges set to FALSE")
            ranges <- FALSE
        }
    }
    timescales <- NULL
    data(timescales, envir = environment())
    timescale <- timescales[[vers]]
    if (quat.rm == TRUE) {
        timescale[(timescale[, "Midpoint"] < 3), "Name"] <- NA
    }
    tscale.data <- matrix(ncol = 3, nrow = 6)
    colnames(tscale.data) <- c("srt", "Depth", "size")
    rownames(tscale.data) <- c("Eon", "Era", "Period", "Epoch",
        "Age", "User")
    if (direction == "upwards") {
        tscale.data[, "srt"] <- c(90, 90, 90, 0, 0, 0)
    }
    else tscale.data[, "srt"] <- c(0, 0, 0, 90, 90, 90)
    tscale.data[, "Depth"] <- c(1, 1, 1, 2, 3.5, 3.5)
    tscale.data[, "size"] <- c(1, 1, 1, 0.8, 0.8, 0.8)
    if (!missing(erotate) && !is.numeric(erotate)) {
        return(cat("\n value for protate must be numeric."))
    }
    if (!missing(arotate) && !is.numeric(arotate)) {
        return(cat("\n value for arotate must be numeric."))
    }
    if (!missing(urotate) && !is.numeric(urotate)) {
        return(cat("\n value for urotate must be numeric."))
    }
    if (!missing(erotate)) {
        tscale.data["Epoch", "srt"] <- erotate
    }
    if (!missing(arotate)) {
        tscale.data["Age", "srt"] <- arotate
    }
    if (!missing(urotate)) {
        tscale.data["User", "srt"] <- urotate
    }
    units <- rownames(tscale.data)[sort(match(units, rownames(tscale.data)),
        decreasing = T)]
    if (!missing(x.lim)) {
        x.lim <- sort(root.age - x.lim)
    }
    else if (ranges == TRUE && !missing(ages) && missing(x.lim)) {
        x.lim <- (root.age - min(ages)) + diff(range(ages)) *
            0.05
    }
    else {
        x.lim <- NULL
    }
    timescale <- timescale[order(timescale[, 1], decreasing = T),
        ]
    timescale.rescaled <- timescale
    timescale.rescaled[, c("Start", "End", "Midpoint")] <- root.age -
        timescale[, c("Start", "End", "Midpoint")]
    first_la <- tree$root.time - dist.nodes(tree)[1, Ntip(tree) +
        1]
    if (ranges == TRUE && missing(ages) == FALSE) {
        offset <- array(dim = length(tree$tip.label), data = 1)
        offset.correction <- diff(range(ages)) * 0.01
        taxon.ranges <- root.age - ages[, c("FAD", "LAD")]
        if (first_la != ages[1, "LAD"]) {
            if (!missing(label.offset)) {
                offset <- array(dim = length(ages[, "FAD"]),
                  data = (ages[, "FAD"] - ages[, "LAD"]) + label.offset)
            }
            else {
                offset <- array(dim = length(ages[, "FAD"]),
                  data = (ages[, "FAD"] - ages[, "LAD"]) + offset.correction)
            }
        }
    }
    else if (!missing(label.offset)) {
        offset <- label.offset
    }
    else {
        offset = 1
    }
    if (tick.scale != "n" | tick.scale != "no") {
        if (tick.scale == "myr" | is.numeric(tick.scale)) {
            scale.ticks = 1
            if (is.numeric(tick.scale)) {
                scale.ages <- tick.scale
            }
            else {
                scale.ages = 10
            }
            tick.position <- root.age - seq(0, 4600, scale.ticks)
            age.name <- seq(0, 4600, scale.ages)
            age.position <- root.age - age.name
            lwd <- c(1, 0.5, 0.5, 0.5, 0.5, 0.7, 0.5, 0.5, 0.5,
                0.5)
            col <- c("black", "grey", "grey", "grey", "grey")
        }
        if (tick.scale != "myr" & is.numeric(tick.scale) == FALSE) {
            age.name <- subset(timescale, timescale[, "Type"] ==
                tick.scale & timescale[, "Source"] == "ICS")
            age.name <- sort(unique(c(age.name[, "Start"], age.name[,
                "End"])))
            age.position <- tick.position <- root.age - age.name
            lwd = 1
            col = "black"
        }
        if (tick.scale == "User") {
            age.name <- sort(unique(c(user.scale[, "Start"],
                user.scale[, "End"])))
            age.position <- tick.position <- root.age - age.name
            lwd = 1
            col = "black"
        }
    }
    if (direction == "upwards") {
        par(fig = c(0, ts.width, 0, 1))
        par(mar = c(3, 1, 2, 5))
        plot.phylo(tree, plot = FALSE, no.margin = T, y.lim = x.lim,
            direction = "upwards", cex = cex.tip, ...)
        timescale.rescaled.names <- timescale.rescaled
        timescale.rescaled.names <- timescale.rescaled.names[timescale.rescaled.names[,
            "End"] > par()$usr[3], ]
        timescale.rescaled.names[timescale.rescaled.names[, "Start"] <
            par()$usr[3], "Start"] <- par()$usr[3]
        if (min(timescale.rescaled.names[, "End"]) < par()$usr[4]) {
            timescale.rescaled.names <- timescale.rescaled.names[timescale.rescaled.names[,
                "Start"] < par()$usr[4], ]
            timescale.rescaled.names[timescale.rescaled.names[,
                "End"] > par()$usr[4], "End"] <- par()$usr[4]
        }
        timescale.rescaled.names[, "Midpoint"] <- (timescale.rescaled.names[,
            "Start"] + timescale.rescaled.names[, "End"])/2
        unit.depths <- tscale.data[units, "Depth"]
        if (tick.scale == "n" | tick.scale == "no") {
            unit.depths <- c(unit.depths, 0.5)
        }
        else if (length(units) <= 3) {
            unit.depths <- c(unit.depths, 2)
        }
        else if (length(units) > 3) {
            unit.depths <- c(unit.depths, 2)
        }
        unit.depths <- cumsum(unit.depths/sum(unit.depths))
        unit.depths <- c(par()$usr[2], par()$usr[2] - (unit.depths *
            (par()$usr[2] - par()$usr[1])))
        depth <- unit.depths[length(unit.depths) - 1] - unit.depths[length(unit.depths)]
        if (tick.scale != "n" && tick.scale != "no") {
            text((unit.depths[length(unit.depths)] + depth *
                0.3), age.position, age.name, cex = cex.age,
                srt = 0)
            segments((unit.depths[length(unit.depths) - 1]),
                tick.position, (unit.depths[length(unit.depths)] +
                  depth * 0.75), tick.position, lwd = lwd, col = col)
        }
        for (t in 1:length(units)) {
            if (units[t] == "User") {
                tscale <- user.scale
                tscale[, c("Start", "End", "Midpoint")] <- root.age -
                  tscale[, c("Start", "End", "Midpoint")]
                tscale.names <- tscale
            }
            else {
                tscale <- subset(timescale.rescaled, timescale.rescaled[,
                  "Type"] == units[t])
                tscale.names <- subset(timescale.rescaled.names,
                  timescale.rescaled.names[, "Type"] == units[t])
            }
            if (ts.col == TRUE & units[t] != "User") {
                rect(unit.depths[t], tscale[, "Start"], unit.depths[t +
                  1], tscale[, "End"], col = rgb(tscale[, "Col_R"],
                  tscale[, "Col_G"], tscale[, "Col_B"], maxColorValue = 255))
            }
            else rect(unit.depths[t], tscale[, "Start"], unit.depths[t +
                1], tscale[, "End"], col = "white")
            text((unit.depths[t] + unit.depths[t + 1])/2, tscale.names[,
                "Midpoint"], tscale.names[, "Name"], cex = cex.ts *
                tscale.data[match(units[t], rownames(tscale.data)),
                  "size"], srt = tscale.data[match(units[t],
                rownames(tscale.data)), "srt"])
        }
        par(fig = c(ts.width, 1, 0, 1), new = T)
        par(mar = c(3, 0, 2, 2))
        plot.phylo(tree, plot = FALSE, no.margin = T, y.lim = x.lim,
            direction = "upwards", cex = cex.tip, ...)
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        if (!missing(boxes) && boxes != "no" && boxes != "n") {
            if (boxes == "User") {
                tscale <- root.age - user.scale[, c("Start",
                  "End")]
            }
            else {
                tscale <- subset(timescale.rescaled, timescale.rescaled[,
                  "Type"] == boxes)
            }
            rect(par()$usr[3], tscale[, "Start"], par()$usr[4],
                tscale[, "End"], col = c("grey90", "white"),
                border = NA)
        }
        par(fig = c(ts.width, 1, 0, 1), new = T)
        par(mar = c(3, 0, 2, 2))
        plot.phylo(tree, label.offset = offset, edge.width = width,
            no.margin = T, y.lim = x.lim, cex = cex.tip, direction = "upwards",
            ...)
        if (ranges == TRUE) {
            par(lend = 1)
            segments(lastPP$xx[c(1:length(tree$tip.label))],
                taxon.ranges[, "FAD"], lastPP$xx[c(1:length(tree$tip.label))],
                taxon.ranges[, "LAD"], col = "black", lwd = width *
                  2)
        }
        if (units[1] == "User") {
            segments(par()$usr[1], min(root.age - user.scale[,
                "Start"]), par()$usr[1], max(root.age - user.scale[,
                "End"]))
        }
        else {
            segments(par()$usr[1], min(timescale.rescaled[timescale.rescaled[,
                "Type"] == units[1], "Start"]), par()$usr[1],
                max(timescale.rescaled[timescale.rescaled[, "Type"] ==
                  units[1], "End"]))
        }
    }
    else {
        par(fig = c(0, 1, 0, ts.width))
        par(mar = c(1, 3, 0, 2))
        plot.phylo(tree, plot = FALSE, no.margin = T, x.lim = x.lim,
            direction = "rightwards", cex = cex.tip, ...)
        timescale.rescaled.names <- timescale.rescaled
        timescale.rescaled.names <- timescale.rescaled.names[timescale.rescaled.names[,
            "End"] > par()$usr[1], ]
        timescale.rescaled.names[timescale.rescaled.names[, "Start"] <
            par()$usr[1], "Start"] <- par()$usr[1]
        if (min(timescale.rescaled.names[, "End"]) < par()$usr[2]) {
            timescale.rescaled.names <- timescale.rescaled.names[timescale.rescaled.names[,
                "Start"] < par()$usr[2], ]
            timescale.rescaled.names[timescale.rescaled.names[,
                "End"] > par()$usr[2], "End"] <- par()$usr[2]
        }
        timescale.rescaled.names[, "Midpoint"] <- (timescale.rescaled.names[,
            "Start"] + timescale.rescaled.names[, "End"])/2
        unit.depths <- tscale.data[units, "Depth"]
        if (tick.scale == "n" & tick.scale == "no") {
            unit.depths <- c(unit.depths, 0.5)
        }
        else if (length(units) <= 3) {
            unit.depths <- c(unit.depths, 2)
        }
        else if (length(units) > 3) {
            unit.depths <- c(unit.depths, 2)
        }
        unit.depths <- cumsum(unit.depths/sum(unit.depths))
        unit.depths <- c(par()$usr[4], par()$usr[4] - (unit.depths *
            (par()$usr[4] - par()$usr[3])))
        depth <- unit.depths[length(unit.depths) - 1] - unit.depths[length(unit.depths)]
        if (tick.scale != "n" && tick.scale != "no") {
            text(age.position, (unit.depths[length(unit.depths)] +
                depth * 0.3), age.name, cex = cex.age, srt = 90)
            segments(tick.position, (unit.depths[length(unit.depths) -
                1]), tick.position, (unit.depths[length(unit.depths)] +
                depth * 0.6), lwd = lwd, col = col)
        }
        for (t in 1:length(units)) {
            if (units[t] == "User") {
                tscale <- user.scale
                tscale[, c("Start", "End", "Midpoint")] <- root.age -
                  tscale[, c("Start", "End", "Midpoint")]
                tscale.names <- tscale
            }
            else {
                tscale <- subset(timescale.rescaled, timescale.rescaled[,
                  "Type"] == units[t])
                tscale.names <- subset(timescale.rescaled.names,
                  timescale.rescaled.names[, "Type"] == units[t])
            }
            if (ts.col == TRUE & units[t] != "User") {
                rect(tscale[, "Start"], unit.depths[t], tscale[,
                  "End"], unit.depths[t + 1], col = rgb(tscale[,
                  "Col_R"], tscale[, "Col_G"], tscale[, "Col_B"],
                  maxColorValue = 255))
            }
            else rect(tscale[, "Start"], unit.depths[t], tscale[,
                "End"], unit.depths[t + 1], col = "white")
            text(tscale.names[, "Midpoint"], (unit.depths[t] +
                unit.depths[t + 1])/2, tscale.names[, "Name"],
                cex = cex.ts * tscale.data[match(units[t], rownames(tscale.data)),
                  "size"], srt = tscale.data[match(units[t],
                  rownames(tscale.data)), "srt"])
        }
        par(fig = c(0, 1, ts.width, 1), new = T)
        par(mar = c(0, 3, 2, 2))
        plot.phylo(tree, plot = FALSE, no.margin = T, x.lim = x.lim,
            cex = cex.tip, ...)
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        if (!missing(boxes) && boxes != "no" && boxes != "n") {
            if (boxes == "User") {
                tscale <- root.age - user.scale[, c("Start",
                  "End")]
            }
            else {
                tscale <- subset(timescale.rescaled, timescale.rescaled[,
                  "Type"] == boxes)
            }
            rect(tscale[, "Start"], par()$usr[3], tscale[, "End"],
                par()$usr[4], col = c("grey90", "white"), border = NA)
        }
        par(fig = c(0, 1, ts.width, 1), new = T)
        par(mar = c(0, 3, 2, 2))
        plot.phylo(tree, label.offset = offset, edge.width = width,
            no.margin = T, x.lim = x.lim, cex = cex.tip, ...)
        if (ranges == TRUE) {
            par(lend = 1)
            segments(taxon.ranges[, "FAD"], lastPP$yy[c(1:length(tree$tip.label))],
                taxon.ranges[, "LAD"], lastPP$yy[c(1:length(tree$tip.label))],
                col = "black", lwd = width * 2)
        }
        if (units[1] == "User") {
            segments(min(root.age - user.scale[, "Start"]), par()$usr[3],
                max(root.age - user.scale[, "End"]), par()$usr[3])
        }
        else {
            segments(min(timescale.rescaled[timescale.rescaled[,
                "Type"] == units[1], "Start"]), par()$usr[3],
                max(timescale.rescaled[timescale.rescaled[, "Type"] ==
                  units[1], "End"]), par()$usr[3])
        }
    }
}

# Plot and save phylogeny with geological timescale (to use with Keynote)
pdf(
  "figures/my_geoscalePhylo_plot.pdf",
  width = unit(3.5, units = "pt"), height = unit(7, units = "pt")
)
my_geoscalePhylo(
  ladderize(Schoeneae_tree3@phylo, right = FALSE), root.edge = TRUE,
  units = c("Period", "Epoch"),
  boxes = "User", user.scale = data.frame(
    Start = c(80, 70, 60, 50, 40, 30, 20, 10),
    End   = c(70, 60, 50, 40, 30, 20, 10,  0),
    Name  = rep(" ", times = 8)
  ),
  x.lim   = c(80, 0),
  ts.col  = FALSE,
  cex.age = 0.75,
  cex.ts  = 0.75,
  quat.rm = TRUE
)
dev.off()

# Adjust height-related node-data by the tree-height
# to get HPDs to plot correctly (not backwards)
tree_height <- max(nodeHeights(Schoeneae_tree@phylo))
Schoeneae_tree@data <- Schoeneae_tree@data %>%
  mutate(
    height          = tree_height - height,
    height_median   = tree_height - height_median,
    height_0.95_HPD = purrr::map(height_0.95_HPD,
                        ~purrr::map_dbl(.,
                          ~tree_height - .
                        )
                      )
  )

subtribes <- subtribes %>%
  mutate(taxa = str_remove(paste(genus, species), " NA"))

# Identify nodes for clades to highlight ---------------------------------------

# .... Other Schoeneae subtribes -----------------------------------------------

subtribe_names <- unique(na.exclude(subtribes$subtribe))
subtribe_names <- subtribe_names[subtribe_names != "Schoeninae"]
clades_to_label <- purrr::map(subtribe_names, ~find_subtribe(
  Schoeneae_tree@phylo,
  subtribe_name = .x,
  subtribes_df = subtribes,
  additional_taxa_to_exclude = ifelse(.x == "Tricostulariinae",
    "Anthelepis paludosa",
    NA
  )
))
names(clades_to_label) <- subtribe_names

# .... Ingroup clades ----------------------------------------------------------

Schoenus_node <-
  find_node(Schoeneae_tree@phylo, "Schoenus")
Clade_A_node <-
  find_node(Schoeneae_tree@phylo, "Schoenus insolitus", "Schoenus sculptus")
Clade_B_node <-
  find_node(Schoeneae_tree@phylo, "Schoenus falcatus",  "Schoenus australis")
Cape_clade_node <-
  find_node(Schoeneae_tree@phylo, "Schoenus dregeanus", "Schoenus australis")

# Plot -------------------------------------------------------------------------

# .... X-axis scaling things ---------------------------------------------------

my_labels <- c(70, 60, 50, 40, 30, 20, 10, 0)
label_positions <- tree_height - my_labels

# .... Make data for grey and white blocks for timescale-background of tree ----

my_panel_grid <- Schoeneae_tree@phylo %>%
  get_tips_in_ape_plot_order() %>%
  map_dfr(~ tibble(
    x       = label_positions - 5,
    species = .x,
    alpha   = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
  )) %>%
  mutate(species = species %>%
    factor(levels = get_tips_in_ape_plot_order(Schoeneae_tree@phylo)) %>%
    as.numeric()
  )

# Plot grid alone and save for use with FigTree/Keynote
my_panel_grid_plot <- ggplot(my_panel_grid) +
  aes(x, species, alpha = alpha) +
  geom_tile(
    data = my_panel_grid,
    fill = "black"
  ) +
  scale_alpha_manual(values = c(0, 0.1), guide = FALSE) +
  theme_void()
ggsave(
  "figures/my_panel_grid_plot.pdf",
  my_panel_grid_plot,
  width = 10, height = 15
)

# .... Main plot ---------------------------------------------------------------

Cyperaceae_tree_plot <-
  ggtree(Schoeneae_tree, ladderize = TRUE) +
  geom_rootedge(rootedge = 5) +
  geom_tile(
    data = my_panel_grid,
    aes(x, species, alpha = alpha),
    fill = "black"
  ) +
  scale_alpha_manual(values = c(0, 0.1), guide = FALSE) +
  geom_tiplab(
    aes(label = paste0('italic(\"', label, '\")')),
    parse = TRUE,
    size = 2.5,
    offset = 2
  ) +
  geom_range(
    "height_0.95_HPD", center = "height_median",
    colour = "darkblue", alpha = 0.5,
    size = 1.5
  ) +
  theme_tree2() +
  scale_x_continuous(name = "Ma",
    limits   = c(-5, 105),  # very wide to make space for annotations (below)
    breaks   = label_positions,
    labels   = my_labels
  ) +
  # Remove empty space above, below tree
  scale_y_continuous(
    limits = c(0, Ntip(Schoeneae_tree@phylo) + 1),
    expand = c(0, 0)
  ) +
  # Remove extra line at left of time axes
  coord_capped_cart(bottom = "right") +
  # Move time axes' titles to the left
  theme(axis.title.x = element_text(hjust = 0.35))

# .... Annotations -------------------------------------------------------------

clade_label_offset <- 25
clade_bar_extension <- 0.2

# Label other Schoeneae subtribes
for (subtribe in subtribe_names) {
  Cyperaceae_tree_plot <- Cyperaceae_tree_plot +
    geom_cladelabel(clades_to_label[[subtribe]],
      label =
        if (subtribe == "Oreobolus") {
          paste0('italic("Oreobolus")~clade')
        } else {
          subtribe
        },
      parse  = TRUE,
      offset = clade_label_offset,
      extend = clade_bar_extension
    )
}

# Label ingroup clades
Cyperaceae_tree_plot <- Cyperaceae_tree_plot +
  geom_cladelabel(Schoenus_node, paste0('italic("Schoenus")'), parse  = TRUE,
    offset = clade_label_offset,
    extend = clade_bar_extension
  ) +
  stat_hilight(
    node = Clade_A_node,
    colour = "darkblue", xmin = 20,
    fill = NA, alpha = NA
  ) +
  stat_hilight(
    node = Clade_B_node,
    colour = "darkblue", xmin = 20,
    fill = NA, alpha = NA
  ) +
  stat_hilight(
    node = Cape_clade_node,
    colour = "darkblue",
    fill = NA, alpha = NA
  ) +
  annotate(
    label = "Clade A", geom = "text",
    x = 25, y = Ntip(Schoeneae_tree@phylo) - 2
  ) +
  annotate(
    label = "Clade B", geom = "text",
    x = 25, y = (0.6 * Ntip(Schoeneae_tree@phylo)) - 2
  ) +
  annotate(
    label = "Cape clade", geom = "text",
    x = 42, y = (0.6 * Ntip(Schoeneae_tree@phylo)) - 2
  )

# Save plot --------------------------------------------------------------------

ggsave(
  "figures/Cyperaceae_tree_plot.pdf",
  Cyperaceae_tree_plot,
  width = 10, height = 15
)

ggsave(
  "figures/Cyperaceae_tree_plot.png",
  Cyperaceae_tree_plot,
  width = 10, height = 15, dpi = 300
)
