## Title:  Schoeneae ancestoral areas
### Name: Tammy Elliott
### Date: June 15, 2020


# Load the required libraries.
#install.packages("rexpokit")
#install.packages("cladoRcpp")
#install.packages("devtools")
#library(devtools)
#devtools::install_github(repo="nmatzke/BioGeoBEARS")

#set working directory
setwd("/Users/tammy/PostDoc/Smuts/Schoeneae/Phylogeny/manuscript/follow-up-analyses")
getwd()

# import your species points
# this creates an object called georef
Cyp.cont<-read.csv("TDWG-continents.csv",header=T)
head(Cyp.cont)
str(Cyp.cont)

#Drop non-Schoeneae taxa
Schoeneae.cont.prelim<-Cyp.cont[grepl("Schoeneae", Cyp.cont$Species),]
#Remove Xyroschoenues and GymnoSchoeneae
remove.species<- c("GymnoSchoeneae_sphaerocephalus", "XyroSchoeneae_hornei")
Schoeneae.cont<-Schoeneae.cont.prelim[!(Schoeneae.cont.prelim$Species %in% remove.species), ]
head(Schoeneae.cont)
str(Schoeneae.cont)

#Extract Schoeneae taxa
Schoeneae.cont.outgroup<-Cyp.cont[grep("Capeobolus brevicaulis|Cyathochaeta avenacea|Cyathocoma hexandra|Evandra aristata|Gahnia tristis|Lepidosperma longitudinale|Machaerina rubiginosa|Neesenbeckia punctoria|Ptilothrix deusta|Chaetospora curvifolia|Ammothryon grandiflorum|Tetraria capillaris|Trianoptiles capensis|Tricostularia pauciflora|Caustis blakei|Oreobolus pectinatus|Tetraria borneensis|XyroSchoeneae_hornei|Chamaedendron fragilis|Costularia leucocarpa|Tetrariopsis octandra|Reedia spathacea|GymnoSchoeneae_sphaerocephalus|Anthelepis paludosa|Morelotia gahniiformis|Mesomelaena tetragona|Tetraria fasciata", Cyp.cont$Species),]
str(Schoeneae.cont.outgroup)

#Create a georeference base file for Schoneaneae
Schoeneae.cont<-rbind(Schoeneae.cont, Schoeneae.cont.outgroup)

#Format to get .txt file that is compatable with text file format required by BioGeoBEARS; This is for Scheneae
Schoeneae.cont.format<-cbind(Schoeneae.cont$Species, Schoeneae.cont$Europe, Schoeneae.cont$Africa, Schoeneae.cont$Asia.temperate,
	Schoeneae.cont$Asia.tropical, Schoeneae.cont$Australasia, Schoeneae.cont$Pacific, Schoeneae.cont$Northern.America, Schoeneae.cont$Southern.America)
head(Schoeneae.cont.format)
colnames(Schoeneae.cont.format)<-c("Species","Europe", "Africa", "AsTemp", "AsTrop", "Aust", "Pac", "NAm", "SAm")
head(Schoeneae.cont.format)
Schoeneae_cont_format<-Schoeneae.cont.format

#Write .txt file 
#write.table(Schoeneae_cont_format, file = "Schoeneae_cont_format.txt", quote=F)

#Format to get .txt file that is compatable with text file format required by BioGeoBEARS; This is for Schoneae
Schoeneae.cont.format<-cbind(Schoeneae.cont$Species, Schoeneae.cont$Europe, Schoeneae.cont$Africa, Schoeneae.cont$Asia.temperate,
	Schoeneae.cont$Asia.tropical, Schoeneae.cont$Australasia, Schoeneae.cont$Pacific, Schoeneae.cont$Northern.America, Schoeneae.cont$Southern.America)
head(Schoeneae.cont.format)
colnames(Schoeneae.cont.format)<-c("Species","Europe", "Africa", "AsTemp", "AsTrop", "Aust", "Pac", "NAm", "SAm")
head(Schoeneae.cont.format)
Schoeneae_cont_format<-Schoeneae.cont.format

#Write .txt file 
#write.table(Schoeneae_cont_format, file = "Schoeneae_cont_format.txt", quote=F)

#The first analysis is using Schoeneae only

#extdata loaad
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir

geo_file_name <- np("Schoeneae_cont_binary.txt")
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=geo_file_name)
tipranges
# Read newick tree file, first drop excess tips
Cyperaceae_tree<-read.tree("Cyperaceae-all-taxa-6calib-st2.newick")
#Only keep Schoeneae species for analysis
Schoeneae_tree_prelim<-drop.tip(Cyperaceae_tree, Cyperaceae_tree$tip.label[!grepl("Schoeneae", Cyperaceae_tree$tip.label)])
Schoeneae_tree<-drop.tip(Schoeneae_tree_prelim, c("XyroSchoeneae_hornei", "GymnoSchoeneae_sphaerocephalus"))
str(Schoeneae_tree)

#See difference in taxa between geo_file name and tree
# extra species will have to be droped from geo_file
(diff_1<-setdiff(Schoeneae.cont$Species, Schoeneae_tree$tip.label))

#Drop data from .text file by hand

#write.tree(Schoeneae_tree, file = "Schoeneae_tree.tre")
#see difference between tip labels and rows in geo_file_name





#Read tree
tree_file_name <- np("Schoeneae_tree.tre")
tree_file_name

tr <- read.tree(tree_file_name)
plot(tr, cex=0.3)
tree<-tr

#Set-up BioGeoBEARS parametres
BioGeoBEARS_run_object <- define_BioGeoBEARS_run()

# Set tree and data set objects
BioGeoBEARS_run_object$trfn <- tree_file_name
BioGeoBEARS_run_object$geogfn <- geo_file_name

# Set parametres
BioGeoBEARS_run_object$max_range_size = 8
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

#Run analysis
results_DEC <- bears_optim_run(BioGeoBEARS_run_object)

#DEC+J Model
#Set up
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

#Location of input files
BioGeoBEARS_run_object$trfn = tree_file_name
BioGeoBEARS_run_object$geogfn = geo_file_name

#Set parametres
BioGeoBEARS_run_object$max_range_size = 8
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

#Set-up model

dstart = results_DEC$outputs@params_table["d","est"]
estart = results_DEC$outputs@params_table["e","es"]

#Set starting values
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

#Add jstart as new free parametre
jstart = 0.0001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

#Run analysis
results_DECJ = bears_optim_run(BioGeoBEARS_run_object)

#Plot ancestral states
#Plot the DEC and DEC+J models

pdffn = "Schoeneae-prac_DEC_vs_DEC+J.pdf"
pdf(pdffn, width=6, height=6)
analysis_titletxt = "DEC on Schoeneae"

scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
res2 = plot_BioGeoBEARS_results(results_DEC, analysis_titletxt, addl_params=list("j"), plotwhat="text",
label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

plot_BioGeoBEARS_results(results_DEC, analysis_titletxt, addl_params=list("j"), plotwhat="pie",
label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

analysis_titletxt ="DEC+J on Schoeneae"
results_object = results_DECJ
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

res1 = plot_BioGeoBEARS_results(results_DECJ, analysis_titletxt, addl_params=list("j"), plotwhat="text",
label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

plot_BioGeoBEARS_results(results_DECJ, analysis_titletxt, addl_params=list("j"), plotwhat="pie",
label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)


#Analysis for Schoeneae
#extdata loaad
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir

geo_file_name <- np("Schoeneae_cont_binary.txt")
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=geo_file_name)
tipranges
# Read newick tree file, first drop excess tips
Cyperaceae_tree<-read.tree("Cyperaceae-all-taxa-6calib-st2.newick")
#Only keep Schoeneae species for analysis
Schoeneae_tree<-drop.tip(Cyperaceae_tree, c("Trianoptiles_capensis", "Carpha_glomerata","Fimbristylis_complanata", "Fimbristylis_ovata",
"Bulbostylis_barbata", "Eleocharis_acicularis", "Eleocharis_palustris", "Fuirena_hirsuta", "Ficinia_indica", "Cyperus_papyrus", "Schoenoplectus_tabernaemontani",
"Dulichium_arundinaceum","Blysmus_compressus","Calliscirpus_criniger","Calliscirpus_brachythrix","Scirpus_pendulus",
"Scirpus_hattorianus","Eriophorum_vaginatum", "Carex_rupestris","Carex_siderosticta","Sumatroscirpus_rupestris","Khaosokia_caricoides","Rhynchospora_rugosa" ,"Rhynchospora_capitellata",      
"Exochogyne_amazonica","Lagenocarpus_alboniger","Cladium_mariscus","Cladium_chinense","Mapania_cuspidata","Hypolytrum_nemorum","Hypolytrum_africanum","Calyptrocarya_glomerulata","Diplacrum_caricinum",
"Scleria_rugosa","Scleria_brownii"))

str(Schoeneae_tree)


#Drop data from .text file by hand

#write.tree(Schoeneae_tree, file = "Schoeneae_tree.tre")
#see difference between tip labels and rows in geo_file_name

#write.tree(Schoeneae_tree, file = "Schoeneae_tree.tre")
#see difference between tip labels and rows in geo_file_name

#Read tree
tree_file_name <- np("Schoeneae_tree.tre")
tree_file_name

tr <- read.tree(tree_file_name)
plot(tr, cex=0.3)
tree<-tr

#Set-up BioGeoBEARS parametres
BioGeoBEARS_run_object <- define_BioGeoBEARS_run()

# Set tree and data set objects
BioGeoBEARS_run_object$trfn <- tree_file_name
BioGeoBEARS_run_object$geogfn <- geo_file_name

# Set parametres
BioGeoBEARS_run_object$max_range_size = 8
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

#Run analysis
results_DEC <- bears_optim_run(BioGeoBEARS_run_object)



#DEC+J Model
#Set up
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

#Location of input files
BioGeoBEARS_run_object$trfn = tree_file_name
BioGeoBEARS_run_object$geogfn = geo_file_name

#Set parametres
BioGeoBEARS_run_object$max_range_size = 8
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

#Set-up model

dstart = results_DEC$outputs@params_table["d","est"]
estart = results_DEC$outputs@params_table["e","es"]

#Set starting values
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

#Add jstart as new free parametre
jstart = 0.0001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

#Run analysis
results_DECJ = bears_optim_run(BioGeoBEARS_run_object)

#Plot ancestral states
#Plot the DEC and DEC+J models

pdffn = "Schoeneae-prac_DEC_vs_DEC+J.pdf"
pdf(pdffn, width=6, height=6)
analysis_titletxt = "DEC on Schoeneae"

scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
res2 = plot_BioGeoBEARS_results(results_DEC, analysis_titletxt, addl_params=list("j"), plotwhat="text",
label.offset=0.45, tipcex=0.4, statecex=0.4, splitcex=0.35, titlecex=0.8, plotsplits=TRUE,
cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

plot_BioGeoBEARS_results(results_DEC, analysis_titletxt, addl_params=list("j"), plotwhat="pie",
label.offset=0.45, tipcex=0.4, statecex=0.4, splitcex=0.35, titlecex=0.8, plotsplits=TRUE,
cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

analysis_titletxt ="DEC+J on Schoeneae"
results_object = results_DECJ
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

res1 = plot_BioGeoBEARS_results(results_DECJ, analysis_titletxt, addl_params=list("j"), plotwhat="text",
label.offset=0.45, tipcex=0.4, statecex=0.4, splitcex=0.35, titlecex=0.8, plotsplits=TRUE,
cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

plot_BioGeoBEARS_results(results_DECJ, analysis_titletxt, addl_params=list("j"), plotwhat="pie",
label.offset=0.45, tipcex=0.4, statecex=0.4, splitcex=0.35, titlecex=0.8, plotsplits=TRUE,
cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)


#Based on the above, use only Schoenus data for model comparison



#Read tree
tree_file_name <- np("Schoenus_tree.tre")
tree_file_name
#Look into the raw Newick file
moref(tree_file_name)


# Look at your phylogeny (plots to a PDF, which avoids issues with multiple graphics in same window):
pdffn = "Schoenus_tree.pdf"
pdf(file=pdffn, width=9, height=12)

tr = read.tree(tree_file_name)
tr
plot(tr)
title("Preliminary Schoenus phylogeny")
axisPhylo() # plots timescale

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)

#Set up geography file
#extdata loaad
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir

geo_file_name <- np("Schoenus_cont_binary.txt")
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=geo_file_name)
tipranges

# Maximum range size observed:
max(rowSums(dfnums_to_numeric(tipranges@df)))

# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
# Schoenus nigricans has 6
max_range_size = 6

#######################################################
#######################################################
# DEC AND DEC+J ANALYSIS
#######################################################
#######################################################

# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = tree_file_name

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geog_file_name

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
# (use more cores to speed it up; this requires
# library(parallel) and/or library(snow). The package "parallel" 
# is now default on Macs in R 3.0+, but apparently still 
# has to be typed on some Windows machines. Note: apparently 
# parallel works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

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
check_BioGeoBEARS_run(BioGeoBEARS_run_object

#Run analysis
results_DEC <- bears_optim_run(BioGeoBEARS_run_object)




#DEC+J Model
#Set up
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

#Location of input files
BioGeoBEARS_run_object$trfn = tree_file_name
BioGeoBEARS_run_object$geogfn = geo_file_name

#Set parametres
BioGeoBEARS_run_object$max_range_size = 6
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE

BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE # Set ancestral states from optim run

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale


#Set-up model (DEC + J)
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = results_DEC$outputs@params_table["d","est"]
estart = results_DEC$outputs@params_table["e","es"]

#Set starting values for d,e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

#Add jstart as new free parametre
jstart = 0.0001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

#Run analysis
results_DECJ = bears_optim_run(BioGeoBEARS_run_object)


#######################################################
#######################################################
# DIVALIKE AND DIVALIKE+J ANALYSIS
#######################################################
#######################################################

# Run DIVALIKE
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = tree_file_name
BioGeoBEARS_run_object$geogfn = geo_file_name
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)


# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DIVALIKE model
# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5


check_BioGeoBEARS_run(BioGeoBEARS_run_object)

#Run analysis
results_DIVALIKE = bears_optim_run(BioGeoBEARS_run_object)



#######################################################
# Run DIVALIKE+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = tree_file_name
BioGeoBEARS_run_object$geogfn = geo_file_name
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    


# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DIVALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

#Run analysis
results_DIVALIKE_J = bears_optim_run(BioGeoBEARS_run_object)



#Plot the DIVA and DIVA+J models

pdffn = "Schoenus-prac_DIVA_vs_DIVA+J.pdf"
pdf(pdffn, width=6, height=6)
analysis_titletxt = "DIVA-like on Schoenus"

scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
res2 = plot_BioGeoBEARS_results(results_DIVALIKE, analysis_titletxt, addl_params=list("j"), plotwhat="text",
label.offset=0.45, tipcex=0.4, statecex=0.4, splitcex=0.35, titlecex=0.8, plotsplits=TRUE,
cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

plot_BioGeoBEARS_results(results_DIVALIKE, analysis_titletxt, addl_params=list("j"), plotwhat="pie",
label.offset=0.45, tipcex=0.4, statecex=0.4, splitcex=0.35, titlecex=0.8, plotsplits=TRUE,
cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

analysis_titletxt ="DIVA-like+J on Schoenus"
results_object = results_DECJ
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

res1 = plot_BioGeoBEARS_results(results_DIVALIKE_J, analysis_titletxt, addl_params=list("j"), plotwhat="text",
label.offset=0.45, tipcex=0.4, statecex=0.4, splitcex=0.35, titlecex=0.8, plotsplits=TRUE,
cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

plot_BioGeoBEARS_results(results_DIVALIKE_J, analysis_titletxt, addl_params=list("j"), plotwhat="pie",
label.offset=0.45, tipcex=0.4, statecex=0.4, splitcex=0.35, titlecex=0.8, plotsplits=TRUE,
cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)


