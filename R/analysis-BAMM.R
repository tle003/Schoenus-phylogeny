# BAMM analysis of Schoenus
# Tammy Elliott & Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(hisse)
library(phytools)
library(BAMMtools)
library(coda)

# ... --------------------------------------------------------------------------

# Read newick tree file, first drop excess tips
#Cyperaceae_tree<-read.tree("Cyperaceae-max-clade-cred-jul29.newick")
#Only keep Schoeneae species for analysis
#Schoenus_tree_prelim<-drop.tip(Cyperaceae_tree, Cyperaceae_tree$tip.label[!grepl("Schoenus", Cyperaceae_tree$tip.label)])
#Schoenus_tree<-drop.tip(Schoenus_tree_prelim, c("XyroSchoeneae_hornei", "GymnoSchoeneae_sphaerocephalus"))
#str(Schoenus_tree)

#phy<-force.ultrametric(Schoenus_tree, method=c("nnls","extend"))
#is.ultrametric(phy)

#Write ultrametric tree
#write.tree(phy, file = "Schoenus_tree_ultrametric.tre", append = FALSE, digits = 10, tree.names = FALSE)
#phy<-read.tree("Schoenus_tree_ultrametric.tre")

# Import data ------------------------------------------------------------------

tree <- read.tree("data/phylogenies/for-BAMM-analysis/Schoenus_tree_ultrametric.tre")

# Make sure tree is acceptable for BAMM's assumptions --------------------------

# Check to make sure that tree makes basic checks
is.ultrametric(tree, option = 2)  # Use option 2 for compatibility reasons
is.binary.tree(tree)

# Count negative branches
sum(tree$edge.length < 0)
# None!

# Count zero length branches
sum(tree$edge.length == 0)
# None!




#Set priors; run this to see how to adapt block
setBAMMpriors(tree)

#Generate control file
#the block will have to be altered for the final runs
#I am using a sampling fraction of .7 since we have sampled about 70% of Schoenus species
#Assumes species sampled at random
generateControlFile("divcontrol.txt",
  type = "diversification",
  params = list(
    treefile                     = "Schoenus_tree_ultrametric.tre",
    runInfoFilename              = "run_info.txt",
    runMCMC                      = 1,
    numberOfGenerations          = 3500000,
    overwrite                    = "1",
    lambdaInitPrior              = "2.65133706894448",
    lambdaShiftPrior             = "0.0216229837749123",
    muInitPrior                  = "2.65133706894448",
    expectedNumberOfShifts       = "1",
    useGlobalSamplingProbability = 1,
    globalSamplingFraction       = 0.70,
    mcmcWriteFreq                = 10000,
    eventDataWriteFreq           = 10000,
    printFreq                    = 1000,
    acceptanceResetFreq          = 10000
  )
)

# Run BAMM with terminal
system("C:\\Users\\USER\\Downloads\\bamm-2.5.0-Windows\\bamm-2.5.0-Windows\\bamm-2.5.0-Windows\\bamm.exe -c C:\\Users\\USER\\Desktop\\Schoenus-phylogeny\\divcontrol.txt --seed 1234")

# Read BAMMdata
edata <- getEventData(tree, eventdata = "event_data.txt", burnin = 0.1)

# Access MCMC convergence
mcmcout <- read.csv("mcmc_out.txt")
plot(mcmcout$logLik ~ mcmcout$generation)

burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

#Access effective sample sizes; this should be ESS of >200
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

#Also can access BAMM convergence by doing different runs and analyzing 
#branch-specific marginal rate shift probabilities (marginalShiftProbsTree)


#Number of rate shifts
post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)

#Plot barplot
barplot(post_probs, ylab="Shift posterior distribution", xlab="Number of rate shifts", ylim=c(0,1))



#Compute the posterior odds ration for two models
post_probs['0'] / post_probs['1']

#Alternative method: summarize the posterior distribution of the number of shifts using the summary method
# Gives posteior probabilities of eac rate shift count
(shift_probs <- summary(edata))



#Prior distribution; compute Bayes factors
#postfile <- "post_mcmc_out.txt"
postfile <- "mcmc_out.txt"
(bfmat <- computeBayesFactors(postfile, expectedNumberOfShifts=1, burnin=0.1))
#BF >12 are considered by some to have at least some effect on 'significance'

#Visualize the prior and posterior simultaneously
plotPrior(postfile, expectedNumberofShifts=1)

#Mean phylorate plot
summary(edata)
plot.bammdata(edata,lwd=2)

#Add interactive legend
plot.bammdata(edata, lwd=2, legend=T)


#Plot phylograte plot
#This shows example of different options for mapping colors to rates in phylorate plots
plot.bammdata(edata, lwd=3, method="polar", pal="temperature")

#Plot x samplle from the posterior 
index <- 25
e2 <- subsetEventData(edata, index = index)
plot.bammdata(e2, lwd=2)
addBAMMshifts(e2, cex=2)

#Plot credible set of macroeveolutionary rate configurations using the marginal odds ratio
#95% credible set - account for 95% of the probability of the data
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=3, set.limit = 0.95)
plot(css)

#Number of distinct shift configurations in the data
css$number.distinct

summary(css)

#Generate phylorate plots for each of the shift configurations with the highest posterior probabilities
plot.credibleshiftset(css)

#Maximum a posteriori probability (MAP) shift configuration; this is the shift configuation with the highest posterior probability (the one that was sampled most)
#This represents the best shift configuation
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1, threshold=5)
plot.bammdata(best, lwd=2)
addBAMMshifts(best, cex=2.5)

#Access list of vector of the most-probable shift configurations
#This gives samples assigned to the most probable shift configuation
css$indices[[1]]

#Plot any specific sample from any shift configuation
#The code below alogs us to sample the 5th sample assigned to the most-probable shift configuration
index <- css$indices[[1]][5]
rsample <- subsetEventData(edata, index=index)
plot.bammdata(rsample)
addBAMMshifts(rsample, cex=2)


#Random shift configurations
dsc <- distinctShiftConfigurations(edata, expectedNumberOfShifts=1, threshold=5)
# Here is one random sample with the BEST shift configuration
plot.bammshifts(dsc, edata, rank=1, legend=F)
# Here is another (read the label text):
plot.bammshifts(dsc, edata, rank=1, legend=F)


#Plotting rate shifts using plot.phylo
mysample <-55

#Get total number of rate regimes on a tree of this sample
#1 = no rate shift
nrow(edata$eventData[[ mysample ]])

#Get node numbers
shiftnodes <- getShiftNodesFromIndex(edata, index = mysample)

#Plot these node on tree
plot.phylo(phy, cex=0.3)
nodelabels(node = shiftnodes, pch=21, col="red", cex=1.5)


#Calculate marginal shift probabilities
#Marginal probability that each branch contains one or more shift events
marg_probs <- marginalShiftProbsTree(edata)
plot.phylo(marg_probs)


#Distinct shift configurations
#distinct shift configurations” within a given dataset as well as the posterior probability of each configuration.
#Each distinct shift configuration may have been sampled multiple times during simulation of the posterior.

cset <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=3)
plot.credibleshiftset(cset, lwd=2.5)


#Prior probabilities of rate shifts on branches; branch shift priors are now equal to branch lengths in the phy object
branch_priors <- getBranchShiftPriors(edata, expectedNumberOfShifts = 1)
plot(branch_priors)

#Marginal Odds Ratio; gives relatives odds that a shift occurred on a specific branch given a shift occurred at all
mo <- marginalOddsRatioBranches(edata, branch_priors)

#Rate variation through time:Color density and grey scale
plot.new()
st <- max(branching.times(phy))
plotRateThroughTime(edata, intervalCol="darkgreen", avgCol="darkgreen", start.time=st, ylim=c(0,1), cex.axis=2)
text(x=30, y= 0.8, label="Schoenus", font=4, cex=2.0, pos=4)

plot.new()
plotRateThroughTime(edata, avgCol="black", start.time=st, ylim=c(0,1), cex.axis=2, intervalCol='gray80', intervals=c(0.05, 0.95), opacity=1)
text(x=30, y= 0.8, label="Schoenus", font=4, cex=2.0, pos=4)








#three side-by-side; shows different colour ramp options
par(mfrow=c(1,3), mar=c(1, 0.5, 0.5, 0.5), xpd=TRUE)

q <- plot.bammdata(edata, tau=0.001, breaksmethod='linear', lwd=2)
addBAMMshifts(edata, par.reset=FALSE, cex=2)
title(sub='linear',cex.sub=2, line=-1)
addBAMMlegend(q, location=c(0, 1, 140, 220))

q <- plot.bammdata(edata, tau=0.001, breaksmethod='linear', color.interval=c(NA,0.12), lwd=2)
addBAMMshifts(edata, par.reset=FALSE, cex=2)
title(sub='linear - color.interval',cex.sub=2, line=-1)
addBAMMlegend(q, location=c(0, 1, 140, 220))

q <- plot.bammdata(edata, tau=0.001, breaksmethod='jenks', lwd=2)
addBAMMshifts(edata, par.reset=FALSE, cex=2)
title(sub='jenks',cex.sub=2, line=-1)
addBAMMlegend(q, location=c(0, 1, 140, 220))


#Phylorate plots for three different evolutionary rate configurations
#Columns show the 10th, 20th, and 30th samples from the posterior distribution
#Every sample from the posterior contains a potentially unique configuration of shift locations, 
        #which are sampled in proportion to their posterior probability. 

ixx <- rep(c(10, 30, 40), 2);
plot.new()
par(mfrow=c(2,3));
colschemes <- list();
colschemes[1:3] <- 'temperature'
colschemes[4:6] <- list(c('blue', 'gray', 'red'))

for (i in 1:length(ixx)) {
        par(mar=c(0,0,0,0))
        index <- ixx[i]
        eventsub <- subsetEventData(edata, index=index);
        plot.bammdata(eventsub, method='polar', pal= colschemes[[i]], par.reset=FALSE, lwd=3)
        addBAMMshifts(eventsub, method='polar', index=1, col='white', bg='black', cex=5, par.reset=FALSE)
}




#do over posterior distribution
# Read newick tree file, first drop excess tips
Cyperaceae_trees<-read.tree("Cyperaceae-max-clade-cred-jul29.newick")
#Only keep Schoeneae species for analysis
#Schoenus_tree_prelim<-drop.tip(Cyperaceae_tree, Cyperaceae_tree$tip.label[!grepl("Schoenus", Cyperaceae_tree$tip.label)])
#Schoenus_tree<-drop.tip(Schoenus_tree_prelim, c("XyroSchoeneae_hornei", "GymnoSchoeneae_sphaerocephalus"))
#str(Schoenus_tree)




#Over 100 trees

#Make sure tree is acceptable
trees <- read.tree("Schoenus.trees.100.trees")
#Check to make sure that tree makes basic checks
lapply(trees, function(x) is.ultrametric(x))
#make ultrametric
phys<-lapply(trees, function(x) force.ultrametric(x, method=c("nnls","extend")))
lapply(phys, function(x) is.ultrametric(x))

#Check to see if binary
lapply(phys, function(x) is.binary.tree(x))
# count negative branches
lapply(phys, function(x) sum(x$edge.length < 0))

# count zero length branches
lapply(phys, function(x) sum(x$edge.length == 0))

#Tree 36 has a branch length of 0

"multiPhylo"<-class(phys)

#Set priors; run this to see how to adapt block; this needs to be iterated over 100 trees and will make 100 different prior.files


setBAMMpriors(phy = phys)





#Now we need 100 different divcontrol files with different lambdaInitPriors, lambdaShiftPriors, mininitPriors
#Generate control file
#the block will have to be altered for the final runs
#I am using a sampling fraction of .7 since we have sampled about 70% of Schoenus species
#Assumes species sampled at random
generateControlFile('divcontrol.txt', type = 'diversification', params = list(
        treefile = "Schoenus_tree_ultrametric.tre",
        runInfoFilename = "run_info.txt",
        runMCMC = 1,
        numberOfGenerations = 3500000,
        overwrite = '1',
        lambdaInitPrior = '2.65133706894448',
        lambdaShiftPrior = '0.0216229837749123',
        muInitPrior = '2.65133706894448',
        expectedNumberOfShifts = '1',
        useGlobalSamplingProbability = 1,
        globalSamplingFraction = 0.70,
        mcmcWriteFreq = 10000,
                eventDataWriteFreq = 10000,
                printFreq = 1000,
                acceptanceResetFreq = 10000))


#Run BAMM with terminal; run for all 100 divcontrol.txt files
#bamm -c divcontrol.txt --seed 1234


#Read Bammdata: 100 times
edata <- getEventData(tree, eventdata = "event_data.txt", burnin=0.1)

#Access MCMC convergence:100 times
mcmcout <- read.csv("mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

#Access effective sample sizes; this should be ESS of >200: 100 times
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

#Also can access BAMM convergence by doing different runs and analyzing 
#branch-specific marginal rate shift probabilities (marginalShiftProbsTree)


#Number of rate shifts; this is for 1 run; is it possible to do for 100 runs
post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)

#Plot barplot: i fun; can you do for 100 runs
barplot(post_probs, ylab="Shift posterior distribution", xlab="Number of rate shifts", ylim=c(0,1))





#Compute the posterior odds ration for two models;can you do for 100 runs
post_probs['0'] / post_probs['1']

#Alternative method: summarize the posterior distribution of the number of shifts using the summary method
# Gives posteior probabilities of eac rate shift count;can you do for 100 runs
(shift_probs <- summary(edata))



#Prior distribution; compute Bayes factors;can you do for 100 runs
#postfile <- "post_mcmc_out.txt"
postfile <- "mcmc_out.txt"
(bfmat <- computeBayesFactors(postfile, expectedNumberOfShifts=1, burnin=0.1))
#BF >12 are considered by some to have at least some effect on 'significance'

#Visualize the prior and posterior simultaneously;can you do for 100 runs
plotPrior(postfile, expectedNumberofShifts=1)



#Plot credible set of macroeveolutionary rate configurations using the marginal odds ratio;can you do for 100 runs
#95% credible set - account for 95% of the probability of the data
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=3, set.limit = 0.95)
plot(css)
