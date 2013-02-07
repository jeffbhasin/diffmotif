###############################################################################
## Usage of Rmotif.R functions for the differentially methylated sites
##
###############################################################################

###############################################
## Load Rmotif functions
###############################################
source("Rmotif.R")

###############################################
## Settings
###############################################
nCores <- 15

###############################################
## Packages
###############################################
library(doMC)
registerDoMC(nCores)

###############################################
## Main Body
###############################################
#Sequences in which to find motifs (FASTA)
seq.path <- "../RunMEME/ColonSeqGR2012/CIMPHyperMe.fasta"
genome.path <- "../RunMEME/ColonSeqGR2012/hg18.fa"

#filter of q values from FIMO output
q.cutoff <- 0.05

