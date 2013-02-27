# #############################################################################
# Test code for implementing propensity score matching
# Author: Jeffrey Bhasin <jmb85@case.edu>
# Created: 2013-02-27
# #############################################################################

# =============================================================================
# Packages and Globals

nCores <- 15

library(ggplot2)
library(gridExtra)
library(doMC)
registerDoMC(nCores)

source("Rmotif.R")

# =============================================================================

# =============================================================================
# Local Functions

# -----------------------------------------------------------------------------
# Calculate GC content of a DNAString
# Input: DNAStringSet object
# Output: Vector of GC contents
getGC <- function(seq)
{
	g <- alphabetFrequency(seq)[,3]
	c <- alphabetFrequency(seq)[,2]
	a <- alphabetFrequency(seq)[,1]
	t <- alphabetFrequency(seq)[,4]

	gc <- (g+c) / (a+t+g+c)
	gc
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Draws background from a source using propensity score matching
# Input:
# Output: 
library(Matching)
drawBackgroundSetPropensity <- function()
{
	# TODO - draw a background set using propensity score matching

	# estimate propensity model first
	#
}
# -----------------------------------------------------------------------------

# =============================================================================

# =============================================================================
# Analysis Code

# -----------------------------------------------------------------------------
# Section 1

# load source set of sequences
seq1.path <- "../RunMEME/ColonSeqGR2012/CIMPHyperMe.fasta"
seq1 <- readDNAStringSet(seq1.path)
seq1.nSeqs <- length(seq1)

# load background "pool" we want to draw from and make match the distribution of variable from the target
seq2.path <- "ignore/ns.111.fasta"
seq2 <- readDNAStringSet(seq2.path)
seq2.nSeqs <- length(seq2)

# make dataframe of variables for each sequence - GC and size for now, could add distance from TSS, etc in the future
seq1.meta <- data.frame(name=names(seq1),size=width(seq1),gc=getGC(seq1))
seq2.meta <- data.frame(name=names(seq2),size=width(seq2),gc=getGC(seq2))

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 2

# -----------------------------------------------------------------------------

# =============================================================================
