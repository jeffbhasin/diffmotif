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
library(Matching)
library(rms)
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

drawBackgroundSetPropensity <- function(target.meta, pool.meta)
{
	# TODO - draw a background set using propensity score matching
	# TODO - make function take a variable of which model/covars to use

	# estimate propensity model first

	# instead of "treatment" the variable is 1 if in target distribution and 0 if in background pool distribution
	# the linear regression is then of this categorical assignment versus a model of the other covariates we want

	# create "treatment" variable

	# run logistic regression of treatment ~ covar model

	# match based on propensity score

	# test quality of matching

}
# -----------------------------------------------------------------------------

# =============================================================================

# =============================================================================
# Analysis Code

# -----------------------------------------------------------------------------
# Load seqs and run via function

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
# Testing of Matching package options and linear models

target.meta <- seq1.meta
pool.meta <- seq2.meta

target.meta$treat <- 0
pool.meta$treat <- 1
all.meta <- rbind(target.meta, pool.meta)

# GC
lrm.out.gc <- lrm(treat ~ gc, data=all.meta)
#png(filename="output/propscore.gc.lrm.png",width=1000,height=800,res=120)
lrm.out.gc.fitted <- predict.lrm(lrm.out.gc,type="fitted")
#qplot(x=gc,y=treat,data=all.meta) + geom_line(aes(x, y), data.frame(x=all.meta$gc,y=lrm.out.gc.fitted)) + ggplot.clean()
#dev.off()



# Size
#glm.out <- glm(treat ~ size, family=binomial(link=logit), data=all.meta)
lrm.out.size <- lrm(treat ~ size, data=all.meta)
#png(filename="output/propscore.size.lrm.png",width=1000,height=800,res=120)
lrm.out.size.fitted <- predict.lrm(lrm.out.size,type="fitted")
#qplot(x=size,y=treat,data=all.meta) + geom_line(aes(x, y), data.frame(x=all.meta$size,y=lrm.out.size.fitted)) + ggplot.clean()
#dev.off()


# Size + GC
lrm.out.sizegc <- lrm(treat ~ size + gc, data=all.meta)
lrm.out.sizegc.fitted <- predict.lrm(lrm.out.sizegc,type="fitted")
# can't plot without separating dimmensions
#png(filename="output/propscore.size-gc.lrm.png",width=1000,height=800,res=120)
#qplot(x=size+gc,y=treat,data=all.meta) + geom_point(aes(x, y), data.frame(x=all.meta$size+all.meta$gc,y=1-glm.out$fitted))
#dev.off()


#Do Matching
#rr.gc  <- Match(Y=NULL, Tr=all.meta$treat, X=lrm.out.gc.fitted, M=1)
#rr.size  <- Match(Y=NULL, Tr=all.meta$treat, X=lrm.out.size.fitted, M=1, version="fast")
#rr.sizegc  <- Match(Y=NULL, Tr=all.meta$treat, X=lrm.out.size.fitted, M=1, version="fast")

#rr.test <- Match(Y=NULL, Tr=all.meta$treat), X=head(lrm.out.gc.fitted), M=1)

#summary(rr.size)

#mb  <- MatchBalance(treat ~ size + gc, data=all.meta, match.out=rr, nboots=10)


# -----------------------------------------------------------------------------

# =============================================================================
