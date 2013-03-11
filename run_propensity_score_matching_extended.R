# #############################################################################
# Test code for implementing propensity score matching
# Author: Jeffrey Bhasin <jmb85@case.edu>
# Created: 2013-02-27
# #############################################################################

# This version adds additional models beyond just size and GC content

# =============================================================================
# Packages and Globals

source("Rmotif.R")

nCores <- 15
registerDoMC(nCores)

# =============================================================================

# =============================================================================
# Local Functions

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
# Load seqs

# load source set of sequences
seq1.path <- "../RunMEME/ColonSeqGR2012/CIMPHyperMe.fasta"
seq1 <- readDNAStringSet(seq1.path)
seq1.nSeqs <- length(seq1)

# load background "pool" we want to draw from and make match the distribution of variable from the target
seq2.path <- "ignore/ns.111.fasta"
seq2.unfiltered <- readDNAStringSet(seq2.path)

# make dataframe of variables for each sequence - GC and size for now, could add distance from TSS, etc in the future
seq1.meta <- data.frame(name=names(seq1),size=width(seq1),gc=getGC(seq1))
seq2.meta.unfiltered <- data.frame(name=names(seq2.unfiltered),size=width(seq2.unfiltered),gc=getGC(seq2.unfiltered))

# filter out all sequences longer than or shorter than the target from the background pool
seq2.meta <- seq2.meta.unfiltered[(seq2.meta.unfiltered$size<=max(seq1.meta$size))&(seq2.meta.unfiltered$size>=min(seq1.meta$size)),]

seq2 <- seq2.unfiltered[(seq2.meta.unfiltered$size<=max(seq1.meta$size))&(seq2.meta.unfiltered$size>=min(seq1.meta$size))]

seq2.nSeqs <- length(seq2)

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Add extra variables from annotation

# repeatPer

# distTSS

# distTSE

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Testing of Matching package options and linear models

target.meta <- seq1.meta
pool.meta <- seq2.meta

# setting binary value for group assignment
target.meta$treat <- 1
pool.meta$treat <- 0
all.meta <- rbind(target.meta, pool.meta)

# randomize sort order - order can bias when Match(..., replace=FALSE)
index.random <- sample(seq(1:nrow(all.meta)),nrow(all.meta), replace=FALSE)
all.meta.shuffle <- all.meta[index.random,]

# GC
# run logistic model
lrm.out.gc <- lrm(treat ~ gc, data=all.meta.shuffle)
# obtain values
lrm.out.gc.fitted <- predict.lrm(lrm.out.gc,type="fitted")
# plot
#png(filename="output/propscore.gc.lrm.png",width=1000,height=800,res=120)
#qplot(x=gc,y=treat,data=all.meta.shuffle) + geom_line(aes(x, y), data.frame(x=all.meta.shuffle$gc,y=lrm.out.gc.fitted)) + ggplot.clean()
#dev.off()
# match
rr.gc  <- Match(Y=NULL, Tr=all.meta.shuffle$treat, X=lrm.out.gc.fitted, M=1, version="fast", replace=FALSE)
# check match
#summary(rr.size)
#mb.gc  <- MatchBalance(treat ~ gc, data=all.meta.shuffle, match.out=rr.gc, nboots=10)
# make new sequence set
matched.meta.gc <- all.meta.shuffle[rr.gc$index.control,]
m <- match(as.character(matched.meta.gc$name),names(seq2))
seq.resamp.gc <- seq2[m]

# plot histograms of new sequence set versus old
png(filename="output/dists.filter.resamp.propen.gc.png",width=1000,height=800,res=120)
plotSequenceHistograms(seq1,seq.resamp.gc)
dev.off()

png(filename="output/dists.filter.resamp.propen.gc.overlap.png",width=1400,height=800,res=120)
plotSequenceHistogramsResampled(seq1,seq2,seq.resamp.gc)
dev.off()


# run motif analysis using the resampled background


# Size
#glm.out <- glm(treat ~ size, family=binomial(link=logit), data=all.meta)
lrm.out.size <- lrm(treat ~ size, data=all.meta.shuffle)
#png(filename="output/propscore.size.lrm.png",width=1000,height=800,res=120)
lrm.out.size.fitted <- predict.lrm(lrm.out.size,type="fitted")
#qplot(x=size,y=treat,data=all.meta) + geom_line(aes(x, y), data.frame(x=all.meta$size,y=lrm.out.size.fitted)) + ggplot.clean()
#dev.off()
# match
rr.size  <- Match(Y=NULL, Tr=all.meta.shuffle$treat, X=lrm.out.size.fitted, M=1, version="fast", replace=FALSE)
# check match
#summary(rr.size)
#mb.gc  <- MatchBalance(treat ~ gc, data=all.meta.shuffle, match.out=rr.gc, nboots=10)
# make new sequence set
matched.meta.size <- all.meta.shuffle[rr.size$index.control,]
m <- match(as.character(matched.meta.size$name),names(seq2))
seq.resamp.size <- seq2[m]

# plot histograms of new sequence set versus old
png(filename="output/dists.filter.resamp.propen.size.png",width=1000,height=800,res=120)
plotSequenceHistograms(seq1,seq.resamp.size)
dev.off()

png(filename="output/dists.filter.resamp.propen.size.overlap.png",width=1400,height=800,res=120)
plotSequenceHistogramsResampled(seq1,seq2,seq.resamp.size)
dev.off()

#png(filename="output/dists.resamp.propen.qq.png",width=1400,height=800,res=120)
#plotSequenceQQ4(seq1,seq2,seq.resamp.gc,seq.resamp.size,c("Original","GC","Size"))
#dev.off()

# Size + GC
lrm.out.sizegc <- lrm(treat ~ size + gc, data=all.meta.shuffle)
lrm.out.sizegc.fitted <- predict.lrm(lrm.out.sizegc,type="fitted")
# can't plot without separating dimmensions
#png(filename="output/propscore.size-gc.lrm.png",width=1000,height=800,res=120)
#qplot(x=size+gc,y=treat,data=all.meta) + geom_point(aes(x, y), data.frame(x=all.meta$size+all.meta$gc,y=1-glm.out$fitted))
#dev.off()
# match
rr.sizegc  <- Match(Y=NULL, Tr=all.meta.shuffle$treat, X=lrm.out.sizegc.fitted, M=1, version="fast", replace=FALSE)
# check match
#summary(rr.size)
#mb.gc  <- MatchBalance(treat ~ gc, data=all.meta.shuffle, match.out=rr.gc, nboots=10)
# make new sequence set
matched.meta.sizegc <- all.meta.shuffle[rr.sizegc$index.control,]
m <- match(as.character(matched.meta.sizegc$name),names(seq2))
seq.resamp.sizegc <- seq2[m]

# plot histograms of new sequence set versus old
png(filename="output/dists.filter.resamp.propen.sizegc.png",width=1000,height=800,res=120)
plotSequenceHistograms(seq1,seq.resamp.sizegc)
dev.off()

png(filename="output/dists.filter.resamp.propen.sizegc.overlap.png",width=1400,height=800,res=120)
plotSequenceHistogramsResampled(seq1,seq2,seq.resamp.sizegc)
dev.off()

# Plots


png(filename="output/dists.filter.resamp.propen.qqr2.png",width=1400,height=800,res=120)
plotSequenceQQR5(seq1,seq2,seq.resamp.gc,seq.resamp.size,seq.resamp.sizegc,c("Original","GC","Size","Size+GC"))
dev.off()

save.image(file="ignore/propensity.Rd")




# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Run binomial test on the matches we made

motif.path <- "../TRANSFAC/ConvertToMEME/TRAJAS.Human.named.meme"
q.cutoff <- 0.05

# read in FIMO run on target
run.name <- "CIMPHyperMe"
nSeqs <- length(seq)
#runFIMO(run.name,seq.path,motif.path)
fimo.out <- readFIMO(paste("output/fimo_out_",run.name,"/fimo.txt",sep=""))
fimo.out.counts <- calcMotifCounts(fimo.out,q.cutoff)


# must have resampled seq set at this point
seq.target <- seq1
seq.background <- seq.resamp.sizegc

sim.nSeqs <- length(seq.background)

unlink("output/simseq.fasta")
writeXStringSet(seq.background, "output/simseq.fasta", append=FALSE, format="fasta")

#unlink("output/fimo_out_sim")
#runFIMO("sim","output/simseq.fasta",motif.path)
fimo.out.sim <- readFIMO("output/fimo_out_sim/fimo.txt")
fimo.out.sim.counts <- calcMotifCounts(fimo.out.sim,q.cutoff)
results <- calcEnrichmentBinom(fimo.out.counts,nSeqs,fimo.out.sim.counts,sim.nSeqs)
results.sort <- results[order(results$pvalue, decreasing=FALSE),]
head(results.sort)
tail(results.sort)
write.csv(results.sort,file="output/diffmotif.bkg.propensity.sizegc.csv",row.names=FALSE)

# -----------------------------------------------------------------------------

# =============================================================================
