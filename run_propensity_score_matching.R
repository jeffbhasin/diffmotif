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

mySeqMeta <- function(myseq)
{
	# name, size and gc - easy
	name <- names(myseq)
	size <- width(myseq)
	gc <- getGC(myseq)

	# add chr/start/end fields parsed from sequence name
	name.split <- data.frame(do.call('rbind', strsplit(as.character(name),'-',fixed=TRUE)))
	names(name.split) <- c("chr","start","end")
	chr <- name.split$chr
	start <- as.numeric(as.character(name.split$start))
	end <- as.numeric(as.character(name.split$end))

	# Add extra variables from annotation

	seq.ranges <- GRanges(seqnames=chr,ranges=IRanges(start=start,end=end))

	repeatPer <- getRepeatPercentFast(seq.ranges,rmsk)
	sizeLog <- log10(size)
	distTSS <- getDistTSS(seq.ranges,ann)
	distTSE <- getDistTSE(seq.ranges,ann)
	distTSSCenter <- getDistTSSCenter(seq.ranges,ann)
	distTSECenter <- getDistTSECenter(seq.ranges,ann)
	distTSSCenterLogX1 <- log10(getDistTSSCenter(seq.ranges,ann)+1)
	distTSECenterLogX1 <- log10(getDistTSECenter(seq.ranges,ann)+1)
	freqCpG <- getFreqCpG(myseq)

	# combine into our covariate dataframe
	#seq.meta <- data.frame(name, chr, start, end, size, sizeLog, gc, freqCpG, repeatPer, distTSS, distTSSCenter, distTSSCenterLogX1, distTSE, distTSECenter, distTSECenterLogX1)
	seq.meta <- data.frame(name, chr, start, end, size, sizeLog, gc, freqCpG, repeatPer, distTSSCenterLogX1, distTSECenterLogX1)
}

myMatching <- function(name, formula)
{
	seq.ref <- drawBackgroundSetPropensity(seq1,seq1.meta,seq2,seq2.meta,formula)

	seq.ref.meta <- mySeqMeta(seq.ref)

	# plot hists
	#png(filename=paste("output/hists.covars.",name,".auto.png",sep=""),width=1400,height=800,res=120)
	#print(plotCovarHistogramsOverlap(seq1.meta,seq.ref.meta,cols,plot.ncols=4))
	#dev.off()

	list(seq=seq.ref,meta=seq.ref.meta)
}

myTesting <- function(name,ref)
{
	# run FIMO
	sim.nSeqs <- length(ref$seq)

	fasta.path <- paste("output/",name,".fasta",sep="")
	unlink(fasta.path)
	writeXStringSet(seq.ref, fasta.path, append=FALSE, format="fasta")

	out.path <- paste("output/fimo_out_",name,sep="")
	unlink(out.path)
	runFIMO(name,fasta.path,motif.path)
	fimo.out.sim <- readFIMO(paste(out.path,"/fimo.txt",sep=""))
	fimo.out.sim.counts <- calcMotifCounts(fimo.out.sim,q.cutoff)

	# binomial test
	results <- calcEnrichmentBinom(fimo.out.counts,seq1.nSeqs,fimo.out.sim.counts,sim.nSeqs)
	#results.sort <- results[order(results$pvalue, decreasing=FALSE),]
	#head(results.sort)
	#tail(results.sort)


	# adjusted p-value	

	results$neglogP <- -log10(results$pvalue)

	csv.path <- paste("output/binom.covars.",name,".csv",sep="")
	write.csv(results,file=csv.path,row.names=FALSE)

	results
}

# =============================================================================

# =============================================================================
# Analysis Code

# -----------------------------------------------------------------------------
# Load seqs

# load annotation for environment
ann <- readUCSCAnnotation(genome="hg18",path="ignore/ucsc/")
rmsk <- readRepeatMasker("hg18","ignore/ucsc/")

# load target set of sequences
seq1.path <- "../RunMEME/ColonSeqGR2012/CIMPHyperMe.fasta"
seq1 <- readDNAStringSet(seq1.path)
seq1.nSeqs <- length(seq1)

# calculate covars of target
seq1.meta <- mySeqMeta(seq1)

# load background "pool" we want to draw from and make match the distribution of variable from the target
seq2.path <- "ignore/ns.111.fasta"
seq2.unfiltered <- readDNAStringSet(seq2.path)

# filter out all sequences longer than or shorter than the target from the background pool
seq2.unfiltered.size <- width(seq2.unfiltered)
seq2 <- seq2.unfiltered[(seq2.unfiltered.size<=max(seq1.meta$size))&(seq2.unfiltered.size>=min(seq1.meta$size))]
seq2.nSeqs <- length(seq2)

# calculate covars of reference pool
seq2.meta <- mySeqMeta(seq2)

# set which cols have covars in them
cols <- 5:ncol(seq1.meta)

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# FIMO on target seqs

motif.path <- "../TRANSFAC/ConvertToMEME/TRAJAS.Human.named.meme"
q.cutoff <- 0.05

# read in FIMO run on target
run.name <- "CIMPHyperMe"
#runFIMO(run.name,seq.path,motif.path)
fimo.out <- readFIMO(paste("output/fimo_out_",run.name,"/fimo.txt",sep=""))
fimo.out.light <- with(fimo.out,data.frame(X.pattern.name,q.value,sequence.name))
fimo.out.counts <- calcMotifCounts(fimo.out,q.cutoff)

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# entire pool

# plot overlapping hists of initial covars
#png(filename="output/hists.covars.orig.png",width=1400,height=800,res=120)
#plotCovarHistogramsOverlap(seq1.meta,seq2.meta,cols,plot.ncols=4)
#dev.off()
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Single Run Example
#name <- "size_gc"
#formula <- as.formula("treat ~ size + gc")

#ref.sizegc <- myMatching(name, formula)
#ref.sizegc$results <- myTesting(ref.sizegc)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Run many formulas in a parallel loop

# size, sizeLog, gc, freqCpG, repeatPer, distTSSCenterLogX1, distTSECenterLogX1

# vector of formulas to use
mine <- c("size", "sizeLog", "sizeLog + gc","freqCpG", "freqCpG + repeatPer", "freqCpG + repeatPer + distTSSCenterLogX1", "freqCpG + repeatPer + distTSECenterLogX1", "sizeLog + gc + freqCpG + repeatPer + distTSSCenterLogX1", "sizeLog + gc + freqCpG + repeatPer + distTSECenterLogX1")

ref.names <- gsub(" ","", mine)
ref.names <- gsub("+","_", ref.names,fixed=TRUE)
mine <- paste("treat ~ ",mine,sep="")

ref <- foreach(i=1:length(mine)) %dopar%
{
	print(paste("Doing ",mine[i],sep=""))
	formula <- as.formula(mine[i])	
	name <- gsub(" ","", as.character(formula)[3])
	name <- gsub("+","_", name,fixed=TRUE)

	myref <- myMatching(name, formula)
	#myref$results <- myTesting(myref.sizegcrep)

	myref
}

ref.results <- foreach(i=1:length(mine)) %dopar%
{
	print(paste("Doing ",mine[i],sep=""))
	#myTesting(ref.names[i], ref[[i]])
}

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Output: Multipage PDF of plots

pdf(file="output/Rmotif.plots.pdf", width=10.5, height=8, paper="USr")
# covar plot
mymeta <- list()
mymeta$pool <- seq2.meta
for(i in 1:length(mine))
{
	mymeta[[ref.names[i]]] <- ref[[i]]$meta
}
plotCovarDistance(seq1.meta, mymeta, cols)
plotCovarQQ(seq1.meta, mymeta, cols, plot.ncols=4)
# overlap histograms

print(plotCovarHistogramsOverlap(seq1.meta,seq2.meta, cols, plot.ncols=4, main="pool"))
for(i in 1:length(mine))
{
	print(plotCovarHistogramsOverlap(seq1.meta, ref[[i]]$meta, cols, plot.ncols=4, main=mine[i]))
}
dev.off()

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Output: Single CSV of test results

#mycsv <- data.frame()
#for(i in 1:length(mine))
#{
#	cbind(mycsv,ref[[i]]$results)
#}
#write.csv(file="",row.names=FALSE)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Output: Save R Workspace

save.image(file="output/run_propensity_score_matching.R")

# -----------------------------------------------------------------------------


# =============================================================================
