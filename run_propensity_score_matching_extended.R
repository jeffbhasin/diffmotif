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

# add chr/start/end fields parsed from sequence name
new <- data.frame(do.call('rbind', strsplit(as.character(seq1.meta$name),'-',fixed=TRUE)))
names(new) <- c("chr","start","end")
seq1.meta <- cbind(seq1.meta,new)
seq1.meta$start <- as.numeric(as.character(seq1.meta$start))
seq1.meta$end <- as.numeric(as.character(seq1.meta$end))


new2 <- data.frame(do.call('rbind', strsplit(as.character(seq2.meta$name),'-',fixed=TRUE)))
names(new2) <- c("chr","start","end")
seq2.meta <- cbind(seq2.meta,new2)
seq2.meta$start <- as.numeric(as.character(seq2.meta$start))
seq2.meta$end <- as.numeric(as.character(seq2.meta$end))

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# FIMO on target seqs

motif.path <- "../TRANSFAC/ConvertToMEME/TRAJAS.Human.named.meme"
q.cutoff <- 0.05

# read in FIMO run on target
run.name <- "CIMPHyperMe"
#runFIMO(run.name,seq.path,motif.path)
fimo.out <- readFIMO(paste("output/fimo_out_",run.name,"/fimo.txt",sep=""))
fimo.out.counts <- calcMotifCounts(fimo.out,q.cutoff)

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Add extra variables from annotation

seq1.ranges <- with(seq1.meta,GRanges(seqnames=chr,ranges=IRanges(start=start,end=end)))
seq2.ranges <- with(seq2.meta,GRanges(seqnames=chr,ranges=IRanges(start=start,end=end)))

# read annotation from disk
ann <- readUCSCAnnotation(genome="hg18",path="ignore/ucsc/")
rmsk <- readRepeatMasker("hg18","ignore/ucsc/")

# repeatPer
#seq1.repeatPer <- getRepeatPercentFast(seq1.ranges,rmsk)
#seq2.repeatPer <- getRepeatPercentFast(seq2.ranges,rmsk)
#save(seq1.repeatPer, seq2.repeatPer, file="ignore/propenex.Rd")
load("ignore/propenex.Rd")

# distTSS
seq1.distTSS <- getDistTSS(seq1.ranges,ann)
seq2.distTSS <- getDistTSS(seq2.ranges,ann)

# distTSE
seq1.distTSE <- getDistTSE(seq1.ranges,ann)
seq2.distTSE <- getDistTSE(seq2.ranges,ann)

# distTSSCenter
seq1.distTSSCenter <- getDistTSSCenter(seq1.ranges,ann)
seq2.distTSSCenter <- getDistTSSCenter(seq2.ranges,ann)

# distTSECenter
seq1.distTSECenter <- getDistTSECenter(seq1.ranges,ann)
seq2.distTSECenter <- getDistTSECenter(seq2.ranges,ann)

# freqCpG
seq1.freqCpG <- getFreqCpG(seq1)
seq2.freqCpG <- getFreqCpG(seq2)

# combine into our covariate dataframe
seq1.meta <- data.frame(seq1.meta, repeatPer=seq1.repeatPer, distTSS=seq1.distTSS, distTSSCenter=seq1.distTSSCenter, distTSE=seq1.distTSE, distTSECenter=seq1.distTSECenter, freqCpG=seq1.freqCpG)
seq2.meta <- data.frame(seq2.meta, repeatPer=seq2.repeatPer, distTSS=seq2.distTSS, distTSSCenter=seq2.distTSSCenter, distTSE=seq2.distTSE, distTSECenter=seq2.distTSECenter, freqCpG=seq2.freqCpG)

# plot all variables
png(filename="output/hists.covars.orig.png",width=1000,height=800,res=120)
plotCovarHistograms(seq1.meta,cols)
dev.off()

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# treat ~ gc + length

# draw seqs
formula <- as.formula("treat ~ size + gc")
seq.ref <- drawBackgroundSetPropensity(seq1,seq1.meta,seq2,seq2.meta,formula)

# calculate covars for drawn seqs
seq.ref.meta <- data.frame(name=names(seq.ref),size=width(seq.ref),gc=getGC(seq.ref))
new <- data.frame(do.call('rbind', strsplit(as.character(seq.ref.meta$name),'-',fixed=TRUE)))
names(new) <- c("chr","start","end")
seq.ref.meta <- cbind(seq.ref.meta,new)
seq.ref.meta$start <- as.numeric(as.character(seq.ref.meta$start))
seq.ref.meta$end <- as.numeric(as.character(seq.ref.meta$end))
seq.ref.ranges <- with(seq.ref.meta,GRanges(seqnames=chr,ranges=IRanges(start=start,end=end)))
seq.ref.repeatPer <- getRepeatPercentFast(seq.ref.ranges,rmsk)
seq.ref.distTSS <- getDistTSS(seq.ref.ranges,ann)
seq.ref.distTSE <- getDistTSE(seq.ref.ranges,ann)
seq.ref.distTSSCenter <- getDistTSSCenter(seq.ref.ranges,ann)
seq.ref.distTSECenter <- getDistTSECenter(seq.ref.ranges,ann)
seq.ref.freqCpG <- getFreqCpG(seq.ref)

# combine into our covariate dataframe
seq.ref.meta.all <- data.frame(seq.ref.meta, repeatPer=seq.ref.repeatPer, distTSS=seq.ref.distTSS, distTSSCenter=seq.ref.distTSSCenter, distTSE=seq.ref.distTSE, distTSECenter=seq.ref.distTSECenter, freqCpG=seq.ref.freqCpG)

# plot all variables
cols <- c(2,3,7,8,9,10,11,12)
png(filename="output/hists.covars.sizegc.png",width=1000,height=800,res=120)
plotCovarHistograms(seq.ref.meta.all,cols)
dev.off()

# run FIMO
sim.nSeqs <- length(seq.ref)
unlink("output/simseq.fasta")
writeXStringSet(seq.ref, "output/simseq.fasta", append=FALSE, format="fasta")
unlink("output/fimo_out_sim")
runFIMO("sim","output/simseq.fasta",motif.path)
fimo.out.sim <- readFIMO("output/fimo_out_sim/fimo.txt")
fimo.out.sim.counts <- calcMotifCounts(fimo.out.sim,q.cutoff)

# binomial test
results <- calcEnrichmentBinom(fimo.out.counts,seq1.nSeqs,fimo.out.sim.counts,sim.nSeqs)
results.sort <- results[order(results$pvalue, decreasing=FALSE),]
head(results.sort)
tail(results.sort)
write.csv(results.sort,file="output/diffmotif.bkg.propensity.sizegc.usingfunction.csv",row.names=FALSE)

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# treat ~ gc + length + distTSS

# draw seqs
formula <- as.formula("treat ~ size + gc + distTSS")
seq.ref <- drawBackgroundSetPropensity(seq1,seq1.meta,seq2,seq2.meta,formula)

# run FIMO
sim.nSeqs <- length(seq.ref)
unlink("output/simseq.fasta")
writeXStringSet(seq.ref, "output/simseq.fasta", append=FALSE, format="fasta")
unlink("output/fimo_out_sim")
runFIMO("sim","output/simseq.fasta",motif.path)
fimo.out.sim <- readFIMO("output/fimo_out_sim/fimo.txt")
fimo.out.sim.counts <- calcMotifCounts(fimo.out.sim,q.cutoff)

# binomial test
results <- calcEnrichmentBinom(fimo.out.counts,seq1.nSeqs,fimo.out.sim.counts,sim.nSeqs)
results.sort <- results[order(results$pvalue, decreasing=FALSE),]
head(results.sort)
tail(results.sort)
write.csv(results.sort,file="output/diffmotif.bkg.propensity.sizegc.usingfunction.csv",row.names=FALSE)

# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# treat ~ gc + length + repeatPer

# draw seqs
formula <- as.formula("treat ~ size + gc + repeatPer")
seq.ref <- drawBackgroundSetPropensity(seq1,seq1.meta,seq2,seq2.meta,formula)

# run FIMO
sim.nSeqs <- length(seq.ref)
unlink("output/simseq.fasta")
writeXStringSet(seq.ref, "output/simseq.fasta", append=FALSE, format="fasta")
unlink("output/fimo_out_sim")
runFIMO("sim","output/simseq.fasta",motif.path)
fimo.out.sim <- readFIMO("output/fimo_out_sim/fimo.txt")
fimo.out.sim.counts <- calcMotifCounts(fimo.out.sim,q.cutoff)

# binomial test
results <- calcEnrichmentBinom(fimo.out.counts,seq1.nSeqs,fimo.out.sim.counts,sim.nSeqs)
results.sort <- results[order(results$pvalue, decreasing=FALSE),]
head(results.sort)
tail(results.sort)
write.csv(results.sort,file="output/diffmotif.bkg.propensity.size-gc-repeat.csv",row.names=FALSE)


# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# treat ~ gc + length + repeatPer + distTSS

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# treat ~ gc + length + repeatPer + distTSE

# -----------------------------------------------------------------------------

# =============================================================================
