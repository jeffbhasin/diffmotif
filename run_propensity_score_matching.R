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

myMatching <- function(name, formula, index.random)
{
	seq.ref <- drawBackgroundSetPropensity(seq1,seq1.meta,seq2,seq2.meta,formula,start.order=index.random)

	seq.ref.meta <- mySeqMeta(seq.ref)

	# plot hists
	#png(filename=paste("output/hists.covars.",name,".auto.png",sep=""),width=1400,height=800,res=120)
	#print(plotCovarHistogramsOverlap(seq1.meta,seq.ref.meta,cols,plot.ncols=4))
	#dev.off()

	list(seq=seq.ref,meta=seq.ref.meta)
}

myMatchingGen <- function(name, formula)
{
	seq.ref <- drawBackgroundSetPropensityGenMatch(seq1,seq1.meta,seq2,seq2.meta,formula)

	seq.ref.meta <- mySeqMeta(seq.ref)

	# plot hists
	#png(filename=paste("output/hists.covars.",name,".auto.png",sep=""),width=1400,height=800,res=120)
	#print(plotCovarHistogramsOverlap(seq1.meta,seq.ref.meta,cols,plot.ncols=4))
	#dev.off()

	list(seq=seq.ref,meta=seq.ref.meta)
}

myTesting <- function(prefix,name,ref)
{
	# run FIMO
	sim.nSeqs <- length(ref$seq)

	fasta.path <- paste(prefix,name,".fasta",sep="")
	unlink(fasta.path)
	writeXStringSet(ref$seq, fasta.path, append=FALSE, format="fasta")

	out.path <- paste(prefix,"fimo_out_",name,sep="")
	unlink(out.path)
	runFIMO(out.path,fasta.path,motif.path)
	fimo.out.sim <- readFIMO(paste(out.path,"/fimo.txt",sep=""))
	fimo.out.sim.counts <- calcMotifCounts(fimo.out.sim,q.cutoff)

	# binomial test
	results <- calcEnrichmentBinom(fimo.out.counts,seq1.nSeqs,fimo.out.sim.counts,sim.nSeqs)
	#results.sort <- results[order(results$pvalue, decreasing=FALSE),]
	#head(results.sort)
	#tail(results.sort)


	# adjusted p-value	

	results$p.adj <- p.adjust(results$pvalue, method="fdr")
	results$neglogP <- -log10(results$p.adj)

	#csv.path <- paste("output/binom.covars.",name,".csv",sep="")
	#write.csv(results,file=csv.path,row.names=FALSE)

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

#mine <- c("size", "gc")

ref.names <- gsub(" ","", mine)
ref.names <- gsub("+","_", ref.names,fixed=TRUE)
mine <- paste("treat ~ ",mine,sep="")

# take one common draw of the starting list order of the pool now so it is used for all matching operations in the loop below
#nTotalSeqs <- nrow(seq1.meta) + nrow(seq2.meta)
#index.random <- sample(seq(1:nTotalSeqs),nTotalSeqs, replace=FALSE)
#save(index.random,.Random.seed,file="output/index.random.Rd")

# load in so we can replicate the plots
load("output/propensity/1start/index.random.Rd")

ref <- foreach(i=1:length(mine)) %dopar%
{
	print(paste("Matching ",mine[i],sep=""))
	formula <- as.formula(mine[i])	
	name <- gsub(" ","", as.character(formula)[3])
	name <- gsub("+","_", name,fixed=TRUE)

	myref <- myMatching(name, formula, index.random)
	#myref$results <- myTesting(myref.sizegcrep)

	myref
}

# add on a resampled draw (no matching)
mydraws <- sample.int(length(seq2), size=length(seq1),replace=FALSE)
myseq <- seq2[mydraws]
mymeta <- mySeqMeta(myseq)
myref <- list(seq=myseq,meta=mymeta)
ref[[length(ref)+1]] <- myref
ref.names <- c(ref.names, "poolsubsample")

pdf(file="output/propensity/1start/ModelPerformance.pdf", width=10.5, height=8, paper="USr")
# covar plot
mymeta <- list()
mymeta$pool <- seq2.meta
for(i in 1:length(ref.names))
{
	mymeta[[ref.names[i]]] <- ref[[i]]$meta
}
plotCovarDistance(seq1.meta, mymeta, cols)
plotCovarQQ(seq1.meta, mymeta, cols, plot.ncols=4)
# overlap histograms

print(plotCovarHistogramsOverlap(seq1.meta,seq2.meta, cols, plot.ncols=4, main="pool"))
for(i in 1:length(ref.names))
{
	print(plotCovarHistogramsOverlap(seq1.meta, ref[[i]]$meta, cols, plot.ncols=4, main=ref.names[i]))
}
dev.off()

# build adjusted pvalue matrix
ref.results <- foreach(i=1:length(ref.names)) %dopar%
{
	print(paste("Testing ",mine[i],sep=""))
	myTesting(prefix="output/propensity/1start/fimo/",ref.names[i], ref[[i]])
}

#HEATMAP
val <- as.numeric(ref.results[[1]]$p.adj)
val2 <- -log10(val)
val2[val2=="Inf"] <- max(val2[val2!="Inf"])

mat <- matrix(val2)
for(i in 2:length(ref.names))
{
	val <- as.numeric(ref.results[[i]]$p.adj)
	val2 <- -log10(val)
	val2[val2=="Inf"] <- max(val2[val2!="Inf"])

	mat <- cbind(mat,val2)
}
rownames(mat) <- ref.results[[1]]$motif
colnames(mat) <- ref.names

# order rows by mean
mydis <- dist(mat, method="euclidean")
myden <- as.dendrogram(hclust(mydis))
#order.dendrogram(myden)
myden.order <- reorder(myden, -1*rowMeans(mat[order.dendrogram(myden),]))
#order.dendrogram(myden.order)

pdf(file="output/propensity/1start/LogLikHeatmap.pdf", width=10.5, height=8, paper="USr")
#png(filename="output/tmp.png",width=1400,height=800,res=120)

#binary by sig:
mybreaks <- c(min(mat),seq(from=1.3,to=max(mat), length.out=20))
mycols <- c("#FFFFFF",colorRampPalette(brewer.pal(9,"Blues")[6:9])(length(mybreaks)-2))

#white=NS, ascending blue=more sig
#mybreaks <- c(min(mat),1.3,max(mat))
#mycols <- colorRampPalette(brewer.pal(9,"Blues"))


hm <- heatmap.2(mat, Rowv = myden.order, key=F, keysize=1.5, trace="none", scale="none", cexCol=1.2, labRow=NA, breaks=mybreaks, main="Heatmap of -log10(adj p)",labCol=seq(1,ncol(mat)), col=mycols)
#col=colorRampPalette(brewer.pal(9,"Blues"))
write.table(data.frame(column=seq(1,ncol(mat)),model=colnames(mat)), file="output/propensity/1start/LogLikHeatmapColKey.txt", row.names=FALSE)
dev.off()

pdf(file="output/propensity/1start/LogLikHeatmap.binary.pdf", width=10.5, height=8, paper="USr")
#png(filename="output/tmp.png",width=1400,height=800,res=120)

#binary by sig:
#mybreaks <- c(min(mat),seq(from=1.3,to=max(mat), length.out=20))
#mycols <- c("#FFFFFF",colorRampPalette(brewer.pal(9,"Blues")[6:9])(length(mybreaks)-2))

#white=NS, ascending blue=more sig
mybreaks <- c(min(mat),1.3,max(mat))
mycols <- colorRampPalette(brewer.pal(9,"Blues"))


hm <- heatmap.2(mat, Rowv = myden.order, key=F, keysize=1.5, trace="none", scale="none", cexCol=1.2, labRow=NA, breaks=mybreaks, main="Heatmap of -log10(adj p)",labCol=seq(1,ncol(mat)), col=mycols)
#col=colorRampPalette(brewer.pal(9,"Blues"))
write.table(data.frame(column=seq(1,ncol(mat)),model=colnames(mat)), file="output/propensity/1start/LogLikHeatmapColKey.txt", row.names=FALSE)
dev.off()

#output one CSV
mycsv <- with(ref.results[[1]],data.frame(motif, percent_seqs, p.adj))
names(mycsv) <- c("motif", "percent_seqs", ref.names[1])
for(i in 2:length(ref.results))
{
	my.new <- data.frame(ref.results[[i]]$p.adj)
	names(my.new) <- c(ref.names[i])
	mycsv <- cbind(mycsv,my.new)
}
write.csv(mycsv,file="output/propensity/1start/DiffMotifs.csv",row.names=FALSE)

# save the final data so we don't have to re-run fimo, etc
save(ref,ref.results,mine,ref.names,file="output/propensity/1start/ref.Rd")

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Multiple Draws of Same Formula (sensitivity testing)
# vector of formulas to use
mine <- rep("sizeLog + gc + freqCpG + repeatPer + distTSSCenterLogX1", 30)

#mine <- c("size", "gc")

ref.names <- gsub(" ","", mine)
ref.names <- gsub("+","_", ref.names,fixed=TRUE)
#ref.names <- paste(ref.names,seq(1,length(ref.names)),sep="-perm")
ref.names <- paste("allTSS",seq(1,length(ref.names)),sep="-perm")
mine <- paste("treat ~ ",mine,sep="")

ref <- foreach(i=1:length(mine)) %dopar%
{
	print(paste("Matching ",mine[i],sep=""))
	formula <- as.formula(mine[i])	
	name <- gsub(" ","", as.character(formula)[3])
	name <- gsub("+","_", name,fixed=TRUE)
	name <- ref.names[i]
	myref <- myMatching(name, formula)
	#myref$results <- myTesting(myref.sizegcrep)

	myref
}

pdf(file="output/Rmotif.plots2.formulaperm2.pdf", width=10.5, height=8, paper="USr")
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

# FIMO RUNS: build adjusted pvalue matrix
ref.results <- foreach(i=1:length(mine)) %dopar%
{
	print(paste("Testing ",mine[i],sep=""))
	myTesting(ref.names[i], ref[[i]])
}


#HEATMAP
library(gplots)

val <- as.numeric(ref.results[[1]]$p.adj)
val2 <- -log10(val)
val2[val2=="Inf"] <- max(val2[val2!="Inf"])

mat <- matrix(val2)
for(i in 2:length(mine))
{
	val <- as.numeric(ref.results[[i]]$p.adj)
	val2 <- -log10(val)
	val2[val2=="Inf"] <- max(val2[val2!="Inf"])

	mat <- cbind(mat,val2)
}
rownames(mat) <- ref.results[[1]]$motif
colnames(mat) <- ref.names

pdf(file="output/heattest.formulaperm2.pdf", width=10.5, height=8, paper="USr")
#png(filename="output/tmp.png",width=1400,height=800,res=120)
hm <- heatmap.2(mat, key=T, keysize=1.5, trace="none",cexCol=1.2, labRow=NA, col=colorRampPalette(brewer.pal(9,"Blues")), main="Heatmap of -log10(adj p)",labCol=seq(1,ncol(mat)))
data.frame(column=seq(1,ncol(mat)),model=colnames(mat))
dev.off()
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Multiple random subsamples with nDraw==nTarget from the reference pool
# vector of formulas to use
nDraws <- 30

mine <- rep("subsamples", nDraws)

#mine <- c("size", "gc")

ref.names <- gsub(" ","", mine)
ref.names <- gsub("+","_", ref.names,fixed=TRUE)
#ref.names <- paste(ref.names,seq(1,length(ref.names)),sep="-perm")
ref.names <- paste("poolsubsamples",seq(1,length(ref.names)),sep="-perm")
#mine <- paste("treat ~ ",mine,sep="")

ref <- foreach(i=1:length(ref.names)) %dopar%
{
	print(paste("Sampling ",ref.names[i],sep=""))
	name <- ref.names[i]

	#myref <- myMatching(name, formula)
	#myref$results <- myTesting(myref.sizegcrep)

	mydraws <- sample.int(length(seq2), size=length(seq1),replace=FALSE)
	myseq <- seq2[mydraws]
	mymeta <- mySeqMeta(myseq)

	myref <- list(seq=myseq,meta=mymeta)

	myref
}

pdf(file="output/Rmotif.plots.poolsubsamples.pdf", width=10.5, height=8, paper="USr")
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

# FIMO RUNS: build adjusted pvalue matrix
ref.results <- foreach(i=1:length(mine)) %dopar%
{
	print(paste("Testing ",mine[i],sep=""))
	myTesting(ref.names[i], ref[[i]])
}


#HEATMAP
library(gplots)

val <- as.numeric(ref.results[[1]]$p.adj)
val2 <- -log10(val)
val2[val2=="Inf"] <- max(val2[val2!="Inf"])

mat <- matrix(val2)
for(i in 2:length(mine))
{
	val <- as.numeric(ref.results[[i]]$p.adj)
	val2 <- -log10(val)
	val2[val2=="Inf"] <- max(val2[val2!="Inf"])

	mat <- cbind(mat,val2)
}
rownames(mat) <- ref.results[[1]]$motif
colnames(mat) <- ref.names

pdf(file="output/heattest.poolsubsamples.pdf", width=10.5, height=8, paper="USr")
#png(filename="output/tmp.png",width=1400,height=800,res=120)
hm <- heatmap.2(mat, key=T, keysize=1.5, trace="none",cexCol=1.2, labRow=NA, col=colorRampPalette(brewer.pal(9,"Blues")), main="Heatmap of -log10(adj p)",labCol=seq(1,ncol(mat)))
data.frame(column=seq(1,ncol(mat)),model=colnames(mat))
dev.off()
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Multiple Draws of Same Formula (sensitivity testing) - GENETIC MATCHING
# vector of formulas to use
mine <- rep("sizeLog + gc + freqCpG + repeatPer + distTSSCenterLogX1", 30)

#mine <- c("size", "gc")

ref.names <- gsub(" ","", mine)
ref.names <- gsub("+","_", ref.names,fixed=TRUE)
#ref.names <- paste(ref.names,seq(1,length(ref.names)),sep="-perm")
ref.names <- paste("allTSSgen",seq(1,length(ref.names)),sep="-perm")
mine <- paste("treat ~ ",mine,sep="")

ref <- foreach(i=1:length(mine)) %dopar%
{
	print(paste("Matching ",mine[i],sep=""))
	formula <- as.formula(mine[i])	
	name <- gsub(" ","", as.character(formula)[3])
	name <- gsub("+","_", name,fixed=TRUE)
	name <- ref.names[i]
	myref <- myMatchingGen(name, formula)
	#myref$results <- myTesting(myref.sizegcrep)

	myref
}

pdf(file="output/Rmotif.plots2.formulaperm2.gen.pdf", width=10.5, height=8, paper="USr")
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

# FIMO RUNS: build adjusted pvalue matrix
ref.results <- foreach(i=1:length(mine)) %dopar%
{
	print(paste("Testing ",mine[i],sep=""))
	myTesting(ref.names[i], ref[[i]])
}


#HEATMAP
library(gplots)

val <- as.numeric(ref.results[[1]]$p.adj)
val2 <- -log10(val)
val2[val2=="Inf"] <- max(val2[val2!="Inf"])

mat <- matrix(val2)
for(i in 2:length(mine))
{
	val <- as.numeric(ref.results[[i]]$p.adj)
	val2 <- -log10(val)
	val2[val2=="Inf"] <- max(val2[val2!="Inf"])

	mat <- cbind(mat,val2)
}
rownames(mat) <- ref.results[[1]]$motif
colnames(mat) <- ref.names

pdf(file="output/heattest.formulaperm2.gen.pdf", width=10.5, height=8, paper="USr")
#png(filename="output/tmp.png",width=1400,height=800,res=120)
hm <- heatmap.2(mat, key=T, keysize=1.5, trace="none",cexCol=1.2, labRow=NA, col=colorRampPalette(brewer.pal(9,"Blues")), main="Heatmap of -log10(adj p)",labCol=seq(1,ncol(mat)))
data.frame(column=seq(1,ncol(mat)),model=colnames(mat))
dev.off()
# -----------------------------------------------------------------------------


# =============================================================================
