# #############################################################################
# Test code for implementing propensity score matching
# Author: Jeffrey Bhasin <jmb85@case.edu>
# Created: 2013-02-27
# #############################################################################

# This version runs multiple trials of models to look at the variation induced by different runs of Match()

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
#seq2.meta <- mySeqMeta(seq2)
#save(seq2.meta,file="output/propensity/seq2.meta.Rd")
load("output/propensity/seq2.meta.Rd")

# set which cols have covars in them
cols <- 5:ncol(seq1.meta)

#save(seq1.meta, cols, file="output/propensity/seq1.meta.Rd")

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
# Draw Random Starting Order - Use same start order for all runs

#nTotalSeqs <- nrow(seq1.meta) + nrow(seq2.meta)
#index.random <- sample(seq(1:nTotalSeqs),nTotalSeqs, replace=FALSE)
#save(index.random,file="output/propensity/30starts/index.random.Rd")

# load in so we can replicate the plots
load("output/propensity/30starts/index.random.Rd")
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Multiple Draws of Random Subsets

nDraws <- 30

sub.mine <- rep("subsamples", nDraws)

#mine <- c("size", "gc")

ref.sub.names <- gsub(" ","", mine)
ref.sub.names <- gsub("+","_", ref.sub.names,fixed=TRUE)
#ref.names <- paste(ref.names,seq(1,length(ref.names)),sep="-perm")
ref.sub.names <- paste("poolsubsamples",seq(1,length(ref.sub.names)),sep="-perm")
#mine <- paste("treat ~ ",mine,sep="")

ref.sub <- foreach(i=1:length(ref.sub.names)) %dopar%
{
	print(paste("Sampling ",ref.sub.names[i],sep=""))
	name <- ref.sub.names[i]

	#myref <- myMatching(name, formula)
	#myref$results <- myTesting(myref.sizegcrep)

	mydraws <- sample.int(length(seq2), size=length(seq1),replace=FALSE)
	myseq <- seq2[mydraws]
	mymeta <- mySeqMeta(myseq)

	myref <- list(seq=myseq,meta=mymeta)

	myref
}

pdf(file="output/propensity/30starts/ModelPerformance.poolsubsamples.pdf", width=10.5, height=8, paper="USr")
# covar plot
mymeta <- list()
mymeta$pool <- seq2.meta
for(i in 1:length(sub.mine))
{
	mymeta[[ref.sub.names[i]]] <- ref.sub[[i]]$meta
}
plotCovarDistance(seq1.meta, mymeta, cols)
plotCovarQQ(seq1.meta, mymeta, cols, plot.ncols=4)
# overlap histograms

print(plotCovarHistogramsOverlap(seq1.meta,seq2.meta, cols, plot.ncols=4, main="pool"))
for(i in 1:length(sub.mine))
{
	print(plotCovarHistogramsOverlap(seq1.meta, ref.sub[[i]]$meta, cols, plot.ncols=4, main=ref.sub.names[i]))
}
dev.off()

# FIMO RUNS: build adjusted pvalue matrix
ref.sub.results <- foreach(i=1:length(sub.mine)) %dopar%
{
	print(paste("Testing ",ref.sub.names[i],sep=""))
	myTesting(prefix="output/propensity/30starts/fimo/",ref.sub.names[i], ref.sub[[i]])
}

save(ref.sub,ref.sub.results,sub.mine,ref.sub.names,file="output/propensity/30starts/ref.sub.Rd")

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Multiple Draws of sizeLog

sizelog.mine <- rep("sizeLog", 30)

ref.sizelog.names <- gsub(" ","", sizelog.mine)
ref.sizelog.names <- gsub("+","_", ref.sizelog.names,fixed=TRUE)
ref.sizelog.names <- paste("sizeLog",seq(1,length(ref.sizelog.names)),sep="-perm")
sizelog.mine <- paste("treat ~ ",sizelog.mine,sep="")

ref.sizelog <- foreach(i=1:length(sizelog.mine)) %dopar%
{
	print(paste("Matching ",ref.sizelog.names[i],sep=""))
	formula <- as.formula(sizelog.mine[i])	
	name <- gsub(" ","", as.character(formula)[3])
	name <- gsub("+","_", name,fixed=TRUE)
	name <- ref.sizelog.names[i]
	myref <- myMatching(name, formula)
	#myref$results <- myTesting(myref.sizelog.sizegcrep)

	myref
}

pdf(file="output/propensity/30starts/ModelPerformance.sizelog.pdf", width=10.5, height=8, paper="USr")
# covar plot
mymeta <- list()
mymeta$pool <- seq2.meta
for(i in 1:length(sizelog.mine))
{
	mymeta[[ref.sizelog.names[i]]] <- ref.sizelog[[i]]$meta
}
plotCovarDistance(seq1.meta, mymeta, cols)
plotCovarQQ(seq1.meta, mymeta, cols, plot.ncols=4)
# overlap histograms

print(plotCovarHistogramsOverlap(seq1.meta,seq2.meta, cols, plot.ncols=4, main="pool"))
for(i in 1:length(sizelog.mine))
{
	print(plotCovarHistogramsOverlap(seq1.meta, ref.sizelog[[i]]$meta, cols, plot.ncols=4, main=ref.sizelog.names[i]))
}
dev.off()

# FIMO RUNS: build adjusted pvalue matrix
ref.sizelog.results <- foreach(i=1:length(sizelog.mine)) %dopar%
{
	print(paste("Testing ",sizelog.mine[i],sep=""))
	myTesting(prefix="output/propensity/30starts/fimo/",ref.sizelog.names[i], ref.sizelog[[i]])
}

save(ref.sizelog,ref.sizelog.results,sizelog.mine,ref.sizelog.names,file="output/propensity/30starts/ref.sizelog.Rd")

# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Multiple Draws of all variables (and distTSS or distTSE)

alltss.mine <- rep("sizeLog + gc + freqCpG + repeatPer + distTSSCenterLogX1", 30)

ref.alltss.names <- gsub(" ","", alltss.mine)
ref.alltss.names <- gsub("+","_", ref.alltss.names,fixed=TRUE)
ref.alltss.names <- paste("alltss",seq(1,length(ref.alltss.names)),sep="-perm")
alltss.mine <- paste("treat ~ ",alltss.mine,sep="")

ref.alltss <- foreach(i=1:length(alltss.mine)) %dopar%
{
	print(paste("Matching ",ref.alltss.names[i],sep=""))
	formula <- as.formula(alltss.mine[i])	
	name <- gsub(" ","", as.character(formula)[3])
	name <- gsub("+","_", name,fixed=TRUE)
	name <- ref.alltss.names[i]
	myref <- myMatching(name, formula)
	#myref$results <- myTesting(myref.alltss.sizegcrep)

	myref
}

pdf(file="output/propensity/30starts/ModelPerformance.alltss.pdf", width=10.5, height=8, paper="USr")
# covar plot
mymeta <- list()
mymeta$pool <- seq2.meta
for(i in 1:length(alltss.mine))
{
	mymeta[[ref.alltss.names[i]]] <- ref.alltss[[i]]$meta
}
plotCovarDistance(seq1.meta, mymeta, cols)
plotCovarQQ(seq1.meta, mymeta, cols, plot.ncols=4)
# overlap histograms

print(plotCovarHistogramsOverlap(seq1.meta,seq2.meta, cols, plot.ncols=4, main="pool"))
for(i in 1:length(alltss.mine))
{
	print(plotCovarHistogramsOverlap(seq1.meta, ref.alltss[[i]]$meta, cols, plot.ncols=4, main=ref.alltss.names[i]))
}
dev.off()

# FIMO RUNS: build adjusted pvalue matrix
ref.alltss.results <- foreach(i=1:length(alltss.mine)) %dopar%
{
	print(paste("Testing ",alltss.mine[i],sep=""))
	myTesting(prefix="output/propensity/30starts/fimo/",ref.alltss.names[i], ref.alltss[[i]])
}

save(ref.alltss,ref.alltss.results,alltss.mine,ref.alltss.names,file="output/propensity/30starts/ref.alltss.Rd")

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Multiple Draws of all variables (and distTSS or distTSE)

alltse.mine <- rep("sizeLog + gc + freqCpG + repeatPer + distTSECenterLogX1", 30)

ref.alltse.names <- gsub(" ","", alltse.mine)
ref.alltse.names <- gsub("+","_", ref.alltse.names,fixed=TRUE)
ref.alltse.names <- paste("alltse",seq(1,length(ref.alltse.names)),sep="-perm")
alltse.mine <- paste("treat ~ ",alltse.mine,sep="")

ref.alltse <- foreach(i=1:length(alltse.mine)) %dopar%
{
	print(paste("Matching ",ref.alltse.names[i],sep=""))
	formula <- as.formula(alltse.mine[i])	
	name <- gsub(" ","", as.character(formula)[3])
	name <- gsub("+","_", name,fixed=TRUE)
	name <- ref.alltse.names[i]
	myref <- myMatching(name, formula)
	#myref$results <- myTesting(myref.alltse.sizegcrep)

	myref
}

pdf(file="output/propensity/30starts/ModelPerformance.alltse.pdf", width=10.5, height=8, paper="USr")
# covar plot
mymeta <- list()
mymeta$pool <- seq2.meta
for(i in 1:length(alltse.mine))
{
	mymeta[[ref.alltse.names[i]]] <- ref.alltse[[i]]$meta
}
plotCovarDistance(seq1.meta, mymeta, cols)
plotCovarQQ(seq1.meta, mymeta, cols, plot.ncols=4)
# overlap histograms

print(plotCovarHistogramsOverlap(seq1.meta,seq2.meta, cols, plot.ncols=4, main="pool"))
for(i in 1:length(alltse.mine))
{
	print(plotCovarHistogramsOverlap(seq1.meta, ref.alltse[[i]]$meta, cols, plot.ncols=4, main=ref.alltse.names[i]))
}
dev.off()

# FIMO RUNS: build adjusted pvalue matrix
ref.alltse.results <- foreach(i=1:length(alltse.mine)) %dopar%
{
	print(paste("Testing ",alltse.mine[i],sep=""))
	myTesting(prefix="output/propensity/30starts/fimo/",ref.alltse.names[i], ref.alltse[[i]])
}

save(ref.alltse,ref.alltse.results,alltse.mine,ref.alltse.names,file="output/propensity/30starts/ref.alltse.Rd")

# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Plot Combined Heatmap with custom clustering orders

# Want 3 Runs: Random Draws from Pool, sizeLog only, and All Varibles+TSS

#HEATMAP

#merge together our cols:

#subsamps
val <- as.numeric(ref.sub.results[[1]]$p.adj)
val2 <- -log10(val)
val2[val2=="Inf"] <- max(val2[val2!="Inf"])

mat <- matrix(val2)
for(i in 2:length(ref.sub.names))
{
	val <- as.numeric(ref.sub.results[[i]]$p.adj)
	val2 <- -log10(val)
	val2[val2=="Inf"] <- max(val2[val2!="Inf"])

	mat <- cbind(mat,val2)
}
rownames(mat) <- ref.sub.results[[1]]$motif
colnames(mat) <- ref.sub.names
mat1 <- mat

#sizelog
val <- as.numeric(ref.sizelog.results[[1]]$p.adj)
val2 <- -log10(val)
val2[val2=="Inf"] <- max(val2[val2!="Inf"])

mat <- matrix(val2)
for(i in 2:length(ref.sizelog.names))
{
	val <- as.numeric(ref.sizelog.results[[i]]$p.adj)
	val2 <- -log10(val)
	val2[val2=="Inf"] <- max(val2[val2!="Inf"])

	mat <- cbind(mat,val2)
}
rownames(mat) <- ref.sizelog.results[[1]]$motif
colnames(mat) <- ref.sizelog.names
mat2 <- mat

#alltss
val <- as.numeric(ref.alltss.results[[1]]$p.adj)
val2 <- -log10(val)
val2[val2=="Inf"] <- max(val2[val2!="Inf"])

mat <- matrix(val2)
for(i in 2:length(ref.alltss.names))
{
	val <- as.numeric(ref.alltss.results[[i]]$p.adj)
	val2 <- -log10(val)
	val2[val2=="Inf"] <- max(val2[val2!="Inf"])

	mat <- cbind(mat,val2)
}
rownames(mat) <- ref.alltss.results[[1]]$motif
colnames(mat) <- ref.alltss.names
mat3 <- mat

mat1.un <- mat1
mat2.un <- mat2
mat3.un <- mat3

#cluster the columns only within their own groups
mydis <- dist(t(mat1), method="euclidean")
myden <- as.dendrogram(hclust(mydis))
#order.dendrogram(myden)
myden.order <- reorder(myden, -1*rowMeans(t(mat1)[order.dendrogram(myden),]))
mat1 <- mat1[,order.dendrogram(myden.order)]

mydis <- dist(t(mat2), method="euclidean")
myden <- as.dendrogram(hclust(mydis))
#order.dendrogram(myden)
myden.order <- reorder(myden, -1*rowMeans(t(mat2)[order.dendrogram(myden),]))
mat2 <- mat2[,order.dendrogram(myden.order)]

mydis <- dist(t(mat3), method="euclidean")
myden <- as.dendrogram(hclust(mydis))
#order.dendrogram(myden)
myden.order <- reorder(myden, -1*rowMeans(t(mat3)[order.dendrogram(myden),]))
mat3 <- mat3[,order.dendrogram(myden.order)]


#combine into one big matrix
mat <- cbind(mat1,mat2,mat3)


#order clustered by only the alltss stuff
# order rows by mean
mydis <- dist(mat3, method="euclidean")
myden <- as.dendrogram(hclust(mydis))
#order.dendrogram(myden)
myden.order <- reorder(myden, -1*rowMeans(mat3[order.dendrogram(myden),]))
mat <- mat[order.dendrogram(myden.order),]




pdf(file="output/propensity/30starts/LogLikHeatmap.pdf", width=10.5, height=8, paper="USr")
#png(filename="output/tmp.png",width=1400,height=800,res=120)

#binary by sig:
mybreaks <- c(min(mat),seq(from=1.3,to=max(mat), length.out=20))
mycols <- c("#FFFFFF",colorRampPalette(brewer.pal(9,"Blues")[6:9])(length(mybreaks)-2))

#white=NS, ascending blue=more sig
#mybreaks <- c(min(mat),1.3,max(mat))
#mycols <- colorRampPalette(brewer.pal(9,"Blues"))

col.labels <- rep("",ncol(mat))
col.labels[1] <- "sub"
col.labels[31] <- "sizeLog"
col.labels[61] <- "alltss"

hm <- heatmap.2(mat, Rowv = FALSE, Colv=FALSE, dendrogram="none", key=F, keysize=1.5, trace="none", scale="none", cexCol=1.2, labRow=NA, breaks=mybreaks, main="Heatmap of -log10(adj p)",labCol=col.labels, col=mycols)
#col=colorRampPalette(brewer.pal(9,"Blues"))
write.table(data.frame(column=seq(1,ncol(mat)),model=colnames(mat)), file="output/propensity/30starts/LogLikHeatmapColKey.txt", row.names=FALSE)
dev.off()


#BINARY
pdf(file="output/propensity/30starts/LogLikHeatmap.binary.pdf", width=10.5, height=8, paper="USr")
#png(filename="output/tmp.png",width=1400,height=800,res=120)

#binary by sig:
#mybreaks <- c(min(mat),seq(from=1.3,to=max(mat), length.out=20))
#mycols <- c("#FFFFFF",colorRampPalette(brewer.pal(9,"Blues")[6:9])(length(mybreaks)-2))

#white=NS, ascending blue=more sig
mybreaks <- c(min(mat),1.3,max(mat))
mycols <- colorRampPalette(brewer.pal(9,"Blues"))


hm <- heatmap.2(mat, Rowv = FALSE, Colv=FALSE, dendrogram="none", key=F, keysize=1.5, trace="none", scale="none", cexCol=1.2, labRow=NA, breaks=mybreaks, main="Heatmap of -log10(adj p)",labCol=seq(1,ncol(mat)), col=mycols)
#col=colorRampPalette(brewer.pal(9,"Blues"))
write.table(data.frame(column=seq(1,ncol(mat)),model=colnames(mat)), file="output/propensity/30starts/LogLikHeatmapColKey.binary.txt", row.names=FALSE)
dev.off()

#output one CSV
mycsv <- with(ref.sub.results[[1]],data.frame(motif, percent_seqs, p.adj))
names(mycsv) <- c("motif", "percent_seqs", ref.sub.names[1])
for(i in 2:length(ref.sub.results))
{
	my.new <- data.frame(ref.sub.results[[i]]$p.adj)
	names(my.new) <- c(ref.sub.names[i])
	mycsv <- cbind(mycsv,my.new)
}
mycsv1 <- mycsv

mycsv <- with(ref.sizelog.results[[1]],data.frame(motif, percent_seqs, p.adj))
names(mycsv) <- c("motif", "percent_seqs", ref.sizelog.names[1])
for(i in 2:length(ref.sizelog.results))
{
	my.new <- data.frame(ref.sizelog.results[[i]]$p.adj)
	names(my.new) <- c(ref.sizelog.names[i])
	mycsv <- cbind(mycsv,my.new)
}
mycsv2 <- mycsv

mycsv <- with(ref.alltss.results[[1]],data.frame(motif, percent_seqs, p.adj))
names(mycsv) <- c("motif", "percent_seqs", ref.alltss.names[1])
for(i in 2:length(ref.alltss.results))
{
	my.new <- data.frame(ref.alltss.results[[i]]$p.adj)
	names(my.new) <- c(ref.alltss.names[i])
	mycsv <- cbind(mycsv,my.new)
}
mycsv3 <- mycsv

mycsv <- cbind(mycsv1,mycsv2,mycsv3)

write.csv(mycsv,file="output/propensity/30starts/DiffMotifs.csv",row.names=FALSE)

#VENN DIAGRAMS

#shared sig factors between the 3 groups
#convert matrix to binary if sig
makeBinary <- function(x)
{
	if(x >= -log10(0.05))
	{
		0
	} else {
		1
	}
}
mat1.bin <- apply(mat1.un,MARGIN=c(1,2),FUN=makeBinary)
mat2.bin <- apply(mat2.un,MARGIN=c(1,2),FUN=makeBinary)
mat3.bin <- apply(mat3.un,MARGIN=c(1,2),FUN=makeBinary)

mat1.sigs <- rowSums(mat1.bin)
mat2.sigs <- rowSums(mat2.bin)
mat3.sigs <- rowSums(mat3.bin)

sigs <- cbind(mat1.sigs,mat2.sigs,mat3.sigs)
sigs <- data.frame(sigs)
names(sigs) <- c("subsamps","sizeLog","allTSS")
write.csv(sigs,file="output/propensity/30starts/DiffMotifs-counts.csv",row.names=TRUE)

names(sigs) <- c("mat1","mat2","mat3")

n1 <- nrow(sigs[(sigs$mat1 > 0),]) 
n2 <- nrow(sigs[(sigs$mat2 > 0),]) 
n3 <-nrow(sigs[(sigs$mat3 > 0),]) 
n12 <- nrow(sigs[(sigs$mat1 > 0) & (sigs$mat2 > 0),])
n23 <- nrow(sigs[(sigs$mat2 > 0) & (sigs$mat3 > 0),])
n13 <- nrow(sigs[(sigs$mat1 > 0) & (sigs$mat3 > 0),])
n123 <- nrow(sigs[(sigs$mat1 > 0) & (sigs$mat2 > 0) & (sigs$mat3 > 0),])

pdf(file="output/propensity/30starts/Venn.pdf", width=10.5, height=8, paper="USr")
venn.plot <- draw.triple.venn(
    area1 = n1,
    area2 = n2,
    area3 = n3,
    n12 = n12,
    n23 = n23,
    n13 = n13,
    n123 = n123,
    category = c("Subsamples", "sizeLog", "allTSS"),
    fill = c("blue", "red", "green"),
    lty = "blank",cex = 2,
cat.cex = 2,
cat.col = c("blue", "red", "green")
);
dev.off()

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Plot Combined Heatmap with custom clustering orders - TSE

# Want 3 Runs: Random Draws from Pool, sizeLog only, and All Varibles+TSE

#HEATMAP

#merge together our cols:

#subsamps
val <- as.numeric(ref.sub.results[[1]]$p.adj)
val2 <- -log10(val)
val2[val2=="Inf"] <- max(val2[val2!="Inf"])

mat <- matrix(val2)
for(i in 2:length(ref.sub.names))
{
	val <- as.numeric(ref.sub.results[[i]]$p.adj)
	val2 <- -log10(val)
	val2[val2=="Inf"] <- max(val2[val2!="Inf"])

	mat <- cbind(mat,val2)
}
rownames(mat) <- ref.sub.results[[1]]$motif
colnames(mat) <- ref.sub.names
mat1 <- mat

#sizelog
val <- as.numeric(ref.sizelog.results[[1]]$p.adj)
val2 <- -log10(val)
val2[val2=="Inf"] <- max(val2[val2!="Inf"])

mat <- matrix(val2)
for(i in 2:length(ref.sizelog.names))
{
	val <- as.numeric(ref.sizelog.results[[i]]$p.adj)
	val2 <- -log10(val)
	val2[val2=="Inf"] <- max(val2[val2!="Inf"])

	mat <- cbind(mat,val2)
}
rownames(mat) <- ref.sizelog.results[[1]]$motif
colnames(mat) <- ref.sizelog.names
mat2 <- mat

#alltse
val <- as.numeric(ref.alltse.results[[1]]$p.adj)
val2 <- -log10(val)
val2[val2=="Inf"] <- max(val2[val2!="Inf"])

mat <- matrix(val2)
for(i in 2:length(ref.alltse.names))
{
	val <- as.numeric(ref.alltse.results[[i]]$p.adj)
	val2 <- -log10(val)
	val2[val2=="Inf"] <- max(val2[val2!="Inf"])

	mat <- cbind(mat,val2)
}
rownames(mat) <- ref.alltse.results[[1]]$motif
colnames(mat) <- ref.alltse.names
mat3 <- mat

mat1.un <- mat1
mat2.un <- mat2
mat3.un <- mat3

#cluster the columns only within their own groups
mydis <- dist(t(mat1), method="euclidean")
myden <- as.dendrogram(hclust(mydis))
#order.dendrogram(myden)
myden.order <- reorder(myden, -1*rowMeans(t(mat1)[order.dendrogram(myden),]))
mat1 <- mat1[,order.dendrogram(myden.order)]

mydis <- dist(t(mat2), method="euclidean")
myden <- as.dendrogram(hclust(mydis))
#order.dendrogram(myden)
myden.order <- reorder(myden, -1*rowMeans(t(mat2)[order.dendrogram(myden),]))
mat2 <- mat2[,order.dendrogram(myden.order)]

mydis <- dist(t(mat3), method="euclidean")
myden <- as.dendrogram(hclust(mydis))
#order.dendrogram(myden)
myden.order <- reorder(myden, -1*rowMeans(t(mat3)[order.dendrogram(myden),]))
mat3 <- mat3[,order.dendrogram(myden.order)]


#combine into one big matrix
mat <- cbind(mat1,mat2,mat3)


#order clustered by only the alltse stuff
# order rows by mean
mydis <- dist(mat3, method="euclidean")
myden <- as.dendrogram(hclust(mydis))
#order.dendrogram(myden)
myden.order <- reorder(myden, -1*rowMeans(mat3[order.dendrogram(myden),]))
mat <- mat[order.dendrogram(myden.order),]




pdf(file="output/propensity/30startsTSE/LogLikHeatmap.pdf", width=10.5, height=8, paper="USr")
#png(filename="output/tmp.png",width=1400,height=800,res=120)

#binary by sig:
mybreaks <- c(min(mat),seq(from=1.3,to=max(mat), length.out=20))
mycols <- c("#FFFFFF",colorRampPalette(brewer.pal(9,"Blues")[6:9])(length(mybreaks)-2))

#white=NS, ascending blue=more sig
#mybreaks <- c(min(mat),1.3,max(mat))
#mycols <- colorRampPalette(brewer.pal(9,"Blues"))

col.labels <- rep("",ncol(mat))
col.labels[1] <- "sub"
col.labels[31] <- "sizeLog"
col.labels[61] <- "alltse"

hm <- heatmap.2(mat, Rowv = FALSE, Colv=FALSE, dendrogram="none", key=F, keysize=1.5, trace="none", scale="none", cexCol=1.2, labRow=NA, breaks=mybreaks, main="Heatmap of -log10(adj p)",labCol=col.labels, col=mycols)
#col=colorRampPalette(brewer.pal(9,"Blues"))
write.table(data.frame(column=seq(1,ncol(mat)),model=colnames(mat)), file="output/propensity/30startsTSE/LogLikHeatmapColKey.txt", row.names=FALSE)
dev.off()


#BINARY
pdf(file="output/propensity/30startsTSE/LogLikHeatmap.binary.pdf", width=10.5, height=8, paper="USr")
#png(filename="output/tmp.png",width=1400,height=800,res=120)

#binary by sig:
#mybreaks <- c(min(mat),seq(from=1.3,to=max(mat), length.out=20))
#mycols <- c("#FFFFFF",colorRampPalette(brewer.pal(9,"Blues")[6:9])(length(mybreaks)-2))

#white=NS, ascending blue=more sig
mybreaks <- c(min(mat),1.3,max(mat))
mycols <- colorRampPalette(brewer.pal(9,"Blues"))


hm <- heatmap.2(mat, Rowv = FALSE, Colv=FALSE, dendrogram="none", key=F, keysize=1.5, trace="none", scale="none", cexCol=1.2, labRow=NA, breaks=mybreaks, main="Heatmap of -log10(adj p)",labCol=seq(1,ncol(mat)), col=mycols)
#col=colorRampPalette(brewer.pal(9,"Blues"))
write.table(data.frame(column=seq(1,ncol(mat)),model=colnames(mat)), file="output/propensity/30startsTSE/LogLikHeatmapColKey.binary.txt", row.names=FALSE)
dev.off()

#output one CSV
mycsv <- with(ref.sub.results[[1]],data.frame(motif, percent_seqs, p.adj))
names(mycsv) <- c("motif", "percent_seqs", ref.sub.names[1])
for(i in 2:length(ref.sub.results))
{
	my.new <- data.frame(ref.sub.results[[i]]$p.adj)
	names(my.new) <- c(ref.sub.names[i])
	mycsv <- cbind(mycsv,my.new)
}
mycsv1 <- mycsv

mycsv <- with(ref.sizelog.results[[1]],data.frame(motif, percent_seqs, p.adj))
names(mycsv) <- c("motif", "percent_seqs", ref.sizelog.names[1])
for(i in 2:length(ref.sizelog.results))
{
	my.new <- data.frame(ref.sizelog.results[[i]]$p.adj)
	names(my.new) <- c(ref.sizelog.names[i])
	mycsv <- cbind(mycsv,my.new)
}
mycsv2 <- mycsv

mycsv <- with(ref.alltse.results[[1]],data.frame(motif, percent_seqs, p.adj))
names(mycsv) <- c("motif", "percent_seqs", ref.alltse.names[1])
for(i in 2:length(ref.alltse.results))
{
	my.new <- data.frame(ref.alltse.results[[i]]$p.adj)
	names(my.new) <- c(ref.alltse.names[i])
	mycsv <- cbind(mycsv,my.new)
}
mycsv3 <- mycsv

mycsv <- cbind(mycsv1,mycsv2,mycsv3)

write.csv(mycsv,file="output/propensity/30startsTSE/DiffMotifs.csv",row.names=FALSE)

#VENN DIAGRAMS

#shared sig factors between the 3 groups
#convert matrix to binary if sig
makeBinary <- function(x)
{
	if(x >= -log10(0.05))
	{
		0
	} else {
		1
	}
}
mat1.bin <- apply(mat1.un,MARGIN=c(1,2),FUN=makeBinary)
mat2.bin <- apply(mat2.un,MARGIN=c(1,2),FUN=makeBinary)
mat3.bin <- apply(mat3.un,MARGIN=c(1,2),FUN=makeBinary)

mat1.sigs <- rowSums(mat1.bin)
mat2.sigs <- rowSums(mat2.bin)
mat3.sigs <- rowSums(mat3.bin)

sigs <- cbind(mat1.sigs,mat2.sigs,mat3.sigs)
sigs <- data.frame(sigs)
names(sigs) <- c("subsamps","sizeLog","allTSE")
write.csv(sigs,file="output/propensity/30startsTSE/DiffMotifs-counts.csv",row.names=TRUE)

names(sigs) <- c("mat1","mat2","mat3")

n1 <- nrow(sigs[(sigs$mat1 > 0),]) 
n2 <- nrow(sigs[(sigs$mat2 > 0),]) 
n3 <-nrow(sigs[(sigs$mat3 > 0),]) 
n12 <- nrow(sigs[(sigs$mat1 > 0) & (sigs$mat2 > 0),])
n23 <- nrow(sigs[(sigs$mat2 > 0) & (sigs$mat3 > 0),])
n13 <- nrow(sigs[(sigs$mat1 > 0) & (sigs$mat3 > 0),])
n123 <- nrow(sigs[(sigs$mat1 > 0) & (sigs$mat2 > 0) & (sigs$mat3 > 0),])

pdf(file="output/propensity/30startsTSE/Venn.pdf", width=10.5, height=8, paper="USr")
venn.plot <- draw.triple.venn(
    area1 = n1,
    area2 = n2,
    area3 = n3,
    n12 = n12,
    n23 = n23,
    n13 = n13,
    n123 = n123,
    category = c("Subsamples", "sizeLog", "allTSE"),
    fill = c("blue", "red", "green"),
    lty = "blank",cex = 2,
cat.cex = 2,
cat.col = c("blue", "red", "green")
);
dev.off()
# -----------------------------------------------------------------------------


# =============================================================================
