###############################################################################
## Finds differentially enriched motifs in a FIMO output
##
###############################################################################

###############################################
## Packages
###############################################
library(ShortRead)
library(plyr)
library(reshape)
library(foreach)
library(stringr)
library(ggplot2)
library(gridExtra)
library(Matching)
library(rms)
library(doMC)

###############################################
## Local Dependencies
###############################################
source("GenomeInfo.R")

###############################################
## Utility
###############################################
ggplot.clean <- function()
{
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.key.size = unit(0.8, "lines"), axis.line = element_line(colour = "grey50"))
}

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

###############################################
## MEME Suite Interaction
###############################################
runFIMO <- function(runname,fasta.path,motifs.path)
{
	system(paste("fimo -oc output/fimo_out_",runname," ",motifs.path," ",fasta.path,sep=""))
}

readFIMO <- function(fimo.out.path)
{
	read.table(file=fimo.out.path,header=TRUE,comment.char="",sep="\t")
}


###############################################
## Background Sequence Generation by Drawing
###############################################
drawBackgroundSet <- function(seq,nSimSeqs=10000,windowSize=50)
{
	#input: sequence set of originals read in from FASTA as a DNAStringSet
	#output: drawn background sequence as DNAStringSet

	#TODO: fix bug if chr freq = 0

	genome.path <- "../RunMEME/ColonSeqGR2012/hg18.fa"

	#create distribution of sequences sizes for simulated set by drawing from real set
	sim.sizes <- sample(width(seq) + windowSize, nSimSeqs, replace=TRUE)

	#draw chr assignments for sim seqs

	#parse out chr names and positions from string set names
	#TODO: should take an annotation data frame with chr/start/end rather than parsing from the DNAStringSet names (this should be done by the user in the main body code depending on how they have it formatted)
	seq.chr.names <- str_split_fixed(names(seq),"-",n=3)[,1]
	seq.chr.starts <- str_split_fixed(names(seq),"-",n=3)[,2]
	seq.chr.ends <- str_split_fixed(names(seq),"-",n=3)[,3]
	seq.annot <- data.frame(chr=seq.chr.names,start=seq.chr.starts,end=seq.chr.ends,stringsAsFactors=FALSE)

	#find max end point for any seq on each chr in set
	seq.chr.range.max <- by(seq.annot, seq.annot$chr, function(x) max(x$end))
	#count number of seqs from each chr in set
	seq.chr.table <- table(seq.annot$chr)

	#make random draws of chr names for the simulated seq set using frequencies from the actual set
	sim.chrs <- sample(names(seq.chr.table), nSimSeqs, replace=TRUE, prob=as.vector(seq.chr.table)/sum(seq.chr.table))

	#data frame of chr frequencies in simulated set we just drew and max end points for each chr
	sim.chr.meta <- data.frame(cbind(freq=as.vector(table(sim.chrs)), max=seq.chr.range.max))

	#make random draw of chr positions - for each chr, draw a position for each seq from all integers from 1 to the max endpoint seen in the real data
	sim.pos <- apply(sim.chr.meta, 1, function(x) sample.int(x[2], x[1]))

	sim.chr.names <- row.names(sim.chr.meta)
	sim.chr.freqs <- as.vector(sim.chr.meta$freq)

	#generate the corresponding chr names for each simulated position (repeat the chr name by the freq in the simulated draw)
	sim.pos.chrs <- sapply(seq(sim.chr.freqs), function(x, sim.chr.names, sim.chr.freqs) rep(sim.chr.names[x], sim.chr.freqs[x]), sim.chr.names=sim.chr.names, sim.chr.freqs=sim.chr.freqs)

	#build data frame with position information for the simulated seq draws
	sim.annot <- data.frame(chr=unlist(sim.pos.chrs), end=unlist(sim.pos), stringsAsFactors=FALSE)
	
	#create start points by subtracting the random draws of sizes from the end points
	sim.annot <- data.frame(sim.annot, start=sim.annot$end-sim.sizes, stringsAsFactors=FALSE)
	sim.annot <- data.frame(chr=sim.annot$chr, start=sim.annot$start, end=sim.annot$end, stringsAsFactors=FALSE)

	#add a col with an ID number
	sim.annot <- cbind(sim.annot, id=seq(nrow(sim.annot)))

	#write simulated sequences to FASTA

	#make bed file of simulated draw positions
	sim.bed <- data.frame(chr=sim.annot$chr,start=sim.annot$start,end=sim.annot$end)
	#sim.bed[sim.bed$start<0,]$start <- 0
	write.table(sim.bed,file="output/simseq.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
	system(paste("fastaFromBed -fi ",genome.path," -bed output/simseq.bed -fo output/simseq.fasta",sep=""))
	system("sed -i 's/:/-/' output/simseq.fasta")

	sim.seq <- readDNAStringSet("output/simseq.fasta")
	sim.seq
}


drawBackgroundSetFromRegions <- function(seq,regions,nSimSeqs=10000,windowSize=50)
{
	#input: sequence set to generate background for, regions dataframe with regions to draw background from
	#output: background set of nearest matching size sequences from supplied regions list
	#TODO

	sim.sizes <- sample(width(seq) + windowSize, nSimSeqs, replace=TRUE)

	sim.sizes.table <- table(sim.sizes)

	background.out <- data.frame()
	for (size in names(sim.sizes.table))
	{
		cat(size)
		sub1 <- subset(regions, subset=(end-start+windowSize)==as.integer(size))
		print(nrow(sub1))

		if(nrow(sub1) <10) sub1 <- subset(regions, subset=(end-start)%in%(as.integer(size)+seq(-100,100, by=50)))
		background.out <- rbind(background.out, sub1[sample.int(nrow(sub1), sim.sizes.table[size]),])
		cat("\tDONE!\n")
	}

	#need to match them back to their chrs so we can pull the seq from genome
	#pull FASTAS
	#read back and return DNAStringSet
	#make bed file of simulated draw positions
	sim.bed <- data.frame(chr=background.out$chr,start=background.out$start,end=background.out$end)
	#sim.bed[sim.bed$start<0,]$start <- 0
	write.table(sim.bed,file="output/simseq.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
	system(paste("fastaFromBed -fi ",genome.path," -bed output/simseq.bed -fo output/simseq.fasta",sep=""))
	system("sed -i 's/:/-/' output/simseq.fasta")

	#TODO: do we need to add some unique ID to any duplicates so FIMO isn't confused or skips them?

	sim.seq <- readDNAStringSet("output/simseq.fasta")
	sim.seq
}

# -----------------------------------------------------------------------------
# Use propensity score matching to create a background set
# Input: formula to match by
# Output: DNAStringSet of propensity matched background sequences
drawBackgroundSetPropensity <- function(target.seq, target.meta, pool.seq, pool.meta, formula)
{
	# setting binary value for group assignment
	target.meta$treat <- 1
	pool.meta$treat <- 0
	all.meta <- rbind(target.meta, pool.meta)

	# randomize sort order - order can bias when Match(..., replace=FALSE)
	index.random <- sample(seq(1:nrow(all.meta)),nrow(all.meta), replace=FALSE)
	all.meta.shuffle <- all.meta[index.random,]

	# run logistic model
	lrm.out <- lrm(formula, data=all.meta.shuffle)

	# obtain values
	lrm.out.fitted <- predict.lrm(lrm.out,type="fitted")

	# match
	rr <- Match(Y=NULL, Tr=all.meta.shuffle$treat, X=lrm.out.fitted, M=1, version="fast", replace=FALSE)
	#summary(rr)

	# make new sequence set
	matched.meta <- all.meta.shuffle[rr$index.control,]
	m <- match(as.character(matched.meta$name),names(pool.seq))
	seq.resamp <- pool.seq[m]
}
# -----------------------------------------------------------------------------

###############################################
## Background Sequence Generators by Shuffles
###############################################
randomBackgroundSet <- function(seq)
{
	#create background set with same base frequencies as input set

	doRandom <- function(x)
	{
		freq <- alphabetFrequency(x)[1:4]
		rand <- paste(sample(c("A", "C", "G", "T"), length(x), replace=TRUE, prob=freq), collapse="")
		#DNAString(rand)
		rand
	}
	new <- unlist(lapply(seq,FUN=doRandom))
	DNAStringSet(new)
}

shuffleBackgroundSet <- function(seq)
{
	#TODO - shuffle of seqs to make the background set
	#input: DNSStringSet of query seqs to shuffle
	#output: DNAStringSet of shuffled seqs
	
	#preserves nucleotide ratios exactly

	doShuf <- function(x)
	{
		x <- DNAString(as.character(x))
		shuf <- x[sample(length(x))]
		as.character(shuf)
	}
	new <- unlist(lapply(seq,FUN=doShuf))
	DNAStringSet(new)
}

shuffleDinucleotides <- function(seq)
{
	new <- unlist(lapply(seq,FUN=uShuffle,klet=2))
	DNAStringSet(new)
}

uShuffle <- function(string, klet)
{
	dyn.load("uShuffle/ushuffle.so")
	
	#input string
	s <- as.character(string)
	#output string
	t <- as.character(string)
	#length of string
	l <- as.integer(nchar(string))
	#k-let size to use
	k <- as.integer(klet)

	cdata <- .C("rshuffle", s=s, t=t, l=l, k=k)
	cdata$t
}

plotSequenceHistograms <- function(seq1,seq2)
{
	breaks <- seq(0,max(seq1.meta$size,seq2.meta$size)+100,by=100)
	breaks.gc <- seq(0,1,by=0.05)

	seq1.meta <- data.frame(name=names(seq1),size=width(seq1),gc=getGC(seq1))
	seq2.meta <- data.frame(name=names(seq2),size=width(seq2),gc=getGC(seq2))

	par(mfrow=c(2,3))
	hist(seq1.meta$size, breaks=breaks,freq=FALSE)
	hist(seq2.meta$size, breaks=breaks,freq=FALSE)
	qqplot(seq1.meta$size,seq2.meta$size)
	hist(seq1.meta$gc,breaks=breaks.gc,freq=FALSE)
	hist(seq2.meta$gc,breaks=breaks.gc,freq=FALSE)
	qqplot(seq1.meta$gc,seq2.meta$gc)
}

plotSequenceHistogramsResampled <- function(seq1,seq2,seq3)
{
	seq1.meta <- data.frame(name=names(seq1),size=width(seq1),gc=getGC(seq1))
	seq2.meta <- data.frame(name=names(seq2),size=width(seq2),gc=getGC(seq2))
	resamp <- data.frame(name=names(seq3),size=width(seq3),gc=getGC(seq3))

	#filter anything in seq2 bigger than the max in seq1
	seq2.meta <- seq2.meta[seq2.meta$size<=max(seq1.meta$size),]

	breaks <- seq(0,max(seq1.meta$size,seq2.meta$size,resamp$size)+100,by=100)
	breaks.gc <- seq(0,1,by=0.05)

	h1 <- hist(seq1.meta$size,breaks=breaks,plot=FALSE)
	h2 <- hist(seq2.meta$size,breaks=breaks,plot=FALSE)
	h2.1 <- hist(resamp$size,breaks=breaks,plot=FALSE)

	h3 <- hist(seq1.meta$gc,breaks=breaks.gc,plot=FALSE)
	h4 <- hist(seq2.meta$gc,breaks=breaks.gc,plot=FALSE)
	h4.1 <- hist(resamp$gc,breaks=breaks.gc,plot=FALSE)

	ggplot.hist.overlap <- function(h1,h2,name1="Hist 1",name2="Hist 2",xlab="",ylab="",main="")
	{
		plot.data <- data.frame(bin=h1$mids,freq=h1$count/sum(h1$count),seq=name1)
		plot.data <- rbind(plot.data,data.frame(bin=h2$mids,freq=h2$count/sum(h2$count),seq=name2))
		ggplot(plot.data, aes(x=bin,y=freq,fill=seq)) + geom_bar(data=subset(plot.data,seq == name1), stat="identity",alpha=0.5) + geom_bar(data=subset(plot.data,seq == name2), stat="identity",alpha=0.5) + scale_fill_manual("", values = c("#007FFF","#FF007F")) + labs(x=xlab, y=ylab, title = main) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.key.size = unit(0.8, "lines"), axis.line = element_line(colour = "grey50"))
	}

	p1 <- ggplot.hist.overlap(h1,h2,"Target","Background","Bin Mid","Prob.","Size: Original")
	p2 <- ggplot.hist.overlap(h1,h2.1,"Target","Background","Bin Mid","Prob.","Size: Resampled") + theme(legend.position = "none")
	p3 <- ggplot.hist.overlap(h3,h4,"Target","Background","Bin Mid","Prob.","GC: Original") + theme(legend.position = "none")
	p4 <- ggplot.hist.overlap(h3,h4.1,"Target","Background","Bin Mid","Prob.","GC: Resampled") + theme(legend.position="none")
	#p5 <- ggplot.hist.overlap(h1,h2.2,"Target","Background","Bin Mid","Prob.","Size: \"IS\" Resamp") + theme(legend.position = "none")
	#p6 <- ggplot.hist.overlap(h3,h4.2,"Target","Background","Bin Mid","Prob.","GC: \"IS\" Resamp") + theme(legend.position="none")

	#g1 <- arrangeGrob(p1 + theme(legend.position = "none"),p2,p5,p3,p4,p6,ncol=3)
	g1 <- arrangeGrob(p1 + theme(legend.position = "none"),p2,p3,p4,ncol=2)

	#plot overlapping hists

	tmp <- ggplot_gtable(ggplot_build(p1))
	leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
	legend <- tmp$grobs[[leg]]


	grid.arrange(g1, legend, widths=c(6/7,1/7), nrow=1)
}

plotSequenceQQ4 <- function(seq1,seq2,seq3,seq4,labels)
{
	seq1.meta <- data.frame(name=names(seq1),size=width(seq1),gc=getGC(seq1))
	seq2.meta <- data.frame(name=names(seq2),size=width(seq2),gc=getGC(seq2))
	resamp <- data.frame(name=names(seq3),size=width(seq3),gc=getGC(seq3))
	resamp2 <- data.frame(name=names(seq4),size=width(seq4),gc=getGC(seq4))

	#filter anything in seq2 bigger than the max in seq1
	seq2.meta <- seq2.meta[seq2.meta$size<=max(seq1.meta$size),]

	breaks <- seq(0,max(seq1.meta$size)+100,by=100)
	breaks.gc <- seq(0,1,by=0.05)

	h1 <- hist(seq1.meta$size,breaks=breaks,plot=FALSE)
	h2 <- hist(seq2.meta$size,breaks=breaks,plot=FALSE)
	h2.1 <- hist(resamp$size,breaks=breaks,plot=FALSE)

	h3 <- hist(seq1.meta$gc,breaks=breaks.gc,plot=FALSE)
	h4 <- hist(seq2.meta$gc,breaks=breaks.gc,plot=FALSE)
	h4.1 <- hist(resamp$gc,breaks=breaks.gc,plot=FALSE)

	h2.2 <- hist(resamp2$size,breaks=breaks,plot=FALSE)
	h4.2 <- hist(resamp2$gc,breaks=breaks.gc,plot=FALSE)

	####QQ Plots + Distance Metric Data Table
	#equations to find one difference measure between two histograms
	distanceEuclidean <- function(h1,h2)
	{
		d1 <- h1$counts/sum(h1$counts)
		d2 <- h2$counts/sum(h2$counts)
		sqrt(sum(abs(d1-d2)^2))
	}

	distanceBhattacharyya <- function(h1,h2)
	{
		d1 <- h1$counts/sum(h1$counts)
		d2 <- h2$counts/sum(h2$counts)
		#-log(sum(sqrt(d1*d2)))
		sum(sqrt(d1*d2)) #value of 1 should indicate a perfect match (finding square of angle between two vecotrs of positions to cosine of angle between them is 1 if identical) - should not be effected by bin widths
	}
	dist1 <- round(c(distanceEuclidean(h1,h2), distanceEuclidean(h1,h2.1), distanceEuclidean(h1,h2.2)),digits=3)
	dist2 <- round(c(distanceEuclidean(h3,h4), distanceEuclidean(h3,h4.1), distanceEuclidean(h3,h4.2)),digits=3)
	dist3 <- round(c(distanceBhattacharyya(h1,h2), distanceBhattacharyya(h1,h2.1), distanceBhattacharyya(h1,h2.2)),digits=3)
	dist4 <- round(c(distanceBhattacharyya(h3,h4), distanceBhattacharyya(h3,h4.1) , distanceBhattacharyya(h3,h4.2)),digits=3)

	distance.stats <- data.frame(size.dist.euclid=dist1,gc.dist.euclid=dist2,size.dist.bhatt=dist3,gc.dist.bhatt=dist4,row.names=labels)


	ggplot.qqplot.geom <- function(data1,data2)
	{
		d <- as.data.frame(qqplot(data1, data2, plot.it=FALSE))
		ggplot(d) + geom_point(aes(x=x, y=y),stat = "identity", position = "identity")
	}

	d <- data.frame(as.data.frame(qqplot(seq1.meta$size, seq2.meta$size, plot.it=FALSE)),Sequences=labels[1])
	d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$size, resamp$size, plot.it=FALSE)),Sequences=labels[2]))
	d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$size, resamp2$size, plot.it=FALSE)),Sequences=labels[3]))
	p1 <- ggplot(d) + geom_point(aes(x=x, y=y,color=Sequences),stat = "identity", position = "identity", ) + ggplot.clean() + labs(x="Target", y="Background",title="Size") + geom_abline(slope = 1, intercept=0)
	d <- data.frame(as.data.frame(qqplot(seq1.meta$gc, seq2.meta$gc, plot.it=FALSE)),Sequences=labels[1])
	d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$gc, resamp$gc, plot.it=FALSE)),Sequences=labels[2]))
	d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$gc, resamp2$gc, plot.it=FALSE)),Sequences=labels[3]))
	p2 <- ggplot(d) + geom_point(aes(x=x, y=y,color=Sequences),stat = "identity", position = "identity", ) + ggplot.clean() + labs(x="Target", y="Background",title="GC") + geom_abline(slope = 1, intercept=0)
	g1 <- arrangeGrob(p1,p2,ncol=2)

	distance.stats.grob <- arrangeGrob(tableGrob(distance.stats))
	grid.arrange(g1, distance.stats.grob, ncol=1,heights=c(4/5,1/5))
}

plotSequenceQQ5 <- function(seq1,seq2,seq3,seq4,seq5,labels)
{
	seq1.meta <- data.frame(name=names(seq1),size=width(seq1),gc=getGC(seq1))
	seq2.meta <- data.frame(name=names(seq2),size=width(seq2),gc=getGC(seq2))
	resamp <- data.frame(name=names(seq3),size=width(seq3),gc=getGC(seq3))
	resamp2 <- data.frame(name=names(seq4),size=width(seq4),gc=getGC(seq4))
	resamp3 <- data.frame(name=names(seq5),size=width(seq5),gc=getGC(seq5))

	#filter anything in seq2 bigger than the max in seq1
	seq2.meta <- seq2.meta[seq2.meta$size<=max(seq1.meta$size),]

	breaks <- seq(0,max(seq1.meta$size,resamp3$size)+100,by=100)
	breaks.gc <- seq(0,1,by=0.05)

	h1 <- hist(seq1.meta$size,breaks=breaks,plot=FALSE)
	h2 <- hist(seq2.meta$size,breaks=breaks,plot=FALSE)
	h2.1 <- hist(resamp$size,breaks=breaks,plot=FALSE)

	h3 <- hist(seq1.meta$gc,breaks=breaks.gc,plot=FALSE)
	h4 <- hist(seq2.meta$gc,breaks=breaks.gc,plot=FALSE)
	h4.1 <- hist(resamp$gc,breaks=breaks.gc,plot=FALSE)

	h2.2 <- hist(resamp2$size,breaks=breaks,plot=FALSE)
	h4.2 <- hist(resamp2$gc,breaks=breaks.gc,plot=FALSE)

	h2.3 <- hist(resamp3$size,breaks=breaks,plot=FALSE)
	h4.3 <- hist(resamp3$gc,breaks=breaks.gc,plot=FALSE)

	####QQ Plots + Distance Metric Data Table
	#equations to find one difference measure between two histograms
	distanceEuclidean <- function(h1,h2)
	{
		d1 <- h1$counts/sum(h1$counts)
		d2 <- h2$counts/sum(h2$counts)
		sqrt(sum(abs(d1-d2)^2))
	}

	distanceBhattacharyya <- function(h1,h2)
	{
		d1 <- h1$counts/sum(h1$counts)
		d2 <- h2$counts/sum(h2$counts)
		#-log(sum(sqrt(d1*d2)))
		sum(sqrt(d1*d2)) #value of 1 should indicate a perfect match (finding square of angle between two vecotrs of positions to cosine of angle between them is 1 if identical) - should not be effected by bin widths
	}
	dist1 <- round(c(distanceEuclidean(h1,h2), distanceEuclidean(h1,h2.1), distanceEuclidean(h1,h2.2), distanceEuclidean(h1,h2.3)),digits=3)
	dist2 <- round(c(distanceEuclidean(h3,h4), distanceEuclidean(h3,h4.1), distanceEuclidean(h3,h4.2), distanceEuclidean(h3,h4.3)),digits=3)
	dist3 <- round(c(distanceBhattacharyya(h1,h2), distanceBhattacharyya(h1,h2.1), distanceBhattacharyya(h1,h2.2), distanceBhattacharyya(h1,h2.3)),digits=3)
	dist4 <- round(c(distanceBhattacharyya(h3,h4), distanceBhattacharyya(h3,h4.1) , distanceBhattacharyya(h3,h4.2), distanceBhattacharyya(h3,h4.3)),digits=3)

	distance.stats <- data.frame(size.dist.euclid=dist1,gc.dist.euclid=dist2,size.dist.bhatt=dist3,gc.dist.bhatt=dist4,row.names=labels)


	ggplot.qqplot.geom <- function(data1,data2)
	{
		d <- as.data.frame(qqplot(data1, data2, plot.it=FALSE))
		ggplot(d) + geom_point(aes(x=x, y=y),stat = "identity", position = "identity")
	}

	d <- data.frame(as.data.frame(qqplot(seq1.meta$size, seq2.meta$size, plot.it=FALSE)),Sequences=labels[1])
	d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$size, resamp$size, plot.it=FALSE)),Sequences=labels[2]))
	d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$size, resamp2$size, plot.it=FALSE)),Sequences=labels[3]))
	d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$size, resamp3$size, plot.it=FALSE)),Sequences=labels[4]))
	p1 <- ggplot(d) + geom_point(aes(x=x, y=y,color=Sequences),stat = "identity", position = "identity", ) + ggplot.clean() + labs(x="Target", y="Background",title="Size") + geom_abline(slope = 1, intercept=0)
	d <- data.frame(as.data.frame(qqplot(seq1.meta$gc, seq2.meta$gc, plot.it=FALSE)),Sequences=labels[1])
	d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$gc, resamp$gc, plot.it=FALSE)),Sequences=labels[2]))
	d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$gc, resamp2$gc, plot.it=FALSE)),Sequences=labels[3]))
	d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$gc, resamp3$gc, plot.it=FALSE)),Sequences=labels[4]))
	p2 <- ggplot(d) + geom_point(aes(x=x, y=y,color=Sequences),stat = "identity", position = "identity", ) + ggplot.clean() + labs(x="Target", y="Background",title="GC") + geom_abline(slope = 1, intercept=0)
	g1 <- arrangeGrob(p1,p2,ncol=2)

	distance.stats.grob <- arrangeGrob(tableGrob(distance.stats))
	grid.arrange(g1, distance.stats.grob, ncol=1,heights=c(4/5,1/5))
}

plotSequenceQQR5 <- function(seq1,seq2,seq3,seq4,seq5,labels)
{
	seq1.meta <- data.frame(name=names(seq1),size=width(seq1),gc=getGC(seq1))
	seq2.meta <- data.frame(name=names(seq2),size=width(seq2),gc=getGC(seq2))
	resamp <- data.frame(name=names(seq3),size=width(seq3),gc=getGC(seq3))
	resamp2 <- data.frame(name=names(seq4),size=width(seq4),gc=getGC(seq4))
	resamp3 <- data.frame(name=names(seq5),size=width(seq5),gc=getGC(seq5))

	#filter anything in seq2 bigger than the max in seq1
	seq2.meta <- seq2.meta[seq2.meta$size<=max(seq1.meta$size),]

	breaks <- seq(0,max(seq1.meta$size,resamp3$size)+100,by=100)
	breaks.gc <- seq(0,1,by=0.05)

	h1 <- hist(seq1.meta$size,breaks=breaks,plot=FALSE)
	h2 <- hist(seq2.meta$size,breaks=breaks,plot=FALSE)
	h2.1 <- hist(resamp$size,breaks=breaks,plot=FALSE)

	h3 <- hist(seq1.meta$gc,breaks=breaks.gc,plot=FALSE)
	h4 <- hist(seq2.meta$gc,breaks=breaks.gc,plot=FALSE)
	h4.1 <- hist(resamp$gc,breaks=breaks.gc,plot=FALSE)

	h2.2 <- hist(resamp2$size,breaks=breaks,plot=FALSE)
	h4.2 <- hist(resamp2$gc,breaks=breaks.gc,plot=FALSE)

	h2.3 <- hist(resamp3$size,breaks=breaks,plot=FALSE)
	h4.3 <- hist(resamp3$gc,breaks=breaks.gc,plot=FALSE)

	####QQ Plots + Distance Metric Data Table
	#get R^2 values on the QQ plots
	qqcor <- function(c1,c2)
	{
		qqd <- qqplot(c1, c2, plot.it=FALSE)
		cor(qqd$x,qqd$y)
	}

#	qqcorlm <- function(c1,c2)
#	{
#		lm1 <- lm(0 + qqd$y ~ 0 + qqd$x)
#		s <- summary(lm1)
#		s$r.squared
#png(filename="tmp.png",res=120)
#plot(qqd$x,qqd$y)
#abline(0, coef(lm1))
#dev.off()

#	}

	dist1 <- round(c(qqcor(seq1.meta$size,seq2.meta$size), qqcor(seq1.meta$size,resamp$size), qqcor(seq1.meta$size,resamp2$size), qqcor(seq1.meta$size,resamp3$size)),digits=3)
	dist2 <- round(c(qqcor(seq1.meta$gc,seq2.meta$gc), qqcor(seq1.meta$gc,resamp$gc), qqcor(seq1.meta$gc,resamp2$gc), qqcor(seq1.meta$gc,resamp3$gc)),digits=3)
	distance.stats <- data.frame(size.r=dist1,gc.r=dist2,row.names=labels)


	ggplot.qqplot.geom <- function(data1,data2)
	{
		d <- as.data.frame(qqplot(data1, data2, plot.it=FALSE))
		ggplot(d) + geom_point(aes(x=x, y=y),stat = "identity", position = "identity")
	}

	d <- data.frame(as.data.frame(qqplot(seq1.meta$size, seq2.meta$size, plot.it=FALSE)),Sequences=labels[1])
	d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$size, resamp$size, plot.it=FALSE)),Sequences=labels[2]))
	d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$size, resamp2$size, plot.it=FALSE)),Sequences=labels[3]))
	d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$size, resamp3$size, plot.it=FALSE)),Sequences=labels[4]))
	p1 <- ggplot(d) + geom_point(aes(x=x, y=y,color=Sequences),stat = "identity", position = "identity", ) + ggplot.clean() + labs(x="Target", y="Background",title="Size") + geom_abline(slope = 1, intercept=0)
	d <- data.frame(as.data.frame(qqplot(seq1.meta$gc, seq2.meta$gc, plot.it=FALSE)),Sequences=labels[1])
	d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$gc, resamp$gc, plot.it=FALSE)),Sequences=labels[2]))
	d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$gc, resamp2$gc, plot.it=FALSE)),Sequences=labels[3]))
	d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$gc, resamp3$gc, plot.it=FALSE)),Sequences=labels[4]))
	p2 <- ggplot(d) + geom_point(aes(x=x, y=y,color=Sequences),stat = "identity", position = "identity", ) + ggplot.clean() + labs(x="Target", y="Background",title="GC") + geom_abline(slope = 1, intercept=0)
	g1 <- arrangeGrob(p1,p2,ncol=2)

	distance.stats.grob <- arrangeGrob(tableGrob(distance.stats))
	grid.arrange(g1, distance.stats.grob, ncol=1,heights=c(4/5,1/5))
}


###############################################
## Enrichment Testing
###############################################
calcEnrichmentBinom <- function(seq1.counts,seq1.nSeqs,seq2.counts,seq2.nSeqs)
{
	#input: count matrices (sequences vs motifs) for two runs of FIMO, number of seqs for each set
	#output: enrichment p-values for each motif

	#calculate frequencies for each motif in each set of sequences
	#consider each occurance of motif only once: convert counts >1 to be 1 so we can easily sum each row to get counts
	makeBinary <- function(value)
	{
		if(value>1)
		{
			value <- 1
		}
		value
	}
	counts1.bin <- apply(seq1.counts,MARGIN=c(1,2),FUN=makeBinary)	
	counts2.bin <- apply(seq2.counts,MARGIN=c(1,2),FUN=makeBinary)

	pvalues <- foreach(i=1:nrow(counts1.bin),.combine=rbind) %dopar%
	{
		#Using binom.test:
		#x = num successes = total count of sequences with one or more instance of motif found by FIMO
		myX <- sum(counts1.bin[i,])
		#n = num trials = number of sequences
		myN <- seq1.nSeqs
		#p = prob of success = count of seqs with >=1 instance in background divided by total seqs in background	
		if(rownames(counts1.bin)[1] %in% rownames(counts2.bin))
		{
			#if motif shows up in background set
			#extract row with matching motif name and calculate frequency
			myP <- sum(counts2.bin[rownames(counts2.bin) %in% rownames(counts1.bin)[i],])/seq2.nSeqs
		} else
		{
			#if motif did not show up in background set
			#print("Instance w/o same motif in both")
			myP <- 0
		}

		pValue <- binom.test(x=myX, n=myN, p=myP, alternative="greater")$p.value

		data.frame(motif=rownames(counts1.bin)[i],pvalue=pValue,percent_seqs=((myX/myN)*100))
	}

	pvalues
}

calcMotifCounts <- function(fimo.out, q.cutoff)
{
	#input: dataframe of fimo's text output, q value cutoff
	#output: matrix of frequencies with a row for every motif and column for every sequence

	fimo.out.filtered <- fimo.out[fimo.out$q.value<q.cutoff,]

	#count of occurances for every combination of sequence and motif
	pairs <- ddply(fimo.out.filtered,.(sequence.name,X.pattern.name),nrow,.drop=TRUE,.parallel=TRUE)

	#reshape into a pairwise table
	table <- cast(pairs,X.pattern.name~sequence.name,value="V1")

	#convert to matrix
	mat <- as.matrix(table)

	#replace the NAs with 0s
	replaceNA <- function(value)
	{
		if(is.na(value))
		{
			value <- 0
		}
		value
	}

	mat <- apply(mat,MARGIN=c(1,2),FUN=replaceNA)

	mat
}

