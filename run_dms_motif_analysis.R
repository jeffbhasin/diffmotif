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
library(ggplot2)
library(gridExtra)
library(doMC)
registerDoMC(nCores)

###############################################
## Analysis Functions
###############################################

###############################################
## Main Body
###############################################

#Initial run of FIMO on sequence set
run.name <- "CIMPHyperMe"
seq.path <- "../RunMEME/ColonSeqGR2012/CIMPHyperMe.fasta"
motif.path <- "../TRANSFAC/ConvertToMEME/TRAJAS.Human.named.meme"
runFIMO(run.name,seq.path,motif.path)
fimo.out <- readFIMO(paste("output/fimo_out_",run.name,"/fimo.txt",sep=""))
q.cutoff <- 0.05
fimo.out.counts <- calcMotifCounts(fimo.out,q.cutoff)
seq <- readDNAStringSet(seq.path)
nSeqs <- length(seq)

###############################################
#Find differential motifs from CIMP HyperMe versus background drawn from whole genome

sim.seq <- drawBackgroundSet(seq,10000)
sim.nSeqs <- length(sim.seq)
runFIMO("sim","output/simseq.fasta",motif.path)
fimo.out.sim <- readFIMO("output/fimo_out_sim/fimo.txt")
fimo.out.sim.counts <- calcMotifCounts(fimo.out.sim,q.cutoff)
results <- calcEnrichmentBinom(fimo.out.counts,nSeqs,fimo.out.sim.counts,sim.nSeqs)
results.sort <- results[order(results$pvalue, decreasing=FALSE),]
head(results.sort)
tail(results.sort)
write.csv(results.sort,file="output/diffmotif.bkg.genomic.csv",row.names=FALSE)

###############################################
#Find differential motifs from CIMP HyperMe versus background drawn from non-sig methylated regions (replicate of GR2012 methods)
#TODO: replace with actual NS DMS sites once data obtained
regions.all <- data.frame(chr=c("chr1","chr2"),start=c(10,300),end=c(50,500),p=c(0.6,1.4))
#filter to only give NS regions with p>0.5
regions <- regions.all[regions.all$p>0.5,]
sim.seq <- drawBackgroundSetFromRegions(seq,regions,10000)

sim.nSeqs <- length(sim.seq)
runFIMO("sim","output/simseq.fasta",motif.path)
fimo.out.sim <- readFIMO("output/fimo_out_sim/fimo.txt")
fimo.out.sim.counts <- calcMotifCounts(fimo.out.sim,q.cutoff)
results <- calcEnrichmentBinom(fimo.out.counts,nSeqs,fimo.out.sim.counts,sim.nSeqs)
results.sort <- results[order(results$pvalue, decreasing=FALSE),]
head(results.sort)
tail(results.sort)

###############################################
#Find differential motifs from CIMP HyperMe versus background of dinucleotide shuffled source
sim.seq <- shuffleDinucleotides(seq)

sim.nSeqs <- length(sim.seq)

writeXStringSet(sim.seq, "output/simseq.fasta", append=FALSE, format="fasta")

runFIMO("sim","output/simseq.fasta",motif.path)
fimo.out.sim <- readFIMO("output/fimo_out_sim/fimo.txt")
fimo.out.sim.counts <- calcMotifCounts(fimo.out.sim,q.cutoff)
results <- calcEnrichmentBinom(fimo.out.counts,nSeqs,fimo.out.sim.counts,sim.nSeqs)
results.sort <- results[order(results$pvalue, decreasing=FALSE),]
head(results.sort)
tail(results.sort)
write.csv(results.sort,file="output/diffmotif.bkg.dinucleotide.csv",row.names=FALSE)

###############################################
#Find differential motifs from CIMP HyperMe versus background of shuffled source
sim.seq <- shuffleBackgroundSet(seq)

sim.nSeqs <- length(sim.seq)

writeXStringSet(sim.seq, "output/simseq.fasta", append=FALSE, format="fasta")

runFIMO("sim","output/simseq.fasta",motif.path)
fimo.out.sim <- readFIMO("output/fimo_out_sim/fimo.txt")
fimo.out.sim.counts <- calcMotifCounts(fimo.out.sim,q.cutoff)
results <- calcEnrichmentBinom(fimo.out.counts,nSeqs,fimo.out.sim.counts,sim.nSeqs)
results.sort <- results[order(results$pvalue, decreasing=FALSE),]
head(results.sort)
tail(results.sort)
write.csv(results.sort,file="output/diffmotif.bkg.shuffle.csv",row.names=FALSE)

###############################################
#Find differential motifs from CIMP HyperMe versus background of random draw source
sim.seq <- randomBackgroundSet(seq)

sim.nSeqs <- length(sim.seq)

writeXStringSet(sim.seq, "output/simseq.fasta", append=FALSE, format="fasta")

runFIMO("sim","output/simseq.fasta",motif.path)
fimo.out.sim <- readFIMO("output/fimo_out_sim/fimo.txt")
fimo.out.sim.counts <- calcMotifCounts(fimo.out.sim,q.cutoff)
results <- calcEnrichmentBinom(fimo.out.counts,nSeqs,fimo.out.sim.counts,sim.nSeqs)
results.sort <- results[order(results$pvalue, decreasing=FALSE),]
head(results.sort)
tail(results.sort)
write.csv(results.sort,file="output/diffmotif.bkg.random.csv",row.names=FALSE)



###############################################
#Find differential motifs from CIMP HyperMe versus resampled background set - made up "nearest" distance method

#want subset of seq2 that best matches length and gc distribution of seq1

seq1 <- seq

seq2.path <- "../RunMEME/ColonSeqGR2012/AllHyperMe.fasta"
seq2 <- readDNAStringSet(seq2.path)

getGC <- function(seq)
{
	g <- alphabetFrequency(seq)[,3]
	c <- alphabetFrequency(seq)[,2]
	a <- alphabetFrequency(seq)[,1]
	t <- alphabetFrequency(seq)[,4]

	gc <- (g+c) / (a+t+g+c)
	gc
}

seq1.meta <- data.frame(name=names(seq1),size=width(seq1),gc=getGC(seq1))
seq2.meta <- data.frame(name=names(seq2),size=width(seq2),gc=getGC(seq2))

#filter anything in seq2 bigger than the max of seq1
seq2.meta <- seq2.meta[seq2.meta$size<=max(seq1.meta$size),]

breaks <- seq(0,max(seq1.meta$size)+100,by=100)
breaks.gc <- seq(0,1,by=0.05)
#filter out 

#freq=TRUE
#make plots of the starting distributions
png(filename="output/dists.png",width=1000,height=800,res=120)
par(mfrow=c(2,3))
hist(seq1.meta$size,breaks=breaks,freq=FALSE)
hist(seq2.meta$size,breaks=breaks,freq=FALSE)
qqplot(seq1.meta$size,seq2.meta$size)
hist(seq1.meta$gc,breaks=breaks.gc,freq=FALSE)
hist(seq2.meta$gc,breaks=breaks.gc,freq=FALSE)
qqplot(seq1.meta$gc,seq2.meta$gc)
dev.off()

#method to find and draw best matches, using each only once to create new distribution (will have same number of seqs as distribution being matched)

pool <- seq2.meta
#pool$used <- 0

resamp <- foreach(seq=iter(seq1.meta, by='row'),.combine=rbind) %do%
{
	#find closest other sizes
	nearests <- which(abs(pool$size-seq$size)==min(abs(pool$size-seq$size)))

	#now look at their GC and take the one with the closest GC
	gc <- pool[nearests,]$gc
	mine <- pool[nearests,]
	mine2 <- mine[which(abs(gc-seq$gc)==min(abs(gc-seq$gc))),]

	#if there are more than one pick one at random	
	use <- sample(mine2$name,1)

	#mark as used - just going to remove it to make it easier
	mine.row <- pool[as.character(pool$name)==as.character(use),]
	pool <- pool[as.character(pool$name)!=as.character(use),]

	#return as member of new distribution
	mine.row
}

#subset original DNAStringSet to just get the seqs we want
seq.resamp <- seq2[names(seq2) %in% resamp$name]

#replot distributions to see how well we did
png(filename="output/dists.resamp.png",width=1000,height=800,res=120)
par(mfrow=c(2,3))
hist(seq1.meta$size,breaks=breaks,freq=FALSE)
hist(resamp$size,breaks=breaks,freq=FALSE)
qqplot(seq1.meta$size,resamp$size)
hist(seq1.meta$gc,breaks=breaks.gc,freq=FALSE)
hist(resamp$gc,breaks=breaks.gc,freq=FALSE)
qqplot(seq1.meta$gc,resamp$gc)
dev.off()

#pdf difference function tests
h1 <- hist(seq1.meta$size,breaks=breaks,plot=FALSE)
h2 <- hist(seq2.meta$size,breaks=breaks,plot=FALSE)
h2.1 <- hist(resamp$size,breaks=breaks,plot=FALSE)

h3 <- hist(seq1.meta$gc,breaks=breaks.gc,plot=FALSE)
h4 <- hist(seq2.meta$gc,breaks=breaks.gc,plot=FALSE)
h4.1 <- hist(resamp$gc,breaks=breaks.gc,plot=FALSE)




#test run of FIMO
sim.seq <- seq.resamp

sim.nSeqs <- length(sim.seq)

writeXStringSet(sim.seq, "output/simseq.fasta", append=FALSE, format="fasta")

runFIMO("sim","output/simseq.fasta",motif.path)
fimo.out.sim <- readFIMO("output/fimo_out_sim/fimo.txt")
fimo.out.sim.counts <- calcMotifCounts(fimo.out.sim,q.cutoff)
results <- calcEnrichmentBinom(fimo.out.counts,nSeqs,fimo.out.sim.counts,sim.nSeqs)
results.sort <- results[order(results$pvalue, decreasing=FALSE),]
head(results.sort)
tail(results.sort)
write.csv(results.sort,file="output/diffmotif.bkg.resamp-allhyperme.csv",row.names=FALSE)

###############################################
#Find differential motifs from CIMP HyperMe versus resampled background set - importance sampling

##Try to do real importance sampling - Boot R Package
#not working
#dum <- function(a)
#{
#	a
#}

#library(boot)

#b <- boot(seq2.meta$gc,data.frame,R=1)

#b <- boot(seq2.meta$gc,sum,R=1)

#wt <- imp.weights(b)

#boot(seq2.meta$gc,sum,R=1,weights=wt)


###############################################
#Find differential motifs from CIMP HyperMe versus resampled background set - attempt at importance sampling using R functions

##weight each sample by frequency of corresponding bin in other histogram
##take a random sample permutation using those frequencies as the weights

#try for GC
freqs <- h3$counts/sum(h3$counts)
freqs <- data.frame(bin=1:length(freqs),freq=freqs)

#compute bin membership for each data point in data we're drawing from
binnums <- cut(seq2.meta$gc,breaks.gc,labels=FALSE)

weights <- freqs[binnums,]$freq

drawn.index <- sample(1:nrow(seq2.meta),1000000,replace=TRUE,prob=weights)

resamp2 <- seq2.meta[drawn.index,]

#replot distributions to see how well we did
png(filename="output/dists.myimpsamp.png",width=1000,height=800,res=120)
par(mfrow=c(2,3))
hist(seq1.meta$size,breaks=breaks,freq=FALSE)
hist(resamp2$size,breaks=breaks,freq=FALSE)
qqplot(seq1.meta$size,resamp2$size)
hist(seq1.meta$gc,breaks=breaks.gc,freq=FALSE)
hist(resamp2$gc,breaks=breaks.gc,freq=FALSE)
qqplot(seq1.meta$gc,resamp2$gc)
dev.off()

#try for size

freqs <- h1$counts/sum(h1$counts)
freqs <- data.frame(bin=1:length(freqs),freq=freqs)

#compute bin membership for each data point in data we're drawing from
binnums <- cut(seq2.meta$size,breaks,labels=FALSE)

weights <- freqs[binnums,]$freq

drawn.index <- sample(1:nrow(seq2.meta),1000000,replace=TRUE,prob=weights)

resamp2 <- seq2.meta[drawn.index,]

#replot distributions to see how well we did
png(filename="output/dists.myimpsamp.sizeonly.png",width=1000,height=800,res=120)
par(mfrow=c(2,3))
hist(seq1.meta$size,breaks=breaks,freq=FALSE)
hist(resamp2$size,breaks=breaks,freq=FALSE)
qqplot(seq1.meta$size,resamp2$size)
hist(seq1.meta$gc,breaks=breaks.gc,freq=FALSE)
hist(resamp2$gc,breaks=breaks.gc,freq=FALSE)
qqplot(seq1.meta$gc,resamp2$gc)
dev.off()

h2.2 <- hist(resamp2$size,breaks=breaks,freq=FALSE)
h4.2 <- hist(resamp2$gc,breaks=breaks.gc,freq=FALSE)

#now using ggplot for histogram plots
#plot.hist.overlap <- function(hist.out1,hist.out2,...)
#{
#	plot( hist.out1, col=rgb(0,0,1,1/4),...)  # first histogram
#	plot( hist.out2, col=rgb(1,0,0,1/4), add=T,...)  # second
#}



#plot.hist.overlap <- function(hist.out1,hist.out2,...)
#{
#	barplot( hist.out1$counts,names.arg=hist.out1$mids, col=rgb(0,0,1,1/4),...)  # first histogram
#	barplot( hist.out1$counts,names.arg=hist.out1$mids, col=rgb(1,0,0,1/4), add=T,...)  # second
#}

#png(filename="output/test.png",width=800,height=800,res=120)
#par(mfrow=c(2,3))
#plot.hist.overlap(h2,h1,main="Size Overlap")
#plot.hist.overlap(h4,h3,main="GC Overlap2")
#dev.off()


#plot.data <- data.frame(name=seq1.meta$name,size=seq1.meta$size,gc=seq1.meta$gc,seq="seq1")
#plot.data <- rbind(plot.data,data.frame(name=seq2.meta$name,size=seq2.meta$size,gc=seq2.meta$gc,seq="seq2"))
#ggplot(plot.data, aes(x=gc, fill=seq)) + geom_histogram(binwidth=0.05, alpha=.5, position="identity")

#h1 <- hist(rnorm(500,4))
#h2 <- hist(rnorm(500,6))



##summary plots/stats of all resampling methods


#plot QQ plots - now doing a combined plot of all
#ggplot.qqplot <- function(data1,data2)
#{
#	d <- as.data.frame(qqplot(data1, data2, plot.it=FALSE))
#	ggplot(d) + geom_point(aes(x=x, y=y),stat = "identity", position = "identity")
#}

#ggplot.clean <- function()
#{
#	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.key.size = unit(0.8, "lines"), axis.line = element_line(colour = "grey50"))
#}

#p1 <- ggplot.qqplot(seq1.meta$size,seq2.meta$size) + ggplot.clean() + labs(x="CIMP Hyper Me", y="All Tumors Hyper Me",title="Size")
#p2 <- ggplot.qqplot(seq1.meta$size,resamp$size) + ggplot.clean() + labs(x="CIMP Hyper Me", y="All Tumors Hyper Me",title="Size: Nearest Resamp")
#p3 <- ggplot.qqplot(seq1.meta$size,resamp2$size) + ggplot.clean() + labs(x="CIMP Hyper Me", y="All Tumors Hyper Me",title="Size: IS Resamp")

#p4 <- ggplot.qqplot(seq1.meta$gc,seq2.meta$gc) + ggplot.clean() + labs(x="CIMP Hyper Me", y="All Tumors Hyper Me",title="GC")
#p5 <- ggplot.qqplot(seq1.meta$gc,resamp$gc) + ggplot.clean() + labs(x="CIMP Hyper Me", y="All Tumors Hyper Me",title="GC: Nearest Resamp")
#p6 <- ggplot.qqplot(seq1.meta$gc,resamp2$gc) + ggplot.clean() + labs(x="CIMP Hyper Me", y="All Tumors Hyper Me",title="GC: IS Resamp")

#g1 <- arrangeGrob(p1,p2,p3,p4,p5,p6,ncol=3)


####Overlap Plots of All Resampling Methods:
ggplot.hist.overlap <- function(h1,h2,name1="Hist 1",name2="Hist 2",xlab="",ylab="",main="")
{
	plot.data <- data.frame(bin=h1$mids,freq=h1$count/sum(h1$count),seq=name1)
	plot.data <- rbind(plot.data,data.frame(bin=h2$mids,freq=h2$count/sum(h2$count),seq=name2))
	ggplot(plot.data, aes(x=bin,y=freq,fill=seq)) + geom_bar(data=subset(plot.data,seq == name1), stat="identity",alpha=0.5) + geom_bar(data=subset(plot.data,seq == name2), stat="identity",alpha=0.5) + scale_fill_manual("", values = c("#007FFF","#FF007F")) + labs(x=xlab, y=ylab, title = main) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.key.size = unit(0.8, "lines"), axis.line = element_line(colour = "grey50"))
}

p1 <- ggplot.hist.overlap(h1,h2,"CIMP HyperMe","All HyperMe","Bin Mid","Prob.","Size: Original")
p2 <- ggplot.hist.overlap(h1,h2.1,"CIMP HyperMe","All HyperMe","Bin Mid","Prob.","Size: \"Nearest\" Resamp") + theme(legend.position = "none")
p3 <- ggplot.hist.overlap(h3,h4,"CIMP HyperMe","All HyperMe","Bin Mid","Prob.","GC: Original") + theme(legend.position = "none")
p4 <- ggplot.hist.overlap(h3,h4.1,"CIMP HyperMe","All HyperMe","Bin Mid","Prob.","GC: \"Nearest\" Resamp") + theme(legend.position="none")
p5 <- ggplot.hist.overlap(h1,h2.2,"CIMP HyperMe","All HyperMe","Bin Mid","Prob.","Size: \"IS\" Resamp") + theme(legend.position = "none")
p6 <- ggplot.hist.overlap(h3,h4.2,"CIMP HyperMe","All HyperMe","Bin Mid","Prob.","GC: \"IS\" Resamp") + theme(legend.position="none")
g1 <- arrangeGrob(p1 + theme(legend.position = "none"),p2,p5,p3,p4,p6,ncol=3)

#plot overlapping hists

tmp <- ggplot_gtable(ggplot_build(p1))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]


png(filename="output/resamp.hists.png",width=1400,height=800,res=120)
grid.arrange(g1, legend, widths=c(6/7,1/7), nrow=1)
dev.off()

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

distance.stats <- data.frame(size.dist.euclid=dist1,gc.dist.euclid=dist2,size.dist.bhatt=dist3,gc.dist.bhatt=dist4,row.names=c("Original","Nearest","IS"))


ggplot.qqplot.geom <- function(data1,data2)
{
	d <- as.data.frame(qqplot(data1, data2, plot.it=FALSE))
	ggplot(d) + geom_point(aes(x=x, y=y),stat = "identity", position = "identity")
}

d <- data.frame(as.data.frame(qqplot(seq1.meta$size, seq2.meta$size, plot.it=FALSE)),Sequences="Originals")
d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$size, resamp$size, plot.it=FALSE)),Sequences="Nearest"))
d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$size, resamp2$size, plot.it=FALSE)),Sequences="IS"))
p1 <- ggplot(d) + geom_point(aes(x=x, y=y,color=Sequences),stat = "identity", position = "identity", ) + ggplot.clean() + labs(x="CIMP Hyper Me", y="All Tumors Hyper Me",title="Size") + geom_abline(slope = 1, intercept=0)
d <- data.frame(as.data.frame(qqplot(seq1.meta$gc, seq2.meta$gc, plot.it=FALSE)),Sequences="Originals")
d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$gc, resamp$gc, plot.it=FALSE)),Sequences="Nearest"))
d <- rbind(d,data.frame(as.data.frame(qqplot(seq1.meta$gc, resamp2$gc, plot.it=FALSE)),Sequences="IS"))
p2 <- ggplot(d) + geom_point(aes(x=x, y=y,color=Sequences),stat = "identity", position = "identity", ) + ggplot.clean() + labs(x="CIMP Hyper Me", y="All Tumors Hyper Me",title="GC") + geom_abline(slope = 1, intercept=0)
g1 <- arrangeGrob(p1,p2,ncol=2)

png(filename="output/resamp.qq.png",width=1400,height=800,res=120)
distance.stats.grob <- arrangeGrob(tableGrob(distance.stats))
grid.arrange(g1, distance.stats.grob, ncol=1,heights=c(4/5,1/5))
dev.off()


###############################################
#Find differential motifs from CIMP HyperMe versus All Tumors HyperMe
#size and GC match farily well based on QQ plots (after filtering out long seqs from the all group)
sim.seq <- seq2[names(seq2) %in% seq2.meta$name]

sim.nSeqs <- length(sim.seq)

writeXStringSet(sim.seq, "output/simseq.fasta", append=FALSE, format="fasta")

runFIMO("sim","output/simseq.fasta",motif.path)
fimo.out.sim <- readFIMO("output/fimo_out_sim/fimo.txt")
fimo.out.sim.counts <- calcMotifCounts(fimo.out.sim,q.cutoff)
results <- calcEnrichmentBinom(fimo.out.counts,nSeqs,fimo.out.sim.counts,sim.nSeqs)
results.sort <- results[order(results$pvalue, decreasing=FALSE),]
head(results.sort)
tail(results.sort)
write.csv(results.sort,file="output/diffmotif.bkg.allhyperme.csv",row.names=FALSE)

