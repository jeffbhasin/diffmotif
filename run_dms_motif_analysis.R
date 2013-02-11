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
## Analysis Functions
###############################################

###############################################
## Main Body
###############################################

###############################################
#Find differential motifs from CIMP HyperMe versus background drawn from whole genome
run.name <- "CIMPHyperMe"
seq.path <- "../RunMEME/ColonSeqGR2012/CIMPHyperMe.fasta"
motif.path <- "../TRANSFAC/ConvertToMEME/TRAJAS.Human.named.meme"
runFIMO(run.name,seq.path,motif.path)
fimo.out <- readFIMO(paste("output/fimo_out_",run.name,"/fimo.txt",sep=""))
q.cutoff <- 0.05
fimo.out.counts <- calcMotifCounts(fimo.out,q.cutoff)
seq <- readDNAStringSet(seq.path)
nSeqs <- length(seq)
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
#Find differential motifs from CIMP HyperMe versus resampled background set (importance based resampling)

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
hist(seq1.meta$size,breaks=breaks)
hist(seq2.meta$size,breaks=breaks)
qqplot(seq1.meta$size,seq2.meta$size)
hist(seq1.meta$gc,breaks=breaks.gc)
hist(seq2.meta$gc,breaks=breaks.gc)
qqplot(seq1.meta$gc,seq2.meta$gc)
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

