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


#Find differential motifs from CIMP HyperMe versus background of dinucleotide shuffled source

#Find differential motifs from CIMP HyperMe versus resampled background set (importance based resampling)
