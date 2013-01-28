###############################################################################
## Finds differentially enriched motifs in a FIMO output
##
###############################################################################

###############################################
## Settings
###############################################

#Sequences in which to find motifs (FASTA)
seq.path <- "../RunMEME/ColonSeqGR2012/CIMPHyperMe.fasta"

#filter of q values from FIMO output
q.cutoff <- 0.01

###############################################
## Libraries
###############################################
library(ShortRead)
library(plyr)
library(reshape)
library(foreach)
library(doMC)
registerDoMC(15)

###############################################
## Functions
###############################################
runFIMO <- function(runname,fasta.path,motifs.path)
{
	system(paste("fimo -oc fimo_out_",runname," ",motifs.path," ",fasta.path," &> fimo_log_",runname,".txt"),sep="")
}

readFIMO <- function(fimo.out.path)
{
	read.table(file=fimo.out.path,header=TRUE,comment.char="",sep="\t")
}

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

		data.frame(motif=rownames(counts1.bin)[i],pvalue=pValue)
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

drawBackgroundSet <- function()
{
	#TODO
}

shuffleBackgroundSet <- function()
{
	#TODO
}

###############################################
## Run FIMO
###############################################
seq <- readDNAStringSet(seq.path)
nSeqs <- length(seq)

#TODO: system() call wrapped with R function to run FIMO directly

fimo.out.path <- "../RunMEME/FIMORun/fimo_out_CIMPHyperMe/fimo.txt"

###############################################
## Load FIMO Output
###############################################

#load FIMO text output file
fimo.out <- read.table(file=fimo.out.path,header=TRUE,comment.char="",sep="\t")

#build matrix of motif occurance counts for each sequence
fimo.out.counts <- calcMotifCounts(fimo.out,q.cutoff)

###############################################
## Generate Background Sequence Set
###############################################
#calcMotifFreqs(fimo.out.sim)
#nSimSeqs <- 10000


###############################################
## Run FIMO on Background Sequences
###############################################
#system("fimo ")

###############################################
## Perform Enrichment Test
###############################################

#just going to load in some other seq set for now
fimo.out.sim.path <- "../RunMEME/FIMORun/fimo_out_AllHyperMe/fimo.txt"
fimo.out.sim <- readFIMO(fimo.out.sim.path)
fimo.out.sim.counts <- calcMotifCounts(fimo.out.sim,q.cutoff)
seq.path2 <- "../RunMEME/ColonSeqGR2012/AllHyperMe.fasta"
seq2 <- readDNAStringSet(seq.path2)
nSeqs2 <- length(seq2)

results <- calcEnrichmentBinom(fimo.out.counts,nSeqs,fimo.out.sim.counts,nSeqs2)
results <- results[order(results$pvalue, decreasing=FALSE),]
