###############################################################################
## Finds differentially enriched motifs in a FIMO output
##
###############################################################################

###############################################
## Settings
###############################################

#Sequences in which to find motifs (FASTA)
seq.path <- "../RunMEME/ColonSeqGR2012/CIMPHyperMe.fasta"
genome.path <- "../RunMEME/ColonSeqGR2012/hg18.fa"

#filter of q values from FIMO output
q.cutoff <- 0.05

###############################################
## Libraries
###############################################
library(ShortRead)
library(plyr)
library(reshape)
library(foreach)
library(doMC)
registerDoMC(15)
library(stringr)

###############################################
## Functions
###############################################
runFIMO <- function(runname,fasta.path,motifs.path)
{
	system(paste("fimo -oc output/fimo_out_",runname," ",motifs.path," ",fasta.path," &> output/fimo_log_",runname,".txt",sep=""))
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

drawBackgroundSet <- function(seq,nSimSeqs=10000,windowSize=50)
{
	#input: sequence set of originals read in from FASTA as a DNAStringSet
	#output: drawn background sequence as DNAStringSet

	#TODO: fix bug if chr freq = 0

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
	sim.bed[sim.bed$start<0,]$start <- 0
	write.table(sim.bed,file="output/simseq.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
	system(paste("fastaFromBed -fi ",genome.path," -bed output/simseq.bed -fo output/simseq.fasta",sep=""))
	system("sed -i 's/:/-/' output/simseq.fasta")

	sim.seq <- readDNAStringSet("output/simseq.fasta")
	sim.seq
}


drawBackgroundSetFromRegions <- function()
{
	#TODO - draw a background set only from another sequence database (FASTA)
}


drawBackgroundSetResampled <- function()
{
	#TODO - draw a background set using importance based resampling to match by both length and GC content
}

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
	#TODO make it apply over whole DNAStringSet object
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
fimo.out <- readFIMO(fimo.out.path)

#build matrix of motif occurance counts for each sequence
fimo.out.counts <- calcMotifCounts(fimo.out,q.cutoff)

###############################################
## Generate Background Sequence Set
###############################################
sim.seq <- drawBackgroundSet(seq,10000)
sim.nSeqs <- length(sim.seq)

###############################################
## Run FIMO on Background Sequences
###############################################
runFIMO("sim","output/simseq.fasta","../TRANSFAC/ConvertToMEME/TRAJAS.Human.named.meme")
fimo.out.sim <- readFIMO("output/fimo_out_sim/fimo.txt")

###############################################
## Perform Enrichment Test
###############################################

#just going to load in some other seq set for now
#fimo.out.sim.path <- "../RunMEME/FIMORun/fimo_out_AllHyperMe/fimo.txt"
#fimo.out.sim <- readFIMO(fimo.out.sim.path)

fimo.out.sim.counts <- calcMotifCounts(fimo.out.sim,q.cutoff)
#seq.path2 <- "../RunMEME/ColonSeqGR2012/AllHyperMe.fasta"
#seq2 <- readDNAStringSet(seq.path2)

results <- calcEnrichmentBinom(fimo.out.counts,nSeqs,fimo.out.sim.counts,sim.nSeqs)
results.sort <- results[order(results$pvalue, decreasing=FALSE),]
head(results.sort)
tail(results.sort)
