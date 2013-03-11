# #############################################################################
# Usage of annotation function to annotate a list of regions
# Author: Jeffrey Bhasin <jmb85@case.edu>
# Created: 2013-03-11
# #############################################################################

# =============================================================================
# Packages and Globals
source("Rmotif.R")
nCores <- 15
registerDoMC(nCores)

# =============================================================================

# =============================================================================
# Local Functions

# -----------------------------------------------------------------------------
# Table reading function
# Input: path to flatfile table
# Output: dataframe of table
read.table.ucsc <- function(path)
{
	#one function to read in all tables exported from UCSC the same way
	read.table(file=path, comment.char="", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="")
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Table reading function for big tables - autodetects types from first 5 rows
# Input: path to flatfile table
# Output: dataframe of table
read.table.ucsc.big <- function(path)
{
	#detects colClasses on first 5 rows to speed up a big read in
	tab5rows <- read.table(file=path, comment.char="", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="", nrows = 5)
	classes <- sapply(tab5rows, class)
	tabAll <- read.table(file=path, comment.char="", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="", colClasses = classes)
	tabAll
}
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Read UCSC table files from disk and join all related tables
# Input: genome id, path with trailing "/"
# Output: joined annotation table dataframe (one row per gene isoform)
readUCSCAnnotation <- function(genome="hg19",path="")
{
	#load downloaded tables
	knownGene.path <- paste(path,genome,".knownGene.txt",sep="")
	knownCanonical.path <- paste(path,genome,".knownCanonical.txt",sep="")
	kgXref.path <- paste(path,genome,".kgXref.txt",sep="")
	cytoBand.path <- paste(path,genome,".cytoBand.txt",sep="")

	knownGene <- read.table.ucsc(knownGene.path)
	knownCanonical <- read.table.ucsc(knownCanonical.path)
	kgXref <- read.table.ucsc(kgXref.path)
	cytoBand <- read.table.ucsc.big(cytoBand.path)

	#join in gene symbols
	kg.sub <- knownGene
	names(kg.sub)[1] <- "name" 
	kgXref.sub <- kgXref
	names(kgXref.sub)[1] <- "name"

	kg.ann.1 <- join(kg.sub,kgXref.sub,by="name",type="left")

	#join in canonical annotation (largest coding seq from nearest cluster)
	kgCanon.sub <- data.frame(name=knownCanonical$transcript,canonical="1",stringsAsFactors=FALSE)
	
	kg.ann.2 <- join(kg.ann.1,kgCanon.sub,by="name",type="left")

	kg.ann <- kg.ann.2

	#join in cytological bands based on gene position
	kg.ranges <- GRanges(seqnames=kg.ann$chrom,ranges=IRanges(start=kg.ann$txStart,end=kg.ann$txEnd),strand=kg.ann$strand,id=kg.ann$name)
	cb.ranges <- with(cytoBand, GRanges(seqnames=X.chrom,ranges=IRanges(start=chromStart,end=chromEnd),band=name,gstain=gieStain))
	
	#band genes start in
	kg.ranges$CytoBandStartIndex <- findOverlaps(kg.ranges,cb.ranges,select="first")
	#band genes end in (in case they span more than one)
	kg.ranges$CytoBandEndIndex <- findOverlaps(kg.ranges,cb.ranges,select="last")
	#total bands spanned
	kg.ranges$CytoBandsSpan <- countOverlaps(kg.ranges,cb.ranges)
	#create DF from GR metadata
	cyto.map <- data.frame(name=kg.ranges$id,count=kg.ranges$CytoBandsSpan,start=kg.ranges$CytoBandStartIndex,end=kg.ranges$CytoBandEndIndex,stringsAsFactors=FALSE)

	cyto.map$startname <- "NA"
	cyto.map[is.na(cyto.map$start)==FALSE,]$startname <- cb.ranges[cyto.map[is.na(cyto.map$start)==FALSE,]$start]$band
	cyto.map$endname <- "NA"
	cyto.map[is.na(cyto.map$start)==FALSE,]$endname <- cb.ranges[cyto.map[is.na(cyto.map$start)==FALSE,]$end]$band
	cyto.map <- data.frame(name=cyto.map$name,cytoBandSpan=cyto.map$count,cytoBandStart=cyto.map$startname,cytoBandEnd=cyto.map$endname,stringsAsFactors=FALSE)

	#join cyto map data back into main annotation by id
	kg.ann <- join(kg.ann,cyto.map,by="name",type="left")

	#add 1-based start coordinate column to remind viewer about this aspect of UCSC data
	#all ucsc start coords are 0 based and this code maintains that
	kg.ann$txStart.1based <- kg.ann$txStart + 1

	kg.ann$cdsStart.1based <- kg.ann$cdsStart + 1

	kg.ann
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Counts of overlapping gene isoforms
# Input: genomic ranges object
# Output: vector of counts of how many gene isoforms overlap with region
getGenicOverlap <- function(regions.ranges, ann)
{
	ann.ranges <- with(ann, GRanges(seqnames=chrom,ranges=IRanges(start=txStart.1based,end=txEnd)))
	overlap <- countOverlaps(regions.ranges,ann.ranges)
	#overlap[!is.na(overlap)] <- 1
	#overlap[is.na(overlap)] <- 0
	overlap
}

# -----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# List of gene names overlapped by each range
# Input:
# Output: one column with a list of overlapping gene symbols for each row in the query regions
getGenicOverlapGenes <- function(regions.ranges, ann)
{
	ann.ranges <- with(ann, GRanges(seqnames=chrom,ranges=IRanges(start=txStart.1based,end=txEnd)))
	overlap <- findOverlaps(regions.ranges,ann.ranges)
	#overlap[!is.na(overlap)] <- 1
	#overlap[is.na(overlap)] <- 0

	overlap <- as.data.frame(overlap)
	overlap$name <- ann[overlap$subjectHits,]$geneSymbol

	out <- foreach(i=1:length(regions.ranges),.verbose=FALSE,.combine="c") %do%
	{
		hits <- overlap[overlap$queryHits==i,]
		genes <- unique(hits$name)
		paste(genes,collapse=", ")
	}

	out
}

# -----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Find overlaps with upstream regions of genes
# Input: before = number bps upstream of TSS, after = number bps downstream of TSS
# Output:
getUpstreamOverlap <- function(regions.ranges, ann, before=1000, after=500)
{
	# add offsets, accounting for strandedness
	ups <- with(ann,data.frame(chr=chrom,start=txStart.1based,end=txEnd,strand=strand))
	ups$start.us <- NA
	ups$end.us <- NA

	ups[ups$strand=="+",]$start.us <- ups[ups$strand=="+",]$start - before
	ups[ups$strand=="+",]$end.us <- ups[ups$strand=="+",]$start + after


	ups[ups$strand=="-",]$start.us <- ups[ups$strand=="-",]$end - after
	ups[ups$strand=="-",]$end.us <- ups[ups$strand=="-",]$end + before

	ann.ranges <- with(ups, GRanges(seqnames=chr,ranges=IRanges(start=start.us,end=end.us)))
	overlap <- countOverlaps(regions.ranges,ann.ranges)
	#overlap[!is.na(overlap)] <- 1
	#overlap[is.na(overlap)] <- 0
	overlap
}
# -----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Find overlaps with downstream regions of genes
# Input: before = number bps upstream of txEnd, after = number bps downstream of txEnd
# Output:
getDownstreamOverlap <- function(regions.ranges, ann, before=500, after=1000)
{
	# add offsets, accounting for strandedness
	downs <- with(ann,data.frame(chr=chrom,start=txStart.1based,end=txEnd,strand=strand))
	downs$start.ds <- NA
	downs$end.ds <- NA

	downs[downs$strand=="+",]$start.ds <- downs[downs$strand=="+",]$end - before
	downs[downs$strand=="+",]$end.ds <- downs[downs$strand=="+",]$end + after


	downs[downs$strand=="-",]$start.ds <- downs[downs$strand=="-",]$start - after
	downs[downs$strand=="-",]$end.ds <- downs[downs$strand=="-",]$start + before

	ann.ranges <- with(downs, GRanges(seqnames=chr,ranges=IRanges(start=start.ds,end=end.ds)))
	overlap <- countOverlaps(regions.ranges,ann.ranges)
	#overlap[!is.na(overlap)] <- 1
	#overlap[is.na(overlap)] <- 0
	overlap
}
# -----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# 3' UTR Overlaps (gaps between cdsEnd and txEnd)
# Input:
# Output:
get3primeUTROverlap <- function(regions.ranges, ann)
{
	# filter for genes with a 3' UTR, accounting for strandedness
	utr <- with(ann,data.frame(chr=chrom, txStart=txStart.1based, txEnd=txEnd, strand=strand, cdsStart=cdsStart.1based, cdsEnd=cdsEnd))

	# filter non-coding transcripts which UCSC codes as cdsStart==cdsEnd
	utr <- utr[utr$cdsStart!=(utr$cdsEnd+1),]

	# filter out if cdsEnd == txEnd for (+) strand
	utr.p <- utr[(utr$strand=="+")&(utr$cdsEnd!=utr$txEnd),]

	# filter out if cdsStart == txStart for (-) strand because these are really the ends
	utr.m <- utr[(utr$strand=="-")&(utr$cdsStart!=utr$txStart),]

	# build list of utr regions
	utr.p$start.utr <- utr.p$cdsEnd
	utr.p$end.utr <- utr.p$txEnd

	utr.m$start.utr <- utr.m$txStart
	utr.m$end.utr <- utr.m$cdsStart

	reg <- rbind(utr.p,utr.m)

	ann.ranges <- with(reg, GRanges(seqnames=chr,ranges=IRanges(start=start.utr,end=end.utr)))
	overlap <- countOverlaps(regions.ranges,ann.ranges)
	#overlap[!is.na(overlap)] <- 1
	#overlap[is.na(overlap)] <- 0
	overlap
}
# -----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Parse exon lists
# Input: data.frame from annotation reader
# Output: data.frame with each row representing an exon range
parseExons <- function(ann)
{
	#split by chrs to make it go faster
	chrs <- unique(ann$chrom)

	#my.ann <- ann[ann$chrom=="chr1",]	

	# parse out exon lists to give regions list where each individual exon is a range

	#ex <- foreach(j=1:length(chrs))
	ex <- foreach(j=1:length(chrs),.verbose=TRUE,.combine="rbind") %dopar%
	{
		my.ann <- ann[ann$chrom==chrs[j],]
		foreach(i=1:nrow(my.ann),.verbose=FALSE,.combine="rbind") %do%
		{
			chr <- my.ann[i,]$chrom
			starts <- as.numeric(unlist(strsplit(my.ann[i,]$exonStarts,",")))
			# correct for UCSC's 0-based system
			starts <- starts + 1
			ends <- as.numeric(unlist(strsplit(my.ann[i,]$exonEnds,",")))

			data.frame(chr=chr,start=starts,end=ends)
		}
	}


	ex
}

# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Overlaps with exons
# Input: exon table from parseExons
# Output:
getExonOverlap <- function(regions.ranges, ann.ex)
{
	ann.ranges <- with(ann.ex, GRanges(seqnames=chr,ranges=IRanges(start=start,end=end)))
	overlap <- countOverlaps(regions.ranges,ann.ranges)
	#overlap[!is.na(overlap)] <- 1
	#overlap[is.na(overlap)] <- 0
	overlap
}
# -----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# Function 2
# Input:
# Output:

# -----------------------------------------------------------------------------

# =============================================================================

# =============================================================================
# Analysis Code

# -----------------------------------------------------------------------------
# Produce annotation table and summary table for a given set of regions
ann <- readUCSCAnnotation(genome="hg18",path="ignore/ucsc/")
ann.ex <- parseExons(ann)

regions.path <- "ignore/DiffMe.GR2012ST1.csv"
regions <- read.csv(file=regions.path,header=TRUE,stringsAsFactors=FALSE)

# want to use 0-based coords internally as UCSC does
# regions$start <- regions$start-1
# seems like they are already 0-based, going to make them 1-based so widths come out to be multiples of 50
# since GenomicRanges seems to want 1-based coords, we'll use those in all range operations

regions.ranges <- with(regions,GRanges(seqnames=chr,ranges=IRanges(start=Start+1,end=End)))

genic <- getGenicOverlap(regions.ranges,ann)
genic.genes <- getGenicOverlapGenes(regions.ranges,ann)
upstream <- getUpstreamOverlap(regions.ranges,ann,1000,500)
downstream <- getDownstreamOverlap(regions.ranges,ann,500,1000)
utr <- get3primeUTROverlap(regions.ranges,ann)
exons <- getExonOverlap(regions.ranges,ann.ex)

regions.annotated <- cbind(regions,genic,upstream,downstream,utr,exons,genic.genes)
write.csv(regions.annotated,file="output/annotation/overlap_dms.csv",row.names=FALSE)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Produce summary of counts/percents for each type from contextual annotation above
dms <- regions.annotated
dms$type <- ""
dms[dms$Pattern=="11",]$type <- "AllHyperMe"
dms[dms$Pattern=="100",]$type <- "AllHypoMe"
dms[dms$Pattern=="1",]$type <- "CIMPHyperMe"
dms[dms$Pattern=="110",]$type <- "CIMPHypoMe"
dms[dms$Pattern=="10",]$type <- "NonCIMPHyperMe"
dms[dms$Pattern=="101",]$type <- "NonCIMPHypoMe"

#counts by type to double check
ddply(dms,.(type),nrow)

#how many of each type have regions in categories
#overlap with at least one of a region type

#convert counts to binaries
dms[dms$genic>0,]$genic <- 1
dms[dms$exons>0,]$exons <- 1
dms[dms$utr>0,]$utr <- 1
dms[dms$upstream>0,]$upstream <- 1
dms[dms$downstream>0,]$downstream <- 1

count_genic <- ddply(dms,.(type,genic),nrow,.drop=FALSE)
count_exons <- ddply(dms,.(type,exons),nrow,.drop=FALSE)
count_utr <- ddply(dms,.(type,utr),nrow,.drop=FALSE)
count_upstream <- ddply(dms,.(type,upstream),nrow,.drop=FALSE)
count_downstream <- ddply(dms,.(type,downstream),nrow,.drop=FALSE)

counts <- data.frame(type=count_genic$type,present=count_genic$genic,genic=count_genic$V1,exons=count_exons$V1,utr=count_utr$V1,upstream=count_upstream$V1,downstream=count_downstream$V1)

#counts$total <- rowSums(total[,3:ncol(total)])
#sumrow <- data.frame("Total","",t(colSums(counts[,3:ncol(counts)])))
#names(sumrow) <- names(counts)
#counts <- rbind(counts,sumrow)

#make the same thing just with percents
#what % of cimp hyper me sites are genic?

percents <- function(mp,ma)
{
	m <- mp/(mp+ma)
	m <- m*100
	m <- round(m,digits=0)
	m
}

#matrix of presents
mp <- counts[counts$present==1,3:ncol(counts)]

#matrix of absents
ma <- counts[counts$present==0,3:ncol(counts)]

#feed not present row as v1 and present row as v2 to get back vector of percents
counts_percents <- percents(mp,ma)
counts_percents <- cbind(type=unique(counts[,1]),counts_percents)


counts_overlap <- mp
counts_overlap <- cbind(type=unique(counts[,1]),counts_overlap)

#write CSV summary
write.csv(counts_overlap,file="output/annotation/overlap_counts.csv",row.names=FALSE)
write.csv(counts_percents,file="output/annotation/overlap_percents.csv",row.names=FALSE)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Produce covariate annotation for use with propensity score matching

seq.list <- dms[dms$type=="CIMPHyperMe",]
seq.meta <- with(seq.list,data.frame(chr=chr,start=Start,end=End))
regions.ranges <- with(seq.meta,GRanges(seqnames=chr,ranges=IRanges(start=start+1,end=end)))

# repeatPer
readRepeatMasker <- function(genome,path)
{
	rmsk.path <- paste(path,genome,".rmsk.txt",sep="")
	rmsk <- read.table.ucsc.big(rmsk.path)
	rmsk
}

getRepeatPercent <- function(regions.ranges,rmsk)
{
	#intersect with rmsk table, then collapse into nonoverlapping regions covering all repeats - want the % coverage of these regions versus the entire sequence length

	rmsk.ranges <- with(rmsk,GRanges(seqnames=genoName,ranges=IRanges(start=genoStart+1,end=genoEnd)))

	#per <- foreach(i=1:length(regions.ranges),.combine="c") %do%
	per <- foreach(i=1:length(regions.ranges),.combine="c",.verbose=TRUE) %dopar%
	{
		cov <- intersect(regions.ranges[i],rmsk.ranges)
		cov <- reduce(cov)
		#after reduction, the sum of the widths is the total coverage and we just need to divide this by the original width of the whole sequence
		rpt <- sum(width(cov)) / width(regions.ranges[i])
		rpt
	}

	per
}

rmsk <- readRepeatMasker("hg18","ignore/ucsc/")

repeatPer <- getRepeatPercent(regions.ranges)

# distTSS


# distTSE

# -----------------------------------------------------------------------------

# =============================================================================
