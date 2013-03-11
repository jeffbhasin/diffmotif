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
rmsk <- readRepeatMasker("hg18","ignore/ucsc/")

repeatPer <- getRepeatPercent(regions.ranges)

# distTSS
distTSS <- getDistTSS(regions.ranges,ann)

# distTSE
distTSE <- getDistTSE(regions.ranges,ann)

head(data.frame(seq.list,repeatPer,distTSS,distTSE))

#write a file with these variables included
#regions.annotated2 <- cbind(regions,repeatPer,genic,upstream,downstream,utr,exons,genic.genes)
#write.csv(regions.annotated2,file="output/annotation/overlap_dms_covars.csv",row.names=FALSE)

# -----------------------------------------------------------------------------

# =============================================================================
