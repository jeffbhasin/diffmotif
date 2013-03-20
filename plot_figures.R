# #############################################################################
# Create Figures for Manuscript Submission
# Author: Jeffrey Bhasin <jmb85@case.edu>
# Created: 2013-03-20
# #############################################################################

# =============================================================================
# Packages and Globals

#source("Rmotif.R")

#nCores <- 15
#registerDoMC(nCores)

# =============================================================================

# =============================================================================
# Local Functions

# -----------------------------------------------------------------------------
# Custom version of covar balance plotting
# Input: 
# Output: 
myPlotCovarDistance <- function(orig.meta,list.meta,cols)
{
	stddist <- function(d1, d2)
	{
		#(mean(d1)-mean(d2))/(sd(c(d1,d2))/2)
		#(100*abs(mean(d1)-mean(d2)))/sqrt(((sd(d1)^2)+(sd(d2)^2))/2)
		(100*(mean(d1)-mean(d2)))/sqrt(((sd(d1)^2)+(sd(d2)^2))/2)
	}

	# calculate distances for each variable
	dists <- foreach(i=1:length(list.meta),.combine="rbind") %do%
	{
		vec <- foreach(u=1:length(cols),.combine="c") %do%
		{		
			stddist(orig.meta[,cols[u]],list.meta[[i]][,cols[u]])
		}
	}
	colnames(dists) <- names(orig.meta)[cols]

	#rownames(dists) <- names(list.meta)
	rownames(dists) <- c("Reference Pool", rep("Random Sample",30), rep("log(Length Model)",30), rep("Full Model",30))

	plot.data <- melt(dists)
	names(plot.data) <- c("matching","variable","stddist")

	plot.data.sub <- plot.data[plot.data$matching==rownames(dists)[1],]
	plot.data.sub$variable <- as.character(plot.data.sub$variable)
	mylevs <- plot.data.sub[sort(plot.data.sub$stddist, index.return=TRUE)$ix,]$variable

	#mylevs <- levels(reorder(x=plot.data[plot.data$matching=="pool",]$variable,X=plot.data[plot.data$matching=="pool",]$stddist, order=FALSE))

	plot.data$matching <- factor(plot.data$matching,levels=unique(rownames(dists)))

	plot.data$variable <- factor(plot.data$variable,levels=mylevs)

	ggplot(plot.data, aes(x=variable,y=stddist,col=matching)) + geom_boxplot(data=plot.data, outlier.shape=NA) + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(linetype=3, colour="grey50"), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.key.size = unit(0.8, "lines"), axis.line = element_line(colour = "grey50"), axis.text=element_text(colour="black")) + geom_abline(intercept=0,slope=0,col="grey50") + coord_flip() + labs(main="Covariate Balance",y="Standardized Distance") + scale_colour_manual(values = brewer.pal(9,"Set1"))
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Function 2
# Input:
# Output:

# -----------------------------------------------------------------------------

# =============================================================================

# =============================================================================
# Analysis Code

# -----------------------------------------------------------------------------
# Load previously run data

#load("output/propensity/seq1.meta.Rd") # also contains other variables (annotation, metadata cols, etc)
#load("output/propensity/seq2.meta.Rd")
#load("output/propensity/30starts/ref.alltse.Rd")
#load("output/propensity/30starts/ref.alltss.Rd")
#load("output/propensity/30starts/ref.sizelog.Rd")
#load("output/propensity/30starts/ref.sub.Rd")

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Figure: Covariate Balance

pdf(file="output/propensity/figures/CovariateBalance.pdf", width=10.5, height=8, paper="USr")
# covar plot
mymeta <- list()
mymeta$pool <- seq2.meta
for(i in 1:length(sub.mine))
{
	mymeta[[ref.sub.names[i]]] <- ref.sub[[i]]$meta
}
for(i in 1:length(sizelog.mine))
{
	mymeta[[ref.sizelog.names[i]]] <- ref.sizelog[[i]]$meta
}
for(i in 1:length(alltss.mine))
{
	mymeta[[ref.alltss.names[i]]] <- ref.alltss[[i]]$meta
}

#customize what fields get plotted
cols.new <- c(6,7,8,9,10)
seq1.meta.new <- seq1.meta
names(seq1.meta.new)[cols.new] <- c("log(Length)", "GC Content", "CpG Islands", "% Repeats", "log(Distance to TSS + 1)")

myPlotCovarDistance(seq1.meta.new, mymeta, cols.new)
dev.off()

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Figure: Heatmap

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

postscript(file="output/propensity/figures/Dendro.Subs.ps", paper="letter")
plot(myden)
dev.off()


mydis <- dist(t(mat2), method="euclidean")
myden <- as.dendrogram(hclust(mydis))
#order.dendrogram(myden)
myden.order <- reorder(myden, -1*rowMeans(t(mat2)[order.dendrogram(myden),]))
mat2 <- mat2[,order.dendrogram(myden.order)]

postscript(file="output/propensity/figures/Dendro.Length.ps", paper="letter")
plot(myden)
dev.off()

mydis <- dist(t(mat3), method="euclidean")
myden <- as.dendrogram(hclust(mydis))
#order.dendrogram(myden)
myden.order <- reorder(myden, -1*rowMeans(t(mat3)[order.dendrogram(myden),]))
mat3 <- mat3[,order.dendrogram(myden.order)]

postscript(file="output/propensity/figures/Dendro.All.ps", paper="letter")
plot(myden)
dev.off()

#combine into one big matrix
mat <- cbind(mat1,mat2,mat3)


#order clustered by only the alltss stuff
# order rows by mean
mydis <- dist(mat3, method="euclidean")
myden <- as.dendrogram(hclust(mydis))
#order.dendrogram(myden)
myden.order <- reorder(myden, -1*rowMeans(mat3[order.dendrogram(myden),]))
mat <- mat[order.dendrogram(myden.order),]

postscript(file="output/propensity/figures/Dendro.RowsAllOnly.ps", paper="letter")
plot(myden)
dev.off()

# make the scale
mybreaks <- c(min(mat),seq(from=1.3,to=max(mat), length.out=20))
mycols <- c("#FFFFFF",colorRampPalette(brewer.pal(9,"Blues")[6:9])(length(mybreaks)-2))

postscript(file="output/propensity/figures/Dendro.RowsAllOnly.ps", paper="letter")
plot(myden)
dev.off()

# plot the heatmap
postscript(file="output/propensity/figures/LogLikHeatmap.ps", paper="letter")
#png(filename="output/tmp.png",width=1400,height=800,res=120)

#binary by sig:

#white=NS, ascending blue=more sig
#mybreaks <- c(min(mat),1.3,max(mat))
#mycols <- colorRampPalette(brewer.pal(9,"Blues"))

col.labels <- rep("",ncol(mat))
col.labels[1] <- "sub"
col.labels[31] <- "sizeLog"
col.labels[61] <- "alltss"

hm <- heatmap.2(mat, Rowv = FALSE, Colv=FALSE, dendrogram="none", key=F, keysize=1.5, trace="none", scale="none", cexCol=1.2, labRow=NA, breaks=mybreaks, main="Heatmap of -log10(adj p)",labCol=col.labels, col=mycols)
#col=colorRampPalette(brewer.pal(9,"Blues"))
write.table(data.frame(column=seq(1,ncol(mat)),model=colnames(mat)), file="output/propensity/figures/LogLikHeatmapColKey.txt", row.names=FALSE)
dev.off()

# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Figure: Venn diagram

#VENN DIAGRAMS

#shared sig factors between the 3 groups
#convert matrix to binary if sig
makeBinary <- function(x)
{
	if(x >= -log10(0.05))
	{
		1
	} else {
		0
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
write.csv(sigs,file="output/propensity/figures/DiffMotifs-counts.csv",row.names=TRUE)

names(sigs) <- c("mat1","mat2","mat3")

n1 <- nrow(sigs[(sigs$mat1 > 0),]) 
n2 <- nrow(sigs[(sigs$mat2 > 0),]) 
n3 <-nrow(sigs[(sigs$mat3 > 0),]) 
n12 <- nrow(sigs[(sigs$mat1 > 0) & (sigs$mat2 > 0),])
n23 <- nrow(sigs[(sigs$mat2 > 0) & (sigs$mat3 > 0),])
n13 <- nrow(sigs[(sigs$mat1 > 0) & (sigs$mat3 > 0),])
n123 <- nrow(sigs[(sigs$mat1 > 0) & (sigs$mat2 > 0) & (sigs$mat3 > 0),])

pdf(file="output/propensity/figures/Venn.pdf", width=10.5, height=8, paper="USr")
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
# Figure: Workflow

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Supp Table: Differential Motifs CSV

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

write.csv(mycsv,file="output/propensity/figures/DiffMotifs.csv",row.names=FALSE)


# -----------------------------------------------------------------------------


# =============================================================================
