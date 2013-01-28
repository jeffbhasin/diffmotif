## Filename:
## By: Yaomin Xu
## Date:
## Notes:
## =============================================================================

order.factor.chr <- function(X, var) {
  if(missing(var)) {
    stopifnot("chr" %in% names(X))
    var <- 'chr'
  }
  this.fac <- X[,var]
  ##fac.lev <- levels(this.fac)
  ##fac.lab <- levels(this.fac)
  fac.ord <- sub("chrX", "chr23", as.character(this.fac))
  fac.ord <- sub("chrY", "chr24", as.character(fac.ord))
  fac.ord <- as.numeric(sub("chr","",fac.ord))
  this.fac.reorder <- reorder(this.fac, fac.ord, mean)
  X[,var] <- this.fac.reorder
  X
}

order.factor.pattern <- function(X,var, digit=3) {
  if(missing(var)) {
    var <- 'pattern'
  }
  this.fac <- X[,var]
  if(class(this.fac) != "integer") this.fac <- as.integer(as.character(this.fac))
  fac.new <- as.factor(substr(as.character(10^digit+this.fac), 2, digit+1))
  X[,var] <- fac.new
  X
}

factor.chr.names <- function(names) {
 
  chrstr <- unlist(lapply(strsplit(names, "-"), function(x) x[2]))
  chrstr
}

site.extend <- function(S001, S011) {

  R1 <- as.matrix(S001)
  R2 <- as.matrix(S011)

  delta1 <- R1[,2]-R2[,1]
  delta2 <- R1[,1]-R2[,2]
  olp <- (delta1*delta2)<0

  off1 <- R1[,1]-R2[,1]
  off2 <- R1[,2]-R2[,2]

  inside <- (off1*off2)<0

  ext1 <- -1*off1*(off1<0)
  ext2 <- off2*(off2>0)

  size1 <- R1[,2]- R1[,1]
  size2 <- R2[,2]-R2[,1]

  ext.inside <- ext1 + ext2
  ext.inside[!olp] <- size1[!olp]
  if.ext <- !(inside & (off1 >0))
  data.frame(if.ext=if.ext, nonCIMP=size2, CIMP = size2+ext.inside, ext=ext.inside*if.ext)
}

region.pattern.size <- function(sites.report,
                                patt.normal=c("110","101","100"),
                                patt.nonCIMP=c("011","110","010"),
                                patt.CIMP=c("011","001","101")) {
  R.ids <- unique(sites.report$R.id)
  output <- NULL
  for (rid in R.ids) {
    sites.rid <- subset(sites.report, subset=R.id == rid)
    R.patt <- paste(sites.rid$pattern, collapse="_", sep="")
    size.normal <- size.nonCIMP <- size.CIMP <- 0
    size.normal <- sum(subset(sites.rid, subset=pattern%in%patt.normal, select="size",drop=T))
    size.nonCIMP <- sum(subset(sites.rid, subset=pattern%in%patt.nonCIMP, select="size",drop=T))
    size.CIMP <- sum(subset(sites.rid, subset=pattern%in%patt.CIMP, select="size",drop=T))
    output <- rbind(output, data.frame(R.id=rid,
                                       normal=size.normal,
                                       nonCIMP=size.nonCIMP,
                                       CIMP=size.CIMP,
                                       R.patt.long=R.patt,
                                       R.patt=as.character(R.patt)))
  }
  output
}

makegroup.region.patt.size <- function(region.pattsize,
                                       ext.pos.grp=c(
                                         '011_011_001',
                                         '011_001_011_001',
                                         '011_001_011',
                                         '011_001_001_011',
                                         '011_001',
                                         '001_011_011',
                                         '001_011_001_011',
                                         '001_011_001',
                                         '001_011',
                                         '011_011_001',
                                         '011_001_011_011',
                                         '011_011_001_011'),
                                       ext.neg.grp=c(
                                         '110_011',
                                         '011_010_011',
                                         '011_010',
                                         '010_011'),
                                       indv.nonCIMP=c(
                                         '110_010',
                                         '110',
                                         '010_110',
                                         '010'
                                         ),
                                       indv.both=c(
                                         '011_011',
                                         '011'
                                         ),
                                       indv.CIMP=c(
                                         '101',
                                         '001_101',
                                         '001'
                                         )) {

  ##region.pattsize.sel <- subset(region.pattsize, subset=R.patt!='100')
  ext.region <- region.pattsize$R.patt %in% c(ext.pos.grp,ext.neg.grp)
  ext.dir <- rep(0, nrow(region.pattsize));
  ext.dir[region.pattsize$R.patt %in% ext.pos.grp] <- '1'
  ext.dir[region.pattsize$R.patt %in% ext.neg.grp] <- '-1'
  indv.type <- rep("", nrow(region.pattsize))
  indv.type[region.pattsize$R.patt %in% indv.nonCIMP] <- "nonCIMP"
  indv.type[region.pattsize$R.patt %in% indv.CIMP] <- "CIMP"
  indv.type[region.pattsize$R.patt %in% indv.both] <- "both"
  indv.size <- region.pattsize$CIMP
  indv.size[region.pattsize$R.patt %in% indv.nonCIMP] <- region.pattsize$nonCIMP[region.pattsize$R.patt %in% indv.nonCIMP]
  data.frame(region.pattsize,
            ext.region=ext.region,
            ext.dir=ext.dir,
            indv.type=indv.type,
            indv.size=indv.size)
}

prepare.site.genomepos <- function(sites) {
  require(quantsmooth)
  CHR <- sub("chr", "", sites$chr)
  MapInfo <- with(sites, (start + end)/2)
  data.frame(CHR, MapInfo)
}

my.prepareGenomePlot <- function (chrompos, cols = "grey50", plot=T, paintCytobands = FALSE, 
                               bleach = 0, topspace = 1, organism, sexChromosomes = FALSE, 
                               units = c("bases", "cM", "ISCN"), ...) 
{
  cytobandWidth <- 0.075
  par(mar = c(1, 4, 2, 3) + 0.1)
  if (!missing(organism)) {
    units <- match.arg(units)
    organism <- match.arg(organism, c("hsa", "mmu", "rno"))
    chrom.n <- switch(organism, hsa = 22, mmu = 19, rno = 20)
    chrs2 <- factor(numericCHR(chrompos[, "CHR"]), levels = c(1:chrom.n, if (sexChromosomes) c(98, 99) else NULL))
    if (organism == "hsa") 
      lens <- lengthChromosome(levels(chrs2), units = units)
    else lens <- sapply(split(chrompos[, "MapInfo"], chrs2), function(x) max(c(0, x)))
    names(lens) <- characterCHR(names(lens))
    cols <- rep(cols, length.out = length(lens))
    names(cols) <- names(lens)
    dwidth <- NULL
    for (i in 1:(chrom.n%/%2)) dwidth[i] <- lens[i] + lens[chrom.n + 1 - i]
    if (chrom.n%%2 == 1) dwidth <- c(dwidth, lens[chrom.n%/%2 + 1])
    if (sexChromosomes)  dwidth <- c(dwidth, lens["X"] + lens["Y"])
    maxdwidth <- max(dwidth) * 1.05
    
    ##leftrow <- c(if (sexChromosomes) "X" else NULL, ((chrom.n + 1)%/%2):1)
    leftrow <- c(if (sexChromosomes) 98 else NULL, ((chrom.n + 1)%/%2):1)
    ##rightrow <- c(if (sexChromosomes) "Y" else NULL,
    rightrow <- c(if (sexChromosomes) 99 else NULL,
                  if (chrom.n%%2 ==1) "" else NULL,
                  ((chrom.n + 1)%/%2 + 1):chrom.n)
    if(plot) {
      plot(c(0, maxdwidth), c(0.5, 0.5 + length(dwidth) + topspace), 
           type = "n", ylab = "Chromosome", xlab = "", axes = FALSE, 
           las = 2, ...)
      axis(2, c(1:length(dwidth)), characterCHR(leftrow), las = 2)
      axis(4, c(1:length(dwidth)), characterCHR(rightrow),las = 2)
      if (paintCytobands && organism == "hsa") {
        for (i in 1:length(dwidth)) {
          if (lens[characterCHR(leftrow[i])] > 0) 
            paintCytobands(leftrow[i],
                           c(0, i + cytobandWidth/2), 
                           "bases",
                           width = cytobandWidth,
                           length.out = lens[characterCHR(leftrow[i])], 
                           legend = FALSE,
                           bleach = bleach)
          if (rightrow[i] != "" && lens[characterCHR(rightrow[i])] > 0)  
            paintCytobands(rightrow[i],
                           c(maxdwidth - lens[characterCHR(rightrow[i])], i + cytobandWidth/2), "bases",
                           width = cytobandWidth,
                           length.out = lens[characterCHR(rightrow[i])],
                           legend = FALSE,
                           bleach = bleach)
        }
      }
      else {
        for (i in 1:length(dwidth)) {
          lines(c(0, lens[characterCHR(leftrow[i])]), c(i, i), col = cols[leftrow[i]], 
                lwd = 2)
          if (rightrow[i] != "") 
            lines(c(maxdwidth - lens[characterCHR(rightrow[i])], maxdwidth), 
                  c(i, i), col = cols[rightrow[i]], lwd = 2)
        }
      }
    }
    dchrompos <- matrix(0, nrow = nrow(chrompos), ncol = 2, 
                        dimnames = list(rownames(chrompos), c("CHR", "MapInfo")))
    for (i in 1:length(rightrow)) if (rightrow[i] != "") {
      probes <- numericCHR(chrompos[, "CHR"]) == rightrow[i]
      dchrompos[probes, 2] <- chrompos[probes, "MapInfo"] + maxdwidth - lens[characterCHR(rightrow[i])]
      dchrompos[probes, 1] <- i
    }
    for (i in 1:length(leftrow)) {
      probes <- numericCHR(chrompos[, "CHR"]) == leftrow[i]
      dchrompos[probes, 2] <- chrompos[probes, "MapInfo"]
      dchrompos[probes, 1] <- i
    }
  }
  else {
    chrs2 <- factor(numericCHR(chrompos[, "CHR"]))
    lens <- sapply(split(chrompos[, "MapInfo"], chrs2), max)
    m <- length(lens)
    cols <- rep(cols, length.out = m)
    maxdwidth <- max(lens)
    if(plot) {
      plot(c(0, maxdwidth), c(0.5, m + 0.5 + topspace), type = "n", 
           ylab = "Chromosome", xlab = "", axes = FALSE, las = 2, 
           ...)
      axis(2, c(m:1), characterCHR(names(lens)), las = 2)
      for (i in 1:m) lines(c(0, lens[i]), c(m + 1 - i, m + 
                                            1 - i), col = cols[as.numeric(names(lens))], lwd = 2)
    }
    dchrompos <- chrompos
    dchrompos[, 1] <- m + 1 - as.numeric(chrs2)
  }
  dchrompos
}

reformat.no111 <- function(sites,offset=50) {
  sites.no111 <- subset(sites, subset=pattern!='111')
  sites.no111 <- order.factor.pattern(sites.no111)
  sites.no111 <- order.factor.chr(sites.no111)
  sites.no111 <- with(sites.no111, transform(sites.no111, size=end-start+offset))
  sites.no111
}

###... Check the reasons for non-sharp patterns
collapse.1.pattern <- function(x){
  fullp <- unlist(strsplit(x, "_"))
  .patt.1 <- fullp[-length(fullp)]
  .patt.2 <- fullp[-1]
  idx <- c(which(.patt.1!=.patt.2), length(fullp))
  paste(fullp[idx], collapse="_", sep="")
}

collapse.pattern <- function(x) {
  sapply(x, collapse.1.pattern, simplify=T)
}

extended.pattern <- function(x) {
  extp.idx <- sapply(x, function(y) length(unlist(strsplit(y, "_"))) > 1, simplify=T)
  x[extp.idx]
}


      
  

