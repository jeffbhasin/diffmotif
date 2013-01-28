## Filename:
## By: Yaomin Xu
## Date:
## Notes: Data process for motif analyses
## =============================================================================
setwd("DPT_ws")
library(session)
restore.session("../assets/biology.data.RSession")
source("../motif.fun.R")
source("../fun.R")
source("src/DiffSeq.access.R")
source("experiment.config.R")

ncore <- 10
winsize <- 50
chrs <- get.ws.chrs('full')

TRANSFORM.MOTIF.DATABASES2MEME <- F
WRITE.COMBINE.TRANSFAC.JASPAR.MEME <- T
WRITE.SITE.SEQUENCE.FILES <- T
WRITE.SIM.SITE.SEQUENCE.FILES <- T

refgenome <- read.refgenome(chrs, path='../../Sess.Analysis.MS/Sess.GABP/run.GABP-50bp/DPT_ws/FASTA/')
if(TRANSFORM.MOTIF.DATABASES2MEME) {
  jaspar2meme.command <- paste('jaspar2meme',
                               '-pfm',
                               '-logodds',
                               '~yxu/lib/jaspar/jaspar_CORE/non_redundant/by_tax_group/vertebrates/FlatFileDir',
                               '>',
                               "../jaspar.meme")
  system(jaspar2meme.command)
  
  transfac2meme.command <- paste('transfac2meme',
                                 '-use_acc',
                                 '-species human',
                                 '-logodds',
                                 '~yxu/lib/transfac/TRANSFAC/dat/matrix.dat',
                                 '>',
                               "../transfac.meme")
  system(transfac2meme.command)
}
###... P011
sites.P011 <- subset(sites.report, subset=contrast=="011"&pattern=="011")
sites.P011$siteID <- row.names(sites.P011)

###... P001
sites.P001 <- subset(sites.report, subset=contrast=="001"&pattern=="001")
sites.P001$siteID <- row.names(sites.P001)

###... Extended regions
sites.ext <- subset(sites.report, subset=sites.report$R.id%in%names(regions.extended))
sites.ext$siteID <- row.names(sites.ext)

###... P001 uni-directional
R.ext.P001.uni <- names(regions.extended)[(regions.extended%in%c("001_011","001_011_011","011_001","011_011_001"))]
sites.ext.P001.uni <- subset(subset(sites.report,
                                    subset=sites.report$R.id%in% R.ext.P001.uni),
                             pattern=="001")
sites.ext.P001.uni$siteID <- row.names(sites.ext.P001.uni)
###... P001 bi-directional
R.ext.P001.bi <- names(regions.extended)[(regions.extended%in%c("001_011_001"))]
sites.ext.P001.bi <- subset(subset(sites.report,
                                   subset=sites.report$R.id%in% R.ext.P001.bi),
                            pattern=="001")
sites.ext.P001.bi$siteID <- row.names(sites.ext.P001.bi)

if(WRITE.SITE.SEQUENCE.FILES) {

  siteStrings <- c("sites.P011",
                   "sites.P001",
                   "sites.ext",
                   "sites.ext.P001.uni",
                   "sites.ext.P001.bi")

  for(sitestr in siteStrings) {
    cat(sitestr, "-","\n")
    process.site.sequence <- function(i) {
      chrom <- chrs[i]
      cat(chrom)
      sites.chr <- subset(get(sitestr), subset=chr==chrom)
      if(nrow(sites.chr)>0) {
        site.seq.chr <- get.sequence.chr(sites.chr, refgenome[[chrom]],winsize,ext=150)
      } else {
        site.seq.chr <- NULL
      }
      cat("\tDONE!\n")
      site.seq.chr
    }
    cat("query reference genome: \n")
    sites.seqs <- sapply(seq(along=chrs), process.site.sequence)
    names(sites.seqs) <- chrs
    cat("write site sequence file: \n")
    sites.seqs.filepath <- paste("..", paste(sitestr,"fa", sep="."), sep="/")
    write.sites.seqs(sites.seqs[sapply(sites.seqs, length)>0], sites.seqs.filepath)
  }

}  
### DISCARD THESE CODE
### Commands are right. Running in R does not work.
## fimo.cut <- 1e-4
## fimo.command.transfac <- paste('fimo',
##                                '--verbosity 1',
##                                '--output-qthresh 0.05',
##                                '--o ../fimo_out_transfac',
##                                '--oc ../fimo_out_transfac',
##                                '../transfac.meme',
##                                sites.seqs.filepath,
##                                sep=" ")
## system(fimo.command.transfac)
### DISCARD

###. Motif analysis
##____________________________________________________________________

###.. TRANSFAC and JASPAR Databases
##--------------------------------------

###... Combine databases

transfac.memelist <- read.meme2list("../transfac.meme")
jaspar.memelist <- read.meme2list("../jaspar.meme")
both.memelist <- combine.memelist(transfac.memelist, jaspar.memelist)

if(WRITE.COMBINE.TRANSFAC.JASPAR.MEME) {
  write.memelist2meme(both.memelist, "../both.meme")
}
###... Motif meta data
transfac.matrices <- read.transfacMatrix2list("/home/yxu/lib/transfac/TRANSFAC/dat/matrix.dat")
transfac.meta <- data.frame(get.info.trmatrix(transfac.matrices),
                            cname=toupper(get.info.trmatrix(transfac.matrices)$name))

jaspar.meta <- data.frame(ac=get.ids.memelist(jaspar.memelist),
                          name=get.names.memelist(jaspar.memelist),
                          cname=toupper(get.names.memelist(jaspar.memelist)))
row.names(jaspar.meta) <- jaspar.meta$ac

## combine both
motif.meta <-rbind(transfac.meta[,c('ac','name','cname')],
                   jaspar.meta[,c('ac','name','cname')])
## cleanup motif.meta
motif.meta$cname[motif.meta$cname=="TFAP2A"] <- "AP-2ALPHA"
motif.meta$cname[motif.meta$cname=="E2A"] <- "E47"
motif.meta$cname[motif.meta$cname=="E2F-1"] <- "E2F"
motif.meta$cname[motif.meta$cname=="KLF4"] <- "GKLF"
motif.meta$cname[motif.meta$cname=="NHLH1"] <- "HEN1"
motif.meta$cname[motif.meta$cname=="MTF1"] <- "MTF-1"
motif.meta$cname[motif.meta$cname=="NF-Y"] <- "NFYA"
motif.meta$cname[motif.meta$cname=="REST"] <- "NRSF"
motif.meta$cname[motif.meta$cname=="PAX4"] <- "PAX-4"
motif.meta$cname[motif.meta$cname=="PAX5"] <- "PAX-5"
motif.meta$cname[motif.meta$cname=="RREB-1"] <- "RREB1"
motif.meta$cname[motif.meta$cname=="SREBP1"] <- "SREBP"
motif.meta$cname[motif.meta$cname=="TATA"] <- "TBP"
motif.meta$cname[motif.meta$cname=="NRSE"] <- "NRSF"
motif.meta$cname[motif.meta$cname=="EGR"] <- "KROX"

###.. FIMO
##--------------------------------------

##COMMAND LINE FIMO Runs-------------
systerm("fimo --verbosity 1 --output-qthresh 0.05 --o fimo_out_ext_both --oc fimo_out_ext_both both.meme sites.ext.fa",
        wait=F)
system("fimo --verbosity 1 --output-qthresh 0.05 --o fimo_out_P001_both --oc fimo_out_P001_both both.meme sites.P001.fa",
       wait=F)
system("fimo --verbosity 1 --output-qthresh 0.05 --o fimo_out_P011_both --oc fimo_out_P011_both both.meme sites.P011.fa",
       wait=F)
system("fimo --verbosity 1 --output-qthresh 0.05 --o fimo_out_ext_uni_both --oc fimo_out_ext_uni_both both.meme sites.ext.P001.uni.fa",
       wait=F)
system("fimo --verbosity 1 --output-qthresh 0.05 --o fimo_out_ext_bi_both --oc fimo_out_ext_bi_both both.meme sites.ext.P001.bi.fa",
       wait=F)

##--------------------------------------
### Wait until all the FIMO processes done!!!
topmotif.q.cutoff <- 0.01

###... Readin fimo results
fimo.uni <- read.table("../fimo_out_ext_uni_both/fimo.txt")
fimo.bi <- read.table("../fimo_out_ext_bi_both/fimo.txt")
fimo.ext <- read.table("../fimo_out_ext_both/fimo.txt")
fimo.P001 <- read.table("../fimo_out_P001_both/fimo.txt")
fimo.P011 <- read.table("../fimo_out_P011_both/fimo.txt")

###... Summarize the fimo results
uni.summary <- process.site.freq.summary(fimo.uni, sites.ext.P001.uni, motif.meta, topmotif.q.cutoff)
bi.summary <- process.site.freq.summary(fimo.bi, sites.ext.P001.bi, motif.meta, topmotif.q.cutoff)
ext.summary <- process.site.freq.summary(fimo.ext, sites.ext, motif.meta, topmotif.q.cutoff)
P001.summary <- process.site.freq.summary(fimo.P001, sites.P001, motif.meta, topmotif.q.cutoff)
P011.summary <- process.site.freq.summary(fimo.P011, sites.P011, motif.meta, topmotif.q.cutoff)

###.. Motif logos 
##--------------------------------------
pct.cutoff <- 0
motifs.found.uni <- row.names(uni.summary$hits)[uni.summary$site.freq$pct>pct.cutoff]
motifs.found.bi <- row.names(bi.summary$hits)[bi.summary$site.freq$pct>pct.cutoff]
motifs.found.ext <- row.names(ext.summary$hits)[ext.summary$site.freq$pct>pct.cutoff]
motifs.found.P001 <- row.names(P001.summary$hits)[P001.summary$site.freq$pct>pct.cutoff]
motifs.found.P011 <- row.names(P011.summary$hits)[P011.summary$site.freq$pct>pct.cutoff]

motifs.found <- unique(c(motifs.found.uni,
                         motifs.found.bi,
                         motifs.found.ext,
                         motifs.found.P001,
                         motifs.found.P011))

write.sub.memelist <- function(fullml, ids, file) {
  ml.sub <- select.memelist(fullml, ids)
  write.memelist2meme(ml.sub, file)
}

write.sub.memelist(both.memelist, motifs.found.uni, "../motifs.found.uni.meme")
write.sub.memelist(both.memelist, motifs.found.bi, "../motifs.found.bi.meme")
write.sub.memelist(both.memelist, motifs.found.ext, "../motifs.found.ext.meme")
write.sub.memelist(both.memelist, motifs.found.P001, "../motifs.found.P001.meme")
write.sub.memelist(both.memelist, motifs.found.P011, "../motifs.found.P011.meme")
write.sub.memelist(both.memelist, motifs.found, "../motifs.found.meme")

logo.n <- length(found.motifs)
for(i in seq(logo.n)) {
  cat(found.motifs[i],"\t")
  inputstr <- paste(paste("-i",i,sep=""), "../motifs.found.meme", sep=" ")
  outstr <- paste("-o", paste("../motiflogos.found/",paste(found.motifs[i],"eps",sep="."),sep=""))
  system(paste("ceqlogo -f EPS",
               inputstr,
               outstr)
         )
  cat("DONE\n")
}

###.. MAST analysis
##--------------------------------------
system("mast -oc ../mast_output_uni ../motifs.found.uni.meme ../sites.ext.P001.uni.fa", wait=F)
system("mast -oc ../mast_output_bi ../motifs.found.bi.meme ../sites.ext.P001.bi.fa", wait=F)
system("mast -oc ../mast_output_ext ../motifs.found.ext.meme ../sites.ext.fa", wait=F)
system("mast -oc ../mast_output_P001 ../motifs.found.P001.meme ../sites.P001.fa", wait=F)
system("mast -oc ../mast_output_P011 ../motifs.found.P011.meme ../sites.P011.fa", wait=F)

###. Simulation-based enrichment test 
##____________________________________________________________________

ref.dat <- sites.P011
sample.n <- 10000

sizes <- with(ref.dat, sample(end - start + winsize, sample.n, replace=T))

##chr.range.min <- by(ref.dat, ref.dat$chr, function(x) min(x$start))
chr.range.max <- by(ref.dat, ref.dat$chr, function(x) max(x$end))

chr.table <- table(ref.dat$chr)

chrs <- sample(names(chr.table), sample.n, replace=T, prob=as.vector(chr.table)/sum(chr.table))


df.1 <- data.frame(cbind(freq=as.vector(table(chrs)), max=chr.range.max))

pois <- apply(df.1,
              1,
              function(x) sample.int(x[2], x[1]))

nms <- row.names(df.1)
fqs <- df.1$freq

chr.pois <- sapply(seq(df.1$freq), function(x, nms, fqs) rep(nms[x], fqs[x]), nms=nms, fqs=fqs)

o.dat.1 <- data.frame(chr=unlist(chr.pois), end=unlist(pois), stringsAsFactors=F)
o.dat <- data.frame(o.dat.1, start=o.dat.1$end-sizes, stringsAsFactors=F)

sites.sim <- cbind(o.dat, siteID=seq(nrow(o.dat)))

##WRITE.SIM.SITE.SEQUENCE.FILES <- T
if(WRITE.SIM.SITE.SEQUENCE.FILES) {

  sitestrs <- c("sites.sim")
  chrstr <- unique(sites.sim$chr)

  for(sitestr in sitestrs) {
    cat(sitestr, "-","\n")
    process.site.sequence <- function(i) {
      chrom <- chrstr[i]
      cat(chrom)
      sites.chr <- subset(get(sitestr), subset=chr==chrom)
      if(nrow(sites.chr)>0) {
        site.seq.chr <- get.sequence.chr(sites.chr, refgenome[[chrom]],ext=0)
      } else {
        site.seq.chr <- NULL
      }
      cat("\tDONE!\n")
      site.seq.chr
    }    
    cat("query reference genome: \n")
    sites.seqs <- sapply(seq(along=chrstr), process.site.sequence)
    names(sites.seqs) <- chrstr
    cat("write site sequence file: \n")
    sites.seqs.filepath <- paste("..", paste(sitestr,"fa", sep="."), sep="/")
    write.sites.seqs(sites.seqs[sapply(sites.seqs, length)>0], sites.seqs.filepath)
  }

}  

###... fimo on simulated data
system("fimo --verbosity 1 --output-qthresh 0.05 --o ../fimo_out_sim_both --oc ../fimo_out_sim_both ../both.meme ../sites.sim.fa",wait=F)
fimo.sim <- read.table("../fimo_out_sim_both/fimo.txt")
sim.summary <- process.site.freq.summary(fimo.sim, sites.sim, motif.meta, topmotif.q.cutoff, round.digit=10)

###... Enrichment test

fimo.enrichment.test <- function(summary.tb, s.n, freq.summ.tb,sim.n=10000) {
  
  s.freq <- summary.tb$site.freq
  ref.freq <- freq.summ.tb$site.freq
  pvalues.out <- rep(NA, nrow(s.freq))
  for (i in seq(nrow(s.freq))) {
    s.freq.i <- s.freq[i,]
    ref.freq.p <- if(row.names(s.freq)[i] %in% row.names(ref.freq)) ref.freq[row.names(s.freq)[i],'pct']/100 else 1/sim.n 
    pvalues.out[i] <- binom.test(s.freq.i$count, s.n, p=ref.freq.p, alternative="greater")$p.value
  }
  pvalues.out
}

pvalues.P011 <- fimo.enrichment.test(P011.summary, nrow(sites.P011), sim.summary)
pvalues.P001 <- fimo.enrichment.test(P001.summary, nrow(sites.P001), sim.summary)
pvalues.ext <- fimo.enrichment.test(ext.summary, nrow(sites.ext), sim.summary)
pvalues.bi <- fimo.enrichment.test(bi.summary, nrow(sites.ext.P001.bi), sim.summary)
pvalues.uni <- fimo.enrichment.test(uni.summary, nrow(sites.ext.P001.uni), sim.summary)

print.P011 <- cbind(P011.summary$site.freq, p.value=pvalues.P011)
print.P001 <- cbind(P001.summary$site.freq, p.value=pvalues.P001)
print.ext <- cbind(ext.summary$site.freq, p.value=pvalues.ext)
print.bi <- cbind(bi.summary$site.freq, p.value=pvalues.bi)
print.uni <- cbind(uni.summary$site.freq, p.value=pvalues.uni)

###.. Subset of null sites 
##--------------------------------------
notsigs <- subset(test.report.flat, subset=p.value.1>0.5)

size.tb <- table(sizes)

notsigs.out <- NULL
for (s.i in names(size.tb)) {
  cat(s.i)
  notsigs.i <- subset(notsigs, subset=(end-start+winsize)==as.integer(s.i))
  if(nrow(notsigs.i) <10) notsigs.i <- subset(notsigs, subset=(end-start)%in%(as.integer(s.i)+seq(-100,100, by=50)))
  notsigs.out <- rbind(notsigs.out, notsigs.i[sample.int(nrow(notsigs.i), size.tb[s.i]),])
  cat("\tDONE!\n")
}
add.gw.siteID <- function (tst) 
{
  sid <- paste(tst$chr, tst$pattern, tst$ID, sep = "-")
  out <- transform(s.id = sid, siteID=sid,tst)
  row.names(out) <- seq(nrow(out))
  out
}
sites.notsigs <- add.gw.siteID(notsigs.out)

WRITE.NULL.SITE.SEQUENCE.FILE <- T
if(WRITE.NULL.SITE.SEQUENCE.FILE) {
  sitestrs <- c("sites.notsigs")
  chrstr <- as.character(unique(get(sitestrs)$chr))

  for(sitestr in sitestrs) {
    cat(sitestr, "-","\n")
    process.site.sequence <- function(i) {
      chrom <- chrstr[i]
      cat(chrom)
      sites.chr <- subset(get(sitestr), subset=chr==chrom)
      if(nrow(sites.chr)>0) {
        site.seq.chr <- get.sequence.chr(sites.chr, refgenome[[chrom]],ext=0)
      } else {
        site.seq.chr <- NULL
      }
      cat("\tDONE!\n")
      site.seq.chr
    }    
    cat("query reference genome: \n")
    sites.seqs <- sapply(seq(along=chrstr), process.site.sequence)
    names(sites.seqs) <- chrstr
    cat("write site sequence file: \n")
    sites.seqs.filepath <- paste("..", paste(sitestr,"fa", sep="."), sep="/")
    write.sites.seqs(sites.seqs[sapply(sites.seqs, length)>0], sites.seqs.filepath)
  }
}
  
###... fimo on simulated data
system("fimo --verbosity 1 --output-qthresh 0.05 --o ../fimo_out_notsigs_both --oc ../fimo_out_notsigs_both ../both.meme ../sites.notsigs.fa", wait=F)
fimo.notsigs <- read.table("../fimo_out_notsigs_both/fimo.txt")
notsigs.summary <- process.site.freq.summary(fimo.notsigs, sites.notsigs, motif.meta, topmotif.q.cutoff, round.digit=10)

###... Enrichment test

fimo.enrichment.test <- function(summary.tb, s.n, freq.summ.tb,sim.n=10000) {
  
  s.freq <- summary.tb$site.freq
  ref.freq <- freq.summ.tb$site.freq
  pvalues.out <- rep(NA, nrow(s.freq))
  for (i in seq(nrow(s.freq))) {
    s.freq.i <- s.freq[i,]
    ref.freq.p <- if(row.names(s.freq)[i] %in% row.names(ref.freq)) ref.freq[row.names(s.freq)[i],'pct']/100 else 1/sim.n 
    pvalues.out[i] <- binom.test(s.freq.i$count, s.n, p=ref.freq.p, alternative="greater")$p.value
  }
  pvalues.out
}

pvalues.P011 <- fimo.enrichment.test(P011.summary, nrow(sites.P011), notsigs.summary)
pvalues.P001 <- fimo.enrichment.test(P001.summary, nrow(sites.P001), notsigs.summary)
pvalues.ext <- fimo.enrichment.test(ext.summary, nrow(sites.ext), notsigs.summary)
pvalues.bi <- fimo.enrichment.test(bi.summary, nrow(sites.ext.P001.bi), notsigs.summary)
pvalues.uni <- fimo.enrichment.test(uni.summary, nrow(sites.ext.P001.uni), notsigs.summary)

print.P011 <- cbind(P011.summary$site.freq[,c(1,2,4)], p.value=pvalues.P011)
print.P001 <- cbind(P001.summary$site.freq[,c(1,2,4)], p.value=pvalues.P001)
print.ext <- cbind(ext.summary$site.freq[,c(1,2,4)], p.value=pvalues.ext)
print.bi <- cbind(bi.summary$site.freq[,c(1,2,4)], p.value=pvalues.bi)
print.uni <- cbind(uni.summary$site.freq[,c(1,2,4)], p.value=pvalues.uni)

###.. Clean up motifs
##--------------------------------------
## motifs.by.cname <- by(motif.meta, motif.meta$cname, function(x) x$ac)
## genes.2cleanup <- motifs.by.cname[sapply(motifs.by.cname, length)>1]
###.. GC content
##--------------------------------------
sitestrs <- c("sites.P011",
              "sites.P001",
              "sites.ext",
              "sites.ext.P001.uni",
              "sites.ext.P001.bi",
              "sites.sim",
              "sites.notsigs")
chrstr <- as.character(unique(get(sitestrs)$chr))

for(sitestr in sitestrs) {
  cat(sitestr, "-","\n")
  process.site.sequence <- function(i) {
    chrom <- chrstr[i]
    cat(chrom)
    sites.chr <- subset(get(sitestr), subset=chr==chrom)
    if(nrow(sites.chr)>0) {
      site.seq.chr <- get.sequence.chr(sites.chr, refgenome[[chrom]],ext=0)
    } else {
      site.seq.chr <- NULL
    }
    cat("\tDONE!\n")
    site.seq.chr
  }    
  cat("query reference genome: \n")
  sites.seqs <- sapply(seq(along=chrstr), process.site.sequence)
  names(sites.seqs) <- chrstr
  cat("calculate alphabet frequencies: \n")
  j <- 0
  assign(paste(sitestr, "alphabetFreq",sep="."),
         sapply(sites.seqs, function(x) {
           j <<- j+1;
           if(is.null(x)) NULL else alphabetFrequency(x)
         }
                )
         )
  cat(" Done on",j,"sites.","\n",sep=" ")  
}

ss.set <- c("sim", "notsigs", "P011", "P001", "ext", "ext.P001.uni", "ext.P001.bi")
for (ss in ss.set) {
  cat(ss, "--", "\n")
  assign(paste("gc",ss,sep="."),
         sapply(get(paste("sites",ss,"alphabetFreq",sep=".")),
                function(x) {
                  if(is.null(x)) cat("NULL","\n")
                  else {
                    cat(nrow(x),"\n")
                    apply(x,1, function(y) {
                      tol <- sum(y[c('A','C','G','T')]);
                      gc <- sum(y[c('C','G')]);
                      pct <- if(tol>0) gc/tol else NA
                    })}}))
}

pdf(file="../gc.dist.pdf")
ss.colors <- c("gray","black","red","blue","green","brown","purple")
for(ss.i in seq(ss.set)) {
  if(ss.i==1) {
    
    plot(density(na.omit(unlist(get(paste("gc",ss.set[ss.i],sep="."))))),
         col=ss.colors[ss.i],
         ylim=c(0,6.5),
         xlim=c(0,1),
         main="GC content")
  }
  else lines(density(na.omit(unlist(get(paste("gc",ss.set[ss.i],sep="."))))),
             col=ss.colors[ss.i])
}
dev.off()

gc.testsites <- NULL
for (ss in ss.set[3]) {
  cat(ss, "--", "\n")
  gc.testsites <-  c(gc.testsites, unlist(get(paste("gc",ss,sep="."))))
}

gc.breaks <- seq(0,1,0.05)

gc.hist.test <- hist(unlist(gc.testsites), plot=F, breaks=gc.breaks)
gc.hist.null <- hist(unlist(gc.notsigs), plot=F,breaks=gc.breaks)
            

gc.prob.test <- gc.hist.test$counts/sum(gc.hist.test$counts)
gc.prob.null <- gc.hist.null$counts/sum(gc.hist.null$counts)

gc.sample.prob.raw <- gc.prob.test/gc.prob.null
gc.sample.prob.med <- gc.sample.prob.raw<median(gc.sample.prob.raw[!is.na(gc.sample.prob.raw)])

###.. Pool the sites bound by RNF96(TRIM28), CNOT3, and ZFX
##--------------------------------------
prot1 <- "RNF96"
prot2 <- "CNOT3"
prot3 <- "ZFX"

bound.sites <- NULL
for (ifimo in paste("fimo",c("uni","bi","ext","P001","P011"),sep=".")) {
  cat(ifimo, "\n")
  bound.sites <- rbind(bound.sites,
                       sitesBound(c(prot1, prot2, prot3),
                                  get(ifimo),
                                  motif.meta,
                                  0.05))
}
bound.sites <- unique(bound.sites)
