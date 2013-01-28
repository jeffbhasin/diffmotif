## Filename: call.pattRecog_win.R
## By: Yaomin Xu
## Date: Sep 22, 2010
## Notes: call for site pattern recognition.
## =============================================================================
srcPath <- "src"
srcExpConfigPath <- ".."
source(paste(srcPath, "DiffSeq.access.R", sep="/"))
source(paste(srcExpConfigPath, "experiment.config.R", sep="/"))

chr <- get.str.chr()
cat(chr, "\n")
ncore <- get.n.cpus()

##get.ws(ws="mixPoisson", what =c("res", "reads.poi", "pmeans", "pmeans.norm"))
##get.ws(ws="mixPoisson", what =c("res", "reads.poi"))
get.ws(ws="mixPoisson",
       what=c("res",
         "reads.poi",
         "pmeans",
         "pmeans.norm",
         "pmeans.bgcorrect.norm",
         "pmeans.bgnorm.norm"))
if(bg.process == "correct") {
  pmeans.norm <- pmeans.bgcorrect.norm
} else if(bg.process == "norm") pmeans.norm <- pmeans.bgnorm.norm

##win.size <- get.wsoption.para("winsize")
winsize <- get.winsize()

res <- res[sample.select]
##sample.label <- sample.label[sample.select]

##sample.label <- c("N","N","N","T1","T1","T1","T2","T2","T2")

## pattern recognition
##events.res <- cbind(c(0,1,1), c(0,0,1), c(1,0,0), c(1,1,0),
##                    c(0,1,0), c(1,0,1), c(1,1,1), c(0,0,0))
## differential testing pattern
## events.wins <- cbind(events.res[,-which(colSums(events.res) == 0)],
##                      c(NA,NA,1), c(NA,1,NA), c(1,NA,NA))
##events.wins <- events.res[,-which(colSums(events.res) == 0)]

e.res <- events.score(res, events.res, reads.poi, sample.label)

e.cut <- choose.events.cutoff(res,
                              e.res,
                              events.cut,
                              events.cutoff.method)

wins.sel.dup <- win.select(e.res, e.cut)
wins.sel <- identify.win.overlap(wins.sel.dup)


sfInit(parallel=TRUE, cpus=ncore, type="SOCK", nostart=if.parallel(ncore))
sfSource(paste(srcPath, "DiffSeq.access.R", sep="/"))

sites <- search.sites.v2(wins.sel, events.wins, chr, winsize,
                         search.sites.cut, cut.1st=search.sites.cut.2)

sfStop()

###.. Joining, shearing, and cleaning
##--------------------------------------

ranged.sites <- rangedSites(sites, winsize)
gaps.site <- site.gaps(ranged.sites)

sites.disjoined <-disjoin.sites(disjoin.sites(ranged.sites$sites.site,
                                              gaps.site$within,
                                              type="two.ends"),
                                gaps.site$within[width(gaps.site$within)>(shear.gap-winsize)])

### Clean the sites based on their pattern consistancy
e.TF <- ranged.e.TF(e.res, winsize, cutoff=-1)

sfInit(parallel=TRUE, cpus=ncore, type="SOCK", nostart=if.parallel(ncore))
sfLibrary(IRanges)
sfLibrary(plyr)
sfExport(list=c("sites.disjoined","e.TF","patt.summary","pattSummary.2dataframe"))
sites.cleaned <- IRangesList(sfLapply(names(sites.disjoined),
                                      clean.sitePatt,
                                      sites=sites.disjoined,
                                      e.TF=e.TF))
sfStop()

names(sites.cleaned) <- names(sites.disjoined)

wins.full <- ranges(e.TF)[[1]]
sites.js <- format.sites4diffTest(wins.full, sites.cleaned, chr)

###... Post process
if(!is.null(bf.cutoff)) {
  bf.filtered <- bf.filter(lapply(pattern.events(e.res, events.wins),c),
                           data.frame(sites.js, pattern.win=match.pattern.win(wins.sel, sites.js)),
                           winsize,
                           bf.cutoff,
                           reads.poi,
                           pmeans.norm,
                           events.wins,
                           sample.select,
                           which.score=bf.which.score)
  sites.js.ex <- bf.filtered$sites.js.ex
} else {
  bf.filtered <- NULL
  sites.js.ex <- data.frame(sites.js, pattern.win=match.pattern.win(wins.sel, sites.js))
}

sfStop()

##attr(sites.j.ex, 'patterns') <- events.wins
attr(sites.js.ex, 'patterns') <- events.wins

###. Save workspace
##____________________________________________________________________
save.ws(ws="pattRecog",
        what=c("chr","sample.label","winsize",
          "sites", "sites.js.ex", "wins.sel", "e.res", "e.cut", "events.res", "events.wins",
          ##"pmeans", "pmeans.norm", "bf.filtered")
          "bf.filtered"),
        chr=chr
        )
