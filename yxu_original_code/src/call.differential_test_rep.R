## Filename: call.differential_test.R
## By: Yaomin Xu
## Date: Sep 22, 2010
## Notes: Call DiffSeq differetial testing
## =============================================================================
###. Setup parameters
##____________________________________________________________________
srcPath <- "src"
srcExpConfigPath <- ".."
source(paste(srcPath, "DiffSeq.access.R", sep="/"))
source(paste(srcExpConfigPath, "experiment.config.R", sep="/"))

chr <- get.str.chr()
cat(chr, "\n")
ncore <- get.n.cpus()

get.ws(ws="mixPoisson",
       what=c("pmeans",
         "pmeans.norm",
         "pmeans.bgcorrect.norm",
         "pmeans.bgnorm.norm"))

pmeans <- list(mean=pmeans$mean[,sample.select],
               var=pmeans$var[,sample.select])
pmeans.norm <- pmeans.norm[,sample.select]
pmeans.bgcorrect.norm <- pmeans.bgcorrect.norm[,sample.select]
pmeans.bgnorm.norm <- pmeans.bgnorm.norm[,sample.select]
new.names <- paste("s", seq(sum(sample.select)), sep="")
names(pmeans$mean) <- names(pmeans$var) <- new.names
names(pmeans.norm) <- names(pmeans.bgcorrect.norm) <- names(pmeans.bgnorm.norm) <- new.names

##sample.label <- sample.label[sample.select]


if(bg.process == "correct") {
  pmeans.norm <- pmeans.bgcorrect.norm
} else if(bg.process == "norm") pmeans.norm <- pmeans.bgnorm.norm
get.ws("pattRecog",
       c("sites.js.ex", "events.wins"))


###. Computation
##____________________________________________________________________

##-- Parallel start
sfInit(parallel=TRUE, cpus=ncore, nostart=if.parallel(ncore))
##sfLibrary(nlme)
##sfLibrary(MASS)
sfSource(paste(srcPath, "DiffSeq.access.R", sep="/"))
sfSource(paste(srcExpConfigPath, "experiment.config.R", sep="/"))
##sfExport(list=c("sample.label"))
##--
patts <- unique(sites.js.ex$pattern)
events.wins.matched <- events.wins[,match.patt.pattmx(patts, events.wins)]
events.ctr.matched <- events.ctr[match.patt.pattmx(patts, events.wins)]
patt.n <- length(patts)
wins.test <- vector("list", length=patt.n)
for (i in seq(patt.n)) {
  wins.test[[i]] <- list()
  for(j in seq(length(events.ctr.matched[[i]]))) {
    cat(as.character(patts[i]), "\t", j, "\n")
    wins.test[[i]][[j]] <- test.effect.rep(subset(sites.js.ex, subset=pattern==patts[i]),
                                           pmeans.norm, pmeans$var, contr= events.ctr.matched[[i]][[j]],
                                           n.core = ncore,
                                           group.label=sample.label,
                                           rep.label=subject.label)
    attr(wins.test[[i]][[j]],'contrast') <- events.ctr.matched[[i]][[j]]
    attr(wins.test[[i]][[j]],'pattern') <- events.wins.matched[,i]
  }
}

wins.test.report <- vector("list", length=patt.n)
for (i in seq(patt.n)) {
  wins.test.report[[i]] <- list()
  for(j in seq(length(events.ctr.matched[[i]]))) {
    cat(as.character(patts[i]), "\t", j, "\n")
    wins.test.report[[i]][[j]] <- combine.report(subset(sites.js.ex, subset=pattern==patts[i]),
                                                wins.test[[i]][[j]])
    attr(wins.test.report[[i]][[j]],'contrast') <- events.ctr.matched[[i]][[j]]
    attr(wins.test.report[[i]][[j]],'pattern') <- events.wins.matched[,i]
  }
}

##-- Parallel stop
sfStop()
##--
###. Save workspace
##____________________________________________________________________
save.ws(ws="diffTest",
        what=c("chr",
          "pmeans", "pmeans.norm",
          "wins.test", "wins.test.report",
          "events.ctr", "events.wins"
          ),
        chr=chr
        )
