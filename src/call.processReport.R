## Filename: call.processReport.R
## By: Yaomin Xu
## Date: Sep 22, 2010
## Notes: Call report process (this is new, There is an old version)
## =============================================================================

###. Setup parameters
##____________________________________________________________________
srcPath <- "src"
srcExpConfigPath <- ".."
source(paste(srcPath, "DiffSeq.access.R", sep="/"))
source(paste(srcExpConfigPath, "experiment.config.R", sep="/"))
winsize <- get.winsize()
###. Computation
##____________________________________________________________________

if(fragsizefilter.TF) {
  test.report.flat <- subset(get.difftest.flatReport(winsize=winsize),
                             subset=end-start>=fragsizefilter.cutoff)
} else test.report.flat <- get.difftest.flatReport(winsize=winsize)

report.contrast <- make.reportContrast(test.report.flat,
                                       pvalue.var="p.value.1",
                                       effect.var="value.1",
                                       type=2)

test.report.sigs <- select.sig.sites(test.report.flat,
                                     report.contrast,
                                     method=if(exists('diffTest.cutoff.method')) diffTest.cutoff.method else "qvalue",
                                     cutoff.p.sel=if(exists('diffTest.cutoff')) diffTest.cutoff else 0.01,
                                     winsize=winsize,
                                     cutoff.by.winsize=T)
###. save workspace
##____________________________________________________________________
save.ws(ws="report",
        what=c("test.report.flat", "test.report.sigs"))

