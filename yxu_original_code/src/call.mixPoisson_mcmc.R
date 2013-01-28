## Filename: call.mixPoisson_mcmc.R
## By: Yaomin Xu
## Date: Sep 22, 2010
## Notes: Call mixPoisson mcmc
## =============================================================================
###. Setup parameters
##____________________________________________________________________
srcPath <- "src"
srcExpConfigPath <- ".."
source(paste(srcPath, "DiffSeq.access.R", sep="/"))
source(paste(srcExpConfigPath, "experiment.config.R", sep="/"))

chr <- get.str.chr()
cat(chr, "\n")
ncore.in <- get.n.cpus()

count.table <- get.mappedCountTab()

if(initfilter.TF) {
  (init.filter.cutoff <- get.initfilter.cutoff(count.table, 0.5))
} else {
  (init.filter.cutoff <- initfilter.cutoff)
}
get.reads.fromJoinedSamples(get.wsoption.path("joinSamples"), threshold=init.filter.cutoff)

###. Computation
##____________________________________________________________________

###.. MCMC
##--------------------------------------
controls <- list(msize=6000, burnin=1000, thin=1);
priors.e <- list(d=c(0.3, 0.3, 0.3),
                 lambda.prior=list(a=1, b=1/20),
                 signal.prior=list(q=0.1,s=0.1))

ncore <- min(ncore.in, length(s.cols))
sfInit(parallel=TRUE, cpus=ncore, type="SOCK", nostart=if.parallel(ncore))
sfExportAll()
sfLibrary(MCMCpack)

para.mixPois_eWin <- function(i.smp, priors, controls) {
  cat("start sample", i.smp, "\n")
  y1 <- reads[, i.smp+1]
  
  ## Run 2 (w/o initials)
  n <- length(y1)
  p0<- c(sum(y1==0)/n, sum(y1>0&y1<=4)/n, sum(y1>4)/n)
  lm0 <- c(mean(y1[which(y1>0&y1<=4)]), mean(y1[which(y1>4)]))
  signal0 <- y1*(1+1/25)
  
  initials.e <- list(p=p0, lambda=lm0[1], signal=signal0, r=1/25)
  res.i <- mixPois_eWin(y=y1,
                        initials=initials.e,
                        para.priors=priors,
                        controls=controls)
  return(res.i)
}

res <- sfLapply(seq(along=s.cols),
                para.mixPois_eWin,
                priors=priors.e,
                controls=controls)
res <- res[!unlist(lapply(res, is.null))]

sfStop()

###.. Post MCMC
##--------------------------------------
compute.pmeans()
###. Save workspace
##____________________________________________________________________
save.ws(ws="mixPoisson",
        what=c("chr", "res",
          "reads","reads.poi","reads.names",
          "pmeans","pmeans.norm","pmeans.bgnorm.norm","pmeans.bgcorrect.norm"),
        chr=chr
        )
