## Filename: DiffSeq.access.R
## By: Yaomin Xu
## Date: Sep 11, 2010
## Notes: To access the DiffSeq src
## =============================================================================
require(RCurl)
require(MCMCpack)
require(snowfall)
require(nlme)
require(gmodels)
require(MASS)
require(plyr)
require(preprocessCore)
require(qvalue)
require(inline)
require(Rcpp)
require(IRanges)


###... DPT methods
##src <- getURL("https://files.me.com/yaomin/2z6hiy")
## src <- getURLContent("http://dl.dropbox.com/u/10421230/dpt/current/DiffSeq.methods.R")
## eval(parse(text = src))

###... Other includes
##get.sharedSrc("methods.BFCut.R")
##get.src("methods.BFCut.R")
##get.sharedSrc("methods.seqs.R")
##get.src("methods.seqs.R")

diffseq.options.ws <- list(ws.root = "DPT_ws",
                           mixpos.path = "MixPoisson",
                           pattrecog.path = "PattRecog",
                           difftest.path = "DiffTest",
                           joinsample.path = "JoinSamples",
                           log.path = "logs",
                           report.path = "Reports",
                           src.path = "src")

get.src <- function(string, base, dpt.option){
  if(missing(dpt.option)) {
    dpt.option <- list(ws.root = "DPT_ws",
                       mixpos.path = "MixPoisson",
                       pattrecog.path = "PattRecog",
                       difftest.path = "DiffTest",
                       joinsample.path = "JoinSamples",
                       log.path = "logs",
                       report.path = "Reports",
                       src.path = "src")
  }
  attach(dpt.option)
  
  if(!missing(base)) {
    source(file.path(base, string))
  } else {
    if(file.exists(file.path(src.path, string))){
      source(file.path(src.path, string))
    } else if(file.exists(file.path(ws.root,src.path, string))) {
      source(file.path(ws.root, src.path, string))
    } else {stop("Can't allocate the src file")}
  }
}

get.src(string="DiffSeq.methods.R",dpt.option=diffseq.options.ws)

###...functions

