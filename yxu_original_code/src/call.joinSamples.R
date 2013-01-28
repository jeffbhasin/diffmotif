## Filename: joinSamples.R
## By: Yaomin Xu
## Date: 
## Notes: Join samples by extending with zeros
## =======================================================================================
###. Setup parameters
##____________________________________________________________________
srcPath <- "src"
srcExpConfigPath <- ".."
source(paste(srcPath, "DiffSeq.access.R", sep="/"))
source(paste(srcExpConfigPath, "experiment.config.R", sep="/"))
winsize <- get.winsize.input()
print(winsize)

###. Computation
##____________________________________________________________________
save(winsize, file=file.path(get.wsoption.path("joinSamples"),"winsize.RData"))
merge.samples(srcExpConfigPath, mappedFiles)
