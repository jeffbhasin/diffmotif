## Filename:
## By: Yaomin Xu
## Date: Oct 21, 2010
## Notes: Setup the experimental grouping, patterns of interest, and contrasts of interest.
## =============================================================================
mappedFiles <- c('InputData/624N.bow.50bp.txt',
                 'InputData/715N.bow.50bp.txt',
                 'InputData/739N.bow.50bp.txt',
                 'InputData/624T.bow.50bp.txt',
                 'InputData/715T.bow.50bp.txt',
                 'InputData/739T.bow.50bp.txt',
                 'InputData/541T.bow.50bp.txt',
                 'InputData/754T.bow.50bp.txt',
                 'InputData/874T.bow.50bp.txt')

initfilter.TF <- FALSE
initfilter.cutoff <- 1

fragsizefilter.TF <- F
fragsizefilter.cutoff <- 120
###... Samples
sample.label.full <- c("N","N","N","T1","T1","T1","T2","T2","T2")
subject.label.full <- c(624,715,739,624,715,739,541,754,874)
sample.select <- c(T,T,T,T,T,T,T,T,T)
sample.label <- sample.label.full[sample.select]
subject.label <- subject.label.full[sample.select]
###... Patterns of interest
## Probability mapping patterns
events.res <- event.allPossiblePatterns(sample.label, sample.select)

## Patterns to scan and identify
events.wins <- event.patternsToScan(events.res,
                                    exclude=NULL,
                                    noZeroPatt=T)
## Event cutoff method
## Beta-based
events.cutoff.method <- "BF"
events.cut <- -1

##bf.which.score <- "qvalue"
bf.cutoff <- NULL

###... Diff Testing
diffTest.cutoff.method <- "qvalue"
diffTest.cutoff <- 0.05
bg.process <- "norm"

events.ctr <- construct.contrasts(events.wins)
contrast.objs <- construct.contrast.objs(events.wins, events.ctr)

###... Other parameters
join.gap <- 150
shear.gap <- 150
region.gap <- 500

search.sites.cut <- rbind(matrix(rep(c(5e4,2e4), each=ncol(events.wins)-1), ncol=2),
                          c(2500, 2000))
search.sites.cut.2 <- c(5000,2000)
