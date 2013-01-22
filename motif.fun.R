## Filename:
## By: Yaomin Xu
## Date:
## Notes: Functions for motif.R
## =============================================================================
process.site.freq.summary <- function(fimo.out, sites, motif.meta, q.cutoff, round.digit=1) {
  browser()
  here.top.tb <- table(subset(fimo.out, subset=V7<q.cutoff, select="V1"))
  here.top.most <- here.top.tb[which(here.top.tb>0)]

  here.topmost.table <- cbind(motif.meta[names(here.top.most),],
                              freq=here.top.most)
  
  here.motif.xxx <- by(here.topmost.table, factor(here.topmost.table$cname), function(x) as.character(x$ac))
  here.motif.max <- sapply(here.motif.xxx,
                           function(x, tb) {
                             if(length(x)>1) x[which.max(subset(tb, subset=ac %in% x, select="freq", drop=T))]
                             else x
                           },
                           tb=here.topmost.table)
  here.motif.vec <- sapply(here.motif.xxx, paste, collapse=",", sep=",")
  
  here.freq <- sapply(here.motif.xxx,
                     function(x, fimotb) length(unique(subset(fimotb,
                                                              subset=V1%in%as.character(unlist(x)),
                                                              select='V2',
                                                              drop=T))),
                     fimotb=fimo.out)
  
  here.freq.cnt <- data.frame(count=here.freq,
                              pct=round(here.freq/nrow(sites)*100, round.digit),
                              motif=here.motif.vec,
                              motif.maxhit=here.motif.max)
  
  list(hits=here.topmost.table,
       site.freq=here.freq.cnt[order(here.freq.cnt$count, decreasing=T),])
}

read.motif.idname <- function(file) {
  rawstr <- grep("^MOTIF",
                 read.table(file, sep="\t", header=F,stringsAsFactors=F)[,1],
                 v=T)
  idstr <- sapply(rawstr,
                  function(x) unlist(strsplit(x, split=" "))[2])
  namestr <- sapply(rawstr,
                    function(x) paste(unlist(strsplit(x, split=" "))[-(1:2)],
                                      collapse=" "))
  output <- data.frame(id=idstr, name=namestr)
  row.names(output) <- idstr
  output
}

read.meme2list <- function(file) {
  ##browser()
  rawlines <- readLines(file)
  eof.line <- length(rawlines)
  m.lines <- grep("MOTIF ", rawlines)
  m.n <- length(m.lines)
  m.strs <- rawlines[m.lines]
  out.list <- NULL
  if(m.n>0) {
    o.header <- rawlines[1:(m.lines[1]-1)]
    idstr <- as.vector(sapply(m.strs,
                              function(x) unlist(strsplit(x, split=" "))[2]))
                       
    namestr <- as.vector(sapply(m.strs,
                                function(x) paste(unlist(strsplit(x, split=" "))[-(1:2)],
                                                  collapse=" ")))
    o.motifs <- vector("list", m.n)
    for(i in seq(m.n-1)) {
      o.motifs[[i]] <- list(id=idstr[i],
                            name=namestr[i],
                            text=rawlines[m.lines[i]:(m.lines[i+1]-1)])
    }
    o.motifs[[m.n]] <- list(id=idstr[m.n],
                            name=namestr[m.n],
                            text=rawlines[m.lines[m.n]:eof.line])
    ##o.motifs <- sapply(o.motifs, function(x) {names(x) <- c("id","name","text");x})
    out.list <- list(header=o.header,
                     motifs=o.motifs)
  }
  out.list
}

write.memelist2meme <- function(ml, file) {
  con <- file(file, "w")
  writeLines(ml$header, con)
  for (i in seq(along=ml$motifs)) {
    writeLines(ml$motifs[[i]]$text, con)
  }
  close(con)
}

combine.memelist <- function(ml1, ml2, header=1) {
  o.header <- if(header==1) ml1$header else ml2$header
  o.motifs <- c(ml1$motifs, ml2$motifs)
  list(header=o.header,
       motifs=o.motifs)
}

get.ids.memelist <- function(ml) {
  sapply(ml$motifs, function(x) x$id)
}

get.names.memelist <- function(ml) {
  sapply(ml$motifs, function(x) x$name)
}

select.memelist <- function(ml, ids) {
  ml.ids <- get.ids.memelist(ml)
  sel.ids <- match(ids, ml.ids)
  if(sum(is.na(sel.ids))>0) warnings(ids[is.na(sel.ids)], "does not exit")
  list(header=ml$header,
       motifs=ml$motifs[na.omit(sel.ids)])
}

read.transfacMatrix2list <- function(file) {
  ##browser()
  rawlines <- readLines(file)
  eof.line <- length(rawlines)
  m.lines <- grep("//", rawlines)
  m.n <- length(m.lines)
  ##m.strs <- rawlines[m.lines]
  if(m.n>0) {
    o.header <- rawlines[1:(m.lines[1]-1)]
    o.motifs <- vector("list", m.n-1)
    for(i in seq(m.n-1)) {
      i.text <- rawlines[(m.lines[i]+1):(m.lines[i+1]-1)]
      ac.str <- unlist(strsplit(grep("^AC ", i.text, v=T), split=" +"))[2]
      id.str <- unlist(strsplit(grep("^ID ", i.text, v=T), split=" +"))[2]
      na.str <- unlist(strsplit(grep("^NA ", i.text, v=T), split=" +"))[2]
      o.motifs[[i]] <- list(id=id.str,
                            name=na.str,
                            ac=ac.str,
                            text=i.text)
    }
    ## o.motifs[[m.n]] <- list(id=idstr[m.n],
    ##                         name=namestr[m.n],
    ##                         text=rawlines[m.lines[m.n]:eof.line])
    ##o.motifs <- sapply(o.motifs, function(x) {names(x) <- c("id","name","text");x})
    out.list <- list(header=o.header,
                     motifs=o.motifs)
  }
  out.list
}

get.info.trmatrix <- function(trmatrix) {
  ##browser()
  o.df <- data.frame(cbind( ac=unlist(sapply(trmatrix$motifs, function(x) x$ac)),
                            id=unlist(sapply(trmatrix$motifs, function(x) x$id)),
                            name=unlist(sapply(trmatrix$motifs, function(x) x$na))))
  row.names(o.df) <- o.df$ac
  o.df
}

select.trmatrix.list <- function(trmatrix, ids) {
  trm.ids <- get.info.trmatrix(trmatrix)$ac
  sel.ids <- match(ids, trm.ids)
  if(sum(is.na(sel.ids))>0) warnings(ids[is.na(sel.ids)], "does not exit")
  list(header=trmatrix$header,
       motifs=trmatrix$motifs[na.omit(sel.ids)])
}

       
## ACs.by.cname <- function(topmost.tb) {
##   by

find.ac.by.cname <- function(protein, motifmeta) {
  as.character(unlist(subset(motifmeta, subset=toupper(cname)==protein, select=ac)))
}

select.sitesBound <- function(protein,
                              fimo.readin,
                              motif.meta,
                              q.cutoff=0.01) {
  protein.ac <- find.ac.by.cname(protein, motif.meta)
  as.character(unlist(subset(fimo.readin, subset=V1%in%protein.ac&V7<q.cutoff, select="V2")))
}

format.sitesBound <- function(sitesBound.list) {
  .full <- unique(unlist(sitesBound.list))
  cbind(site=as.character(.full),
        as.data.frame(lapply(sitesBound.list, function(x) .full%in%x)))
}

sitesBound <- function(proteins, fimoout, motifmeta, q.cutoff=0.5) {
  
  format.sitesBound(lapply(sapply(proteins,
                                  select.sitesBound,
                                  fimo.readin = fimoout,
                                  motif.meta=motifmeta,
                                  q.cutoff = q.cutoff),
                           function(x) unique(x)))
}

    
    
