#' fisher's exact test and hypergeometric test for gene set
#'
#' This fxn calculates statistics using fisher's exact test and hypergeometric test given a set of genes and gene set
#' 
#' @param PathGenes a vector of genes belonging to a pathway
#' @param interestGenes a vector of genes of interest
#' @param totalGenes a vector of background genes
#' @param comet_exact_test logical for comet exact test of mutual exclusive(requires package cometExactTest)
#' @param mutmatplot plot mutation mat for cometExactTest result
#' @return a list
#' @author guangchun
#' @keywords enrich.test
#' @examples
#' 
#' @rdname enrich.test
#' @export


enrich.test = function(PathGenes, interestGenes, totalGenes,comet_exact_test=F,mutmatplot=T) {
################################################
###           | inPath | not inPath |
###-------------------------------------------
###   Predict |  a     |   b        | nPredict
###notPredict |  c     |   d        |
###-------------------------------------------
###           |  nPath |            | nAll
################################################
  ##remove genes not included in reference gene set
  if(comet_exact_test) {
    if(!ewrefxn::require0('cometExactTest')) stop('error installing/loading cometExactTest,stop')
  }
  PathGenes = intersect(PathGenes,totalGenes)
  interestGenes = intersect(interestGenes, totalGenes)
  
  nPath <- length(PathGenes)
  nPredict <- length(interestGenes)
  a <- length(intersect(PathGenes,interestGenes))
  b <- nPredict - a
  c <- nPath - a
  nAll = length(totalGenes)
  d <- nAll - (a + b + c)
  m <- matrix(c(a,b,c,d), ncol=2,byrow=T)
  colnames(m) <- c('inPath','notInPath');rownames(m) <- c('Predict','notPredict')
  p <- fisher.test(m,alternative='greater')$p.value

  ##using phyper for hypergeometric test
  p1 <- phyper(a-1,nPath, nAll-nPath,nPredict,lower.tail=F)
  p2 <-1 - phyper(a-1,nPath, nAll-nPath,nPredict)
  tmp <- c(a,b,c,d,p,p1,p2)

  
  res = list(mat=m,p=list('fisher.exact'=p,
                          'phyper1'=p1,
                          'phyper2'=p2),
            overlapped=intersect(PathGenes, interestGenes))
  ##mutual exclusive test
  if(comet_exact_test) {
  p.comet = comet_exact_test(c(a,b,c,d),mutmatplot = mutmatplot)
  res$cometTest = p.comet
  }
  res
}

