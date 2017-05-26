
enrich.test = (PathGenes, interestGenes, totalGenes) {
################################################
###           | inPath | not inPath |
###-------------------------------------------
###   Predict |  a     |   b        | nPredict
###notPredict |  c     |   d        |
###-------------------------------------------
###           |  nPath |            | nAll
################################################
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
                          'phyper2'=p2))
  res
}

