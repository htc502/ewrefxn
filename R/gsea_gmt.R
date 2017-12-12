#' read in gmt gene set file from gsea
#'
#' This fxn accepts a gmt file and transform it to a list 
#' 
#' @param x gmt file name
#'
#' @return a list
#' @author guangchun
#' @keywords read.gsea.gmt
#' @examples
#' 
#' @rdname read.gsea.gmt
#' @export

read.gsea.gmt = function(f) {
  if(!file.exists(f)) stop('error finding gmt file:',f,'\n')
  tmp = readLines(f)
  tmp1 = strsplit(tmp, split = '\t') 
  tmp2 = lapply(tmp1, function(e) {
    Url = e[2]
    e=e[-2] ##remove url
    nGene = length(e)
    Name = e[1]
    Genes = e[2:nGene]
    list(Name=Name,
         Genes=Genes,
         url=Url,
         nGene = nGene)
  })
  tmp2
}

gsea.gmt.getName = function(plist) {
  unlist(lapply(plist,function(e) e$Name)) 
}
gsea.gmt.getGenes = function(plist) {
  lapply(plist,function(e) e$Genes) -> tmp
  names(tmp) = gsea.gmt.getName(plist)
  tmp
}
gsea.gmt.getnGene = function(plist) {
  unlist(lapply(plist,function(e) e$nGene)) -> tmp
  names(tmp) = gsea.gmt.getName(plist)
  tmp
}
gsea.gmt.geturl = function(plist) {
  unlist(lapply(plist,function(e) e$url)) 
}


#' re-organize the list returned by read.gsea.gmt
#'
#' This fxn accepts a list from read.gsea.gmt and reorganize it to a list with four elements: Names, Genes, nGene, urls
#' 
#' @param x pathway list
#'
#' @return a list
#' @author guangchun
#' @keywords gsea.gmt.reOrganize
#' @examples
#' 
#' @rdname gsea.gmt.reOrganize
#' @export
gsea.gmt.reOrganize = function(plist) {
  list(Names=gsea.gmt.getName(plist),
       Genes=gsea.gmt.getGenes(plist),
       nGene=gsea.gmt.getnGene(plist),
       urls=gsea.gmt.geturl(plist))
}

#' output with each pathway as a column
#'
#' @param x pathway list
#' @return a matrix
#' @author ghan
#' @export
gsea.gmt.outputMatrix = function(plist) {
  tmplist = gsea.gmt.getGenes(plist)
  names(tmplist) = paste0(gsea.gmt.getName(plist),'(',gsea.gmt.getnGene(plist),')')
  ewrefxn::list2matrix(tmplist)
}
