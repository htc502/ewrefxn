#' combine htseq count result into a matrix
#'
#' @param x htseq count result root dir
#' @param sampleID optional sample names if provided, otherwise will use the subdirectory name
#'
#' @return a data frame
#' @author guangchun
#' @keywords combineHTseq
#' @examples
#' 
#' @rdname combineHTseq
#' @export

combineHTseq <-function(x, sampleID=NA) {
  root_dir = x
  sub_dirs = list.dirs(root_dir, recursive=F)
  if(is.na(sampleID)) sampleID = basename(sub_dirs)
  readscountFile = paste0(basename(sub_dirs),'.readsCount')
  tmp = read.delim(file.path(sub_dirs[1],readscountFile[1]),header=F)
  tmp_stat = tmp[(nrow(tmp)-4):nrow(tmp),]
  tmp = tmp[1:(nrow(tmp)-5),]
  GeneName = unlist(tmp[,1])
  M = matrix(NA,nrow=length(GeneName),ncol=length(sub_dirs));rm(tmp)
  stat = matrix(NA,nrow=5,ncol=length(sub_dirs))
  for(i in seq_along(sub_dirs)) {
    tmp = read.delim(file.path(sub_dirs[i],readscountFile[i]),header=F)
    stati = ( tmp[(nrow(tmp)-4):nrow(tmp),2])
    stat[,i]  =stati
    tmp = tmp[1:(nrow(tmp)-5),]
    M[,i] = (tmp[,2])
  }
  rownames(M) = GeneName; colnames(M) = sampleID
  rownames(stat) = c('__no_feature',
                     '__ambiguous',
                     '__to_low_aQual',
                     '__not_aligned',
                     '__alignment_not_unique')
  colnames(stat) = sampleID
  return(list(M=M,
              stat=stat))
}

