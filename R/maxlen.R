#' inner fxn used in list2matrix
#'
#'return the length of the vector element with max length
#' @param x a list of vectors
#' @keywords length
#' @export
maxlen<-function(x) {
    elemntLen<-vapply(x,length,FUN.VALUE=vector(mode="integer",length=1))
    return(max(elemntLen))
}

