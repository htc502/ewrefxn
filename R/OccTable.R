#' convert a named list to an Ocurrence table(emerged 1, absent 0)
#'
#' this function convert a names list  with length n to a matrix with dim x*n
#' x is the total number of unique elements stored in the list
#' @param lst a deg list
#' @return a matrix with 0s and 1s
#' @export
list2OccTable <- function(lst) {
    if(is.null(attributes(lst)$names))
        stop('input list should have names')
    if(length(lst) <= 1)
        stop('length of input list <= 1')
    elems <- unique(unlist(lst))
    if(length(elems) <= 1)
        stop('number of uniqued elements in the input list <= 1')
    res <- matrix(NA, nrow=length(elems),ncol=length(lst))
    for(i in 1:length(elems)) {
        elem <- elems[i]
        res[ i, ] <- unlist(lapply(lst,
                                   function(x, elem) {
                                       ifelse( elem %in% x, 1,0) },elem= elem))
    }
    rownames(res) <- elems;colnames(res) <- names(lst)
    res
}

#' convert a named list to an Ocurrence table(emerged 1, absent 0), this one is speedy
#'
#' this function convert a names list  with length n to a matrix with dim x*n
#' x is the total number of unique elements stored in the list
#' @param lst a deg list
#' @return a matrix with 0s and 1s
#' @export
list2OccTable1 <- function(lst) {
    if (is.null(attributes(lst)$names))
        stop("input list should have names")
    if (length(lst) <= 1)
        stop("length of input list <= 1")
    elems <- unique(unlist(lst))
    if (length(elems) <= 1)
        stop("number of uniqued elements in the input list <= 1")
    res <- lapply(lst, function(x, elems) {
		   tmp <- rep(0,length=length(elems))
		   tmp[ elems %in% x ] <- 1
		   tmp
	      },elems=elems)
    res1 <- as.matrix(as.data.frame(res))
    rownames(res1) <- elems
    colnames(res1) <- names(lst)
    res1
}
