#' Check if two character vectors has the same element order
#'
#' Usually, we have to check taht the col names(sample name) of
#' an expression matrix is in the same order with the sample names
#' listed in the metainfor data matrix, this function is intended
#' for this purpose
#'
#' @param x vector, must be a character vector with no duplicated elements
#' @param y another vector which is used to be compared against, must be a character vector with no duplicated elements too
#'
#' @return a list with two components: status and position
#' @author guangchun hanguangchunATgmail.com
#' @keywords samOrder
#' @examples
#'
#' x<-sample(letters,6)
#' y<-sample(x, 6)
#'
#' samOrder(x,x)
#' samOrder(x,y)
#'
#' @rdname sameOrder
#' @export

sameOrder <- function(x,y) {
    if(!is.character(x) | !is.character(y) )
        stop('non character inputs found')
    if(anyDuplicated(x) != 0 | anyDuplicated(y) != 0)
        stop('duplications found in the input vectors')
    if(any(sort(x) != sort(y)))
        stop('two input vectors have different elements')
    match(x,y) ->pos
    len <- length(x)
    res <- ifelse(all(pos == seq(1,len,by=1)),T,F)
    return(list(res=res, position=pos))
}
