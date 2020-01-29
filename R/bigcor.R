#' correlation on big matrix
#'
#' This fxn calculate expression on big matrix
#' 
#' @param x same as x in cor
#' @param y same as y in cor
#' @param nblocks number of blocks
#' @author guangchun
#' @keywords cor
#' @examples
#' 
#' @rdname bigcorXY
#' @export


bigcorXY <- function(x,y,nblocks = NULL,verbose = TRUE, ...)
{
    NCOL <- ncol(x);NCOL2 = ncol(y)

    if(NCOL > 1e10 & is.null(nblocks)) stop('specify nblocks when ncol > 1e10')
    if(is.null(nblocks)) nblocks = max(ewrefxn::findfactors(NCOL))
    
    ## test if ncol(x) %% nblocks gives remainder 0
    if (NCOL %% nblocks != 0) stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")
    
    ## preallocate square matrix of dimension
    ## ncol(x) in 'ff' single format
    corMAT <- matrix(NA,nrow = NCOL, ncol = NCOL2)
    
    ## split column numbers into 'nblocks' groups
    SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
    
    ## iterate through each block combination, calculate correlation matrix
    ## between blocks and store them in the preallocated matrix on both
    ## symmetric sides of the diagonal
    for (i in seq_along(SPLIT)) {
        SPLITi = SPLIT[[i]]
        G1 <- SPLITi
        if (verbose) cat("Block", i, "\n")
        COR <- cor(x[, G1], y, ...)
        corMAT[G1, ] <- COR
        COR <- NULL
    }
    
    gc()
    return(corMAT)
}
