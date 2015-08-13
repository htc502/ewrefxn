#' @export
ebi.fetch <- function(...) {
    if(!require0(ArrayExpress))
        stop("error loading arrayexpress package\n")
    ##arrayexpress file name list returned
    ae.flist <- getAE(...)
    ae.flist
}
#' @export
ebi.agilent44k.prep <- function( ae.flist ) {
    if(!("sdrf" %in% names(ae.flist)) )
        stop("sample&data relationship file(sdrf file) not found\n")
    SDRF <- read.delim(ae.flist$sdrf, stringsAsFactors = F,
                       check.names = F)
    cat("array data files we will process are shown bellow, to remove/add samples, manually edit the sdrf file...\n")
    print(SDRF[ ,"Array Data File" ])
    answer <- dual.choice("sure we go ahead?",value1 = "y",value2="n")
    if(answer == "n")
        return(NULL)
    prep.obj <- agilent44k.prep(array.fnames = SDRF[ ,"Array Data File" ])
    prep.obj
}
#' @export
ebi.agilent80k.prep <- function(ae.flist) {
    res <- ebi.agilent44k.prep(ae.flist)
    res
}
#' @export
ebi.affy.st.prep <- function(ae.flist) {
    if(!("sdrf" %in% names(ae.flist)) )
        stop("sample&data relationship file(sdrf file) not found\n")
    SDRF <- read.delim(ae.flist$sdrf, stringsAsFactors = F,
                       check.names = F)
    cat("array data files we will process are shown bellow, to remove/add samples, manually edit the sdrf file...\n")
    print(SDRF[ ,"Array Data File" ])
    answer <- dual.choice("sure we go ahead?",value1 = "y",value2="n")
    if(answer == "n")
        return(NULL)
    prep.obj <- affy.st.prep(array.fnames = SDRF[ ,"Array Data File" ])
    prep.obj

}
#' @export
ebi.affy.ivt.prep <- function(ae.flist) {
    if(!("sdrf" %in% names(ae.flist)) )
        stop("sample&data relationship file(sdrf file) not found\n")
    SDRF <- read.delim(ae.flist$sdrf, stringsAsFactors = F,
                       check.names = F)
    cat("array data files we will process are shown bellow, to remove/add samples, manually edit the sdrf file...\n")
    print(SDRF[ ,"Array Data File" ])
    answer <- dual.choice("sure we go ahead?",value1 = "y",value2="n")
    if(answer == "n")
        return(NULL)
    prep.obj <- affy.ivt.prep(array.fnames = SDRF[ ,"Array Data File" ])
    prep.obj
}
