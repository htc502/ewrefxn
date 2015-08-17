#' fetch gene expression data from EBI arrayexpress
#'
#' this is a wrapper function which is used to download data from EBI arrayexpress
#' @param accession Accession number of the data: E-GEOD-6039 etc
#' @param path the dir where downloaded files will be stored in
#' @param type raw/processed/full for downlaod the raw/processed/both data
#' @param extract default to T, otherwise the downloaded data will not be extracted
#' @param local local files to be read instead of downloading from the server
#' @param sourcedir if local is TRUE, files will be read from this dir
#' @return return
#' @export
ebi.fetch <- function(accession, path = getwd(), type = 'full',
                      extract = T, local = F, sourcedir = path, ...) {
    if(!require0(ArrayExpress))
        stop("error loading arrayexpress package\n")
    ##arrayexpress file name list returned
    ae.flist <- getAE(accession, path = getwd(), type = 'full',
                      extract = T, local = F, sourcedir = path, ...)
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
