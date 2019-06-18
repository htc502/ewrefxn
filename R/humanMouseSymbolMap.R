#' map gene symbols between mouse and human
#'
#' the function convert the symbol of mouse to the corresponding one of human
#' based on the record of homologen database in NCBI
#' @param mousesymbols
#'
#' @return a vector of corresponding human gene symbol
#' @author guangchun
#' @keywords mouse2human
#' @examples
#' m1 <- c('App', 'Apoe')
#' mouse2human(m1)
#' @rdname mouse2human
#' @export
mouse2human <- function(msymb) {
    f = system.file("extdata", "HMD_HumanPhenotype.dat", package = "ewrefxn",
                    mustWork = T)
    humanMouseMap = read.delim(f,row.names = NULL,
                               comment.char = '#')
    if(any(duplicated(msymb))) stop('duplicated values found in input vector')
    df = humanMouseMap
    humanMouseMap1 = aggregate(df$human, list(mouse = df$mouse), paste, collapse = '|')
    pos <- match(msymb, humanMouseMap$mouse)
    humanMouseMap1[ pos, ]$x
}


#' map gene symbols between mouse and human
#'
#' the function convert the symbol of human to the corresponding one of mouse
#' based on the record of homologen database in NCBI
#' @param humansymbols
#'
#' @return a vector of corresponding mouse gene symbol
#' @author guangchun
#' @keywords human2mouse
#' @examples
#' m1 <- c('APP', 'APOE')
#' human2mouse(m1)
#' @rdname human2mouse
#' @export

human2mouse <- function(msymb) {
    f = system.file("extdata", "HMD_HumanPhenotype.dat", package = "ewrefxn",
                    mustWork = T)
    humanMouseMap = read.delim(f,row.names = NULL,
                               comment.char = '#')
    if(any(duplicated(msymb))) stop('duplicated values found in input vector')
    df = humanMouseMap
    humanMouseMap1 = aggregate(df$mouse, list(human = df$human), paste, collapse = '|')
    pos <- match(msymb, humanMouseMap$human)
    humanMouseMap1[ pos, ]$x
}

