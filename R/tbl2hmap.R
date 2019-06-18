#' convert long tibble to pheatmap data input
#'
#' @param tbl.long tibble in long format
#' @param rowVar column name used as row for pheatmap matrix
#' @param colVar column name used as column for pheatmap matrix
#' @param valueVar column name used as the values for pheatmap matrix
#' @param colAnnVars column name(s) used for additonal column annotation
#' @param rowAnnVars column name(s) used for additonal row annotation
#' @return a list with 3 components: mat, rowAnn, colAnn
#' @export

tbl2hmap= function(tbl.long, rowVar, colVar, valueVar,
                   colAnnVars = NULL, rowAnnVars = NULL) {
    Mat0 = dplyr::select(tbl.long, one_of(c(rowVar, colVar,valueVar)))%>%
        tidyr::spread_(key = colVar, value=valueVar)
    Mat = select(Mat0, -one_of(rowVar)) %>% as.matrix()
    rownames(Mat) = unlist(select(Mat0,one_of(rowVar)))

    colAnn = rowAnn= NULL

    if(!is.null(colAnnVars)) {
        colAnn0 = dplyr::select(tbl.long,one_of( c(colVar,colAnnVars))) %>%
            unique()
        colAnn = select(colAnn0, - one_of(colVar)) %>% as.data.frame()
        rownames(colAnn) = unlist(select( colAnn0,one_of(colVar)))
    }

     if(!is.null(rowAnnVars)) {
        rowAnn0 = dplyr::select(tbl.long,one_of(c( rowVar,rowAnnVars))) %>%
            unique()
        rowAnn = select(rowAnn0, -one_of( rowVar)) %>% as.data.frame()
        rownames(rowAnn) = unlist(select(rowAnn0,one_of(rowVar)))
    }

    list(mat = Mat, rowAnn =rowAnn, colAnn = colAnn)
}
