detectP.filter <- function(expr, detectP.mtr, p.cutoff=0.05, detectRate=0.2) {
    ##remove probes with >= detectRate missing values according to p.cutoff
    detectP <- matrix(T, nrow=dim(detectP.mtr)[1], ncol=dim(detectP.mtr)[2])
    detectP[detectP.mtr>p.cutoff]=F
    mode(detectP) <- "logical"
    detectFlag <- apply(detectP, 1, function(x) {(sum(x)/dim(detectP)[2]) >= detectRate})
    expr[detectFlag, ]->res
    res
}

rm.absent<- function(expr, PA.mtr, detectRate=0.2) {
    ##auxilary function used to remove Absent probes defined by mas5call of affy package, used for affymetrix chip processing
    detectFlag <- apply(PA.mtr, 1, function(x) {(sum(x=="P")/dim(PA.mtr)[2]) >= detectRate})
    expr[detectFlag, ]->res
    res
}

expValue.filter <- function(expr, val.cutoff=6, percentage=0.2) {
    ##remove probes with >= percentage missing values according to val.cutoff
    index <- matrix(T, nrow=dim(expr)[1], ncol=dim(expr)[2])
    index[expr<val.cutoff]=F
    mode(index) <- "logical"
    Flag <- apply(index, 1, function(x) {(sum(x)/dim(index)[2]) >= percentage})
    expr[Flag, ]->res
    res
}
#' matrix extractor of an geoquery object
#'
#' this fxn extract the expression matrix from an object returned by GEOquery
#' @param gse.obj the geoquery object returned by GEOquery
#' @param col specify the column number of expression value
#' @return expression matrix with each column correponding to a sample, each row to a gene/probe
#' @export
getMtr <- function(gse.obj, col=2) {
    ##this function extract a specified column from a GSMList and combine these columns to form a data matrix
    if(!require(GEOquery))
        stop("error loading GEOquery package\n")
    probesets <- Table(GPLList(gse.obj)[[1]])$ID
    res <- do.call('cbind', lapply(GSMList(gse.obj), function(x) {
                                              tab <- Table(x)
                                              pos <- match(probesets, tab$ID_REF)
                                              return(tab[pos, col])
                                          } ))
    res <- apply(res, 2 ,function(x) as.numeric(as.character(x) ) )
    rownames(res) <- probesets
    res
}

check.pl.unique <- function(gse.obj) {
    ##check platform uniqueness
    if(!require(GEOquery))
        stop("error loading GEOquery package\n")
    gsmpls=lapply(GSMList(gse.obj), function(x) Meta(x)$platform)
    if(length(unique(unlist(gsmpls) )) != 1){
        cat("multiple platforms detected: \n")
        cat(unique(unlist(gsmpls) ),"\n")
    } else {
          cat(unique(unlist(gsmpls) ),"\n")
          cat("no problem, go ahead...\n")
      }
}

print.datatable <- function(gse.obj) {
    ##print 5 header lines of an expression table
    if(!require(GEOquery))
        stop("error loading GEOquery package\n")
    print(Table(GSMList(gse.obj)[[1]])[1:5,])
}

print.meta <- function(gse.obj) {
    if(!require(GEOquery))
        stop("error loading GEOquery package\n")
    print(str(Meta(GSMList(gse.obj)[[1]])))
}
#' phenotype data extractor of an geoquery object
#'
#' this fxn extract the expression matrix from an object returned by GEOquery
#' @param gse.obj the geoquery object returned by GEOquery
#' @param phe.num specify the column number of expression value
#' @return phenotype information list with each element correponding to a sample
#' @export
getMeta <- function(gse.obj, phe.num) {
    ##this function get Meta information of those samples and return a matrix with each row corespond to a sample
    ##and each column a phenotype
    if(!require(GEOquery))
        stop("error loading GEOquery package\n")
    res <- lapply(GSMList(gse.obj), function(x) {
                                              tmp <- Meta(x)
                                              tmp[[phe.num]]})
    res
}

rm.na1 <- function(expr) {
    ##remove rows filled up with NA
    missing.ind <- which(is.na(expr),2)
    tmp <- table(missing.ind[,2])
    print("missing value information:\n")
    print(tmp)
    gene.names <- rownames(expr)
    if(length(unique(tmp)) == 1) {
        ##remove these lines
        nline <- unique(missing.ind[,1])
        res <- expr[-nline,]
        rownames(res) <- gene.names[-nline]
        return(res)
    } else {
          print("values are not missed in all samples, rm cancelled\n")
          return(expr)
      }
}
impute.na <- function(expr) {
    ##fill value-missing cells with imputation function
    if(!require(impute))
        stop("error loading impute\n")
    res <- impute.knn(expr)$data
    res
}

ilumin.soft.proc <- function(AccNum,detection.p=0.05,detectionRate = 0.5) {
    ##there are datasets in GEO that encapsulating raw illumina microarray data in GEO soft file format,
    ##parse it with this function
    if(!require(GEOquery))
        stop("error loading GEOquery package\n" )
    g.dat <- getGEO(GEO=AccNum, destdir=".", GSEMatrix=F,getGPL=F)
    ##gsm to expressionset according to http://www.bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.pdf
    gpls <- unlist(lapply(GSMList(g.dat), function(x) Meta(x)$platform))
    if(length(unique(gpls)) != 1)
        stop("multiple platforms in the data\n")

    if(!require(lumi))
        stop("error loading lumi package\n")
    ##we need to extract expr, se.expr, beadnum, detection, phenodata to initialize a LumiBatch object
    ##TODO: control probe data...
    ##expr first
    ##check the data table dimension
    gsms <- GSMList(g.dat)
    gsms.nrow <- unlist(lapply(gsms, function(x) dim(Table(x))[1]))
    gsms.ncol <- unlist(lapply(gsms, function(x) dim(Table(x))[2]))
    if( (length(unique(gsms.ncol)) != 1) | (length(unique(gsms.nrow)) != 1))
        stop("demensions of the data table of each samples are not the same, we have to stop\n")
    ## get the probeset ordering
    probesets <- Table(gsms[[1]])$ID_REF
    exprs <- do.call('cbind',lapply(gsms,function(x)
        {tab <- Table(x)
         mymatch <- match(probesets,tab$ID_REF)
         return(tab$VALUE[mymatch])
     }))
    exprs <- apply(exprs,2,function(x) {as.numeric(as.character(x))})

    ##the se.expr
    se.exprs <- do.call('cbind',lapply(gsms,function(x)
        {tab <- Table(x)
         mymatch <- match(probesets,tab$ID_REF)
         return(tab$BEAD_STDERR[mymatch])
     }))
    se.exprs <- apply(se.exprs,2,function(x) {as.numeric(as.character(x))})

    ##nbeads
    beadNum<- do.call('cbind',lapply(gsms,function(x)
        {tab <- Table(x)
         mymatch <- match(probesets,tab$ID_REF)
         return(tab$Avg_NBEADS[mymatch])
     }))
    beadNum<- apply(beadNum,2,function(x) {as.numeric(as.character(x))})

    ##detectionp
    detection<- do.call('cbind',lapply(gsms,function(x)
        {tab <- Table(x)
         mymatch <- match(probesets,tab$ID_REF)
         return(tab[mymatch,5])
     }))
    detection<- apply(detection,2,function(x) {as.numeric(as.character(x))})

    eset <- new("LumiBatch", exprs=exprs, se.exprs=se.exprs, beadNum=beadNum, detection = detection)
    lumiExpresso(eset, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), variance.stabilize = TRUE,
                 varianceStabilize.param = list(), normalize = TRUE, normalize.param = list(method='ssn'), QC.evaluation = TRUE,
                 QC.param = list(), verbose = TRUE)->lumi.N.Q

    expr0 <- exprs(lumi.N.Q)
    featureData0 <- fData(lumi.N.Q)
    detect.count <- detectionCall(lumi.N.Q,Th=detection.p,type="probe")
    detect.rate <- detect.count/ncol(expr0)
    detectP<- detectionCall(lumi.N.Q,Th=detection.p,type="matrix")
    detectP[detectP=="P"]=T
    detectP[detectP=="A"]=F
    mode(detectP) <- "logical"
    expr0[!detectP] <- NA
    ##we plna to use impute.knn to impute NA values, leave those probes detection p rate >= detectionRate
    expr0[detect.rate>=detectionRate,]->expr
    featureData <- featureData0[detect.rate >= detectionRate,]
    if(!require(impute))
        stop("error loading impute package\n")
    expr.impt <- impute.knn(expr)$data
    expr.impt
}
#' expression matrix generator
#'
#' given an accession number(gseXXXX), this fxn download and extract the expression matrix
#' @param AccNum accession number
#' @param filename read data from the local data instead of the remote server
#' @param exp.col column number of expression value
#' @param symbol.col column number of expression value, for arraymatrix data the GPL header is usually SYMBOL;
#'         for agilent44K it is GENE_SYMBOL(10th column of GPL4134)
#' @param entrezID.col column number of entrez id, for affymetrix data it is named Entrez_Gene_ID, for agilent it is
#'         GENE(9th column of GPL4134)
#' @param log2 set to T for logrithm transformation
#' @return a list with expmatrix probeID symbols and entrezID
#' @export
soft2exp <- function(AccNum,filename=NULL,exp.col = 2,symbol.col,entrezID.col ,log2=F) {
    ##GEO soft file processing
    ##use this kind of data is vunerable to data comparison problems(as the data may be lack of data preprocessing)
    if(!require(GEOquery))
        stop("error loading GEOquery package\n" )
    if(!is.null(filename))
        g.dat <- getGEO(filename=filename, destdir=".", GSEMatrix=F,getGPL=F)
    else
        g.dat <- getGEO(GEO=AccNum, destdir=".", GSEMatrix=F,getGPL=F)
    ##gsm to expressionset according to http://www.bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.pdf
    gpls <- unlist(lapply(GSMList(g.dat), function(x) Meta(x)$platform))
    if(length(unique(gpls)) != 1)
        stop("multiple platforms in the data\n")

    gsms <- GSMList(g.dat)
    gsms.nrow <- unlist(lapply(gsms, function(x) dim(Table(x))[1]))
    gsms.ncol <- unlist(lapply(gsms, function(x) dim(Table(x))[2]))
    if( (length(unique(gsms.ncol)) != 1) | (length(unique(gsms.nrow)) != 1))
        stop("demensions of the data table of each samples are not the same, we have to stop\n")

    ## get the probeset ordering
    probesets <- Table(gsms[[1]])$ID_REF
    exprs <- do.call('cbind',lapply(gsms,function(x)
        {tab <- Table(x)
         mymatch <- match(probesets,tab$ID_REF)
         return(tab[mymatch,exp.col])
     }))
    exprs <- apply(exprs,2,function(x) {as.numeric(as.character(x))})

    dat.exp <- as.matrix(exprs )
    mode(dat.exp) <- "numeric"
    gpl.annot <- Table(GPLList(g.dat)[[1]])
    pos <- match(probesets, gpl.annot$ID)
    if(any(is.na(pos)))
        stop("probesets not annoted in gpl table\n")
    symbols <- gpl.annot[pos,symbol.col]
    entrezIDs <- gpl.annot[pos,entrezID.col]
    if(log2)
        dat.exp <- log2(dat.exp)
    res <- list(exp=dat.exp, probeID=probesets, symbols=symbols, entrezIDs=entrezIDs)
    res
}
#' expression matrix generator (with detection pvalue)
#'
#' given an accession number(gseXXXX), this fxn download and extract the expression matrix
#' @param AccNum accession number
#' @param exp.col column number of expression value
#' @param detectp.col column number of detection p value
#' @param detectionRate threshold of sample detection rate for each gene
#' @param log2 set to T for logrithm transformation
#' @return a list with expmatrix probeID symbols and entrezID
#' @export
soft2exp.pval <- function(AccNum,detectp.col = 3,exp.col = 2,detectp=0.05,detectionRate=0.5,log2=F) {
    ##processing soft files with detection pvalue available(do detection p filtration)
    ##use this kind of data is vunerable to data comparison problems(as the data may be lack of data preprocessing)
    if(!require(GEOquery))
        stop("error loading GEOquery package\n" )
    g.dat <- getGEO(GEO=AccNum, destdir=".", GSEMatrix=F,getGPL=T)
    ##gsm to expressionset according to http://www.bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.pdf
    gpls <- unlist(lapply(GSMList(g.dat), function(x) Meta(x)$platform))
    if(length(unique(gpls)) != 1)
        stop("multiple platforms in the data\n")

    gsms <- GSMList(g.dat)
    gsms.nrow <- unlist(lapply(gsms, function(x) dim(Table(x))[1]))
    gsms.ncol <- unlist(lapply(gsms, function(x) dim(Table(x))[2]))
    if( (length(unique(gsms.ncol)) != 1) | (length(unique(gsms.nrow)) != 1))
        stop("demensions of the data table of each samples are not the same, we have to stop\n")

    ## get the probeset ordering
    probesets <- Table(gsms[[1]])$ID_REF
    exprs <- do.call('cbind',lapply(gsms,function(x)
        {tab <- Table(x)
         mymatch <- match(probesets,tab$ID_REF)
         return(tab[mymatch,exp.col])
     }))
    exprs <- apply(exprs,2,function(x) {as.numeric(as.character(x))})

    ##detectionp
    detection<- do.call('cbind',lapply(gsms,function(x)
        {tab <- Table(x)
         mymatch <- match(probesets,tab$ID_REF)
         return(tab[mymatch,detectp.col])
     }))
    detection<- apply(detection,2,function(x) {as.numeric(as.character(x))})

    dat.exp <- as.matrix(exprs )
    mode(dat.exp) <- "numeric"
    dat.detectp <- as.matrix(detection)
    mode(dat.detectp) <- "numeric"
    dat.detectp1 <- matrix(data=T, nrow=nrow(dat.detectp), ncol=ncol(dat.detectp))
    dat.detectp1[dat.detectp > detectp]=F
    mode(dat.detectp1) <- "logical"
    detect.count <- apply(dat.detectp1, 1, sum)
    detect.rate <- detect.count/ncol(dat.detectp1)
    dat.exp[!dat.detectp1] <- NA
    ##we plan to use impute.knn to impute NA values, leave those probes detection p rate >= detectionRate
    dat.exp[detect.rate>=detectionRate,]-> dat.exp1
    probesets[detect.rate >= detectionRate] ->probesets
    gpl.annot <- Table(GPLList(g.dat)[[1]])
    pos <- match(probesets, gpl.annot$ID)
    if(any(is.na(pos)))
        stop("probesets not annoted in gpl table\n")
    symbols <- gpl.annot$Symbol[pos]
    entrezIDs <- gpl.annot$Entrez_Gene_ID[pos]
    if(!require(impute))
        stop("error loading impute package\n")
    expr.impt <- impute.knn(dat.exp1)$data
    if(log2)
        expr.impt <- log2(expr.impt)
    res <- list(exp=expr.impt, probeID=probesets, symbols=symbols, entrezIDs=entrezIDs)
    res
}

##probeID to symbol with geoquery
#' @export
probeID2symbol <- function(id_vec, gpl, symbol.col=NULL) {
    if(!require(GEOquery))
        stop("error loading GEOquery\n")
    id_vec <- as.character(id_vec)
    gpl.dat <- getGEO(gpl, destdir=".")
    annot.df <- gpl.dat@dataTable@table
    pos <- match(id_vec, as.character(annot.df$ID))
    if(is.null(symbol.col))
        res <- as.character(annot.df$GENE_SYMBOL)[pos]
    else
        res <- as.character(annot.df[,symbol.col])[pos]
    names(res) <- id_vec
    res
}

#' @export
geoSupp.fetch <- function(AccNum,
                          destdir = getwd(),local=F,extract=T ) {
    if(!require(GEOquery))
        stop("error loading GEOquery package\n")
    if(!local) {
        geo.flist <- getGEOSuppFiles(AccNum,makeDirectory = F)
        ##by default, extrat the first file with name *RAW*
        rawfile <- grep("\\.*RAW\\.*",rownames(geo.flist),value=T)[1]
    } else {
          geo.flist <- list.files()
          rawfile <- grep("\\.*RAW\\.*",geo.flist,value=T)[1]
      }
    if(extract) {
        untar(rawfile,
              exdir=destdir)
    }
    ##get the file names of the extracted files
    extract.flist <- untar(rawfile, list=T)
}

#' @export
grepGPLtable <- function(query,gpl) {
    ##this function grep query against gpl table generated by GEOquery Table(gpl)
    if(!require(GEOquery))
        stop("error loading GEOquery package\n")
    gpl <- getGEO(gpl, destdir=".")
    tab <- Table(gpl)
    hit.index <- apply(tab, 1, function(x) grepl(query, x) )
    if(sum(hit.index) > 0 & sum(hit.index) <= 10)
        print(tab[hit.index, ])
    else if(sum(hit.index) > 10)
        print("query matches >= 10, try another one\n")
    else
        print("no match found\n")
}
#' agilent 4*44k microchip preprocessing method
#'
#' @param geo.flist geo file list
#' @export
geo.agilent44k.prep <- function(geo.flist) {
     cat("array data files we will process are shown bellow\n")
    print(geo.flist)
    answer <- dual.choice("sure we go ahead?",value1 = "y",value2="n")
    if(answer == "n")
        return(NULL)
    prep.obj <- agilent44k.prep(array.fnames = geo.flist)
    prep.obj
}
#' agilent 4*44k microchip preprocessing method
#'
#' @param geo.flist geo file list
#' @export
geo.agilent80k.prep <- function(geo.flist) {
    res <- geo.agilent44k.prep(geo.flist)
    res
}
