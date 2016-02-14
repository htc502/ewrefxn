#' internal fxn for array data preprocessing
#'
#' this function is used to extract and prepare the data matrix which will be used to do differential expression calculation.
#' if class1 is set to NULL, the returned value will be used for one class DEG analysis.
#' @param dat.mtr an expression data matrix with one column per sample; each gene as a row.
#' @param featureName optional gene name for the expression matrix(vector with length == nrow(dat.mtr)).
#' @param sample.label sample treament(condition) labels for each sample(vector with length == ncol(dat.mtr)).
#' @param class1 class1 and class2 are strings indicating which groups are compared.
#' @param class2 class1 and class2 are strings indicating which groups are compared.
#' @return a list with 3 elements: dat for expression data matrix, cls for the coressponding sample labels, featureName is attached as also.
#' @export
datPrep <- function(dat.mtr, featureName=NULL, sample.label,
                    class1, class2){
    ##set class1 = NULL for one class DEG analysis
    cl <- sample.label
    cl1.index<-rep(F,length=length(cl))
    if(!is.null(class1)) {
        f1<-unique(unlist(strsplit(class1,split="\\+")))
        if(any(f1 == ""))
            stop(class1, "is malformed\n")
        for(i in 1:length(f1) ) {
            cl1.index[cl==f1[i]]<-T
        }
        rm(f1)
    }
    cl2.index<-rep(F,length=length(cl))
    f2<-unique(unlist(strsplit(class2,split="\\+")))
    if(any(f2 == ""))
        stop(class2, "is malformed\n")
    for(i in 1:length(f2) ) {
        cl2.index[cl==f2[i]]<-T
    }
    rm(f2)

    if(sum(cl1.index) == 0 & sum(cl2.index) == 0)
        stop("both ",class1," and ",class2," are not found in cls\n")
    cls <- c( rep(0,length=sum(cl1.index)),
             rep(1,length=sum(cl2.index)) )
    dat.mtr1<- cbind(dat.mtr[,cl1.index], dat.mtr[,cl2.index] )
    res <- list(dat=dat.mtr1,
                cls=cls,
                featureName=featureName)
}
get.arg.value <- function(arg.name,named.arg.list) {
    ##return the arg value if found
    ##else return NULL
    ind <- match(arg.name ,names(named.arg.list))
    res <- ifelse(is.na(ind),NULL, named.arg.list[[ind]])
    res
}

#' internal fxn for array data DEG calc
#'
#' @param dat.obj this is a list returned by dat.prep function.
#' @param diff.method possible values are: RP(for rankprod) limma sam
#' @param ... further parameters passed to RP limma or sam algrithm
#' @return a list with 2 tables: up&down regulated genes expression fold change table
#' @export
diffCalc <- function(dat.obj,diff.method,... ){
    ##return FC p value and fdr without thresholding
    ##leave thresholding to diff.output
    print("We assume that the data input is log2 transformed...")
    additional.args <- list(...)
    if(diff.method == "RP") {
        ##code chunck for RankProd is identical for both
        ##case/control design and single R/G ratio design
        if(!require(RankProd))
            stop("error loading RankProd package,try \n
install.packages(RankProd) to install this package\n")
        rp.out <- RP(data=dat.obj$dat, cl=dat.obj$cls, ...)
        topGene(x=rp.out, cutoff=NULL,
                num.gene = nrow(dat.obj$dat),
                gene.names = dat.obj$featureName,
                logged = T)->topgene
        res <- topgene
        ##steps following makes the output identical among different methods...we define taht a deg result is a list that contains two elements: up regulated(class1<class2) table with FC<1; down regulated table with FC>1, each table contains 4 columns: gene.name FC p.value fdr(pfp)
        res$Table1 <- data.frame(gene=rownames(topgene$Table1),topgene$Table1,check.names=F, stringsAsFactors = F, row.names=NULL)
        res$Table2 <- data.frame(gene=rownames(topgene$Table2),topgene$Table2,check.names=F, stringsAsFactors = F, row.names=NULL)
        res$Table1 <- res$Table1[,c(1,4,6,5)] ##check this when RankProd package updates
        res$Table1 <- res$Table1[ res$Table1[ ,2 ] < 1 , ]
        res$Table2 <- res$Table2[,c(1,4,6,5)] ##check this when RankProd package updates
        res$Table2 <- res$Table2[ res$Table2[ ,2 ] >= 1 , ]
        attributes(res) <- list(method="RP",names=names(res))
        return(res)
    }
    if(diff.method == "limma") {
        ##for case/control analysis, limma generate a
        ##model matrix and a contrast matrix, for R/G ratio
        ##analysis, just omit the design matrix in lmfit
        ##from lmFit funciton help we know that limma assume taht data is log transformed
        if(!require(limma))
            stop("error loading limma package, try \n
install.packages(limma) may help\n")
        if(length(unique(dat.obj$cls) ) == 1) {
            ##for R/G ratio data, should be tested
            ##this is copied from limma user guide page57
            fit <- lmFit(dat.obj$dat)
            fit <- eBayes(fit)
            res <- topTable(fit,number = length(dat.obj$featureName),
                            genelist=dat.obj$featureName)

        } else { ##single channel data
            cl <- factor(dat.obj$cls)
            ##rename cl level names for later use
            levels(cl) <- c("ctrl","case")
            design.mtr <- model.matrix(~0+cl)
            colnames(design.mtr) <- levels(cl)
            contrasts <- paste(levels(cl)[1],"-",levels(cl)[2],sep="")
            contra.mtr <- makeContrasts(contrasts = contrasts , levels = design.mtr)
            fit <- lmFit(dat.obj$dat,design.mtr)
            fit2 <- contrasts.fit(fit, contra.mtr)
            fit2 <- eBayes(fit2)
            res <- topTable(fit2, number = length(dat.obj$featureName),
                            genelist=dat.obj$featureName)
        }
        ##column order for 1,2,8,9: gene, logFC, p, fdr
        Table1 <- res[ res$logFC < 0 , c(1, 2, 5,6) ]
        Table2 <- res[ res$logFC > 0 , c(1, 2, 5,6) ]
        ##to keep in accordance with rankprod, transform logfc to fc
        Table1[ ,2 ] <- 2^Table1[ ,2 ]
        colnames(Table1)[2] <- 'FC'
        Table2[ ,2 ] <- 2^Table2[ ,2 ]
        colnames(Table2)[2] <- 'FC'
        res <- list(Table1=Table1, Table2=Table2)
        attributes(res) <- list(method="limma",names=names(res))
        return(res)
        ##from limma user guide pdf we know that limma is capable of
        ##more complex comparison, amazing, +1 for linea regression model formulae...
    }
    if(diff.method == "sam") {
        if(!require(samr))
            stop("error loading samr package\n")
        ##samr require the cls should be composed of 1 and 2...
        cls <- factor(dat.obj$cls)
        levels(cls) <- c(1,2)
        resp.type <- ifelse(length(unique(dat.obj$cls)) == 1,
                            "One class","Two class unpaired")

        samr.obj <- SAM(x=dat.obj$dat,
                        y=cls,
                        geneid = NULL,
                        genenames = dat.obj$featureName,
                        logged2 = T,
                        random.seed = 123,
                        nperms = 100,
                        fdr.output = 1,
                        resp.type = resp.type)

        ##samr has q-value in the table, just used it as fdr,
        res <- samr.obj$siggenes.table
        Table1 <- res$genes.up
        Table1 <- Table1[ ,c(1,6,7) ]
        Table2 <- res$genes.lo
        Table2 <- Table2[ ,c(1,6,7) ]
        res <- list(Table1=Table1,
                    Table2=Table2)
        ##no pvalue quatity...
        attributes(res) <- list(method="sam",names=names(res))
        res
    }
}
#' internal fxn for thresholding and output DE result
#'
#' @param diff.calc.obj this is a list returned by diff.calc function.
#' @param FC fold change threshold
#' @param fdr false discovery rate threshold
#' @param prefix result file prefix string
#' @return this fxn has no return value
#' @export
diffOutput <- function(diff.calc.obj,FC=1.5,
                       fdr=0.05,prefix = ""){
    if(attributes(diff.calc.obj)$method == "sam") {##sam output is different...
        table1 <- diff.calc.obj$Table1[ diff.calc.obj$Table1[ ,2 ] >= FC , ]
        table1 <- table1[ table1[ ,3 ] <= fdr, ]
        table2 <- diff.calc.obj$Table2[ diff.calc.obj$Table2[ ,2 ] <= 1/FC , ]
        table2 <- table2[ table2[ ,3 ] <= fdr, ]

        tab1Name <- paste("table1.csv",sep="")
        tab2Name <- paste("table2.csv",sep="")
        if(prefix != ""){
            tab1Name <- paste0(prefix, "_", tab1Name)
            tab2Name <- paste0(prefix, "_", tab2Name)
        }
        write.csv(table1,tab1Name,row.names=F)
        write.csv(table2,tab2Name,row.names=F)

    } else if(attributes(diff.calc.obj)$method == "fc" ) {
        table1 <- diff.calc.obj$Table1[ diff.calc.obj$Table1[ ,4 ] <= 1/FC , ]
        table2 <- diff.calc.obj$Table2[ diff.calc.obj$Table2[ ,4 ] >= FC , ]

        tab1Name <- paste("table1.csv",sep="")
        tab2Name <- paste("table2.csv",sep="")
        if(prefix != ""){
            tab1Name <- paste0(prefix, "_", tab1Name)
            tab2Name <- paste0(prefix, "_", tab2Name)
        }
        write.csv(table1,tab1Name,row.names=F)
        write.csv(table2,tab2Name,row.names=F)
    } else {
        ##here we assume that the data structure of the diff.calc.obj is : geneName FC fdr Pvalue
        ##fc is calculated as class1/class2
        ##if we assume that class1 is control grp, then table 1 is up regulated genes, table 2 is down regulated genes
        ##table column: gene fc p.value fdr
        table1 <- diff.calc.obj$Table1[ diff.calc.obj$Table1[ ,2 ] <= 1/FC , ]
        table1 <- table1[ table1[ ,4 ] <= fdr, ]
        table2 <- diff.calc.obj$Table2[ diff.calc.obj$Table2[ ,2 ] >= FC , ]
        table2 <- table2[ table2[ ,4 ] <= fdr, ]

        tab1Name <- paste("table1.csv",sep="")
        tab2Name <- paste("table2.csv",sep="")
        if(prefix != ""){
            tab1Name <- paste0(prefix, "_", tab1Name)
            tab2Name <- paste0(prefix, "_", tab2Name)
        }
        write.csv(table1,tab1Name,row.names=F)
        write.csv(table2,tab2Name,row.names=F)
    }
}

SymbolAggregate<- function(expr, map.db, map.name,get.symbol.only=F) {
    ## first map probe id to symbol using annotaion packages privided by bioconductor;
    ##then aggregate the data matrix according to the symbol
    ##set get.symbol.only=T in case if you only want to convert probeID to symbol
    if(!require(map.db,character.only=T))
        stop("error loading ",map.db,"\n")
    symbol.list=as.list(map.name[rownames(expr)])
    ##remove multiple mapped& non-mapped probes...
    lapply(symbol.list, length) -> tmp
    multi.ind <- unlist(tmp) >= 2
    sum(multi.ind) ->n.multi
    cat(n.multi, " multiple mapped probes\n")
    non.ind <- unlist(lapply(symbol.list, is.na))
    sum(non.ind)->n.nonmapped
    cat(n.nonmapped, " none mapped probes\n")

    probe.rm = multi.ind | non.ind
    expr1 =expr[!probe.rm, ]
    symbol.list=symbol.list[!probe.rm]
    symbol=unlist(symbol.list)
    if(get.symbol.only)
        return(symbol)
    ##aggregate probes to get gene level expression matrix
    aggregate(expr1, by=list(gene=symbol), FUN=median)->aggr.tmp
    rownames(aggr.tmp)=aggr.tmp[,1]
    aggr.tmp=aggr.tmp[,-1]
    exp=aggr.tmp
    exp
}

SymbolAggregate1 <- function(expr, gpl, sym.col, get.symbol.only=F) {
    ##same as SymbolAggregate, but use annotation table provided by GEO gpl table,
    ##downloading gpl table may take minites...
    ##set get.symbol.only=T in case if you only want to convert probeID to symbol
    if(!require(GEOquery))
        stop("error loading GEOquery package\n")
    gpl <- getGEO(gpl, destdir=".")
    tab <- Table(gpl)
    probeid <- rownames(expr)
    pos <- match(probeid, tab$ID)
    missing.ind <- is.na(pos)
    cat(sum(missing.ind), "probes unmapped\n")
    expr <- expr[!missing.ind, ]
    pos <- pos[!missing.ind]
    symbol <- as.character(tab[pos, sym.col])
    if(get.symbol.only)
        return(symbol)
    aggregate(expr, by=list(gene=symbol), FUN=median)->aggr.tmp
    rownames(aggr.tmp)=aggr.tmp[,1]
    aggr.tmp=aggr.tmp[,-1]
    exp=aggr.tmp
    exp
}

remove100 <- function(expr) {
    exp.medi <- apply(expr, 1, median)
    ind <- exp.medi >= 100
    print(paste(sum(!ind)," probes with median expression level < 100 to remove...\n",sep=""))
    res <- expr[ ind, ]
    res
}


##code from microarray_analysis.R in stem_cell dir
##Time function to get time as file name
label.time<-function(){
    Systim<-Sys.time()
    Systim<-as.character(Systim)
    Systim<-strsplit(Systim," ")
    dat<-Systim[[1]][1]
    tim<-Systim[[1]][2]
    tim<-strsplit(tim,":")
    h<-tim[[1]][1]
    m<-tim[[1]][2]
    sec<-tim[[1]][3]
    tim<-paste(h,m,sec,sep="-")
    return(paste(dat,tim,sep="_"))
}

annot.affyProbeID <- function(affy.annot.file, probeIDs, target.col=15) {
    skip=0
    ##read file line by line in R
    con <- file(affy.annot.file, open="r")
    while(1) {
        oneline <- readLines(con, n =1 )
        if(!grepl("^#",oneline))
            break
        skip <- skip + 1
    }
    close(con)
    annot.df <- read.csv(affy.annot.file, skip=skip)
    match(probeIDs, annot.df[,1]) -> pos
    as.character(annot.df[pos, target.col]) -> res
    res
}
#' internal fxn for affymatrix IVT microarray preprocessing
#'
#' @param phenoData a table derived from sdrf file which must have two columns called ArrayDataFile and State
#' @param Preprocessing prep method to use, possible values: mas5, dchip, rma, gcrma
#' @param PMA_cutoff percentage threshold of samples which failed the PMA calling test
#' @param auto.annot set to TURE for automatically annotation of the probes
#' @param ann.file if auto.annot = F, this file will be used to annote the probes
#' @param ann.col if auto.annot = F, this column of the ann.file will be used to annote the probes
#' @param multiple2one method for merge the probes corresponding to the same gene
#' @return an expressionSet object
#' @author dapeng, liang and guangchun, han
#' @export
affy.Microarray.preprocessing = function(phenoData, Preprocessing = "mas5", PMA_cutoff = 0.1, auto.annot=T,ann.file=NULL, ann.col=13, multiple2one = "Variance") {
                                        #Check the parameters
    if(!any(Preprocessing == c("mas5","rma","gcrma","dchip")))
        stop('Preprocessing method should be one of "mas5", "rma", "gcrma", "dchip"')
    if(PMA_cutoff>1 | PMA_cutoff<=0)
        stop("PMA_cutoff should be 0< value <=1")
    if(!any(multiple2one == c("V", "Variance", "M", "Mean", "Aggregation" ,"A")))
        stop('mutiple probsets to one gene should be handled with one of "V","Variance", "M", "Mean", "Aggregation" ,"A"')
    if(!require(affy))
        stop("Please make sure that you have load affy or your chip specific cdf package")
                                        #Reads in the phenoData file
    cat("Reading in " ,phenoData, " file...\n");
    pd = read.AnnotatedDataFrame(phenoData, header=TRUE, sep = "\t", row.names="ArrayDataFile",stringsAsFactors = FALSE)
    if(is.null(rownames(pData(pd))))
        stop( "a column named ArrayDataFile should be there as we use this column as the rownames of phenoData table")

                                        #Reads in the raw data directly from the Affymetrix CEL files and store the sample covaribale
    cat("Reading in raw CEL files and sample covaribale ...\n");
    affydata = ReadAffy(filenames = rownames(pData(pd)),phenoData=pd);

                                        #Processes the data with different methods
    cat("Running", Preprocessing, "and MAS5 calls on the raw data...\n");
    if(Preprocessing == "mas5")
        eset = mas5(affydata);
    if(Preprocessing == "rma")
        eset = rma(affydata);
    if(Preprocessing == "gcrma"){
        if(require(Preprocessing, character.only = TRUE)){
            eset = gcrma(affydata);
        }
        else
            stop("Please make sure that you have load gcrma or corresponding chip-specific probe packages")
    }
    if(Preprocessing == "dchip")
        eset = expresso(affydata, bg.correct = FALSE,normalize.method="invariantset", pmcorrect.method = "pmonly",summary.method="liwong")

                                        #log2 transformation
    if(any(Preprocessing == c("mas5", "dchip")))
        expr_log_matrix = log(exprs(eset), 2)
    else
        expr_log_matrix = exprs(eset)

                                        #PMA_calls Filtering
    data.mas5calls = mas5calls(affydata);
    PMA_calls = exprs(data.mas5calls)
    expr_rm_absent = expr_log_matrix [rowSums(PMA_calls == "P")/ncol(PMA_calls)>=PMA_cutoff, ]

                                        #Control Probe Filtering
    expr_rm_ctrl = expr_rm_absent[-grep("AFFX", rownames(expr_rm_absent)), ]

    if(auto.annot) {
        ##Loads the chip specific annotation libraries
        db.package = paste(annotation(eset), ".db", sep = "", delim = "");
        if(!require(db.package, character.only = TRUE) & !require(annaffy))
            stop("Please make sure that you have load ", db.package,"and ","annaffy")

        ##removing unannoted probes corresponding to unidentified EST
        annoted_probes = getText(aafLocusLink(rownames(expr_rm_ctrl), db.package));
        expr_rm_unannoted = expr_rm_ctrl[annoted_probes != "", ]

        ##probe ids  map to gene identifier
        geneid = getText(aafLocusLink(rownames(expr_rm_unannoted), db.package ));
    } else {
        if(is.null(ann.file) | is.null(ann.col))
            stop("annotation file and target column number is required for mannally annotation mode\n")
        annoted_probes = annot.affyProbeID(ann.file, rownames(expr_rm_ctrl),ann.col);
        expr_rm_unannoted = expr_rm_ctrl[!is.na(annoted_probes ), ]
        geneid =annot.affyProbeID(ann.file, rownames(expr_rm_unannoted),ann.col);
    }

    cat("Get representative probesets for each gene by", multiple2one, "...\n")
                                        #Mutiple ProbSets to the same genes(highest var)
    if(multiple2one == "Variance"|multiple2one == "V"){
        probe_var = apply(expr_rm_unannoted ,1,var)
        get_hVAR = function(duplicate_probe){
            if(length(duplicate_probe) == 1)
                return(duplicate_probe)
            else
                return(duplicate_probe[which.max(probe_var[duplicate_probe])])
        }
        hvar_probe = tapply(rownames(expr_rm_unannoted), geneid, get_hVAR)
        pos.hvar_probe <- match(hvar_probe, rownames(expr_rm_unannoted))
        expr_log_preprocessing = expr_rm_unannoted [pos.hvar_probe, ]

                                        #probe annotation
        geneid = geneid[pos.hvar_probe];
        rownames(expr_log_preprocessing) = geneid
    }

                                        #Mutiple ProbSets to the same genes(highest mean)
    if(multiple2one == "Mean"|multiple2one == "M"){
        probe_mean = apply(expr_rm_unannoted , 1, mean)
        get_hMEAN = function(duplicate_probe){
            if(length(duplicate_probe) == 1)
                return(duplicate_probe)
            else
                return(duplicate_probe[which.max(probe_mean[duplicate_probe])])
        }
        hmean_probe = tapply(rownames(expr_rm_unannoted), geneid, get_hMEAN)
        pos.hmean_probe <- match(hmean_probe, rownames(expr_rm_unannoted))
        expr_log_preprocessing = expr_rm_unannoted [pos.hmean_probe, ]

                                        #probe annotation
        geneid = geneid[pos.hmean_probe];
        rownames(expr_log_preprocessing) = geneid
    }

                                        #Mutiple ProbSets to the same genes(average)
    if(multiple2one == "Aggregation"|multiple2one == "A"){
        expr_log_preprocessing = aggregate(as.data.frame(expr_rm_unannoted), by = list(geneid), mean)
        rownames(expr_log_preprocessing) = expr_log_preprocessing[,1]
        expr_log_preprocessing = as.matrix(expr_log_preprocessing[,-1])
    }

                                        #Retains only expression information for the remaining probes after filtering
    exprs(eset) = expr_log_preprocessing ;

    cat("Save the results\n")
    Tim<-label.time()
    fname<-paste(Tim,"preprocessing.rData",sep="_")
    save(eset,file=fname)

    cat("Returns expression information (expressionset class)\n")
    return(eset);
}
#' internal fxn for illumina microarray preprocessing
#'
#' @param dataFile this is the expression datafile exported from genomeStudio
#' @param controlFile control probe signal file. This file is optional, it is also exported form genomestudio for background correction
#' @param arrayDesignFile this is the `sdrf file`
#' @param detection.p detection pvalue threshold
#' @param detectionRate PMA_cutoff analogue to affymatrix array
#' @param norm.method normalization method
#' @return a list of expressionMatrix, sample metaInfo, and probe annotation info
#' @author guangchun han
#' @export
ilumi.preprocess <- function(dataFile,controlFile=NULL,arrayDesignFile=NULL, detection.p=0.05, detectionRate=0.5,
                             norm.method = 'ssn') {
    ##to be accomplished
    if(!require(lumi))
        stop("error loading lumi package\n")
    lumiR(dataFile) -> lumi.Ds
    ##add control data to lumi.Ds
    if(!is.null(controlFile))
        lumi.Ds <- addControlData2lumi(controlFile,lumi.Ds)
    if(!is.null(arrayDesignFile) ) {
        read.delim(arrayDesignFile)->metaInfo
        tmp <- pData(lumi.Ds)
        pData(lumi.Ds)<-cbind(pData(lumi.Ds),metaInfo[match(tmp$sampleID,metaInfo$arrayID),])
        rm(tmp)
    }

    lumiExpresso(lumi.Ds, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), variance.stabilize = TRUE,
                 varianceStabilize.param = list(), normalize = TRUE, normalize.param = list(method=norm.method), QC.evaluation = TRUE,
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
    list(exp=expr.impt, fd=featureData, pd=pData(lumi.N.Q))->res
    res
}


ask.column <- function(default.value,
                       display.info,
                       max) {
    res <- default.value
    ok <- F
    while(!ok) {
        answer <- readline(display.info)
        if(answer == "")
            return(res)
        answer <- as.numeric(answer)
        if(answer > max) {
            print("illegal input, try again.")
            next
        }
        res <- answer
        ok <- T
    }
    res
}

dual.choice <- function(display.info,value1, value2) {
    while(T) {
        answer <- readline(paste0(display.info,"(",value1,"/",value2,")"))
        if(answer == value1)
            return(answer)
        if(answer == value2)
            return(answer)
        else {
            print(paste0(value1," or ",value2))
            next
        }
    }
}
#' function used to filter expression matrix based on detection pvalue/PMA matrix
#'
#' @export
#' @return a list with express matrix; detection pvalue and an index indicating which probes are filtered out
#' @author guangchun
rm.absent<- function(expr, PA.mtr=NULL, detectRate=0.2,DP.mtr=NULL,dp.cutoff=0.05) {
    ##auxilary function used to remove Absent probes defined by mas5call of affy package, used for affymetrix chip processing
    if(all(c(is.null(PA.mtr),is.null(DP.mtr))) )
        stop('provide at least one of PA.mtr/DP.mtr')
    if(!is.null(PA.mtr)) {
        apply(PA.mtr,2,function(x) {sum(x=='P')/dim(PA.mtr)[1]} )->sm.detectionP
        if(sum(sm.detectionP <= 0.2) > 0)
            print(paste0(paste(which(sm.detectionP <= 0.2),collapse=','),'has less than 50% probes detected, this is to inform you that there may exist bad quality samples...'))
        detectFlag <- apply(PA.mtr, 1, function(x) {(sum(x=="P")/dim(PA.mtr)[2]) >= detectRate})
        expr[detectFlag, ]->expr1
        PA.mtr1 <- PA.mtr[detectFlag, ]
        res <- list(expr=expr1, PAmtr=PA.mtr1,idx=detectFlag)
        return(res)
    } else {
        apply(DP.mtr,2,function(x) {sum(x<=dp.cutoff)/dim(DP.mtr)[1]} )->sm.detectionP
        if(sum(sm.detectionP <= 0.2) > 0)
            print(paste0(paste(which(sm.detectionP <= 0.2),collapse=','),'has less than 50% probes detected, this is to inform you that there may exist bad quality samples...'))
        detectFlag <- apply(DP.mtr, 1, function(x) {(sum(x<=dp.cutoff)/dim(DP.mtr)[2]) >= detectRate})
        expr[detectFlag, ]->expr1
        DP.mtr1 <- DP.mtr[detectFlag, ]
        res <- list(expr=expr1, DPmtr=DP.mtr1,idx=detectFlag)
        return(res)
    }
}

#'agilent 4*44k microarray raw data preprocessing
#' @param array.fnames a character vector of array data file names
#' @param detect.rate PMA_cutoff analogue to affymatrix array
#' @return a list of expressionMatrix and probe annotation info
#' @author guangchun han
#' @export
agilent44k.prep <- function(array.fnames, detect.rate = .2) {
    ##array.fnames is a character vector of array data file names
    ##from limma user guide we know that read.mainage is used for data input for any microarray platform except for illumina and affymetrix
    ##we have to separate annotation from preprocessing, or , at least make
    ##annotation as a modular step in prep...
    if(!require(limma))
        stop("error loading limma package\n")
    elist <- read.maimages(array.fnames,
                           source="agilent", green.only = T)
    ##elist.p : processed elist
    print("applying background corecction and normexp normalization")
    elist.p <- backgroundCorrect(elist, method = "normexp")
    print("between array normalzation with quantile")
    elist.p <- normalizeBetweenArrays(elist.p, method = "quantile")

    ##low expression probe/control probe remove
    neg95 <- apply(elist.p$E[ elist.p$genes$ControlType == -1, ],
                   2, function(x) quantile(x, p=.95) )
    cutoff <- matrix(1.1*neg95, nrow(elist.p), ncol(elist.p),
                     byrow = T)
    isexpr <- rowSums(elist.p$E > cutoff) >= detect.rate * ncol(elist.p)
    print("number of probes left after applying filtration: ")
    print(table(isexpr))
    ##[, has been defined for EList in limma, El[ 1:3, ]is legal...
    elist.res <- elist.p[ elist.p$genes$ControlType == 0 & isexpr, ]
    res <- list(exp=elist.res$E, feature=elist.res$genes$GeneName)
    res
}

#' check if a package has been installed or not
#' @export
is.installed <- function(pkg.name) is.element(pkg.name, installed.packages()[,1])

#' install a package using bioconductor
#' @export
biocInst <- function(pkg.name) {
    source("http://bioconductor.org/biocLite.R")
    biocLite(pkg.name, suppressUpdates=T)
}
#' if hasn't been installed, install and load
#' @export
require0 <- function(pkg.name) {
    if(is.installed(pkg.name)) {
        if(!require(pkg.name, character.only = T, quietly = T))
            return(F)
    } else {
        biocInst(pkg.name)
        if(!require(pkg.name, character.only=T, quietly = T))
            return(F)
    }
    T
}
#' affymatrix ivt microarray preprocessing 'lite' version
#' @param array.fnames a character vector of array data file names
#' @return a list of expressionMatrix and probe annotation info
#' @author guangchun han
#' @export
affy.ivt.prep <- function(array.fnames) {
    ##gcrma for background correction; gcrma also applies quantile normalization by default;it also applies log2 transform  to the data
    print('probe detectionRate 0.2, be careful with this for experiments with unbalanced case-control designs(ie. 100 controls .vs. 10 cases)')
    detection.p.cutoff <- 0.05
    detectionRate <- 0.2
    if(!require(affy))
        stop("error loading affy\n")
    affydata <- ReadAffy(filenames = array.fnames)
    if(!require0("gcrma"))
        stop("error loading gcrma\n")
    eset <- gcrma(affydata)
    pma.calls <- exprs(mas5calls(affydata))
    exp.mtr <- exprs(eset)
    exp.mtr1 <- exp.mtr[ rowSums(pma.calls == "P")/ncol(pma.calls) >= detectionRate, ]
    db.package <- paste(annotation(eset), ".db",sep="")
    if(!require0(db.package))
        stop("error loading ",db.package,"\n")
    if(!require0("annotate"))
        stop("error loading annotation package\n")
    symbols <- getSYMBOL(rownames(exp.mtr1),annotation(eset))
    ##symbols is a vector
    missing <- is.na(symbols)
    exp.mtr1 <- exp.mtr1[!missing, ]
    symbols <- symbols[!missing]
    res <- list(exp=exp.mtr1, feature=symbols)
    res
}
#' affymatrix tilling microarray preprocessing
#' @param array.fnames a character vector of array data file names
#' @return a list of expressionMatrix and probe annotation info
#' @author guangchun han
#' @export
affy.st.prep <- function(array.fnames) {
    ##use RMA normalization method,no PMA filtration
    ##may need further threshold on expression value
    ##TODO: we should remove control probes...
    if(!require(oligo,quietly=T) )
        stop("oligo package is not available\n")
    rawData <- read.celfiles(array.fnames)
    cat("rma normalization on core(gene) level, rmaC0 is an expressionSet object\n")
    rmaC0 <- rma(rawData,target="core")
    featureData(rmaC0) <- getNetAffx(rmaC0, "transcript")
    fd <- fData(rmaC0)
    exp <- exprs(rmaC0)

    res <- list(exp=exp,feature=fd$geneassignment)
    res
}
#' nimblegene array preprocessing(experimentail)
#' @param array.fnames a character vector of array data file names
#' @return a list of expressionMatrix and probe annotation info as well as a sdrf file
#' @author guangchun han
#' @export
nimblegene.prep <- function(array.fnames, detect.rate = .2) {
    ##array.fnames is a character vector of array data file names
    if(!require(limma))
        stop("error loading limma package\n")
    if(!require(Ringo))
        stop("error loading Ringo package\n")
    elist <- read.maimages(array.fnames,
                           source="agilent", green.only = T)
    ##elist.p : processed elist
    print("applying background corecction and normexp normalization")
    elist.p <- backgroundCorrect(elist, method = "normexp")
    print("between array normalzation with quantile")
    elist.p <- normalizeBetweenArrays(elist.p, method = "quantile")

    ##low expression probe/control probe remove
    neg95 <- apply(elist.p$E[ elist.p$genes$ControlType == -1, ],
                   2, function(x) quantile(x, p=.95) )
    cutoff <- matrix(1.1*neg95, nrow(elist.p), ncol(elist.p),
                     byrow = T)
    isexpr <- rowSums(elist.p$E > cutoff) >= detect.rate * ncol(elist.p)
    print("number of probes left after applying filtration: ")
    print(table(isexpr))
    ##[, has been defined for EList in limma, El[ 1:3, ]is legal...
    elist.res <- elist.p[ elist.p$genes$ControlType == 0 & isexpr, ]
    res <- list(exp=elist.res$E, feature=elist.res$gene, sdrf=SDRF)
    res
}

#' genepix array preprocessing(experimentail)
#' @param array.fnames a character vector of array data file names
#' @return a list of expressionMatrix and probe annotation info
#' @author guangchun han
genepix.prep <- function(array.fname) {
    ##prepare function for Phalanx one array
    if(!require(limma))
        stop("error loading limma")
    ##limma user guide
    filter.f <- function(spot) as.numeric(splot$Flag > -99)
    RGlist <- read.maimages(array.fname,
                            source = "genepix",
                            wt.fun = filter.f)
    ##add 50 without considering the data specific condition may
    ##be inadvisible
    print("background correction")
    RGlist <- backgroundCorrect(RGlist,
                                method = "normexp",
                                offset = 50)
    print("within array normalization")
    MA <- normalizeWithinArray(RGlist)
    gene.probe <- MA$genes$Status=="Gene"
    MA <- MA[gene.probe,]

}

##list2table
up.or.down <- function(deg, gene) {
    if(length(intersect(deg$up,deg$dn)) != 0)
        stop("up.or.down error: ",gene," presented in both up and
down regulated gene list\n")
    if(gene %in% deg$up)
        return(1)
    if(gene %in% deg$dn)
        return(-1)
    return(0)
}
#' rearrange degs.list to degs table
#'
#' assume that degs.list is a list of degs. each element is a list of length 2:up regulated genes&down regulated genes
#' @export
list2table <- function(degs.list) {
    ##assume that degs.list is a list of degs
    ##each element is a list of length 2:up regulated genes&down regulated genes
    DEGs <- unique(unlist( ## code below for combine up and down genes
        lapply(degs.list, function(x) unique(unlist(x)))) )
    res <- matrix(0, ncol=length(degs.list), nrow = length(DEGs) )
    rownames(res) <- DEGs
    colnames(res) <- names(degs.list)
    for(i in 1:length(DEGs) ) {
        if(i %% 100 == 0)
            print(paste0("dealing with ",i,"th gene"))
        gene <- DEGs[i]
        tmp <- unlist(lapply(degs.list,up.or.down,
                             gene = gene) )
        res[ i, ] <- tmp
        rm(tmp)
    }
    res
}

#' fc.calc for small sample DEG identification based on FC
#' @param dat.obj object returned by dat.prep function
#' @return list with up&down regulated genes expression tables
#' @export
fc.calc <- function(dat.obj) {
    ##dat.ojb: exp, cls(0,1), genesymbol
    print("the experssion data should be log2 transformed")
    dat <- dat.obj$dat
    cls <- factor(dat.obj$cls)
    levels(cls) <- c("ctrl","case")
    gene <- dat.obj$featureName
    dat.t <- t(dat)
    print("FC calculation")
    aggregate(dat.t, by = list(class=cls),
              FUN = median) ->tmp
    rownames(tmp) <- tmp[,1]
    tmp <- tmp[ ,-1 ]
    dat1 <- t(tmp);rm(tmp)
    FC <- 2^(dat1[,1] - dat1[,2])
    res <- data.frame(dat1, gene = gene, FC = FC,
                      stringsAsFactors = F)
    res.up <- res[ FC < 1, ]
    res.dn <- res[ FC > 1, ]
    res.up <- res.up[ order(res.up$FC), ]
    res.dn <- res.dn[ order(res.dn$FC, decreasing = T), ]
    res <- list(Table1 = res.up, Table2 = res.dn)
    attributes(res) <- list(method="fc",
                            names=names(res) )
    res
}

#' readin arrayexpress sdrf file
#'
#' @export
read.sdrf <- function(fname) {
    res <- read.delim(fname,
                      stringsAsFactors = F,
                      check.names = F)
    res
}

#' 'combine' multiple factors 
#'
#' this function will group factors of the same length into a 'combined' one
#' levels of inidividual factors will be combined together
#' @export
#' @return a new factor with levels combined from each individual factor
#' @param ... factors to be passed in
#' @param sep the separater used for combine the factor levels, with '_' as the default
grpFactor <- function(...) {
    f <- list(...)
    if(any(names(f) == 'sep')) {
        idx <- which(names(f) == 'sep')
        if(length(idx) != 1)
            stop('error: multiple separators found')
        sep = f$sep
        f <- f[-idx]	
    } else 
        sep = '_'

    nf <- length(f)
    newf <- as.factor(do.call(paste,c(f,sep=sep)))
    newf
}

#' calculate DEGs considering additional factors
#'
#' ... for additional factors that is used to group samples into more precisable sub-groups
#' DEG will be calculated in each subgroups
#' @export
#' @return a list with class labels(useful for identify up and down regulated genes based on FC) and DEGs calcualted stored in the list
#' @param dat.preped object returned by datPrep function
#' @param pheno.mtr this is the phenotype matrix which stores all phenotypes
#' @param comp.cls.col DEGS will be calculated based on this class
#' @param ... additional factors that is used to group samples into sub-groups

batchedDEG <- function(dat.preped, pheno.mtr, comp.cls.col,...) {
    addt_cols <- list(...)
    iadd <- length(addt_cols) 
    addlist <- list()
    for(i in 1:iadd) {addlist[[i]] = pheno.mtr[,addt_cols[[i]] ]}
    newcls <- as.factor(do.call(paste,c(addlist,sep='_')))
    lvls <- levels(newcls)
    deg.list <- list()
    for(i in 1:length(lvls) ) {
        ilvl <- lvls[i]
        print(paste0('processing ',ilvl))
        newcls == ilvl -> idx
        expi <- dat.preped$exp[ , idx ]
        label <- pheno.mtr[,comp.cls.col][idx]
        labeluniq <- unique(label)
        names(labeluniq) <- c('class1','class2')
        if(length(labeluniq) != 2)
            stop('there should be two levels in the class')
        datPrep(expi,dat.preped$feature, label, labeluniq[1],labeluniq[2])->dat4deg
        diff.res <- diffCalc(dat4deg,diff.method='RP')
        deg.list[[i]] <- list(diffRes=diff.res,class=labeluniq)
        names(deg.list)[i] <- ilvl
        rm(ilvl);rm(idx);rm(expi);rm(label);rm(dat4deg);rm(diff.res)
    }
    deg.list
}

