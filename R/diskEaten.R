##Du for diskusage status
getDu <- function(dir) {
    cmdGetDiskUsage <- paste0('du -d1 -m ',dir)
    system(cmdGetDiskUsage,intern=T) -> cmdOutput
    list(cmdOutput = cmdOutput, dir=dir)
}
getDuWin <- function(dir) {
    cmdGetDiskUsage <- paste0('du -l 1 -c ',dir)
    system(cmdGetDiskUsage,intern = T) -> cmdOutput
    cmdOutput <- cmdOutput[-1]
    ##tidy work
    gsub('\"','',cmdOutput) -> cmdOutput
    gsub('\\\\','\\/',cmdOutput) -> cmdOutput

    lapply(cmdOutput, function(x) {strsplit(x,',')[[1]][c(6,1)]}) -> cmdOutput
    ##the last row corresponds to the topDir,change it for compatibility
    cmdOutput[[length(cmdOutput)]][2] <- dir

    unlist(lapply(cmdOutput, function(x) paste(x,collapse='\t'))) -> cmdOutput
    list(cmdOutput = cmdOutput, dir=dir)
}
parseDu <- function(getDuRes) {
    topDir <- getDuRes$dir
    getDuRes <- getDuRes$cmdOutput
    tmp <- strsplit(getDuRes, split='\t')
    nBlk <- as.numeric(unlist(lapply(tmp,function(x) x[1])))
    dir <- unlist(lapply(tmp,function(x) x[2]))
    idxTopDir <- which(dir == topDir)
    total <- list(dir = dir[idxTopDir],nBlk=nBlk[idxTopDir])
    list(dir=dir[-idxTopDir],nBlk=nBlk[-idxTopDir]) -> contents
    list(total = total,contents = contents)
}
outLiers <- function(values) {
    qts <- quantile(values,c(.25,.5,.75))
    iqr <- abs(qts[1]-qts[3])
    idxUpper <- values > qts[2];idxLower <- values < qts[2]
    idxUpperOutLier <- (values - qts[3])/iqr >= 2
    idxLowerOutLier <- (values - qts[1])/iqr >= 2
    which(idxUpper & idxUpperOutLier) -> upperOutliers
    which(idxLower & idxLowerOutLier) -> lowerOutliers
    c(upperOutliers, lowerOutliers)
}

findExtreme <- function(parseDuRes) {
    outliers <- outLiers(parseDuRes$contents$nBlk)
    if(length(outliers) == 0) {
        print('there is no extraordinary directories')
        return
    }
    print('total disk usage')
    print(paste0(parseDuRes$total$dir,' : ',parseDuRes$total$nBlk))
    print('--------------------')
    print('extremely large subdirs&files')
    print(paste(parseDuRes$contents$dir[outliers],' : ',parseDuRes$contents$nBlk[outliers],'MB'))
    print('--------------------')
}

#' find the extraordinary large directories under the given dir
#'
#' to use this function under windows, download the du.exe binary from
#' https://technet.microsoft.com/en-us/sysinternals/du.aspx and  put it into PATH.
#' As this function use different underlying system API on win and linux,
#' the disk usage amount reported may not be the same.
#' @param dir the dir under which the contents will be checked
#' @return this funciton return nothing, it prints out the result
#' @export
diskEaten <- function(dir) {
    sysinfo <- Sys.info()
    if(sysinfo['sysname'] == 'Windows')
        duRes <- getDuWin(dir)
    else
        duRes <- getDu(dir)
    parseDu(duRes) -> pDuRes
    findExtreme(pDuRes)
}


