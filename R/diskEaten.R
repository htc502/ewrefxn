##Du for diskusage status
getDu <- function(dir) {
	cmdGetDiskUsage <- paste0('du -d1 ',dir)
		system(cmdGetDiskUsage,intern=T) -> cmdOutput
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
		print('total disk usage')
		print(paste0(parseDuRes$total$dir,' : ',parseDuRes$total$nBlk))
		print('--------------------\n')
		print('extremely large subdirs&files')
		print(paste(parseDuRes$contents$dir[outliers],' : ',parseDuRes$contents$nBlk[outliers]))
		print('--------------------')
}

diskEaten <- function(dir) {
	duRes <- getDu(dir)
		parseDu(duRes) -> pDuRes
		findExtreme(pDuRes)
}


