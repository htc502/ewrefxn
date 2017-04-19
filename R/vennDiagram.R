##VennDiagram with color and font set
#' VennDiagam with color and font set
#' 
#' It's based on VennDiagram package with color, fontfamily set.
#' Color scheme was chosen from color brewer. 
#' It takes the same parameter as venn.diagram
#' You can override the paramters by specifiying them again
#' @export
myVennDiagram <- function(...) {
	Args <- list(...)
	if('x' %in% names(Args)) {
	inputList <- Args$x
	} else {
	inputList <- Args[[1]]
	}
	if(!is.list(inputList))
		stop('input x not found')
	nCat <- length(inputList)
	if(nCat < 3 | nCat > 8) stop('error: # categories should be >= 3 and <= 8')
	if(!require(RColorBrewer))
		stop('error loading RColorBrewer')
	fill <- brewer.pal(nCat,'Set2')
	myArgs <- list(	lwd=0,
			fill=fill,
			main.fontfamily='sans',
			cat.fontfamily='sans',
			fontfamily='sans',
			alpha=rep(.6,length(inputList)))
##remove those args respecified in myargs
	dupArgs <- intersect(names(Args),names(myArgs))
	if(length(dupArgs) != 0) {
		pos <- match(dupArgs, names(myArgs))
			myArgs <- myArgs[-pos]
	}
	newArgs <- c(Args, myArgs)
	if(!require(VennDiagram))
		stop('error loading VennDiagram')
	do.call(venn.diagram,newArgs)
}
