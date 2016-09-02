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
		myArgs <- list(	lwd=0,
				fill=c('#1b9e77','#d95f02','#7570b3'),
				main.fontfamily='sans',
				cat.fontfamily='sans',
				fontfamily='sans',
				alpha=rep(.6,3))
##remove those args respecified in myargs
		dupArgs <- intersect(names(Args),names(myArgs))
		if(length(dupArgs) != 0) {
			pos <- match(dupArgs, names(myArgs))
				myArgs <- myArgs[-pos]
		}
	newArgs <- c(myArgs, Args)
		do.call(venn.diagram,newArgs)
}
