#' Make heatmap based on gsea enrichment results
#'
#' 
#' @param wd the directory which contains gsea result *-pos|neg.xls (actually txt) files.
#' @param showPathway a vector contains the pathways want to show
#' @param myPalette an optional palette function used to color the heatmap. default is brewer.pal RdYlBu.
#' @param colOrder an optional vector specifying the column order
#' @param rowOrder an optional vector specifying the row order
#' @return ggplot2 object for downstream use
#' @export

gseaHeatmap = function(wd = wd, showPathway = NULL, myPalette = NULL, colOrder = NULL,rowOrder = NULL)
{
  if(!require(readr) | ! require(ggplot2)) stop("error loading readr|ggplot2 package")
  gseaFiles = list.files(path = wd, pattern = ".*-(pos|neg).xls",full.names = T)

  gseaMat = c()
  for(f in gseaFiles) {
    fdat = read_tsv(f) %>%
      select(NAME,NES,`FDR q-val`)
    fdat$File = rep(f,nrow(fdat))
    gseaMat = rbind(gseaMat,
                    fdat)
  }
  gseaMat = gseaMat %>%
    mutate(group = gsub("-pos.xls|-neg.xls","",
                        basename(File))) %>%
    dplyr::rename(
             Pathway = NAME,
             q = "FDR q-val"
           )
  gdat = gseaMat %>%
    select(Pathway,
           NES,
           q,
           group)
  if(is.null(showPathway)) showPathway =  unique( filter(gdat,q < 0.05)$Pathway )
  if(is.null(myPalette)) myPalette = colorRampPalette(rev(brewer.pal(11, "RdYlBu")), space="Lab") 
  gdat$q[gdat$q == 0] = 1e-5 ##avoid inf
  gdat = filter(gdat,
                Pathway %in% showPathway) %>%
    mutate(q1 = ifelse(q >= 0.1, 1, q))

  if(!is.null(colOrder)) gdat$group = factor(gdat$group,
                                             levels = colOrder)

  if(!is.null(rowOrder)) gdat$Pathway = factor(gdat$Pathway,
                                               levels = rowOrder)

  gp = ggplot(data = gdat)
  gp = gp + geom_tile(aes(x=group,y=Pathway),fill='grey90',color='black')
  gp = gp + geom_point(aes(x=group,y=Pathway,
                           fill = NES,
                           size=-log10(q1)),
                       shape = 21,
                       color = 'black',
                       stroke = .2)
  gp = gp + scale_fill_gradientn(colours = myPalette(100))
  gp = gp + scale_size(range=c(0,7),breaks = c(0,1,2,4,5))
  gp = gp + theme(axis.text.x = element_text(angle=90,colour='black',hjust = .5,vjust = .5),
                  axis.text.y = element_text(colour='black'),
                  ##                  axis.text.y = element_blank(),
                  panel.background = element_blank(),
                  axis.ticks = element_blank(),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  gp = gp + ylab('') + xlab('')
  gp
}
