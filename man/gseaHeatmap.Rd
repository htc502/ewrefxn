% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsea-heatmap.R
\name{gseaHeatmap}
\alias{gseaHeatmap}
\title{Make heatmap based on gsea enrichment results}
\usage{
gseaHeatmap(
  wd = wd,
  showPathway = NULL,
  myPalette = NULL,
  colOrder = NULL,
  rowOrder = NULL
)
}
\arguments{
\item{wd}{the directory which contains gsea result *-pos|neg.xls (actually txt) files.}

\item{showPathway}{a vector contains the pathways want to show}

\item{myPalette}{an optional palette function used to color the heatmap. default is brewer.pal RdYlBu.}

\item{colOrder}{an optional vector specifying the column order}

\item{rowOrder}{an optional vector specifying the row order}
}
\value{
ggplot2 object for downstream use
}
\description{
Make heatmap based on gsea enrichment results
}
