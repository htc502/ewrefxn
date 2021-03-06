% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tbl2hmap.R
\name{tbl2hmap}
\alias{tbl2hmap}
\title{convert long tibble to pheatmap data input}
\usage{
tbl2hmap(
  tbl.long,
  rowVar,
  colVar,
  valueVar,
  colAnnVars = NULL,
  rowAnnVars = NULL
)
}
\arguments{
\item{tbl.long}{tibble in long format}

\item{rowVar}{column name used as row for pheatmap matrix}

\item{colVar}{column name used as column for pheatmap matrix}

\item{valueVar}{column name used as the values for pheatmap matrix}

\item{colAnnVars}{column name(s) used for additonal column annotation}

\item{rowAnnVars}{column name(s) used for additonal row annotation}
}
\value{
a list with 3 components: mat, rowAnn, colAnn
}
\description{
convert long tibble to pheatmap data input
}
