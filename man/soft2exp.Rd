% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo.R
\name{soft2exp}
\alias{soft2exp}
\title{expression matrix generator}
\usage{
soft2exp(
  AccNum,
  filename = NULL,
  exp.col = 2,
  symbol.col,
  entrezID.col,
  phe.col,
  log2 = F
)
}
\arguments{
\item{AccNum}{accession number}

\item{filename}{read data from the local data instead of the remote server}

\item{exp.col}{column number of expression value}

\item{symbol.col}{column number of expression value, for arraymatrix data the GPL header is usually SYMBOL;
for agilent44K it is GENE_SYMBOL(10th column of GPL4134)}

\item{entrezID.col}{column number of entrez id, for affymetrix data it is named Entrez_Gene_ID, for agilent it is
GENE(9th column of GPL4134)}

\item{phe.col}{specify the phenotype information want to extract from the gsmlist}

\item{log2}{set to T for logrithm transformation}
}
\value{
a list with expmatrix probeID symbols and entrezID
}
\description{
given an accession number(gseXXXX), this fxn download and extract the expression matrix
}
