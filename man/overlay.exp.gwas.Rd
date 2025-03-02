\name{overlay.exp.gwas}
\alias{overlay.exp.gwas}
\title{Get the intersect of genes that are in microarray and in GWAS}
\usage{
  overlay.exp.gwas(exp.meta, gwas)
}
\arguments{
  \item{exp.meta}{is a list returned by the function
  read.expression()}

  \item{gwas}{is dataframe returned by the function
  read.gwas()}
}
\value{
  a list, elementss are GWAS data, GeneID, disease
  expression data, and control expression data
}
\description{
  Get the intersect of genes that are in microarray and in
  GWAS
}

