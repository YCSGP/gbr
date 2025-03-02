\name{read.expression}
\alias{read.expression}
\title{read in the gene expression datasets of disease samples and healthy samples seperately}
\usage{
  read.expression(disease.name, healthy.name)
}
\arguments{
  \item{disease.name}{is a string, the file name of the
  disease gene expression dataset}

  \item{healthy.name}{is a string, the file name of the
  healthy gene expression dataset}
}
\value{
  a list with three components, a vector of gene names in
  microarray, the microarray matrix of disease samples, and
  the microarray matrix of healthy samples
}
\description{
  read in the gene expression datasets of disease samples
  and healthy samples seperately
}

