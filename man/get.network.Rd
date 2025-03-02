\name{get.network}
\alias{get.network}
\title{calculate rewiring information}
\usage{
  get.network(gwas.exp.data)
}
\arguments{
  \item{gwas.exp.data}{is a list, returned by the function
  overlay.exp.gwas()}
}
\value{
  a list of three matrix, the co-expression matrix of
  disease samples, the co-expression matrix of healthy
  samples, and the rewiring matrix between disease and
  healthy subjects
}
\description{
  calculate rewiring information
}

