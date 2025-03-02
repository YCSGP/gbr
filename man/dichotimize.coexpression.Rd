\name{dichotimize.coexpression}
\alias{dichotimize.coexpression}
\title{build the network for Markov random filed, edges are co-expression in disease/healthy condition,
weighted by rewiring degree}
\usage{
  dichotimize.coexpression(gwas.exp.data)
}
\arguments{
  \item{gwas.exp.data}{is a list returned by the function
  overlay.exp.gwas()}
}
\value{
  a list with network information, including rewiring
  matrix, neighbor list, and an edge matrix
}
\description{
  build the network for Markov random filed, edges are
  co-expression in disease/healthy condition, weighted by
  rewiring degree
}

