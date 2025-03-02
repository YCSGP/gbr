\name{ichip.llk}
\alias{ichip.llk}
\title{calculate the likelihood of observed GWAS p-values given a network configuration}
\usage{
  ichip.llk(netwk, nodesLabel, priorPara2 = c(3, 1, 10, 1))
}
\arguments{
  \item{netwk}{is the output of function
  dichotimize.coexpression()}

  \item{nodesLabel}{is a vector of +1/-1, which is a
  configuration of the network}

  \item{priorPara2}{is the hyperparameter of the p-value
  prior distribution, fixed as (3,1,10,1) if not
  otherwisely indicated}
}
\value{
  the log-likelihood of GWAS p-values
}
\description{
  calculate the likelihood of observed GWAS p-values given
  a network configuration
}

