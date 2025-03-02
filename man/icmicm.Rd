\name{icmicm}
\alias{icmicm}
\title{ICM algorithm}
\usage{
  icmicm(netwk, priorPara1 = c(0.01, 0.01),
    priorPara2 = c(3, 1, 10, 1), delta = 0.95, q1 = 0.9)
}
\arguments{
  \item{netwk}{is the output of function
  dichotimize.coexpression()}

  \item{priorPara1}{is the hyper-parameter in Markov Random
  Filed model, prefixed as (0.01,0.01) if not otherwisely
  indicated}

  \item{priorPara2}{is the hyperparameter of the p-value
  prior distribution, fixed as (3,1,10,1)}
}
\value{
  a dataframe, columns are GeneID, GWAS p-values, posterior
  probability of MRF model
}
\description{
  ICM algorithm
}

