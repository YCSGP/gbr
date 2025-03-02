\name{likpw.pval}
\alias{likpw.pval}
\title{calculate likelihood of GWAS data given association state}
\usage{
  likpw.pval(nodesData, nodesLabel, priorPara2 = NULL)
}
\arguments{
  \item{nodesData}{is a vector of "pval" objects, derived
  form pval.class()}

  \item{nodesLabel}{is a vector of association state, 1
  means association, and -1 means}

  \item{priorPara2}{is the hyperparameter of the p-value
  prior distribution}
}
\value{
  a vector of likelihoods of GWAS p-values
}
\description{
  calculate likelihood of GWAS data given association state
}

