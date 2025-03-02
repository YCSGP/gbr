\name{read.gwas}
\alias{read.gwas}
\title{read in GWAS data, and the output is taken as input of the Markov random field approach}
\usage{
  read.gwas(gwas.name, priorPara2 = c(3, 1, 10, 1))
}
\arguments{
  \item{gwas.name}{is string, the name of the file where
  GWAS data is}

  \item{priorPara2}{is the hyperparameter of the p-value
  prior distribution, fixed as (3,1,10,1) if not
  otherwisely indicated}
}
\value{
  a data frame, the columns are gene name, p-value, and the
  LLR
}
\description{
  read in GWAS data, and the output is taken as input of
  the Markov random field approach
}

