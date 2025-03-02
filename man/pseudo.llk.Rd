\name{pseudo.llk}
\alias{pseudo.llk}
\title{calculate the pseudo-likelihood of the Markov Random field model}
\usage{
  pseudo.llk(netwk, nodesLabel, priorPara1 = c(0.01, 0.01),
    priorPara2 = c(3, 1, 10, 1), h1)
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
  a numerical value, which is the pseudo-likelihood of the
  Markov Random field model
}
\description{
  calculate the pseudo-likelihood of the Markov Random
  field model
}

