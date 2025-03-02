\name{updateP}
\alias{updateP}
\title{calculate the posterior probability of Markov random field model}
\usage{
  updateP(netwk, nodesLabel, priorPara1 = c(0.01, 0.01),
    h1)
}
\arguments{
  \item{netwk}{is the output of function
  dichotimize.coexpression()}

  \item{nodesLabel}{is a vector of +1/-1, which is a
  configuration of the network}

  \item{priorPara1}{is the hyper-parameter in Markov Random
  Filed model, prefixed as (0.01,0.01) if not otherwisely
  indicated}

  \item{h1}{is a numeric value, taken from estimate.h1()}
}
\value{
  the posterior probability
}
\description{
  calculate the posterior probability of Markov random
  field model
}

