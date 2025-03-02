\name{graphPrior}
\alias{graphPrior}
\title{calculate the potential of the network configuration}
\usage{
  graphPrior(netwk, nodesLabel, priorPara1 = c(0.01, 0.01),
    delta = 0.95, h1)
}
\arguments{
  \item{nodei}{is a integer, indicating the current node}

  \item{nodesLabel}{is a vector of +1/-1, which is a
  configuration of the network}

  \item{netwk}{is the output of function
  dichotimize.coexpression()}

  \item{h1}{is a numeric value, taken from estimate.h1()}
}
\value{
  the potential of the configuration of the network
}
\description{
  calculate the potential of the network configuration
}

