\name{MRF.condProb}
\alias{MRF.condProb}
\title{calculate the conditional probability of association of a gene given associaiton states and
rewiring of its neighbors}
\usage{
  MRF.condProb(nodei, nodesLabel, netwk, priorPara1, h1,
    delta = 0.95)
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
  the posterior probability of nodei to be associated given
  the configuration of the rest of the network
}
\description{
  calculate the conditional probability of association of a
  gene given associaiton states and rewiring of its
  neighbors
}

