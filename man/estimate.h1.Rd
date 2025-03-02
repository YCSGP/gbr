\name{estimate.h1}
\alias{estimate.h1}
\title{emprically determining the hyper-parameter h1}
\usage{
  estimate.h1(netwk, priorPara1 = c(0.01, 0.01), q1 = 0.9,
    delta = 0.95)
}
\arguments{
  \item{netwk}{is the output of function
  dichotimize.coexpression()}

  \item{priorPara1}{is the hyper-parameter in Markov Random
  Filed model, prefixed as (0.01,0.01) if not otherwisely
  indicated}
}
\value{
  a numeric value, which is the emprical choice of h1
}
\description{
  emprically determining the hyper-parameter h1
}

