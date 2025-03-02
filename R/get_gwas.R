#' transfer GWAS data
#'  
#' @param pv is a vector of p-values from GWAS studies
#' @return an element of class "pval"
pval.class <- function(pv) {
  pv[pv<=.Machine$double.eps] <- .Machine$double.eps
  nodeD <- pv
  class(nodeD) <- "pval"
  attr(nodeD,'p') <- pv
  attr(nodeD,'y') <-  qnorm(1-pv);
  nodeD
}

#' calculate likelihood of GWAS data given association state
#' 
#' @param nodesData is a vector of "pval" objects, derived form pval.class()
#' @param nodesLabel is a vector of association state, 1 means association, and -1 means 
#  non-association
#' @param priorPara2 is the hyperparameter of the p-value prior distribution
#' @return a vector of likelihoods of GWAS p-values
likpw.pval <- function( nodesData, nodesLabel, priorPara2=NULL) {
  y <- attr(nodesData,'y')
  mubar <- priorPara2[1]
  a <- priorPara2[2]
  nu <- priorPara2[3]
  lambda <- priorPara2[4]
  scale0 <- sqrt(lambda*(a+1)/a)
  ifelse(nodesLabel==1, dt((y-mubar)/scale0,nu)/scale0  , dnorm(y,0,1) )
}

#' read in GWAS data, and the output is taken as input of the Markov random field approach
#'
#' @param gwas.name is string, the name of the file where GWAS data is
#' @param priorPara2 is the hyperparameter of the p-value prior distribution, fixed as (3,1,10,1)
#' if not otherwisely indicated
#' @return a data frame, the columns are gene name, p-value, and the LLR
#' @export
read.gwas<-function(gwas.name,priorPara2=c(3,1,10,1)){
  pgwas<-read.table(file=gwas.name,row.names=1,header=T)
  gname<-rownames(pgwas)
  p<-pgwas[,1]
  names(p)<-gname
  p<-p[p<=0.5]
  gname<-names(p)
  obsAsso<-pval.class(p)
  NN<-length(p)
  lik.1<-likpw.pval(obsAsso, rep(1,NN), priorPara2)
  lik.0<-likpw.pval(obsAsso, rep(-1,NN))
  likRatio<-lik.1/lik.0
  gwas<-data.frame(p,likRatio)
  names(gwas)<-c("pvalue","LLR")
  rownames(gwas)<-gname
  return(gwas)
}
