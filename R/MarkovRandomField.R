#' emprically determining the hyper-parameter h1
#'
#' @param netwk is the output of function dichotimize.coexpression()
#' @param priorPara1 is the hyper-parameter in Markov Random Filed model, prefixed
#' as (0.01,0.01) if not otherwisely indicated
#' @return a numeric value, which is the emprical choice of h1
estimate.h1<-function(netwk, priorPara1=c(0.01,0.01),q1=0.9,delta=0.95){
  geneid<-netwk[["GeneID"]]
  gwas<-netwk[["GWAS"]]
  tau11<-priorPara1[1]
  tau22<-priorPara1[2]
  neighbors<-netwk[["NeighborList"]]
  nNodes<-length(geneid)
  genesymbol<-1:nNodes
  nodesLabel<-ifelse(gwas[,1]<=0.05,1,-1)
  pdiff<-netwk[["RewireMatrix"]]
  dpot<-lapply(1:nNodes,function(nodei){
      i.neighbor <- neighbors[[nodei]]
      i.neighbor.1<-is.element(genesymbol,i.neighbor)&nodesLabel==1
      i.neighbor.2<-is.element(genesymbol,i.neighbor)&nodesLabel==-1&pdiff[nodei,]>delta
      logitp<-tau11*sum(pdiff[nodei,]*i.neighbor.1)+tau22*sum(pdiff[nodei,]*i.neighbor.2)
      logitp
  })
  dpot<-unlist(dpot)
  h1<-as.numeric(quantile(dpot,q1)*(-1))
  h1
}

#' calculate the conditional probability of association of a gene given associaiton states and 
#' rewiring of its neighbors
#'
#' @param nodei is a integer, indicating the current node
#' @param nodesLabel is a vector of +1/-1, which is a configuration of the network
#' @param netwk is the output of function dichotimize.coexpression()
#' @param h1 is a numeric value, taken from estimate.h1() 
#' @return the posterior probability of nodei to be associated given the configuration of the 
#' rest of the network
MRF.condProb <- function(nodei,nodesLabel,netwk,priorPara1,h1,delta=0.95)
{
  neighbors<-netwk[["NeighborList"]]
  genesymbol<-1:length(netwk[["GeneID"]])
  pdiff<-netwk[["RewireMatrix"]]
  tau11<-priorPara1[1]
  tau22<-priorPara1[2]
  i.neighbor <- neighbors[[nodei]]
  i.neighbor.1<-is.element(genesymbol,i.neighbor)&nodesLabel==1
  i.neighbor.2<-is.element(genesymbol,i.neighbor)&nodesLabel==-1&pdiff[nodei,]>delta
  logitp<-(h1+tau11*sum(pdiff[nodei,]*i.neighbor.1)+
  tau22*sum(pdiff[nodei,]*i.neighbor.2))
  1/(1+exp(-logitp))
}

#' calculate the potential of the network configuration
#'
#' @param nodei is a integer, indicating the current node
#' @param nodesLabel is a vector of +1/-1, which is a configuration of the network
#' @param netwk is the output of function dichotimize.coexpression()
#' @param h1 is a numeric value, taken from estimate.h1()
#' @return the potential of the configuration of the network
graphPrior <- function(netwk,nodesLabel,priorPara1=c(0.01,0.01),delta=0.95,h1)
{
  tau11<-priorPara1[1]
  tau22<-priorPara1[2]
  edges<-netwk[["EdgeVector"]]
  edge11 <- nodesLabel[edges[,1]] == 1 & nodesLabel[edges[,2]]== 1
  edge22<-nodesLabel[edges[,1]]==-1&nodesLabel[edges[,2]]==-1&edges[,3]>delta
  U<-h1*sum(nodesLabel==1)-tau22*sum(edges[edge22,3])
  U
}

#' calculate the likelihood of observed GWAS p-values given a network configuration
#'
#' @param netwk is the output of function dichotimize.coexpression()
#' @param nodesLabel is a vector of +1/-1, which is a configuration of the network
#' @param priorPara2 is the hyperparameter of the p-value prior distribution, fixed as (3,1,10,1)
#' if not otherwisely indicated
#' @return the log-likelihood of GWAS p-values
ichip.llk <- function(netwk,nodesLabel,priorPara2=c(3,1,10,1))
{
  y<-netwk[["GWAS"]][,2]
  mubar <- priorPara2[1]
  a <- priorPara2[2]
  nu <- priorPara2[3]
  lambda <- priorPara2[4]
  scale0 <- sqrt(lambda*(a+1)/a)
  lk<-ifelse(nodesLabel==1, dt((y-mubar)/scale0,nu)/scale0,dnorm(y,0,1) )
  sum(log(lk))
}

#' calculate the pseudo-likelihood of the Markov Random field model
#'
#' @param netwk is the output of function dichotimize.coexpression()
#' @param nodesLabel is a vector of +1/-1, which is a configuration of the network
#' @param priorPara2 is the hyperparameter of the p-value prior distribution, fixed as (3,1,10,1)
#' if not otherwisely indicated
#' @return a numerical value, which is the pseudo-likelihood of the Markov Random field model
pseudo.llk<-function(netwk,nodesLabel,priorPara1=c(0.01,0.01),priorPara2=c(3,1,10,1),h1)
{
ichip.llk(netwk,nodesLabel)+graphPrior(netwk,nodesLabel,priorPara1,0.95,h1)
}

#' calculate the posterior probability of Markov random field model
#'
#' @param netwk is the output of function dichotimize.coexpression()
#' @param nodesLabel is a vector of +1/-1, which is a configuration of the network
#' @param priorPara1 is the hyper-parameter in Markov Random Filed model, prefixed
#' as (0.01,0.01) if not otherwisely indicated
#' @param h1 is a numeric value, taken from estimate.h1()
#' @return the posterior probability
updateP<-function(netwk,nodesLabel,priorPara1=c(0.01,0.01),h1)
{
  tau11<-priorPara1[1]
  tau22<-priorPara1[2]
  nNodes<-length(netwk[["GeneID"]])
  likRatio<-netwk[["GWAS"]][,2]
  prob<-numeric(nNodes)
  for(nodeInd in 1:nNodes )
  {
    condProb<-MRF.condProb(nodeInd, nodesLabel,netwk,priorPara1,h1)
    prob[nodeInd]<-(likRatio[nodeInd]*condProb)/
    (likRatio[nodeInd]*condProb+(1-condProb))
  }
  prob
}

#' ICM algorithm
#'
#' @param netwk is the output of function dichotimize.coexpression()
#' @param priorPara1 is the hyper-parameter in Markov Random Filed model, prefixed
#' as (0.01,0.01) if not otherwisely indicated
#' @param priorPara2 is the hyperparameter of the p-value prior distribution, fixed as (3,1,10,1)
#' @return a dataframe, columns are GeneID, GWAS p-values, posterior probability of MRF model
#' @export
icmicm<-function(netwk,priorPara1=c(0.01,0.01),priorPara2=c(3,1,10,1),delta=0.95,q1=0.9){
  h1<-estimate.h1(netwk)
  nNodes<-length(netwk[["GeneID"]])
  initLabel<-NULL
  currentNodesLabel <- if(is.null(initLabel)) 2*rbinom(nNodes,1,0.5)-1
  likRatio<-netwk[["GWAS"]][,2]
  prevNodesLabel<-currentNodesLabel
  iteNum<-1;iteT<-50;converged<-F
  randsearch <- function(initNode,currentNodesLabel)
  {
    temp.label<-currentNodesLabel
    for(nodei in sample(nNodes) )
    {
      condProb<-MRF.condProb(nodei,temp.label,netwk,priorPara1,h1)
      prob<-(likRatio[nodei]*condProb)/(likRatio[nodei]*condProb+(1-condProb))
      temp.label[nodei]<-sign(prob-0.5)
    }
    temp.label
  }
  while (iteNum<iteT & !converged)
  {
    prevNodesLabel<-currentNodesLabel
    currentNodesLabel<-randsearch(sample(1:nNodes,1),currentNodesLabel)
    if(all(currentNodesLabel==prevNodesLabel)) {
      converged<-T
    }
    a<-pseudo.llk(netwk,prevNodesLabel,priorPara1,priorPara2,h1)
    b<-pseudo.llk(netwk,currentNodesLabel,priorPara1,priorPara2,h1)
    if(a>b) {
      currentNodesLabel<-prevNodesLabel
      converged<-T
    }
    iteNum <- iteNum+1
  }
  post.gwas<-updateP(netwk,currentNodesLabel,priorPara1,h1)
  outout<-data.frame(netwk[["GeneID"]],netwk[["GWAS"]][,1],1-post.gwas)
  colnames(outout)<-c("GeneID","GWAS","prioritization")
  outout  
}








