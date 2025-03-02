#' calculate rewiring information 
#'
#' @param gwas.exp.data is a list, returned by the function overlay.exp.gwas()
#' @return a list of three matrix, the co-expression matrix of disease samples,
#' the co-expression matrix of healthy samples, and the rewiring matrix
#' between disease and healthy subjects
#' @export
get.network<-function(gwas.exp.data){
  disease.exp<-gwas.exp.data[[3]]
  healthy.exp<-gwas.exp.data[[4]]
  n1<-ncol(disease.exp)
  n2<-ncol(healthy.exp)
  disease.pcc<-cor(t(disease.exp),use="pair"); diag(disease.pcc)<-0
  healthy.pcc<-cor(t(healthy.exp),use="pair"); diag(healthy.pcc)<-0
  FisherTrans<-function(r){log((1+r)/(1-r))*1/2}
  z1<-FisherTrans(disease.pcc)
  z2<-FisherTrans(healthy.pcc)
  z<-(z1-z2)/sqrt(1/(n1-3)+1/(n2-3))
  rewire<-2*pnorm(abs(z))-1
  g<-gwas.exp.data[[2]]
  rownames(disease.pcc)<-g; colnames(disease.pcc)<-g
  rownames(healthy.pcc)<-g; colnames(healthy.pcc)<-g
  rownames(rewire)<-g; colnames(rewire)<-g
  network.data<-list(disease.pcc,healthy.pcc,rewire)
  network.data
}

#' build the network for Markov random filed, edges are co-expression in disease/healthy condition,
#' weighted by rewiring degree
#'
#' @param gwas.exp.data is a list returned by the function overlay.exp.gwas()
#' @return a list with network information, including rewiring matrix, neighbor list, and 
#' an edge matrix 
#' @export
dichotimize.coexpression<-function(gwas.exp.data){
  network.data<-get.network(gwas.exp.data)
  avec<-seq(0.3,0.9,by=0.1)
  gwas<-gwas.exp.data[[1]]
  disease.pcc<-network.data[[1]]
  healthy.pcc<-network.data[[2]]
  get.r2<-function(A){
    degree<-unlist( lapply(1:nrow(healthy.pcc),function(ii){
                           sum(abs(healthy.pcc[ii,])>=A|abs(disease.pcc[ii,])>=A)}
                          ) 
                  )
    degree<-degree[degree>0]
    if(length(degree)<=50) return(c(0,0,0))
    tdeg<-table(degree)
    y<-as.numeric(tdeg)/length(degree)
    x<-as.numeric(names(tdeg))
    lfit<-lm(log(y)~log(x))
    b<-as.numeric(lfit$coefficients[2])
    p<-anova(lfit, test="F")$"Pr(>F)"[1]
    r2<-summary(lfit)$r.squared
    c(b,p,r2)
  }
  fit.res<-matrix(unlist(lapply(avec,"get.r2")),byrow=T,ncol=3)
  A<-avec[min(which(fit.res[,1]<0&fit.res[,2]<=0.01))]
  netwk<-matrix( 0,ncol=ncol(disease.pcc),nrow=nrow(disease.pcc) )
  g<-gwas.exp.data[[2]]
  rownames(netwk)<-g
  colnames(netwk)<-g
  netwk[abs(disease.pcc)>=A|abs(healthy.pcc)>=A]<-1
  degree<-rowSums(netwk)
  gmrf<-g[degree>0]
  rewire<-network.data[[3]][gmrf,gmrf]
  netwk<-netwk[gmrf,gmrf]
  netwk1<-netwk
  netwk1[lower.tri(netwk1)]<-0
  nNodes<-length(gmrf)
  genesymbol<-c(1:nNodes)
  gwas<-gwas[gmrf,]
  edges <- expand.grid(genesymbol,genesymbol)[as.vector(netwk1)==1,]
  NEI<-function(ii){which(netwk[ii,]==1)}
  neighbors<-lapply(1:nNodes,"NEI")
  names(neighbors)<-as.character(1:nNodes)
  edge.weight<-as.vector(rewire)[as.vector(netwk1)==1]
  edges<-cbind(edges,edge.weight)
  network.bundle<-list(gmrf,rewire,edges,neighbors,gwas)
  names(network.bundle)<-c("GeneID","RewireMatrix","EdgeVector","NeighborList","GWAS")
  network.bundle
}  



