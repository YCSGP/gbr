#' read in the gene expression datasets of disease samples and healthy samples seperately
#' 
#' @param disease.name is a string, the file name of the disease gene expression dataset
#' @param healthy.name is a string, the file name of the healthy gene expression dataset
#' @return a list with three components, a vector of gene names in microarray,
#' the microarray matrix of disease samples, and the microarray matrix of healthy samples
#' @export
read.expression<-function(disease.name,healthy.name)
{
  disease.exp<-read.table(file=disease.name,header=T,row.names=1,sep='\t')
  healthy.exp<-read.table(file=healthy.name,header=T,row.names=1,sep='\t')
  g.microarray<-rownames(disease.exp)
  exp.meta<-list(g.microarray,disease.exp,healthy.exp)
}

#' Get the intersect of genes that are in microarray and in GWAS
#'
#' @param exp.meta is a list returned by the function read.expression()
#' @param gwas is dataframe returned by the function read.gwas()
#' @export
overlay.exp.gwas<-function(exp.meta,gwas)
{
  g.microarray<-exp.meta[[1]]
  g.gwas<-rownames(gwas)
  g<-intersect(g.microarray,g.gwas)
  gwas<-gwas[g,]
  disease.exp<-exp.meta[[2]][g,]
  healthy.exp<-exp.meta[[3]][g,]
  gwas.exp.meta<-list(gwas,g,disease.exp,healthy.exp)
  return(gwas.exp.meta)
}
