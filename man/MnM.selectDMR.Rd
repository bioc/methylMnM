\name{MnM.selectDMR}
\alias{MnM.selectDMR}
\title{
Select significants of each comparation.
}
\description{
The function is used to select significants of each comparation.
}
\usage{
MnM.selectDMR(frames = NULL, up =1.45, down = 1/1.45, p.value.MM = 0.01,
 p.value.SAGE = 0.01,q.value = 0.01,cutoff="q-value",  quant= 0.6)
  }
\arguments{
  \item{frames}{
The input dataset and q-values of each bin.
}
  \item{up}{
The ratio of Medip1/Medip2 should be larger than "up" value if we call it significant.
}
  \item{down}{
The ratio of Medip1/Medip2 should be smaller than "down" value if we call it significant.
}
  \item{p.value.MM}{
The p-value of the bin which use MM test should be smaller than "p.value.MM" if we call it significant.
}
  \item{p.value.SAGE}{
The p-value of the bin which use SAGE test should be smaller than "p.value.SAGE" if we call it significant.
}

  \item{q.value}{
The q-value of the bin should be smaller than "q.value" if we call it significant.
}
  \item{cutoff}{
We should use p-value or q-value cutoff to detect DMRs (If cutoff="q-value", then we use q-value to detect DMRs, else we use p-value  ).
}
  \item{quant}{
The rank of absolute value of (Medip1-Medip2) should be larger than "quant" value if we call it significant.
}
}
\value{
The DMRs of the comparation. 
}

\author{
Yan Zhou, Bo Zhang, Nan Lin, BaoXue Zhang and Ting Wang
}
\examples{
  datafile<-system.file("extdata",  package = "methylMnM")
  filepath<-datafile[1]
  file1<-paste(filepath,"/all_CpGsite_chr18.txt",sep="")
  CpGsite<-read.table(file1, header=FALSE,skip=0, nrows=200, as.is=TRUE) 
  winbin<-CpGsite[1:100,1:4]
  winbin[,2]<-seq(0,49500,500)
  winbin[,3]<-winbin[,2]+500
  count<-matrix(rpois(600, lambda=5), nrow=100 )
  count[,6]<-count[,5]
  pvalue<-runif(100, min=0, max=1)
  ts<-rnorm(100, mean=0, sd=1) 
  cpgpq<-cbind(winbin[,1:3],count,pvalue,ts)
  colnames(cpgpq)=c("chr", "chrSt","chrEnd","Medip1","Medip2","MRE1",
  "MRE2","cg","mrecg","pvalue",'Ts')
  f<-MnM.selectDMR(frames=cpgpq, p.value.MM = 0.1, 
  p.value.SAGE = 0.1,cutoff="p-value")
}

