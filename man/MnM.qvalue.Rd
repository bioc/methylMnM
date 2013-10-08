\name{MnM.qvalue}
\alias{MnM.qvalue}
\title{Estimate the q-values for a given set of p-values}
\description{
  The function is used to estimate the q-values for a given set of p-values.  The q-value of a
  test measures the proportion of false positives incurred (called the
  false discovery rate) when that particular test is called significant.
}
\usage{
MnM.qvalue(datafile,writefile=NULL,reportfile=NULL)
}
\arguments{
\item{datafile}{Input data of p-values file (Including all input)}
\item{writefile}{ 
The file path of output result. (If
writefile=NULL,there will return the results back to main program )
}
\item{reportfile}{The path of output results of bin length,
the number of bin, total reads before processing and total reads
after processing.
}
}
\value{
The output file is just add a q-value column to the input file.
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
  pvaluefile<-paste(setwd(getwd()), "/pvalue.bed", sep = "")  
  write.table(cpgpq, pvaluefile,sep="\t", quote=FALSE,row.names =FALSE)
  f<-MnM.qvalue(datafile=pvaluefile)
}
