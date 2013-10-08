\name{CNVnormal}
\alias{CNVnormal}
\title{
Normalize copy number variation (CNV).
}
\description{
The function is used to normalize copy number variation.
}
\usage{
CNVnormal(CNVfile,bincount)
}
\arguments{
  \item{CNVfile}{
 The path of copy number variation file.
}
  \item{bincount}{
Count of all bins.
}
}
\value{
Count of all bins after CNV normalization. 
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
  winbin[,4]<-rpois(100, lambda=5)
  cnv<-winbin[1:5,]
  cnv[,2]<-c(0,10000,20000,30000,40000)
  cnv[,3]<-cnv[,2]+10000
  cnv[,4]<-c(1.2,1.6,1,2,1)
  CNVfile<-paste(setwd(getwd()), "/CNVfile.bed", sep = "")  
  write.table(cnv, CNVfile, quote=FALSE, row.names =FALSE,col.names =FALSE)
  f<-CNVnormal(CNVfile,winbin)
}


