\name{removeblacklist}
\alias{removeblacklist}
\title{
Remove  blacklist.
}
\description{
The function is used to remove  blacklist which we are not interest.
}
\usage{
removeblacklist(file2,cpg)
}
\arguments{
  \item{file2}{
 The path of blacklist site file.
}
  \item{cpg}{
All bins.
}
}
\value{
All bins except blacklist region. 
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
  blacklist<-winbin[1:5,]
  blacklist[,2]<-c(0,10000,20000,30000,40000)
  blacklist[,3]<-blacklist[,2]+1000
  blacklistfile<-paste(setwd(getwd()), "/blacklist.bed", sep = "")  
  write.table(blacklist, blacklistfile, quote=FALSE, 
 row.names =FALSE,col.names =FALSE)
  f<-removeblacklist(blacklistfile,winbin)
}

