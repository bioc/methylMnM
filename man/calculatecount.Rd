\name{calculatecount}
\alias{calculatecount}
\title{
Call C programs to R.
}
\description{
Call C programs to R for calculate MeDIP-seq or CpG count of each bin.
}
\usage{
calculatecount(data2, data3, cpg2, cpg3, datalength, cpglength, 
count=rep(0,cpglength))
}
\arguments{
  \item{data2}{
Start position of each tag. 
}
  \item{data3}{
End position of each tag. 
}
  \item{cpg2}{ 
Start position of each bin. 
}
  \item{cpg3}{
End position of each bin. 
}
  \item{datalength}{
The number of tags 
}
  \item{cpglength}{
The number of bins 
}
  \item{count}{
Read count of each bin.
}
}
\value{
Read count of each bin.
}
\author{
Yan Zhou, Bo Zhang, Nan Lin, BaoXue Zhang and Ting Wang
}
\examples{
  data<-matrix( 1:800, nrow=400 )
  data[,2]<-data[,1]+37
  cpg<-matrix( 1:20, nrow=10)
  cpg[,1]<-seq(0,360,length=10)
  cpg[,2]<-seq(40,400,length=10)
  f <- calculatecount(data[,1], data[,2], cpg[,1], cpg[,2], length(data[,1]),
 length(cpg[,2]), count=rep(0,length(cpg[,2])))
}


