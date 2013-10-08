\name{calculatecount1}
\alias{calculatecount1}
\title{
Call C programs to R.
}
\description{
Call C programs to R for calculate MRE-seq "+" direction count of each bin.
}
\usage{
calculatecount1(data2, data3, cpg2, cpg3, datalength, cpglength,
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
Count of MRE-seq "+" direction of each bin.
}
}
\value{
Count of MRE-seq "+" direction of each bin.
}
\author{
Yan Zhou, Bo Zhang, Nan Lin, BaoXue Zhang and Ting Wang
}
\examples{
  data<-matrix( 1:400, nrow=200 )
  cpg<-matrix( 1:40, nrow=20)
  cpg[,1]<-seq(0,380,length=20)
  cpg[,2]<-seq(20,400,length=20)
  f <- calculatecount1(data[,1], data[,2], cpg[,1], cpg[,2], length(data[,1]),
 length(cpg[,2]), count=rep(0,length(cpg[,2])))
}



