\name{cpgcount}
\alias{cpgcount}
\title{
Call C programs to R.
}
\description{
Call C programs to R for calculate which CpG are contained in MRE-CpG.
}
\usage{
cpgcount(data2, data3, cpg2, cpg3, datalength, cpglength, 
count=rep(0,cpglength))
}
\arguments{
  \item{data2}{
Start position of each MRE-CpG. 
}
  \item{data3}{
End position of each MRE-CpG. 
}
  \item{cpg2}{ 
Start position of each CpG. 
}
  \item{cpg3}{
End position of each CpG. 
}
  \item{datalength}{
The number of MRE-CpG. 
}
  \item{cpglength}{
The number of MRE-CpG. 
}
  \item{count}{
MRE-CpG count of each CpG. 
}
}
\value{
MRE-CpG count of each CpG.
}
\author{
Yan Zhou, Bo Zhang, Nan Lin, BaoXue Zhang and Ting Wang
}
\examples{
  cpg<-matrix( 1:800, nrow=400 )
  cpg[,2]<-cpg[,1]+2
  data<-cpg[3:100,]
  data[,1]<-data[,1]-1
  data[,2]<-data[,2]+1
  f <- cpgcount(data[,1], data[,2], cpg[,1], cpg[,2], length(data[,1]),
 length(cpg[,2]), count=rep(0,length(cpg[,2])))
}

