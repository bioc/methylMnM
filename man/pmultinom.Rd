\name{pmultinom}
\alias{pmultinom}
\title{
Call C programs to R.
}
\description{
Call C programs to R for calculate p-value of each bin with multinomial distribution.
}
\usage{
pmultinom(T, SIZE,length, P1, P2, P3, P4, C1, C2,
 pvalue=rep(0,length(T)))
}
\arguments{
  \item{T}{
Statistic.
}
  \item{SIZE}{
The sum of MeDIP-seq count and MRE-seq count of each bin of two samples.
}
  \item{length}{ 
The number of bins.
}
  \item{P1}{
The probability of MeDIP-seq of sample 1 in multinomial distribution. 
}
  \item{P2}{
The probability of MeDIP-seq of sample 2 in multinomial distribution. 
}
  \item{P3}{
The probability of MRE-seq of sample 1 in multinomial distribution. 
}
  \item{P4}{
The probability of MRE-seq of sample 2 in multinomial distribution. 
}
  \item{C1}{
A constant to balance MeDIP-seq of sample 1 and sample 2. 
}
  \item{C2}{
A constant to balance MRE-seq of sample 1 and sample 2.
}
  \item{pvalue}{
p-values of windows.
}
}
\value{
p-value.
}
\author{
Yan Zhou, Bo Zhang, Nan Lin, BaoXue Zhang and Ting Wang
}
\examples{
T<-4
SIZE<-200
p<-c(0.25,0.25,0.25,0.25)
c1<-1
c2<-1
length<-1
f<-pmultinom(T, SIZE,length, p[1], p[2], p[3], p[4], c1,
 c2, pvalue=rep(0,length(T)))
}
