\name{normpdft1}
\alias{normpdft1}
\title{
Compute p-value with normal distribution.
}
\description{
The function is used to compute p-value with normal distribution.
}
\usage{
normpdft1(t,n,p,c1,c2)
}
\arguments{
  \item{t}{
Statistic.
}
  \item{n}{
The sum of MeDIP-seq count and MRE-seq count of each bin of two samples.
}
  \item{p}{
The probability in multinomial distribution. 
}
  \item{c1}{
A constant to balance MeDIP-seq of sample 1 and sample 2.
}
  \item{c2}{
A constant to balance MRE-seq of sample 1 and sample 2.
}
}
\value{
statistic of a bin.
}

\author{
Yan Zhou, Bo Zhang, Nan Lin, BaoXue Zhang and Ting Wang
}
\examples{
t<-0.1
n<-200
p<-c(0.25,0.25,0.25,0.25)
c1<-1
c2<-1
f<-normpdft1(t,n,p,c1,c2)
}


