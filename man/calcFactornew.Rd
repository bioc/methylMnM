\name{calcFactornew}
\alias{calcFactornew}
\title{
Normalization factor.
}
\description{
Amends of TMM normalization for our methond. 
}
\usage{
calcFactornew(obs, ref, m, k, logratioTrim=.3, sumTrim=0.05,
 doWeighting=TRUE, Acutoff=-1e10) 

}
\arguments{
\item{obs}{Counts of treatment sample.}
\item{ ref}{Counts of control sample.}
\item{m}{ The number of CpG in each bin. }
\item{k}{The number of MRE-CpG in each bin.}
\item{logratioTrim}{amount of trim to use on log-ratios ("M" values)}
\item{sumTrim}{amount of trim to use on the combined absolute levels ("A" values)}
\item{doWeighting}{logical, whether to compute (asymptotic binomial precision) weights}
\item{Acutoff}{cutoff on "A" values to use before trimming}
}
\value{
A real value larger than 0.
}
\author{
Yan Zhou, Bo Zhang, Nan Lin, BaoXue Zhang and Ting Wang
}
\examples{
  d <- matrix( rpois(1000, lambda=5), nrow=200 )
  m<-rep(1,nrow=200 )
  k<-rep(1,nrow=200 )
  f <- calcFactornew(d[,2], d[,1], m, k, logratioTrim=.3, sumTrim=0.05,
 doWeighting=TRUE, Acutoff=-1e10) 
}
