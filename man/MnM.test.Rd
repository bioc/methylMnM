\name{MnM.test}
\alias{MnM.test}
\title{
Compute p-value of each bin.
}
\description{
The function is used to compute p-value of each bin.
}
\usage{
MnM.test(file.dataset=NULL,chrstring=NULL,file.cpgbin=NULL,
file.mrecpgbin=NULL,writefile=NULL,reportfile=NULL,
mreratio=3/7,method="XXYY", psd=2,mkadded=1,a=1e-16,
cut=100,top=500)
}
\arguments{
  \item{file.dataset}{
The files path of sample. (datafile should be
c(datafile1,datafile2,datafile3,datafile4), where datafile1 and
datafile2 are  path of Medip-seq data,   datafile3 and datafile4 are
path of MRE-seq data). 
}
  \item{chrstring}{
The chromosome should be test.
}
  \item{file.cpgbin}{ 
The file path of all CpG number of each bin.
}
  \item{file.mrecpgbin}{ 
 The file path of MRE-CpG number of each
window (If NULL, mrecpgfile will equal to cpgfile).
}
  \item{writefile}{ 
The file path of output result. (If
writefile=NULL, there will return the results back to main program )
}
  \item{reportfile}{ 
 The path of output results of bin length,
the number of bin, total reads before processing and total reads
after processing.
}
  \item{mreratio}{ 
The ratio of total unmethylation level with
total methylation level (Defaulted mreratio is 3/7).
}
  \item{method}{ 
Option different data for the test. 
}
  \item{psd}{ 
The parameters of pseudo count,
which pseudo count added to Medip-seq and MRE-seq count. 
}
  \item{mkadded}{ 
Added
to all CpG and MRE CpG (We set psd=2 and mkadded=1 as defaulted for
robust)
}
  \item{a}{ 
Cut-off for recalculating p-value with multinomial distribution when normal
p-values smaller than a and the sum of observations smaller than
top.
}
  \item{cut}{ 
Cut-off for recalculating p-value with multinomial distribution when the sum of observations smaller than
cut. }
  \item{top}{ 
Cut-off for recalculating p-value with multinomial distribution when normal
p-values smaller than a and the sum of observations smaller than
top.
}
}
\value{
The output file "writefile" will own eleven columns, that is,
"chr", "chrSt", "chrEnd", "Medip1", "Medip2", "MRE1",
"MRE2", "cg", "mrecg", "pvalue" and "plus-minus". We
also output a report file which will include parameters "s1/s2",
"s3/s4", "N1", "N2", "N3", "N4", "c1", "c2",  "Number of
windows" and "Spend time".
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
  winbinfile1<-paste(setwd(getwd()), "/winbinfile1.bed", sep = "")  
  write.table(winbin, winbinfile1,sep="\t", quote=FALSE, row.names =FALSE)
  winbin1<-winbin
  winbin1[,4]<-winbin[,4]+20
  winbinfile2<-paste(setwd(getwd()), "/winbinfile2.bed", sep = "")  
  write.table(winbin1, winbinfile2,sep="\t", quote=FALSE, row.names =FALSE)
  datafile<-c(winbinfile1,winbinfile2)
  cpg<-winbin
  cpg[,4]<-rpois(100, lambda=12)
  cpgfile<-paste(setwd(getwd()), "/cpgfile.bed", sep = "")  
  write.table(cpg, cpgfile, sep="\t", quote=FALSE, row.names =FALSE)
  f<-MnM.test(file.dataset=datafile,file.cpgbin=cpgfile)
}

