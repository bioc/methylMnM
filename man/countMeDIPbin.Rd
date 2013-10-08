\name{countMeDIPbin}
\alias{countMeDIPbin}
\title{
Compute the total MeDIP-seq number of each bin.
}
\description{
The function is used to compute the total MeDIP-seq  number of each bin.
}
\usage{
countMeDIPbin (file.Medipsite,file.blacklist=NULL,file.bin=NULL,
file.CNV=NULL, writefile=NULL, reportfile=NULL, binlength=500)
}
\arguments{
  \item{file.Medipsite}{
 The path of MeDIP-seq site file or sequence
tag file.
}
  \item{file.blacklist}{
 The path of blacklist file (If we do
not use the file, there will be defaulted as NULL).
}
  \item{file.bin}{
 The path of all bins file. For
computing the number of sequence tag of each window, we use the file
as a normalization window position. (If we do not use the file,
there will be defaulted as NULL).
}
  \item{file.CNV}{
 If need, we should input CNV file to normalize
count of each bin.
}
  \item{writefile}{
The path of output results. (If
writefile=NULL, there will return the results back to main program.)
}
  \item{reportfile}{
 The path of output results of bin length,
the number of bin, total reads before processing and total reads
after processing.
}
  \item{binlength}{
 The length of each window.(Defaulted
length is 500).
}
}
\value{
The MeDIP-seq site should include at least three columns "chromosome",
"start position" and "end position". The output file is include four columns, that is "chromosome",
"start position", "end position" and "MeDIP-seq count". Also, the function
output a report for some
parameters.
}
\author{
Yan Zhou, Bo Zhang, Nan Lin, BaoXue Zhang and Ting Wang
}
\examples{
  datafile<-system.file("extdata",  package = "methylMnM") 
  filepath<-datafile[1]
  file.Medipsite<-paste(filepath,"/all_CpGsite_chr18.txt",sep="")
  f<-countMeDIPbin(file.Medipsite, binlength=5000)
}



