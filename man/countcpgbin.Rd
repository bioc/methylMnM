\name{countcpgbin}
\alias{countcpgbin}
\title{
Compute the total CpG number of each bin with each CpG site.
}
\description{
The function is used to compute the total CpG number of each bin with each CpG site.
}
\usage{
countcpgbin(file.cpgsite,file.blacklist=NULL,file.bin=NULL, writefile=NULL,
 reportfile=NULL, binlength=500)
}
\arguments{
  \item{file.cpgsite}{
 The path of cpg site file or sequence
tag file.
}
  \item{file.blacklist}{
 The path of blacklist file (If we do
not use the file, there will be defaulted as NULL).
}
  \item{file.bin}{
 The path of all cpg bin file. For
computing the number of sequence tag of each window, we use the file
as a normalization window position. (If we do not use the file,
there will be defaulted as NULL).
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
The CpG site should include at least three columns "chromosome",
"start position" and "end position". The output file is include four columns, that is "chromosome",
"start position", "end position" and "CpG count". Also, the function
output a report for some
parameters.
}
\author{
Yan Zhou, Bo Zhang, Nan Lin, BaoXue Zhang and Ting Wang
}
\examples{
  datafile<-system.file("extdata",  package = "methylMnM") 
  filepath<-datafile[1]
  file.cpgsite<-paste(filepath,"/all_CpGsite_chr18.txt",sep="")
  f<-countcpgbin(file.cpgsite, binlength=5000)
}



