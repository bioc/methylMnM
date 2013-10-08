\name{countMREbin}
\alias{countMREbin}
\title{
Compute the total MRE-seq number of each bin.
}
\description{
The function is used to compute the total MRE-seq number of each bin.
}
\usage{
countMREbin(file.MREsite,file.blacklist=NULL, file.bin=NULL,
file.CNV=NULL, cutoff=0,writefile=NULL, 
reportfile=NULL, binlength=500)
}
\arguments{
  \item{file.MREsite}{
 The path of MRE-seq sites file.
}
  \item{file.blacklist}{
 The path of blacklist file (If we do
not use the file, there will be defaulted as NULL).
}
  \item{file.bin}{
 The path of all bin file. For
computing the number of sequence tag of each window, we use the file
as a normalization window position. (If we do not use the file,
there will be defaulted as NULL).
}
  \item{file.CNV}{
 If need, we should input CNV file to normalize
count of each bin.
}
  \item{cutoff}{
The critical value of PCR. (If we do not use
the critical value, there will be defaulted as 0.)
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
The MRE-seq sites should include at least three columns "chromosome",
"start position" and "end position". The output file is include four columns, that is "chromosome",
"start position", "end position" and "MRE-seq count". Also, the function
output a report for some
parameters.
}
\author{
Yan Zhou, Bo Zhang, Nan Lin, BaoXue Zhang and Ting Wang
}
\examples{
  datafile<-system.file("extdata",  package = "methylMnM") 
  filepath<-datafile[1]
  file.MREsite<-paste(filepath,"/all_CpGsite_chr18.txt",sep="")
  f<-countMREbin(file.MREsite, binlength=5000)
}



