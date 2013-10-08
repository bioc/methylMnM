\name{countMREcpgbin}
\alias{countMREcpgbin}
\title{
Compute the MRE CpG number of each bin with MRE CpG sites.
}
\description{
The function is used to compute the MRE CpG number of each bin with MRE CpG sites. MRE
CpG is some specific CpGs in genome-wide, such as "CCGG", "GCGC" and
"CCGC". The specific CpG number is directly bound up with each
experiment.
}
\usage{
countMREcpgbin(mrecpg.site,file.allcpgsite,file.bin=NULL,
writefile=NULL, binlength=500)
}
\arguments{
  \item{mrecpg.site}{
 The data of mreCpG site. 
}
  \item{file.allcpgsite}{
 The path of all cpg site file or
sequence tag file.
}
  \item{file.bin}{ 
 The path of all bins file. For
computing the number of sequence tag of each window, we use the file
as a normalize window position. (If we do not use the file, there
will be defaulted as NULL).
}
  \item{writefile}{ 
 The path of output result. (If
writefile=NULL, there will return the results back to main program )
}
  \item{binlength}{ 
The length of each window. (Defaulted length
is 500)
}
}
\value{
The output file
is include four columns, that is "chromosome", "start position",
"end position" and
"MRE CpG count".
}

\author{
Yan Zhou, Bo Zhang, Nan Lin, BaoXue Zhang and Ting Wang
}
\examples{
  datafile<-system.file("extdata",  package = "methylMnM") 
  filepath<-datafile[1]
  file<-paste(filepath,"/three_Mre_CpGsite_chr18.txt",sep="")
  file1<-paste(filepath,"/all_CpGsite_chr18.txt",sep="")
  five_Mre_CpGsite<-read.table(file, header=FALSE, as.is=TRUE)  
  f<-countMREcpgbin(mrecpg.site=five_Mre_CpGsite[1:1000,],
file.allcpgsite=file1,binlength=5000)
}

