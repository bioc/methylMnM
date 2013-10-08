countMREcpgbin <-
function(mrecpg.site,file.allcpgsite,file.bin=NULL,writefile=NULL, binlength=500)
{
      message(paste("Count", "MRE-CpG", "number of each window" ))
      if(is.null(file.bin)){
              cpg2<-NULL
      }else{
              cpg2<-read.table(file.bin, header=TRUE,sep="\t", as.is=TRUE) 
      }
      if(is.null(file.allcpgsite)) stop("Error, No file for all CpG sites") 
      message("Reading data, please wait...") 
      nrec<-length(count.fields(file.allcpgsite))
      #read data
      tag_cpg1<-read.table(file.allcpgsite, header=FALSE, skip=0, nrows=floor(nrec/4),as.is=TRUE)
      tag_cpg2<-read.table(file.allcpgsite, header=FALSE, skip=floor(nrec/4), nrows=floor(nrec/4),as.is=TRUE)
      tag_cpg3<-read.table(file.allcpgsite, header=FALSE, skip=2*floor(nrec/4), nrows=floor(nrec/4),as.is=TRUE)
      tag_cpg4<-read.table(file.allcpgsite, header=FALSE, skip=3*floor(nrec/4), nrows=nrec,as.is=TRUE)
      #####################################
      if(is.null(cpg2)){ 
             chrtype1<-unique(tag_cpg1[,1])
             chrtype2<-unique(tag_cpg2[,1])
             chrtype3<-unique(tag_cpg3[,1])
             chrtype4<-unique(tag_cpg4[,1])
             chrstring<-unique(c(chrtype1,chrtype2,chrtype3,chrtype4))
      }else{
             chrstring<-unique(cpg2[,1])
      }
      chrsizes<-length(chrstring)
      #chrstring<-paste("chr",c(1:22,"X","Y"),sep = "")
      mrecpgcount<-matrix(0,2,4)
      for(i in 1:chrsizes){
                xx141<-tag_cpg1[tag_cpg1[,1]==chrstring[i],]
                xx142<-tag_cpg2[tag_cpg2[,1]==chrstring[i],]
                xx143<-tag_cpg3[tag_cpg3[,1]==chrstring[i],]
                xx144<-tag_cpg4[tag_cpg4[,1]==chrstring[i],]
                tag_cpgchr19<-rbind(xx141,xx142,xx143,xx144)
                tag_cpgchr19<-tag_cpgchr19[order(tag_cpgchr19[,2]),]
                cpg18<-tag_cpgchr19
                mrechr<-mrecpg.site[mrecpg.site[,1]==chrstring[i],]
                mrechr<-mrechr[order(mrechr[,2]),]
                message("Obtaining CpG sites which are included in MRE-CpG sites.")
                cpg18[,4]<-cpgcount(mrechr[,2],mrechr[,3],cpg18[,2],cpg18[,3],length(mrechr[,2]),length(cpg18[,2]),count=rep(0,length(cpg18[,2])))
                mrecpg1<-cpg18[cpg18[,4]!=0,]
                #define all 1KB or 500BP bins 
                if(is.null(cpg2)){
                            n<-ceiling(max(tag_cpgchr19[,3])/binlength)
                            cpg19<-tag_cpgchr19[1:n,1:3]
                            cpg19[,1]<-chrstring[i]
                            cpg19<-cbind(cpg19,cpg19[,2])
                            colnames(cpg19)=c("V1", "V2", "V3", "V4")
                            poststart<-c(0:(n-1))*binlength
                            postend<-c(1:n)*binlength
                            cpg19[,2]<-poststart
                            cpg19[,3]<-postend
                            if(max(tag_cpgchr19[,3])<n*binlength){
                                      cpg19[n,3]<-max(tag_cpgchr19[,3])
                            }  
                            cpg19[,4]<-0
                }else{
                            cpg19<-cpg2[cpg2[,1]==chrstring[i],]
                }
                cpg19<-cpg19[order(cpg19[,2]),]
                mrecpg1<-mrecpg1[order(mrecpg1[,2]),]
                #count CpG of each bin
                message(paste("Count",  "MRE-CpG",  "number of each window of", chrstring[i]))
                cpg19[,4]<-calculatecount(mrecpg1[,2],mrecpg1[,3],cpg19[,2],cpg19[,3],length(mrecpg1[,2]),length(cpg19[,2]),count=rep(0,length(cpg19[,2])))
                mrecpgcount<-rbind(mrecpgcount,cpg19)
                print(i)
      }
      mrecpgcount<-mrecpgcount[-(1:2),]
      mrecpgcount[,2]<-as.integer(mrecpgcount[,2])
      mrecpgcount[,3]<-as.integer(mrecpgcount[,3])
      if(is.null(writefile)){
                return(mrecpgcount)
      }else{
                write.table(mrecpgcount, writefile,sep="\t", quote=FALSE, row.names =FALSE)
      }
     message("The End")
}

