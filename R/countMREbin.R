countMREbin <-function(file.MREsite,file.blacklist=NULL, file.bin=NULL,file.CNV=NULL, cutoff=0,writefile=NULL, reportfile=NULL, binlength=500)
{
   message(paste("Count", " MRE-seq ", "number of each window" ))
   if(is.null(file.bin)){
       allcpg<-NULL
   }else{
       allcpg<-read.table(file.bin, header=TRUE,sep="\t", as.is=TRUE) 
   }

   #read data
   if(is.null(file.MREsite)) stop("Error, No file which we should count for each window") 
   message("Reading data, please wait...") 
   nrec<-length(count.fields(file.MREsite))
   tag_cpg1<-read.table(file.MREsite, header=FALSE, skip=0, nrows=floor(nrec/4),as.is=TRUE)
   tag_cpg2<-read.table(file.MREsite, header=FALSE, skip=floor(nrec/4), nrows=floor(nrec/4),as.is=TRUE)
   tag_cpg3<-read.table(file.MREsite, header=FALSE, skip=2*floor(nrec/4), nrows=floor(nrec/4),as.is=TRUE)
   tag_cpg4<-read.table(file.MREsite, header=FALSE, skip=3*floor(nrec/4), nrows=nrec,as.is=TRUE)
   #####################################
   if(is.null(allcpg)){ 
         chrtype1<-unique(tag_cpg1[,1])
         chrtype2<-unique(tag_cpg2[,1])
         chrtype3<-unique(tag_cpg3[,1])
         chrtype4<-unique(tag_cpg4[,1])
         chrstring<-unique(c(chrtype1,chrtype2,chrtype3,chrtype4))
   }else{
         chrstring<-unique(allcpg[,1])
   }
   chrsizes<-length(chrstring)
   #chrstring<-paste("chr",c(1:22,"X","Y"),sep = "")
   cpg<-matrix(0,2,4)
   for(i in 1:chrsizes){
       #find site of each CpG
       xx141<-tag_cpg1[tag_cpg1[,1]==chrstring[i],]
       xx142<-tag_cpg2[tag_cpg2[,1]==chrstring[i],]
       xx143<-tag_cpg3[tag_cpg3[,1]==chrstring[i],]
       xx144<-tag_cpg4[tag_cpg4[,1]==chrstring[i],]
       tag_cpgchr19<-rbind(xx141,xx142,xx143,xx144)
       #define all 1KB or 500BP bins 
       tag_cpgchr19<-tag_cpgchr19[order(tag_cpgchr19[,2]),]
       #define all 1KB or 500BP bins 
       if(is.null(allcpg)){
             message("Defining length of all windows") 
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
             cpg19<-allcpg[allcpg[,1]==chrstring[i],1:4]
       }
       cpg19<-cpg19[order(cpg19[,2]),]
       #count MRE-seq number of each bin
       #  cpg19[,4]<-calculatecount1(tag_cpgchr19[,2],tag_cpgchr19[,3],cpg19[,2],cpg19[,3],length(tag_cpgchr19[,2]),length(cpg19[,2]),count=rep(0,length(cpg19[,2]))) 
       tag_cpgchr191<-tag_cpgchr19[tag_cpgchr19[,6]=="+",]
       positivecount<- calculatecount1(tag_cpgchr191[, 2], 
       tag_cpgchr191[, 3], cpg19[, 2], cpg19[, 3], length(tag_cpgchr191[, 
                  2]), length(cpg19[, 2]), count = rep(0, length(cpg19[, 
                  2])))       
       tag_cpgchr192<-tag_cpgchr19[tag_cpgchr19[,6]=="-",]
       negativecount<- calculatecountneg(tag_cpgchr192[, 2], 
                tag_cpgchr192[, 3], cpg19[, 2], cpg19[, 3], length(tag_cpgchr192[, 
                  2]), length(cpg19[, 2]), count = rep(0, length(cpg19[, 
                  2])))
       cpg19[, 4] <- positivecount+negativecount
       cpg<-rbind(cpg,cpg19)
       print(i)
   }
   cpg<-cpg[-(1:2),]
   if(is.null(file.blacklist)){
        cpg2<-cpg
   }else{
        cpg2<-removeblacklist(file.blacklist,cpg)
   }
   cpg2<-CNVnormal(file.CNV,cpg2)
   if(cutoff>0 & cutoff<1){
        message("cutoff: Remove PCR effect for MRE-seq")
        ordervalue<-cpg2[cpg2[,4]>0,4]
        n1<-floor(length(ordervalue)*cutoff)
        ordervalue<-ordervalue[order(ordervalue)]
        cutvalue<-ordervalue[length(ordervalue)-n1]
        cpg2[cpg2[,4]>cutvalue,4]<-cutvalue
   }
  
   cpg2[,2]<-as.integer(cpg2[,2])
   cpg2[,3]<-as.integer(cpg2[,3])
   if(!is.null(reportfile)){
        write(paste("binlength:", binlength), reportfile,  append = FALSE)
        write(paste("the number of bin:", length(cpg2[,4])), reportfile,  append = TRUE)
        write(paste("total reads before processing:", nrec), reportfile,  append = TRUE)
        write(paste("total reads after processing:", sum(cpg2[,4])), reportfile,  append = TRUE)        
   }
   if(is.null(writefile)){
        return(cpg2)
   }else{
        write.table(cpg2, writefile,sep="\t", quote=FALSE, row.names =FALSE)
   }
   message("The End")
}

