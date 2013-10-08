CNVnormal<-function(CNVfile,bincount)
{
   bincount1<-bincount
   if(is.null(CNVfile)){
      message("Without CNV normalization")
      return(bincount1)
   }else{
      message("Doing CNV normalization")
      CNV<-read.table(CNVfile, header=FALSE, as.is=TRUE)
      CNV[CNV[,1]=="chr23",1]<-"chrX"
      CNV[CNV[,1]=="chr24",1]<-"chrY"
      #ind<-rep(1,length(bincount1[,1]))
      for(i in 1:length(CNV[,1])) {
          if(CNV[i,4]==0){
               bincount1[bincount[,2]<CNV[i,3] & bincount[,3]>CNV[i,2] & CNV[i,1]==bincount[,1],4]<-0
          }else{
               bincount1[bincount[,2]<CNV[i,3] & bincount[,3]>CNV[i,2] & CNV[i,1]==bincount[,1],4]<-2*bincount1[bincount[,2]<CNV[i,3] & bincount[,3]>CNV[i,2] & CNV[i,1]==bincount[,1],4]/CNV[i,4]
          }
          #ind[bincount[,2]<CNV[i,3] & bincount[,3]>CNV[i,2] & CNV[i,1]==bincount[,1]]<-0
          if(i %% 50==0) print(i)
      }
      #bincount1[ind==1,4]<-bincount1[ind==1,4]
      return(bincount1)
   }
}