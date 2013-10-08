##############################
##############################
#(15)select DMRs with p-value
MnM.selectDMR<-function(frames = NULL, up =1.45, down = 1/1.45, p.value.MM = 0.01, p.value.SAGE = 0.01,q.value = 0.01,cutoff="q-value",  quant= 0.6){
     if(is.null(frames)) stop("Must specify a frames set.")
     # frames[,4]<-frames[,4]*(1e6/sum(frames[,4]))
     # frames[,5]<-frames[,5]*(1e6/sum(frames[,5]))
     # frames[,6]<-frames[,6]*(1e6/sum(frames[,6]))
     # frames[,7]<-frames[,7]*(1e6/sum(frames[,7]))
     D_value<-quantile(abs(frames[,4]-frames[,5]), probs = quant)
     print(paste("Total number of frames: ", length(frames[, 1]), sep = ""), quote = F)
     if(cutoff=="q-value"){
          frames<-frames[frames[,12]<=q.value,]
          print(paste("Remaining number of frames with q.value<=", q.value, ": ", length(frames[, 1]), sep = ""), quote = F)
     }else{
          frames1<-frames[frames[,9]==0,]
          frames1<-frames1[frames1[,10]<=p.value.SAGE,]
          print(paste("Remaining number of frames with MRE-CpG=0 and p.value<=", p.value.SAGE, ": ", length(frames1[, 1]), sep = ""), quote = F)
          frames2<-frames[frames[,9]!=0,]
          frames2<-frames2[frames2[,10]<=p.value.MM,]
          print(paste("Remaining number of frames with MRE-CpG!=0 and p.value<=", p.value.MM, ": ", length(frames2[, 1]), sep = ""), quote = F)
          frames<-rbind(frames1,frames2)
     }  
     frames<-frames[abs(frames[,4]-frames[,5])>= D_value & frames[,8]>0,]
     print(paste("Remaining number of frames with D_value>=", D_value, "& CpG!=0: ", length(frames[, 1]), sep = ""), quote = F)
     ratio<-frames[,4]/frames[,5]
     frames<-frames[ratio<=down|ratio>=up,]
     print(paste("Remaining number of frames with ratio<=",down," or ratio>=",up, ": ", length(frames[, 1]), sep = ""), quote = F)
     if(cutoff=="q-value"){
        print(paste("[Note: There are ", sum(frames[, 10]==0), " frames associated to a q-value==0.]", sep = ""), quote = F)
     }else{
        print(paste("[Note: There are ", sum(frames[, 10]==0), " frames associated to a p-value==0.]", sep = ""), quote = F)
     }
     print(paste("Remaining number of frames called DMR" , ":", length(frames[, 1]), sep = ""), quote = F)
     return(frames)
 }

