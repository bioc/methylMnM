removeblacklist <-
function(file2,cpg)
{
   if(is.null(file2)){
      message("Without removing blacklist")
      return(cpg)
   }else{
      message("Removing blacklist")
      blacklist<-read.table(file2, header=FALSE, as.is=TRUE)
      ind<-rep(1,length(cpg[,1]))
      for(i in 1:length(blacklist[,1])){
         ind[cpg[,2]<blacklist[i,3] & cpg[,3]>blacklist[i,2] & blacklist[i,1]==cpg[,1]]<-0 
         if(i %% 50==0) print(i)
      }
      cpg2<-cpg[ind==1,]
      #cpg2<-cpg2[cpg2[,4]!=0,]
      return(cpg2)
  }
}

