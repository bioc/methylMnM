#(14) q-value
MnM.qvalue<-function(datafile,writefile=NULL,reportfile=NULL){
  message("Estmate q-value for each window" )
  message("Reading data, please wait...") 
  pg<-read.table(datafile, header=TRUE,sep="\t", as.is=TRUE)
  pg<-pg[,1:11]
  pg1<-pg[pg[,9]==0,]
  pg2<-pg[pg[,9]!=0,]
  if(length(pg1[,1])==0 | length(pg2[,1])==0){
        message("Only one group") 
        pg11<-pg[,10]
        l1<-sort(pg11)
        k<-(length(l1))/(1-l1[1])
        for(i in 2:length(l1)){
           lgi<-(length(l1)+1-i)/(1-l1[i])
           if(k<lgi){ break } 
           k<-lgi
        }
        kd<-(floor(lgi)+1)/length(l1)
        pi0<-min(kd,1)
        idex<-qvalue.rank(pg[,10])
        qvalue1<-pg[,10]*(length(pg[,10])/idex)*pi0
        q1<-rep(1,length(qvalue1))
        qvalue<-apply(cbind(qvalue1,q1),1,min)

  }else{
        message("There are two groups.")
        pg11<-pg1[,10]
        pg12<-pg2[,10]
        l1<-sort(pg11)

        k<-(length(l1))/(1-l1[1])
        for(i in 2:length(l1)){
           lgi<-(length(l1)+1-i)/(1-l1[i])
           if(k<lgi){ break } 
           k<-lgi
        }
        kd<-(floor(lgi)+1)/length(l1)
        pi0<-min(kd,1)

        l2<-sort(pg12)
        k1<-(length(l2))/(1-l2[1])
        for(j in 2:length(l2)){
            lgi1<-(length(l2)+1-j)/(1-l2[j])
            if(k1<lgi1){ break } 
            k1<-lgi1
        }

        kd1<-(floor(lgi1)+1)/length(l2)
        pi1<-min(kd1,1)
        ############################
        pi<-(i+j)/length(pg[,1])

        pg11a<-pg11*(pi0/(1-pi0))
        pg12a<-pg12*(pi1/(1-pi1))

        pgs<-pg[,10]
        pgs[pg[,9]==0]<-pg11a
        pgs[pg[,9]!=0]<-pg12a
        idex<-qvalue.rank(pgs)
        qvalue1<-pgs*(length(pgs)/idex)*pi
        q1<-rep(1,length(qvalue1))
        qvalue<-apply(cbind(qvalue1,q1),1,min)
  }
  new<-cbind(pg,qvalue)
    if(!is.null(reportfile)){
            write(paste("Number of windows:", length(new[,1])), reportfile,  append = FALSE)
            write(paste("Number of SAGE windows:", length(pg1[,1])), reportfile,  append = TRUE)
            write(paste("Number of M&M windows:", length(pg2[,1])), reportfile,  append = TRUE)
            write(paste("Pi0 of SAGE:", pi0), reportfile,  append = TRUE) 
            write(paste("Pi0 of M&M:", pi1), reportfile,  append = TRUE)
            write(paste("Pi0 of all windows:",1-pi), reportfile,  append = TRUE) 
         
            write(paste("Number of q<1e-8:", sum(new[,12]<1e-8)), reportfile,  append = TRUE)
            write(paste("Number of q<1e-4:", sum(new[,12]<1e-3)), reportfile,  append = TRUE)     
     }

     new[,4]<-new[,4]*(1e6/sum(new[,4]))
     new[,5]<-new[,5]*(1e6/sum(new[,5]))
     new[,6]<-new[,6]*(1e6/sum(new[,6]))
     new[,7]<-new[,7]*(1e6/sum(new[,7]))

     if(is.null(writefile)){
                return(new)
      }else{
                write.table(new, writefile,sep="\t", quote=FALSE, row.names =FALSE)
      }
      message("The End")
}

