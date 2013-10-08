MnM.test<-
function(file.dataset=NULL,chrstring=NULL,file.cpgbin=NULL,file.mrecpgbin=NULL,writefile=NULL,reportfile=NULL,mreratio=3/7,method="XXYY", psd=2,mkadded=1,a=1e-16,cut=100,top=500)
{
     starttime=proc.time()
     library(edgeR)
     library(statmod)
     #read data
     message("reading data ....")
     if(length(file.dataset)==2) {
            message("The file.dataset is only Medip-seq data")
            data1<-read.table(file.dataset[1], header=TRUE,sep="\t", as.is=TRUE)
            data2<-read.table(file.dataset[2], header=TRUE,sep="\t", as.is=TRUE)
            data3<-data1
            data4<-data2
            data3[,4]<-0
            data4[,4]<-0
     }else{
            if(length(file.dataset)!=2 & length(file.dataset)!=4){stop("The input MeDIP and MRE file should be two or four files")}
            data1<-read.table(file.dataset[1], header=TRUE,sep="\t",as.is=TRUE)
            data2<-read.table(file.dataset[2], header=TRUE,sep="\t", as.is=TRUE)
            data3<-read.table(file.dataset[3], header=TRUE,sep="\t", as.is=TRUE)
            data4<-read.table(file.dataset[4], header=TRUE,sep="\t", as.is=TRUE)
     }
     dataset<-cbind(data1,data2[,4],data3[,4],data4[,4])
     message("Top line of dataset:")
     print(dataset[1,])
     dataset[,4:7] = dataset[,4:7]
     colnames(dataset)=c("V1","V2","V3","V4","V5","V6","V7")
     if(is.null(file.cpgbin)) {
           stop(" No CpG window's file")
     }else{
           cpg<-read.table(file.cpgbin, header=TRUE,sep="\t",as.is=TRUE)
           print(cpg[1,])
     }
     if(is.null(file.mrecpgbin)) {
           mrecpg<-cpg
           message("warning: MRE-CpG window file is same as CpG file")
     }else{
         mrecpg<-read.table(file.mrecpgbin, header=TRUE,sep="\t",as.is=TRUE)
         print(mrecpg[1,])
     }
     cpg[,4] = cpg[,4]
     mrecpg[,4]=mrecpg[,4]
     ########################
     ########################
     if(is.null(mkadded)) mkadded<-1
     cpg[mrecpg[,4]!=0,4]<-cpg[mrecpg[,4]!=0,4]+mkadded
     mrecpg[mrecpg[,4]!=0,4]<-mrecpg[mrecpg[,4]!=0,4]+mkadded
     message("MeDIP-MRE normalizing....")
     newdata<-dataset
     newdata[mrecpg[,4]!=0,6:7]<-newdata[mrecpg[,4]!=0,6:7]*(cpg[mrecpg[,4]!=0,4]/mrecpg[mrecpg[,4]!=0,4])
     NN1<-sum(dataset[,4])
     NN2<-sum(dataset[,5])
     NN3<-sum(dataset[,6])*(sum(cpg[,4])/sum(mrecpg[,4]))
     NN4<-sum(dataset[,7])*(sum(cpg[,4])/sum(mrecpg[,4]))
     if(mreratio>=0)    
     {   if(NN3==0){ cc1<-1 }else{
           cc1<-(NN1/NN3)*mreratio}
         if(NN4==0){ cc2<-1 }else{
           cc2<-(NN2/NN4)*mreratio}
     }else{
          cc1<-1
          cc2<-1
     }
     if(method=="XXYY"){
          newdata[,6]<-round(newdata[,6]*cc1)
          newdata[,7]<-round(newdata[,7]*cc2)
     }else{
          newdata[,6]<-round(newdata[,6]*cc1*(mrecpg[,4]/cpg[,4]))
          newdata[,7]<-round(newdata[,7]*cc2*(mrecpg[,4]/cpg[,4]))
     }
     dataset<-newdata 
     if(is.null(psd)) psd<-2
     dataset[mrecpg[,4]!=0,4:7]<-round(dataset[mrecpg[,4]!=0,4:7]+psd)
     dataset[mrecpg[,4]==0,4:5]<-round(dataset[mrecpg[,4]==0,4:5]+psd)
     Medip<-dataset[,4:5]
     Medip<-as.matrix(Medip)
     fws1s2<-calcNormFactors(Medip,logratioTrim=0.05, sumTrim=0.05)[2]
     mre<-dataset[rowSums(dataset[,6:7])!=2*round(psd) & mrecpg[,4]!=0,6:7]
     if(length(mre[,1])==0){
        fws3s4<-1  
     }else{
        mre<-as.matrix(mre)
        mno2<-cpg[rowSums(dataset[,6:7])!=2*round(psd)& mrecpg[,4]!=0,4]
        kno2<-mrecpg[rowSums(dataset[,6:7])!=2*round(psd)& mrecpg[,4]!=0,4]
        if(method=="XXYY"){
             fws3s4<-calcNormFactors(mre,logratioTrim=0.05, sumTrim=0.05)[2]
        }else{
             fws3s4<-calcFactornew(mre[,2], mre[,1], mno2, kno2,logratioTrim=0.05, sumTrim=0.05)
        }
     }
     #########################
     if(is.null(chrstring)){
            sm<-as.matrix(dataset[,4:7])            
            cpg<-cpg
            kk<-mrecpg
     }else{
            sm<-as.matrix(dataset[dataset[,1]==chrstring,4:7])
            cpg<-cpg[cpg[,1]==chrstring,]
            kk<-mrecpg[mrecpg[,1]==chrstring,]
            if(length(sm[,1])==0) message("No such chrstring data")
     }
     pall<-rep(0,length(sm[,1]))
     if(length(file.dataset)==2) {
            message("Compute all p-values with SAGE.test")
            kk[,4]<-0
      }

     #(1)MRE-CpG=0
     medipct<-sm[kk[,4]==0,1:2]
     cpgS=cpg[kk[,4]==0,4] ##
     message(paste("Computing p-values for only MeDIP-seq windows:", length(medipct[,1])))
     pmedip<-sage.test(medipct[,1], medipct[,2], n1=sum(medipct[,1])/sqrt(fws1s2), n2=sum(medipct[,2])*sqrt(fws1s2))
     pall[kk[,4]==0]<-pmedip
     ts<-cpg[,4]
     ns1=round(sum(medipct[,1])/sqrt(fws1s2))
     ns2=round(sum(medipct[,2])*sqrt(fws1s2))
     probs<-ns1/(ns1+ns2)
     cs<-fws1s2*(sum(medipct[,2])/sum(medipct[,1]))
     tts1<-cpg[kk[,4]==0,4]
     nn<-rowSums(medipct)
     tts1[nn==0]<-0
     
     tts1[nn!=0]<-sqrt(nn[nn!=0])*((cs*medipct[nn!=0,1]/nn[nn!=0]-medipct[nn!=0,2]/nn[nn!=0])-((cs+1)*probs-1))/(sqrt(probs*(1-probs))*(1+cs)) 
     ts[kk[,4]==0]<-tts1

     #(2)MRE-CpG!=0
     if(length(cpg[kk[,4]!=0,1])!=0){
     cpg4<-cpg[kk[,4]!=0,]
     sm1<-sm[kk[,4]!=0,]
     kk1<-kk[kk[,4]!=0,]
     m<-cpg4[,4]
     k<-kk1[,4]
     message(paste("Computing p-values for MeDIP-seq and MRE-seq windows:", length(m)))
     sm1<-round(sm1)
     sm1[sm1[,1]==0,1]<-1
     sm1[sm1[,2]==0,2]<-1
     sm1[sm1[,3]==0,3]<-1
     sm1[sm1[,4]==0,4]<-1
     x1<-sm1[,1]  #numMedip014[,4]
     x2<-sm1[,2]  #numMedip015[,4]
     y1<-sm1[,3]  #nummre10[,4]
     y2<-sm1[,4]  #nummre11[,4]
     z1<-sm1[,3]  #nummre10[,4]
     z2<-sm1[,4]  #nummre11[,4]
     if(method!="XXYY"){
         y1[k!=0]<-y1[k!=0]*(m[k!=0]/k[k!=0])
         y2[k!=0]<-y2[k!=0]*(m[k!=0]/k[k!=0])         
     }
     sm2<-cbind(x1,x2,y1,y2)
     ####################
     ####################
     #estimate rmd21 rmd31 and rmd41
     N1<-sum(sm2[,1])
     N2<-sum(sm2[,2])
     N3<-sum(sm2[,3])
     N4<-sum(sm2[,4])
     #rate of ramda
     rmd21<-fws1s2*(N2/N1)
     rmd43<-fws3s4*(N4/N3)
     #rmd31<-z1/(4*x1)+(N2/N1)*(fws1s2)*(z1/(4*x2))+(N2/N1)*(N3/N4)*(fws1s2)*(z2/(4*x2))/(fws3s4)+(N3/N4)*(z2/(4*x1))/(fws3s4)
     rmd31<-(z1+z2)/(x1+x2)*((1+(N2/N1)*(fws1s2))/(1+(N4/N3)*fws3s4))
     rmd41<-rmd43*rmd31
     rmd<-1+rmd21+rmd31+rmd41
     ##############################
     # estimate probility of each data
     pxy<-matrix(0,length(rmd),4)
     pxy[,1]<-1/rmd
     pxy[,2]<-rmd21/rmd
     pxy[,3]<-rmd31/rmd
     pxy[,4]<-rmd41/rmd
     ######################################
     #calculate p-value.
     T<-abs(rmd21*x1*z2-rmd43*x2*z1)
     size<-round(rowSums(sm1))
     tts3<-(rmd21*x1*z2-rmd43*x2*z1)/size^2  
     tts2<-cpg[kk[,4]!=0,4]
     for(j in 1:length(tts2)){
            if(size[j]==0){tts2[j]<-0}
            else{tts2[j]<-normpdft1(tts3[j],size[j],pxy[j,],rmd21,rmd43)}
     }
     ts[kk[,4]!=0]<-tts2

     pvalue<-rep(2,length(T))
     if(is.null(cut)) cut<-100
     if(sum(size<cut)>0){
           smT<-T[size<cut]
           smsize<-size[size<cut]
           smpxy<-pxy[size<cut,]
           smsm1<-sm1[size<cut,]
           typesm1<-unique(smsm1)
           typermd21<-rmd21
           typermd43<-rmd43
           #typermd31<-typesm1[,3]/(4*typesm1[,1])+(N2/N1)*(fws1s2)*(typesm1[,3]/(4*typesm1[,2]))+(N2/N1)*(N3/N4)*(fws1s2)*(typesm1[,4]/(4*typesm1[,2]))/(fws3s4)+(N3/N4)*(typesm1[,4]/(4*typesm1[,1]))/(fws3s4)
           typermd31<-(typesm1[,3]+typesm1[,4])/(typesm1[,1]+typesm1[,2])*((1+(N2/N1)*(fws1s2))/(1+(N4/N3)*fws3s4))
           typermd41<-typermd43*typermd31
           typermd<-1+typermd21+typermd31+typermd41
           typepxy<-matrix(0,length(typermd),4)
           typepxy[,1]<-1/typermd
           typepxy[,2]<-typermd21/typermd
           typepxy[,3]<-typermd31/typermd
           typepxy[,4]<-typermd41/typermd         
           typeT<-abs(typermd21*typesm1[,1]*typesm1[,4]-typermd43*typesm1[,2]*typesm1[,3])
           typesize<-round(rowSums(typesm1))
           message("Exact compute p-values with multinomial distribution.")
           print(length(typesize))
           if(length(typesize)==1){
               typepvalue<-pmultinom(typeT, typesize, length(typeT), typepxy[1], typepxy[2], typepxy[3], typepxy[4], typermd21, typermd43, pvalue=rep(0,length(typeT)))
           }else{
               typepvalue<-pmultinom(typeT, typesize, length(typeT), typepxy[,1], typepxy[,2], typepxy[,3], typepxy[,4], typermd21, typermd43, pvalue=rep(0,length(typeT)))
           }
          classifypvalue <-function(type1, type2, type3, type4, sm1chring1, sm1chring2, sm1chring3, sm1chring4, p, typelength, sm1chringlength, pvalue=rep(0,sm1chringlength))
          {
                  problity<-.C("pvalueclassify", as.integer(type1), as.integer(type2), as.integer(type3), as.integer   (type4), as.integer(sm1chring1), as.integer(sm1chring2), as.integer(sm1chring3), as.integer(sm1chring4),as.double(p),as.integer(typelength), as.integer(sm1chringlength),as.double(pvalue))
                  return(problity[[12]])
           }
           smpvalue<-classifypvalue(typesm1[,1], typesm1[,2], typesm1[,3], typesm1[,4], smsm1[,1], smsm1[,2],  smsm1[,3], smsm1[,4], typepvalue, length(typesm1[,1]), length(smsm1[,1]), pvalue=rep(0,length(smsm1[,1])))
           pvalue[size<cut]<-smpvalue
     }
     #We use normal distribution to calculate p-value when Rowsum is larger than cut.
     if(sum(size>=cut)>0){
            size1<-size[size>=cut]
            T1<-T[size>=cut]/size1^2
            pxy1<-pxy[size>=cut,]
            pvalue11<-rep(0,length(T1))
            for(j in 1:length(T1)){ 
                    pvalue11[j]<-normpdf(T1[j],size1[j],pxy1[j,],rmd21,rmd43)
            }
            pvalue12<-pvalue11 
            if(is.null(a)) a<-0
            if(is.null(top)) top<-500
            #print(paste("sum(pvalue11>a & pvalue11<0.005 & size1<top):",sum(pvalue11>a & pvalue11<0.005 & size1<top)))
            if(sum(pvalue11>a & pvalue11<0.005 & size1<top)>3000){
                  midpvalue<-pvalue11[pvalue11>a & pvalue11<0.005 & size1<top]
                  T11<-T[size>=cut]
                  T12<-T11[pvalue11>a & pvalue11<0.005 & size1<top]
                  rectsize<-size1[pvalue11>a & pvalue11<0.005 & size1<top]
                  rectpxy<-pxy1[pvalue11>a & pvalue11<0.005 & size1<top,]
                  #sam<-sample(1:length(T12),500)
                  sam<-1:500
                  samT12<-T12[sam]
                  samsize<-rectsize[sam]
                  sampxy<-rectpxy[sam,]
                  sampvalue<-pmultinom(samT12, samsize, length(samT12), sampxy[,1], sampxy[,2], sampxy[,3], sampxy[,4], rmd21, rmd43, pvalue=rep(0,length(samT12)))
                  samnormp<-midpvalue[sam]
                  logsamp<-log(-log10(sampvalue[order(-samnormp)]))
                  logsnop<-log(-log10(samnormp[order(-samnormp)]))
                  lm.D90<- lm(logsamp~ logsnop)
                  est<-lm.D90[1]$coefficients[1]+lm.D90[1]$coefficients[2]*log(-log10(midpvalue))
                  estpvalue<-10^(-exp(est))
                  pvalue12[pvalue11>a & pvalue11<0.005 & size1<top]<-estpvalue
            }else{
                  T11<-T[size>=cut]
                  T12<-T11[pvalue11>a & pvalue11<0.005 & size1<top]
                  rectsize<-size1[pvalue11>a & pvalue11<0.005 & size1<top]
                  rectpxy<-pxy1[pvalue11>a & pvalue11<0.005 & size1<top,]
                  multpvalue<-pmultinom(T12, rectsize, length(T12), rectpxy[,1], rectpxy[,2], rectpxy[,3], rectpxy[,4], rmd21, rmd43, pvalue=rep(0,length(T12))) 
                  pvalue12[pvalue11>a & pvalue11<0.005 & size1<top]<-multpvalue
            }            
            if(sum(pvalue11<a & size1<top)>0  ){                    
                    print(sum(pvalue11<a & size1<top)) 
                    T11<-T[size>=cut]
                    largT<-T11[pvalue11<a & size1<top]
                    largsize<-size1[pvalue11<a & size1<top]
                    largpxy<-pxy1[pvalue11<a & size1<top,]
                    if(sum(pvalue11<a & size1<top)==1){
                        largpvalue<-pmultinom(largT, largsize, length(largT), largpxy[1], largpxy[2], largpxy[3], largpxy[4], rmd21, rmd43, pvalue=rep(0,length(largT)))
                    }else{ 
                        largpvalue<-pmultinom(largT, largsize, length(largT), largpxy[,1], largpxy[,2], largpxy[,3], largpxy[,4], rmd21, rmd43, pvalue=rep(0,length(largT)))
                    }
                    pvalue12[pvalue11<a & size1<top]<-largpvalue 
            }
            pvalue[size>=cut]<-pvalue12   
     }
     pall[kk[,4]!=0]<-pvalue
     }
     endtime=proc.time()
     spendtime=endtime-starttime
     message("Output file or return value.")
     if(!is.null(reportfile)){
            write(paste("method:", method), reportfile,  append = FALSE)
            write(paste("s1/s2:", fws1s2), reportfile,  append = TRUE)
            write(paste("s3/s4:", fws3s4), reportfile,  append = TRUE)
            if(length(cpg[kk[,4]!=0,1])!=0){
              write(paste("N1(sum of MRE-CpG!=0 ):", N1), reportfile,  append = TRUE) 
              write(paste("N2(sum of MRE-CpG!=0 ):", N2), reportfile,  append = TRUE)
              write(paste("N3(sum of MRE-CpG!=0 ):", N3), reportfile,  append = TRUE) 
              write(paste("N4(sum of MRE-CpG!=0 ):", N4), reportfile,  append = TRUE) 
              write(paste("c1:", rmd21), reportfile,  append = TRUE)
              write(paste("c2:", rmd43), reportfile,  append = TRUE) 
            }
            write(paste("Number of windows:", length(pall)), reportfile,  append = TRUE)
            write(paste("Number of p<1e-8:", sum(pall<1e-8)), reportfile,  append = TRUE)
            write(paste("Number of p<1e-3:", sum(pall<1e-3)), reportfile,  append = TRUE)
            write(paste("Spend time:", spendtime), reportfile,  append = TRUE)       
     }

     cpg[,2]<-as.integer(cpg[,2])
     cpg[,3]<-as.integer(cpg[,3])
     ##
     cpgpq<-cbind(cpg[,1:3],sm,cpg[,4],kk[,4],pall,ts)
     colnames(cpgpq)=c("chr", "chrSt","chrEnd","Medip1","Medip2","MRE1","MRE2","cg","mrecg","pvalue",'Ts')
     if(is.null(writefile)){
             message("Return value.")
             return(cpgpq)
     }else{
             message("Output to writefile.")
             write.table(cpgpq, writefile,sep="\t", quote=FALSE,row.names =FALSE)
     }
     message("The End.")
}

