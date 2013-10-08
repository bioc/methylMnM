normpdft1 <-
function(t,n,p,c1,c2)
{ 
    varmutb<-matrix(c(p[1]*(1-p[1]), -p[1]*p[2]   ,  -p[1]*p[3]   ,
                     -p[1]*p[2]   ,p[2]*(1-p[2]),  -p[2]*p[3]   ,
                     -p[1]*p[3]   , -p[2]*p[3]   , p[3]*(1-p[3])),3,3)
    fxy<-c(c1*(1-2*p[1]-p[2]-p[3]),-(c1*p[1]+c2*p[3]),-(c1*p[1]+c2*p[2]))
    sigma<-t(fxy)%*%varmutb%*%fxy
    mu<-c1*p[1]*p[4]-c2*p[2]*p[3]
    t1<-sqrt(n)*(t-mu)/sqrt(sigma)
    #pvalue<-2*(1-pnorm(t1, mean = 0, sd = 1))
    return(t1)
}

