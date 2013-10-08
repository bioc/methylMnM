calcFactornew <-
function(obs, ref, m, k, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10) {
  if( all(obs==ref) )
    return(1)
  m[m==0]<-1
  k[k==0]<-1
  nO <- sum(obs*(m/k))
  nR <- sum(ref*(m/k))
  logR <- log2((obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size
  absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression
  v <- (nO-obs*(m/k))/nO/obs + (nR-ref*(m/k))/nR/ref   # estimated asymptotic variance
  
  # remove infinite values, cutoff based on A
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]
  
  # taken from the original mean() function
  n <- sum(fin)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS
  qvalue.rank <- function(x) {
        idx <- sort.list(x)
        fc <- factor(x)
        nl <- length(levels(fc))
        bin <- as.integer(fc)
        tbl <- tabulate(bin)
        cs <- cumsum(tbl)
        tbl <- rep(cs, tbl)
        tbl[idx] <- tbl
        return(tbl)
  }
  
  keep <- (qvalue.rank(logR) %in% loL:hiL) & (qvalue.rank(absE) %in% loS:hiS)
  if (doWeighting) 
    2^( sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE) )
  else
    2^( mean(logR[keep], na.rm=TRUE) )
}