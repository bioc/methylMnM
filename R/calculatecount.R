calculatecount <-
function(data2, data3, cpg2, cpg3, datalength, cpglength, count=rep(0,cpglength)) 
 {
      bintags<-.C("calculatecount", as.integer(data2),as.integer(data3),as.integer(cpg2),as.integer(cpg3),
            as.integer(datalength),as.integer(cpglength),as.integer(count))
      return(bintags[[7]])
 }

