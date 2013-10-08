pmultinom <-
function(T, SIZE,length, P1, P2, P3, P4, C1, C2, pvalue=rep(0,length(T)))
 {
    problity<-.C("pcalculate", as.double(T),as.integer(SIZE),as.integer(length),as.double(P1),as.double(P2),
            as.double(P3),as.double(P4),as.double(C1),as.double(C2),as.double(pvalue))
    return(problity[[10]])
 }

