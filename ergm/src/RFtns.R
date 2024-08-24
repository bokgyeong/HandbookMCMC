compuhat = function(Sx, Sy, wIS){
  return(Sx - sapply(1:ncol(Sy), function(i) sum(wIS * Sy[,i])))
}

# f_Uhat = function(Sx, Sy){
#   return( as.vector(Sx - colMeans(Sy)) )
# }

compdhat = function(Sx, Sy, wIS){
  eS1 = sum(wIS * Sy[,1])
  eS2 = sum(wIS * Sy[,2])
  eS1S1 = sum(wIS * Sy[,1] * Sy[,1])
  eS2S2 = sum(wIS * Sy[,2] * Sy[,2])
  eS1S2 = sum(wIS * Sy[,1] * Sy[,2])
  
  H11 = - eS1S1 + eS1 * eS1
  H12 = - eS1S2 + eS1 * eS2
  H22 = - eS2S2 + eS2 * eS2
  
  J11 = Sx[1] * Sx[1] -  2 * Sx[1] * eS1 + eS1 * eS1
  J12 = Sx[1] * Sx[2] - Sx[1] * eS2 - Sx[2] * eS1 + eS1 * eS2
  J22 = Sx[2] * Sx[2] - 2 * Sx[2] * eS2 + eS2 * eS2
  
  return( as.vector(c(H11, H12, H22) + c(J11, J12, J22)) )
}


# f_dhat = function(Sx, Sy){
#   eS1 = mean(Sy[,1])
#   eS2 = mean(Sy[,2])
#   eS1S1 = mean(Sy[,1] * Sy[,1])
#   eS2S2 = mean(Sy[,2] * Sy[,2])
#   eS1S2 = mean(Sy[,1] * Sy[,2])
#   
#   H11 = - eS1S1 + eS1 * eS1
#   H12 = - eS1S2 + eS1 * eS2
#   H22 = - eS2S2 + eS2 * eS2
#   
#   J11 = Sx[1] * Sx[1] -  2 * Sx[1] * eS1 + eS1 * eS1
#   J12 = Sx[1] * Sx[2] - Sx[1] * eS2 - Sx[2] * eS1 + eS1 * eS2
#   J22 = Sx[2] * Sx[2] - 2 * Sx[2] * eS2 + eS2 * eS2
#   
#   return( as.vector(c(H11, H12, H22) + c(J11, J12, J22)) )
# }


compVhat = function(d, N){
  n = nrow(d)
  # b = floor(min(n^(1/3), N^(1/3)))
  b = floor(min(n^(1/3), N^(1/2-0.1)))
  a = floor(n/b)
  
  dbarbk = sapply(1:a, function(k) return(colMeans(d[((k - 1) * b + 1):(k * b),])))
  dbar = rowMeans(dbarbk)
  dummy = 0
  for(k in 1:a){
    dummy = dummy + (dbarbk[,k] - dbar) %*% t(dbarbk[,k] - dbar)
  }
  Sigmahat = b * dummy / (a-1)
  
  return( Sigmahat )
}


compACD = function(d, N){
  n = nrow(d)
  dbar = colMeans(d)
  res = n * t(dbar) %*% solve(compVhat(d, N)) %*% dbar
  return(as.vector(res))
}





