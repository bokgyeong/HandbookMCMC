f_uhat = function(Sx, Sy){
  return(Sx - mean(Sy))
}

f_dhat = function(Sx, Sy){
  Sybar = mean(Sy)
  uhat = Sx - Sybar
  return( -mean(Sy^2) + Sybar * Sybar + uhat * uhat)
}

f_Vhat = function(D, N){
  n = length(D)
  b = floor(min(n^(1/3), N^(1/3)))
  a = floor(n/b)

  dbarbk = sapply(1:a, function(k) return(mean(D[((k - 1) * b + 1):(k * b)])))
  dbar = mean(dbarbk)
  sigmahat = b * sum((dbarbk - dbar)^2) / (a-1)

  return( sigmahat )
}

ACD = function(D, N){
  n = length(D)
  dbar = mean(D)
  
  res = n * (dbar^2) / f_Vhat(D, N)
  return(res)
}



