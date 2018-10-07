xp_t.func <- function(t,P=12,
                      delta0=6,delta.sig=0,
                      A0=10,A.sig=0,
                      Pu0=4,Pu.sig=0){
  
  A <- rnorm(length(t),A0,A.sig)
  Pu <- rnorm(length(t),Pu0,Pu.sig)
  delta <- rnorm(length(t),delta0,delta.sig)
  
  xp_t.mat <- matrix(NA,P,length(t))
  
  for(p in 1:P){
    xp_t.mat[p,] <- ((2*A)/(pi*p))*sin(pi*p*Pu/P)*cos(((2*pi*p)/P)*(t-delta))
  }
  return(colSums(xp_t.mat))
  
  
}