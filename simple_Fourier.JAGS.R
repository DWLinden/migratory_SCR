cat("
model {
    
  a0 ~ dunif(-10,10)
  sig ~ dt(0,0.16,1)T(0,)
  tau <- 1/(sig^2)
  
  A.mu ~ dnorm(0,0.001)
  A.sig ~ dunif(0,10)
  #A.sig ~ dt(0,0.16,1)T(0,)
  
  #A ~ dgamma(1, 1)
  
  #eps.sig ~ dunif(0,10)
  
  Pu ~ dunif(4,9)
  #delta ~ dunif(7,8)
  
  #Pu.mu ~ dnorm(0,0.001)
  #Pu.sig ~ dunif(0,10)
  
  for (t in 1:nt){
    #eps[t] ~ dnorm(0,1/(eps.sig^2))
    A[t] ~ dnorm(A.mu,1/(A.sig^2))
    #Pu[t] ~ dnorm(Pu.mu,1/(Pu.sig^2))
  }

  # for (yr in 1:nyrs){
  #   A[yr] ~ dnorm(A.mu,1/(A.sig^2))
  # }
  
  for (i in 1:P){
    for (t in 1:nt){
      xp_t.mat[i,t] <- ((2*A[t])/(pi*i))*sin(pi*i*Pu/P)*cos(((2*pi*i)/P)*(t-delta))
    }}
  
  for (t in 1:nt){
    y[t] ~ dnorm(mu[t],tau)
    mu[t] <- a0 + sum(xp_t.mat[,t]) #+ eps[t]
  }
  
  
    
    }
    ",file="Fourier_model.txt")
