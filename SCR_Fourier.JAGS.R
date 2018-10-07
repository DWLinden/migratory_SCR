cat("
model {
    
  # Fourier priors
  #a0 ~ dnorm(0,.04)
  a0 <- 0
  sig ~ dt(0,0.16,1)T(0,)
  tau <- 1/(sig^2)
  A ~ dgamma(40, 2)
  
  # SCR priors
  pop.sigma ~ dgamma(1, 1)
  sigma ~ dgamma(2, 1)
  #for (t in 1:nT){
    lam0 ~ dunif(-5, 5)
  #}
  for (yr in 1:nyrs){
    beta0[yr] ~ dunif(-10, 10)
  }
  beta1 ~ dunif(-10, 10)
  
  # Fourier likelihood
  for (i in 1:P){
    for (t in 1:nT){
      xp_t.mat[i,t] <- ((2*A)/(pi*i))*sin(pi*i*Pu/P)*cos(((2*pi*i)/P)*(t-delta))
    }}
  
  for (t in 1:nT){
    LAT.mu[t] <- a0 + sum(xp_t.mat[,t])
    pop.LAT[t] ~ dnorm(LAT.mu[t],tau)
    pop.LON[t] ~ dunif(5,15)
  }
  
  # SCR likelihood
  for (t in 1:nT){
    
  for (g in 1:ng) {
    pop.dist[t,g] <- sqrt((ss[g,1]-pop.LON[t])^2 + (ss[g,2]-(pop.LAT[t]+LAT.adj))^2)
    mu[t,g] <- exp(beta0[yr[t]] - (1/(2*pop.sigma^2))*(pop.dist[t,g]^2))*pixArea
    #mu[t,g] <- exp(beta0[yr[t]])*pixArea
    probs[t,g] <- mu[t,g]/EN[t]
  }
  EN[t] <- sum(mu[t,])
  psi[t] <- EN[t]/M
  
  for(i in 1:M) {
    z[i,t] ~ dbern(psi[t])
    g.s[i,t] ~ dcat(probs[t,])
    x.s[i,t] <- ss[g.s[i,t],1]
    y.s[i,t] <- ss[g.s[i,t],2]
    for(j in trap.ind[jT[t],]) {
      dist[i,j,t] <- sqrt((x.s[i,t]-traps[j,1])^2 +
                          (y.s[i,t]-traps[j,2])^2)
      lambda[i,j,t] <- exp(lam0 - (dist[i,j,t]^2)/(2*(sigma^2))) * z[i,t]
      y[i,j,t] ~ dpois(lambda[i,j,t])
    }
  }
  N[t] <- sum(z[,t])
  
  }
  
  
    
}
",file="SCR_Fourier_model.txt")
