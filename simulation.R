#-------------------------------------------------------
# This file sets up a simulated migratory population of
# North Atlantic right whales that move up and down
# the coast following a sinusoidal curve.
#-------------------------------------------------------

library(tidyr); library(plyr); library(dplyr)
library(ggplot2)

# lowest latitude = 30 degrees N; highest latitude = 50 degrees N 
latitudes <- runif(1000,30,50)

# number of years
nyrs <- 5

# period of the cycle (12 months)
P <- 12

# frequency over period P
w0 <- (2*pi)/P

# months in the observed (simulated) data
t <- 1:(12*nyrs)

a0 <- 40

a <- 2
# a <- rnorm(length(t),-10,5)
a <- rev(sort(runif(length(t),-20,-10)))

# b <- 2
# b <- rnorm(length(t),0,.1)
b <- seq(-10,-5,length.out = length(t))


xp_t <- a0 + a*cos(w0*t) + b*sin(w0*t)

plot(t,xp_t,type="l",ylim=c(20,60))

plot(a[c(seq(6,54,by=12))], xp_t[c(seq(6,54,by=12))])



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

xp_t <- a0 + alpha*cos((2*pi/P)*t) + beta*sin((2*pi/P)*t)

eps <- rnorm(length(t),0,1)
xp_t <- a0 + xp_t.func(t,Pu0=6,A.sig=0,Pu.sig=0,delta.sig=0) + eps
  
plot(t,xp_t,type="l")



library(jagsUI)

jags_data <- list(y = scale(xp_t),
             t = t, nt = length(t),
             delta = 6, pi = 3.14159,
             P = 12, Pu = 1)

params <- c(
  "a0", "A", "A.mu", "A.sig", "Pu.mu", "Pu.sig", "sig", "eps", "eps.sig"
)

inits <- NULL

n.chains<-3
n.adapt<-10000  #5000 5e3
n.iter <-20000 #50000 5e4
# new version: 608min for 10000
n.burnin<-0
n.thin<-1


(start<-Sys.time())
out <- jags(jags_data,inits,params,"Fourier_model.txt",
            n.chains,n.adapt,n.iter,n.burnin,n.thin,parallel=T)
(end<-Sys.time()-start)


cat("
model {

a0 ~ dunif(-10,10)
sig ~ dunif(0,10)
tau <- 1/(sig^2)

#A.mu ~ dnorm(0,0.001)
#A.sig ~ dunif(0,10)
A ~ dgamma(1, 1)

eps.sig ~ dunif(0,10)

#Pu.mu ~ dnorm(0,0.001)
#Pu.sig ~ dunif(0,10)

for (t in 1:nt){
  eps[t] ~ dnorm(0,1/(eps.sig^2))
  #A[t] ~ dnorm(A.mu,1/(A.sig^2))
  #Pu[t] ~ dnorm(Pu.mu,1/(Pu.sig^2))
}

for (i in 1:P){
  for (t in 1:nt){
    xp_t.mat[i,t] <- ((2*A)/(pi*i))*sin(pi*i*Pu/P)*cos(((2*pi*i)/P)*(t-delta))
}}

for (t in 1:nt){
  y[t] ~ dnorm(mu[t],tau)
  mu[t] <- a0 + sum(xp_t.mat[,t]) + eps[t]
}



}
",file="Fourier_model.txt")


# examining the right whale sightings data

load("data/NARW_sightings_1900-2017.Rdata")
right_whale_sightings$MONTH <- month(right_whale_sightings$DATE)
right_whale_sightings$YEAR <- year(right_whale_sightings$DATE)

avg_lat <- right_whale_sightings %>% filter(YEAR > 1997) %>% 
  group_by(YEAR,MONTH) %>% summarise(LAT = weighted.mean(LATDD,GROUPSIZE),n=n())
# add ordinal month
avg_lat$ordinal_month <- 1:nrow(avg_lat)

## Add to your plot
ggplot(avg_lat) + 
geom_point(aes(x = ordinal_month, y = LAT), size = 1) +
  geom_line(aes(x = ordinal_month, y = LAT))
#geom_line(data = ps, aes(x = x, y = y))
