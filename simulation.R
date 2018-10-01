#-------------------------------------------------------
# This file sets up a simulated migratory population of
# North Atlantic right whales that move up and down
# the coast following a sinusoidal curve.
#-------------------------------------------------------

library(tidyr); library(plyr); library(dplyr)
library(lubridate)
library(ggplot2)
library(leaflet)
library(mapview)

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

source("simple_Fourier.JAGS.R")

(start<-Sys.time())
out <- jags(jags_data,inits,params,"Fourier_model.txt",
            n.chains,n.adapt,n.iter,n.burnin,n.thin,parallel=T)
(end<-Sys.time()-start)


# examining the right whale sightings data

load("data/NARW_sightings_1900-2017.Rdata")
sightings <- right_whale_sightings
sightings$MONTH <- month(sightings$DATE)
sightings$YEAR <- year(sightings$DATE)

avg_lat <- sightings %>% filter(YEAR > 1997) %>% 
  group_by(YEAR,MONTH) %>% summarise(LAT = weighted.mean(LATDD,GROUPSIZE),n=n())
# add ordinal month
avg_lat$ordinal_month <- 1:nrow(avg_lat)

## Add to your plot
ggplot(avg_lat) + 
geom_point(aes(x = ordinal_month, y = LAT), size = 1) +
  geom_line(aes(x = ordinal_month, y = LAT))
#geom_line(data = ps, aes(x = x, y = y))


# provide a data frame to leaflet()

yr <- 2016

for (mo in 1:12){

m <- sightings %>% filter(YEAR == yr & MONTH == mo) %>%
  leaflet() %>% 
  #addProviderTiles(providers$Esri.WorldStreetMap) %>% 
  addProviderTiles(providers$Esri.OceanBasemap) %>% 
  fitBounds(-60,30,-75,50) %>%
  addCircleMarkers(~LONDD,~LATDD,radius=~GROUPSIZE, color = c("red"))

mapshot(m, file = paste0("sightings_,",yr,"_",mo,".png"))
}
  
  






jags_data <- list(y = scale(avg_lat$LAT),
                  t = avg_lat$ordinal_month,
                  nt = length(avg_lat$ordinal_month),
                  #delta = 6.94, 
                  #Pu = 7.6,
                  pi = 3.14159, P = 12)

params <- c(
  "a0", "Pu","delta","A.mu", "A.sig", "Pu.mu", "Pu.sig", "sig", "eps", "eps.sig","mu","A"
)

inits <- NULL

n.chains<-3
n.adapt<-5000  
n.iter <- 10000 
n.burnin<-0
n.thin<-1

source("simple_Fourier.JAGS.R")

(start<-Sys.time())
out <- jags(jags_data,inits,params,"Fourier_model.txt",
            n.chains,n.adapt,n.iter,n.burnin,n.thin,parallel=T,
            codaOnly = c("mu","A"))
(end<-Sys.time()-start)


plot(out$sims.list$delta,out$sims.list$Pu,pch=16,col=rgb(0,0,0,0.01))


plot(apply(out$sims.list$mu,2,median),jags_data$y)
abline(0,1,col="red")

plot(jags_data$t,apply(out$sims.list$mu,2,median),type="l")
lines(jags_data$t,jags_data$y,lty=2,col="red")
