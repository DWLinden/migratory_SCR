#-------------------------------------------------------
# This file sets up a simulated migratory population of
# North Atlantic right whales that move up and down
# the coast following a sinusoidal curve.
#-------------------------------------------------------

library(tidyr); library(plyr); library(dplyr)
library(lubridate)
library(ggplot2); library(ggthemes)
library(gganimate)
library(maps)
library(leaflet)
library(mapview)
library(jagsUI)

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



source("Fourier_func.R")

xp_t <- a0 + alpha*cos((2*pi/P)*t) + beta*sin((2*pi/P)*t)

eps <- rnorm(length(t),0,1)
xp_t <- a0 + xp_t.func(t,Pu0=6,A.sig=0,Pu.sig=0,delta.sig=0) + eps
  
plot(t,xp_t,type="l")





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
sightings$GROUPSIZE <- as.numeric(as.character(sightings$GROUPSIZE))
# fix weird group size
sightings$GROUPSIZE[sightings$GROUPSIZE==99999] <- 1
# add timing
sightings$MONTH <- month(sightings$DATE)
sightings$YEAR <- year(sightings$DATE)
sightings$YEAR_MONTH <- paste(sightings$YEAR,formatC(sightings$MONTH,width=2,flag="0"),sep="-")

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
for (yr in c(1996,1997,2016,2017)){
for (mo in 1:12){

# yr <- 1996
# mo <- 1
  
m <- sightings %>% filter(YEAR %in% yr & MONTH %in% mo) %>%
  leaflet() %>% 
  #addProviderTiles(providers$Esri.WorldStreetMap) %>% 
  addProviderTiles(providers$Esri.OceanBasemap) %>% 
  fitBounds(-60,30,-75,50) %>%
  addCircleMarkers(~LONDD,~LATDD,radius=~GROUPSIZE, color = c("red"))

mapshot(m, file = paste0("./maps/sightings_",yr,"_",formatC(mo,width=2,flag="0"),".png"))
}}
  
world <- ggplot() + borders("world",colour = "gray85",fill="gray80") + 
  borders("lakes",colour="white", fill="white")+
  borders("state",fill=NA,colour="gray90",lty=2)+
  theme_map()


  

# animation example
yr <- 2006:2016
mo <- 1:12

# gganimate example
plot.sightings <- sightings %>% filter(YEAR %in% yr & MONTH %in% mo)

p <- world + geom_point(data=plot.sightings,aes(x=LONDD,y=LATDD,size=as.numeric(GROUPSIZE))) +
  coord_fixed(xlim=c(-85,-55),ylim=c(25,52),ratio=1.3) +
  theme(legend.position="none",plot.title=element_text(size=rel(2.5))) +
  # transition_time(MONTH) +
  # labs(title = "2016; Month = {frame_time}")
  transition_states(YEAR_MONTH,transition_length=1,state_length=3) +
  labs(title = "{closest_state}",size=2)
animate(p, nframes = length(yr)*length(mo)*3)
#anim_save(filename="2016_sightings.gif")
anim_save(filename="NARW_sightings.gif")
#




jags_data <- list(y = scale(avg_lat$LAT),
                  t = avg_lat$ordinal_month,
                  nt = length(avg_lat$ordinal_month),
                  nyrs = length(unique(avg_lat$YEAR)),
                  YEAR = avg_lat$YEAR-min(avg_lat$YEAR)+1,
                  W.p=diag(2),
                  delta = 7, 
                  Pu = 8,
                  pi = 3.14159, P = 12)

params <- c(
  "a0", "Pu","delta","A.mu", "A.sig", "Pu.mu", "Pu.sig", 
  "rho12.p","sd.p","Sig.p",
  "sig", "eps", "eps.sig","mu","A"
)

inits <- function(){
  list(
    beta0 = c(2.1, 1.96)
    )}

n.chains<-3
n.adapt<-10000  
n.iter <- 20000 
n.burnin<-0
n.thin<-1

source("simple_Fourier.JAGS.R")

(start<-Sys.time())
out <- jags(jags_data,inits,params,"Fourier_model.txt",
            n.chains,n.adapt,n.iter,n.burnin,n.thin,parallel=T,
            codaOnly = c("mu","A","eps"))
(end<-Sys.time()-start)
gc()

#plot(out$sims.list$delta,out$sims.list$Pu,pch=16,col=rgb(0,0,0,0.01))

plot(apply(out$sims.list$mu,2,median),jags_data$y)
abline(0,1,col="red")

plot(jags_data$t,apply(out$sims.list$mu,2,median),type="l",
     ylim=range(jags_data$y))
lines(jags_data$t,jags_data$y,lty=2,col="red")

plot(avg_lat$MONTH,apply(out$sims.list$A,2,median))


avg_lat$A_hat <- apply(out$sims.list$A,2,median)
avg_lat$mu_hat    <- apply(out$sims.list$mu,2,median)

ggplot(avg_lat,aes(x=MONTH,y=mu_hat,col=YEAR,group=YEAR))+geom_line()+
  scale_color_gradient(low="yellow",high="blue")


