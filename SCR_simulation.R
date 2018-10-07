library(oSCR)
library(viridis)
library(ggplot2)
library(jagsUI)
library(raster)

# Simulate SCR data, with moving population center

x.extent <- c(0,20)
y.extent <- c(0,60)

# state space
ss <- expand.grid(x=seq(x.extent[1],x.extent[2],by=2),
                  y=seq(y.extent[1],y.extent[2],by=2))
plot(ss,pch=15,col="gray",asp=1,cex=1)
ss.r <- rasterFromXYZ(data.frame(ss,z=1))
pixArea <- prod(res(ss.r))

traps <- expand.grid(x=seq((10-4.5),(10+4.5),by=1.5),
                     y=c(seq((10-4.5),(10+4.5),by=1.5),
                         #seq((50-3),(50+3)),
                         seq((35-4.5),(35+4.5),by=1.5)))

ng <- nrow(ss)
nj <- nrow(traps)

# population size
N <- 50
# surveys within a month
K <- 5

# years
nyrs <- 5
# period of the cycle
P <- 2  #12 for real
# total time periods (4 seasons within year)
nT <- nyrs * P
t <- 1:nT

# frequency over period P
w0 <- (2*pi)/P

source("Fourier_func.R")

# when P=12 then delta0=7 and Pu0=8 matches the NARW sightings
xp_t <- xp_t.func(t = 1:60, P = 12, delta0 = 7, A0 = 20, Pu0 = 8, A.sig = 0)
ggplot(data.frame(x=1:60,y=xp_t),aes(x=factor(x),y=y)) + geom_line() + geom_point(col="red",size=2,pch=15)


#xp_t <- xp_t.func(t = t, P = 4, delta0 = 3, A0 = 65, Pu0 = 2, A.sig = 0)
xp_t <- xp_t.func(t = t, P = P, delta0 = 0, A0 = 20, Pu0 = 1, A.sig = 0)
LAT.adj <- 22.5
ggplot(data.frame(x=t,y=xp_t),aes(x=x,y=y+LAT.adj)) + geom_line() + geom_point(col="red",size=2,pch=15) +
  coord_cartesian(ylim=y.extent)

pop.LAT <- round(xp_t+LAT.adj,2)
pop.LON <- rep(10,nT)   # could make it stochastic: rnorm(nT,10,1)
pop.xy <- data.frame(x=pop.LON,y=pop.LAT)
pop.sigma <- 1

unique.LAT <- unique(pop.LAT)
unique.LAT.ind <- match(pop.LAT,unique.LAT)

pop.LAT.dist <- e2dist(ss,data.frame(x=pop.LON[1:length(unique.LAT)],y=unique.LAT))

plot(ss,pch=15,cex=.5,
     col=rev(viridis(10))[cut(exp((-1/(2*pop.sigma^2))*pop.LAT.dist[,2]),
                              breaks=c(seq(0,1,by=.1)))],asp=1)

mu.g <- exp(0.675 + (-1/(2*pop.sigma^2))*pop.LAT.dist) * pixArea
sum(mu.g[,1])
# grid cell probabilities for s
pr.g <- mu.g/N

# simulated s
s <- array(NA,dim=c(N,nT,2))
g.s <- matrix(NA,N,nT)
set.seed(123)

for (tt in 1:nT){
  # select the cells with individuals
  g.select <- rmultinom(1,N,prob=pr.g[,unique.LAT.ind[tt]])
  # distribute those cells across the N simulated individuals
  g.s[,tt] <- sample(rep(which(g.select>0),g.select[which(g.select>0)]),N)
  s[,tt,] <- as.matrix(ss[g.s[,tt],])
}
    
plot(ss,pch=15,col="gray90",asp=1)
points(s[,2,],pch=15,col=rgb(1,0,0,.2))


y <- array(NA, c(N, nj, nT)) # capture data
sigma <- 2.5 # scale parameter
lam0 <- -.5 # basal encounter rate
lam <- array(NA, c(N, nj, nT))
for (i in 1:N) {
  for (j in 1:nj) {
    for (tt in 1:nT){
      distSq <- (s[i,tt,1]-traps[j,1])^2 + (s[i,tt,2] - traps[j,2])^2
      lam[i,j,tt] <- exp(lam0 - distSq/(2*sigma^2))
      y[i,j,tt] <- rpois(1,lam[i,j,tt])
    }
  }
}
table(apply(y[,,1],1,sum))
sum(y)
    

M <- 100
n <- nrow(y) # not real n bc all-0s have not been removed

z.known <- apply(y,c(1,3),sum)
z.known[z.known>0] <- 1

z.known.aug <- rbind(z.known,matrix(NA,nrow=M-n,nT))
z.known.aug[z.known.aug==0] <- NA

z.inits <- (z.known+1)/(z.known+1)
z.inits[z.known==1] <- NA

z.inits.aug <- rbind(z.inits,matrix(0,M-n,nT))

y.aug <- array(NA, c(M, nj, nT))
y.aug[1:N,,] <- y



jags_data <- list(y = y.aug,
                  nT = nT,
                  nyrs = nyrs,
                  yr = sort(rep(1:nyrs,P)),
                  
                  ss = ss,
                  ng = ng,
                  traps = traps,
                  nj = nj,
                  jT = unique.LAT.ind,
                  trap.ind = rbind(c(1:49),c(50:98)),
                  pixArea = pixArea,
                  
                  M=M,
                  z=z.known.aug,
                  
                  delta = 0, 
                  Pu = 1,
                  pi = 3.14159, P = 2,
                  LAT.adj = LAT.adj)

params <- c(
  "a0","A", "sig", 
  "pop.sigma", "sigma",
  "N", "beta0", "beta1", "lam0", "psi",
  "pop.LAT","pop.LON"
)

inits <- function(){
  list(
    beta0 = rep(.65,nyrs),
    #beta0 = rep(-3.3,nyrs),
    lam0 = -.5,
    z = z.inits.aug,
    pop.sigma = 1,
    sigma = 2.5,
    #a0 = 22.5,
    A = 20,
    pop.LON = rep(10,nT),
    pop.LAT = rep(c(10,35),nyrs),
    g.s = rbind(g.s,g.s)
    
  )}

n.chains<-3
n.adapt<-1000 #1000  
n.iter <- 5000 #5000 
n.burnin<-0
n.thin<-1

source("SCR_Fourier.JAGS.R")

(start<-Sys.time())
out <- jags(jags_data,inits,params,"SCR_Fourier_model.txt",
            n.chains,n.adapt,n.iter,n.burnin,n.thin,parallel=n.chains>1,
            codaOnly = NULL)
(end<-Sys.time()-start)
gc()
