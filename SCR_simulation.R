# Simulate SCR data, with moving population center


dist <- 0:200
sigma <- 5

plot(dist,.5*exp(-(1/(sigma^2))*dist))
