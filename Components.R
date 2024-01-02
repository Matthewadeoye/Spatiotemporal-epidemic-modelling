# Trend component (r_t)

r <- numeric(time)

kappa_r<- rgamma(1,10,0.001)
sigma_r <- sqrt(1/kappa_r)

# initial values
r[1] <- runif(1)
r[2] <- r[1]

for (t in 3:time) {
  epsilon <- rnorm(1, mean = 0, sd = sigma_r)  
  r[t] <- 2*(r[t - 1]) - r[t - 2] + epsilon
}

#Seasonal component (s_t)

s <- numeric(time)

kappa_s<- rgamma(1,10,0.009)
sigma_s <- sqrt(1/kappa_s)

initial_WN <- rnorm(time, mean = 0, sd = sigma_s)

for (t in 1:time) {
  if (t <= 11) {
    s[t] <- initial_WN[t]
  } else {
    s[t] <- rnorm(1, mean = -sum(s[(t-1):(t-11)]), sd = sigma_s)
  }
}


#Spatial component (u_i)
library(ar.matrix)

Q<- Westmidlands_adjmat
Q<- -1*Q
kappa_u<- rgamma(1,10,0.5)
sigma<- 1/sqrt(kappa_u)

diag(Q)<- -rowSums(Q, na.rm = T)
Q<- (1/sigma^2)*Q

u<- c(sim.GMRF(1,Q))
