library(pracma)

K=200
r=1
alpha=1

verhulstDE <- function(t, N){
  r*N*(1-N/K)
}

verhulstsol <- function(t, N0){
  exp(r*t)*K*N0/(K+N0*(exp(r*t) - 1))
}

vectorfield(verhulstDE, xlim=c(0, 6), ylim=c(0, 200), scale=0.2, col='black')
mytime <- seq(0, 6, 0.01)
for (start in seq(10, 150, 20)){
  trial <- verhulstsol(mytime, start)
  lines(mytime, trial, col='blue')
}

kompertzDE <- function(t, N){
  alpha*N*log(K/N)
}

kompertz <- function(t, N0){
  N0^(exp(-alpha*t)) * K^(1-exp(-alpha*t))
}

vectorfield(kompertzDE, xlim=c(0.1, 6), ylim=c(0, 200), scale=0.2, col='black')

for (start in seq(5, 140, 20)){
  trial <- kompertz(mytime, start)
  lines(mytime, trial, col='blue')
}

r0=1
r1=1

v_nonauto_DE <- function(t, N){
  r0*(1+r1*cos(2*pi*t))*N*(1-N/K)
}

vectorfield(v_nonauto_DE, xlim=c(0.1, 6), ylim=c(0, 200), scale=0.2, col='black')


v_nonauto <- function(t, N0){
  N0 * K * exp(r0*t+r0*r1/2/pi*sin(2*pi*t)) / (K + N0 * (exp(r0*t+r0*r1/2/pi*sin(2*pi*t)) - 1))
}

for (start in seq(5, 140, 20)){
  trial <- v_nonauto(mytime, start)
  lines(mytime, trial, col='blue')
}

