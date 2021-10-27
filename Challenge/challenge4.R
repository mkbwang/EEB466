# challenge 4 code

library(tidyverse)
library(pomp)

# Verhulst-Pearl model
vp <- pomp()
trajectory(
  t0=0, times=seq(from=0,to=10,by=0.01),
  skeleton = vectorfield(
    function (r, K, N, ...) {
      c(N=r*N*(1-N/K))
    }
  ),
  rinit = function (N0, ...) {
    c(N=N0)
  },
  params=c(r=1,K=100,N0=1)
) -> vp

trajectory(...,
  t0 = 0, times = seq(1,30,by=0.1),
  rinit = function (x0, ...) {
    c(x = x0)
  },
  skeleton = vectorfield(
    function (r, e, t, x, ...) {
      c(x=r*(1+e*cos(t))*x)
    }
  ),
  params = c(r=1,e=3,x0=1)
) -> po

# r=1, k=100
pars <- parmat(coef(vp),5)




