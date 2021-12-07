library(pomp)
library(tidyverse)
library(cowplot)


simulate(
  t0=0, times=seq(0, 20, by=1),
  params=c(birth=0.0002, death=0.0002, Beta=4, mu_EI=7, mu_IR=2, N0=1e6, I0=100),
  rprocess = euler(Csnippet("
    double rate[8], dN[8];
    rate[0] = birth; // birth rate
    rate[1] = rate[3] = rate[5] = rate[7] = death; // death rate
    rate[2] = Beta * I / N; // mu_SE
    rate[4] = mu_EI;
    rate[6] = mu_IR;
    dN[0] = rpois(rate[0]*N*dt);
    reulermultinom(2,S,&rate[1],dt,&dN[1]);
    reulermultinom(2,E,&rate[3],dt,&dN[3]);
    reulermultinom(2,I,&rate[5],dt,&dN[5]);
    reulermultinom(1,R,&rate[7],dt,&dN[7]);
    S += dN[0]-dN[1]-dN[2];
    E += dN[2]-dN[3]-dN[4];
    I += dN[4]-dN[5]-dN[6];
    R += dN[6]-dN[7];
    N += dN[0]-dN[1]-dN[3]-dN[5]-dN[7];"), delta.t=0.01),
  rinit = Csnippet("S=N0-I0; E=0; I=I0; R=0; N=N0;"),
  paramnames=c("birth", "death", "Beta", "mu_EI", "mu_IR", "N0", "I0"),
  statenames=c("S", "E", "I", "R", "N"),
  nsim=5
) -> em

em |>
  as.data.frame() ->df

df |>
  ggplot(aes(x=time,y=E,group=.id,color=factor(.id)))+
  geom_line()+
  guides(color="none")+
  expand_limits(y=0)+
  scale_y_sqrt()+
  labs(y=expression(I[t]),x=expression(t))
