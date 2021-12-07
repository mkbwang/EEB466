## ----packages,include=FALSE,cache=FALSE,purl=TRUE-----------------------------
library(pomp)
library(tidyverse)
library(cowplot)
stopifnot(getRversion()>="4.1")
stopifnot(packageVersion("pomp")>="4.0.11")
set.seed(97849913)




## ----ito1---------------------------------------------------------------------
simulate(
  t0=0, times=seq(0,5,by=0.1),
  params=c(Beta=10,gamma=4,N=1e6,I0=10,B=2e4,delta=0.02),
  rprocess = euler(Csnippet("
    double dN[6];
    double lambda = Beta*I/N;
    dN[0] = rnorm(B*dt,sqrt(B*dt));               // births
    dN[1] = rnorm(lambda*S*dt,sqrt(lambda*S*dt)); // infections
    dN[2] = rnorm(delta*S*dt,sqrt(delta*S*dt));   // S deaths
    dN[3] = rnorm(gamma*I*dt,sqrt(gamma*I*dt));   // recoveries
    dN[4] = rnorm(delta*I*dt,sqrt(delta*I*dt));   // I deaths
    dN[5] = rnorm(delta*R*dt,sqrt(delta*R*dt));   // R deaths
    S += dN[0]-dN[1]-dN[2];
    I += dN[1]-dN[3]-dN[4];
    R += dN[3]-dN[5];"), delta.t=0.1),
  rinit = Csnippet("S=nearbyint(N-I0); I = nearbyint(I0); R = 0;"),
  paramnames=c("Beta","gamma","I0","N","B","delta"),
  statenames=c("S","I","R"),
  nsim=20
) -> ito

## ----ito2,dependson="ito1",fig.dim=c(6,3),echo=FALSE--------------------------
ito |>
  as.data.frame() |>
  ggplot(aes(x=time,y=I,group=.id,color=factor(.id)))+
  geom_line()+
  guides(color="none")+
  expand_limits(y=0)+
  scale_y_sqrt()+
  labs(y=expression(I[t]),x=expression(t))


## ----poisson1-----------------------------------------------------------------
ito[[1]] |>
  simulate(
    rprocess = euler(Csnippet("
      double dN[6];
      double lambda = Beta*I/N;
      dN[0] = rpois(B*dt);        // births
      dN[1] = rpois(lambda*S*dt); // infections
      dN[2] = rpois(delta*S*dt);  // S deaths
      dN[3] = rpois(gamma*I*dt);  // recoveries
      dN[4] = rpois(delta*I*dt);  // I deaths
      dN[5] = rpois(delta*R*dt);  // deaths
      S += dN[0]-dN[1]-dN[2];
      I += dN[1]-dN[3]-dN[4];
      R += dN[3]-dN[5];"), delta.t=0.1),
    params=c(Beta=10,gamma=4,N=1000,I0=1,delta=0.02,B=20),
    paramnames=c("Beta","gamma","N","B","delta"),
    statenames=c("S","I","R"),
    nsim=20
  ) -> pois

## ----poisson2,dependson="poisson1",fig.dim=c(6,3),echo=FALSE------------------
pois |>
  as.data.frame() |>
  ggplot(aes(x=time,y=I,group=.id,color=factor(.id)))+
  geom_line()+
  guides(color="none")+
  expand_limits(y=0)+
  scale_y_sqrt()+
  labs(y=expression(I[t]),x=expression(t))


## ----eulermult1---------------------------------------------------------------
ito[[1]] |>
  simulate(
    rprocess = euler(Csnippet("
      double rate[6], dN[6];
      rate[0] = B;
      rate[1] = Beta*I/N;
      rate[2] = rate[4] = rate[5] = delta;
      rate[3] = gamma;
      dN[0] = rpois(rate[0]*dt);
      reulermultinom(2,S,&rate[1],dt,&dN[1]);
      reulermultinom(2,I,&rate[3],dt,&dN[3]);
      reulermultinom(1,R,&rate[5],dt,&dN[5]);
      S += dN[0]-dN[1]-dN[2];
      I += dN[1]-dN[3]-dN[4];
      R += dN[3]-dN[5];"), delta.t=0.1),
    params=c(Beta=10,gamma=4,N=1000,I0=1,delta=0.02,B=20),
    paramnames=c("Beta","gamma","N","B","delta"),
    statenames=c("S","I","R"),
    nsim=20
  ) -> em

## ----eulermult2,dependson="eulermult1",fig.dim=c(6,3),echo=FALSE--------------
em |>
  as.data.frame() |>
  ggplot(aes(x=time,y=I,group=.id,color=factor(.id)))+
  geom_line()+
  guides(color="none")+
  expand_limits(y=0)+
  scale_y_sqrt()+
  labs(y=expression(I[t]),x=expression(t))


## ----comparison1,echo=FALSE---------------------------------------------------
bind_rows(
  deterministic=ito[[1]] |>
    trajectory(
      skeleton=vectorfield(Csnippet("
        DS = B - Beta*S*I/N - delta*S;
        DI = Beta*S*I/N - gamma*I - delta*I;
        DR = gamma*I - delta*I;")),
      params=c(Beta=10,gamma=4,N=100,I0=1,delta=0.02,B=20),
      paramnames=c("Beta","gamma","N","B","delta"),
      statenames=c("S","I","R"),
      format="d"
    ) |>
    mutate(.id=as.character(.id)),
  Ito=ito[[1]] |>
    simulate(
      params=c(Beta=10,gamma=4,N=100,I0=1,delta=0.02,B=20),
      nsim=20,format="d"),
  Poisson=pois[[1]] |>
    simulate(
      params=c(Beta=10,gamma=4,N=100,I0=1,delta=0.02,B=20),
      nsim=20,format="d"),
  Eulermultinomial=em[[1]] |>
    simulate(
      params=c(Beta=10,gamma=4,N=100,I0=1,delta=0.02,B=20),
      nsim=20,format="d"),
  .id="method"
) -> dat1

## ----comparison1B,dependson="comparison1",echo=FALSE,fig.dim=c(6,5)-----------
dat1 |>
  select(method,time,I,.id) |>
  mutate(I=coalesce(I,0)) |>
  ggplot(aes(x=time,y=I,group=.id,color=factor(.id)))+
  geom_line()+
  guides(color="none")+
  facet_grid(
    factor(method,levels=c("deterministic","Ito","Poisson","Eulermultinomial"))~.
  )


## ----comparison2,echo=FALSE---------------------------------------------------
bind_rows(
  Ito=ito[[1]] |>
    simulate(
      times=c(0.5,2,5),
      params=c(Beta=10,gamma=4,N=100,I0=1,delta=0.02,B=20),
      nsim=10000,format="d"),
  Poisson=pois[[1]] |>
    simulate(
      times=c(0.5,2,5),
      params=c(Beta=10,gamma=4,N=100,I0=1,delta=0.02,B=20),
      nsim=10000,format="d"),
  Eulermultinomial=em[[1]] |>
    simulate(
      times=c(0.5,2,5),
      params=c(Beta=10,gamma=4,N=100,I0=1,delta=0.02,B=20),
      nsim=10000,format="d"),
  .id="method"
) -> dat2

## ----comparison2B,dependson="comparison2",echo=FALSE,fig.dim=c(6,4.5)---------
dat2 |>
  select(method,time,I) |>
  mutate(
    I=coalesce(I,-1)
  ) |>
  group_by(method,time) |>
  summarize(
    p=1-sum(I>0)/length(I),
    q=sum(I<0)/length(I)
  ) |>
  ungroup() -> extinct
dat2 |>
  filter(I>0) |>
  ggplot(aes(x=I))+
  geom_histogram(
    aes(y=..density..),
    color=NA,fill="#00274c",
    binwidth=5,boundary=0
  )+
  geom_text(
    data=extinct,
    mapping=aes(label=formatC(signif(p,4),digits=4,flag="#")),
    x=40,y=0.1,color="black"
  )+
  geom_text(
    data=extinct,
    mapping=aes(label=formatC(signif(q,4),digits=4,flag="#")),
    x=40,y=0.05,color="red"
  )+
  labs(x=expression(I[t]~"|"~I[t]>0))+
  facet_grid(
    factor(method,levels=c("Ito","Poisson","Eulermultinomial"))~time,
    labeller=labeller(
      time=label_both,
      method=label_value
    )
  )

