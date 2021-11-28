## ----packages,include=FALSE,cache=FALSE,purl=TRUE-----------------------------
library(pomp)
library(tidyverse)
library(cowplot)
stopifnot(getRversion()>="4.1")
stopifnot(packageVersion("pomp")>="4.0.11")
set.seed(97849913)




## ----wiener-------------------------------------------------------------------
wiener <- function (times, sigma = 1, n = 1) {
  dt <- diff(c(0,times))
  dW <- rnorm(n=n*length(times),mean=0,sd=sigma*sqrt(dt))
  dW |>
    matrix(nrow=length(times)) |>
    apply(2,cumsum) |>
    as.numeric() -> W
  data.frame(
    t=times,
    W=W,
    .id=rep(seq_len(n),each=length(times))
  )
}


## ----wiener_plot,echo=FALSE---------------------------------------------------
t <- c(seq(0,0.1,by=1e-5),seq(0.11,10,by=0.01))
dat <- wiener(t,sigma=1)
dat |>
  filter(t<=0.1) |>
  ggplot(aes(x=t,y=W))+
  geom_line()+
  labs(y=expression(W[t])) -> pl1
dat |>
  ggplot(aes(x=t,y=W))+
  geom_line()+
  labs(y=expression(W[t])) -> pl2
wiener(times=seq(0,10,by=0.01),sigma=1,n=20) |>
  ggplot(aes(x=t,y=W,group=.id,color=factor(.id)))+
  geom_line()+
  guides(color="none")+
  labs(y=expression(W[t])) -> pl3
plot_grid(pl2,pl1,ncol=1,labels=c("A","B")) -> plL
plot_grid(plL,pl3,nrow=1,labels=c(NA,"C"))


## ----gbm1,fig.dim=c(8,3)------------------------------------------------------
wiener(times=seq(0,10,by=0.01),sigma=1,n=20) |>
  mutate(X=exp(W-t/2)) |>
  ggplot(aes(x=t,y=X,group=.id,color=factor(.id)))+
  geom_line()+
  guides(color="none")+
  labs(y=expression(X[t]))


## ----gbm_ex1,include=FALSE,purl=TRUE------------------------------------------
wiener(
  times=exp(seq(from=log(1),to=log(20),length=20)),
  sigma=1,n=1e5
) |>
  mutate(X=1000*exp(W-t/2)) -> dat
dat |>
  group_by(t,gp=.id%%10) |>
  summarize(
    mean=mean(X)
  ) |>
  ungroup() |>
  ggplot(aes(x=t,y=mean,group=gp))+
  geom_point()+
  geom_line()+
  scale_x_log10()
dat |>
  ggplot(aes(x=log10(X),group=t,fill=t,color=t))+
  geom_density(alpha=0.1)


## ----ou-----------------------------------------------------------------------
ou <- function (times, alpha, sigma = 1, x0 = 0, n = 1) {
  t0 <- 0
  x <- rep(x0,each=n)
  n <- length(x)
  X <- array(
    dim=c(length(times),n),
    dimnames=list(time=NULL,rep=seq_len(n))
  )
  for (k in seq_along(times)) {
    t <- times[k]
    M <- x*exp(-alpha*(t-t0))
    V <- (1-exp(-2*alpha*(t-t0)))/2/alpha*sigma*sigma
    x <- X[k,] <- rnorm(n=n,mean=M,sd=sqrt(V))
    t0 <- t
  }
  cbind(t=times,as.data.frame(X)) |>
    pivot_longer(-t,names_to=".id",values_to="X")
}

## ----ou_plot,fig.dim=c(8,3),echo=FALSE----------------------------------------
ou(times=seq(0,10,by=0.1),alpha=1,x0=10,n=10) |>
  ggplot(aes(x=t,y=X,group=.id,color=factor(.id)))+
  geom_line()+
  guides(color="none")+
  labs(y=expression(X[t])) -> pl1
ou(times=seq(0,10,by=0.1),alpha=1,x0=seq(-10,10,by=2),n=5) |>
  ggplot(aes(x=t,y=X,group=.id,color=factor(.id)))+
  geom_line()+
  guides(color="none")+
  labs(y=expression(X[t])) -> pl2
plot_grid(pl1,pl2,nrow=1,labels="AUTO")


## ----ou_em,fig.dim=c(8,3)-----------------------------------------------------
ou_em <- function (delta.t, times, alpha, sigma = 1, x0 = 0, n = 1) {
  t0 <- 0
  x <- rep(x0,each=n)
  n <- length(x)
  X <- array(dim=c(length(times),n))
  for (k in seq_along(times)) {
    t <- times[k]
    nstep <- ceiling((t-t0)/delta.t)
    dt <- (t-t0)/nstep
    for (i in seq_len(nstep)) {
      dw <- rnorm(n=n,mean=0,sd=sqrt(dt))
      x <- x - alpha*x*dt + sigma*dw
    }
    X[k,] <- x
    t0 <- t
  }
  data.frame(
    t=times,
    X=as.numeric(X),
    .id=rep(seq_len(n),each=length(times))
  )
}

## ----<ou_em_plot,fig.dim=c(8,3),echo=FALSE------------------------------------
ou_em(delta.t=0.05,times=seq(0,10,by=0.1),alpha=1,x0=10,n=10) |>
  ggplot(aes(x=t,y=X,group=.id,color=factor(.id)))+
  geom_line()+
  guides(color="none")


## ----ou_em_pomp---------------------------------------------------------------
simulate(
  times=seq(0,10,by=0.1),
  t0=0,
  rinit=Csnippet("X = x0;"),
  rprocess=euler(
    Csnippet("X += -alpha*X*dt + sigma*rnorm(0,sqrt(dt));"),
    delta.t=0.05
  ),
  paramnames=c("alpha","sigma","x0"),
  statenames=c("X"),
  params=c(x0=10,alpha=1,sigma=1)
) -> oup

## ----ou_em_pomp_plot,echo=FALSE,fig.dim=c(7,3)--------------------------------
oup |>
  simulate(
    params={
      p <- parmat(coef(oup),6,names=letters[1:6])
      p["x0",] <- seq(-10,10,length=6)
      p
    },
    nsim=20,
    format="d"
  ) |>
  separate(.id,into=c("x0",".id")) |>
  ggplot(aes(x=time,y=X,group=interaction(.id,x0),color=factor(x0)))+
  guides(color="none")+
  geom_line()


## ----ou_density---------------------------------------------------------------
oup |>
simulate(
  times=c(0,1,2,4,8,16,32),
  params={
    expand_grid(
      x0=runif(n=1e5,min=10,max=20),
      alpha=1,
      sigma=1
    ) %>% t()
  },
  format="data.frame",
  nsim=1
) |>
  ggplot(aes(x=X))+
  geom_histogram(
    aes(y=..density..),color=NA,fill=grey(0.5),
    binwidth=0.1,boundary=20
  )+
  facet_grid(
    time~.,
    scales="free_y",
    labeller=label_bquote(rows=t==.(time))
  )


## ----sl1,fig.dim=c(8,3)-------------------------------------------------------
simulate(
  t0 = 0, times=seq(0,100,by=0.1),
  rinit=Csnippet("X = log(N0);"),
  rprocess=euler(
    Csnippet("X += (1-exp(X)-sigma*sigma/2)*dt+sigma*rnorm(0,sqrt(dt));"),
    delta.t=0.01
  ),
  paramnames=c("sigma","N0"), statenames=c("X"),
  params=c(N0=0.1,sigma=0.1)
) -> sl

sl |>
  as.data.frame() |>
  ggplot(aes(x=time,y=exp(X)))+
  labs(y=expression(N[t]),x=expression(t))+
  geom_line()


## ----stochlog2,echo=FALSE-----------------------------------------------------
expand_grid(
  N0=1,
  sigma=c(0,0.1,0.2,0.4,1,1.5)
) -> params
sl |>
  simulate(params=t(params),format="d") |>
  ggplot(aes(x=time,y=exp(X)))+
  geom_line()+
  labs(y=expression(N[t]),x=expression(t))+
  facet_wrap(~.id,labeller=label_bquote(sigma==.(params$sigma[.id])))

