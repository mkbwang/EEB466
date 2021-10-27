## ----packages-----------------------------------------------------------------
library(tidyverse)
library(pomp)
stopifnot(getRversion()>="4.1")
stopifnot(packageVersion("pomp")>="4.0.4.1")


## ----install_pomp,eval=FALSE--------------------------------------------------
## install.packages("pomp",repos="https://kingaa.github.io/")


## ----vp1----------------------------------------------------------------------
trajectory(
  t0=0, times=seq(from=1,to=10,by=1),
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


## ----vp1_class----------------------------------------------------------------
class(vp)


## ----vp1_plot-----------------------------------------------------------------
plot(vp,type='o')


## ----vp2----------------------------------------------------------------------
pars <- parmat(coef(vp),5)
pars
pars["N0",] <- seq(1,200,length=5)
pars  

vp |>
  trajectory(
    params=pars,
    times=seq(0,20,by=0.1)
  ) -> dat

dat |>
  as.data.frame() |>
  ggplot(aes(x=time,y=N,color=factor(.id),group=.id))+
  geom_line()


## ----lf1----------------------------------------------------------------------
trajectory(
  t0=0, times=seq(0,10,by=0.01),
  skeleton = vectorfield(
    function (x1, x2, a11, a12, a21, a22, ...) {
      c(
        x1 = a11*x1 + a12*x2,
        x2 = a21*x1 + a22*x2
      )
    }
  ),
  rinit=function (x1_0, x2_0, ...) {
    c(x1=x1_0, x2=x2_0)
  },
  params=c(
    a11  = -1, a12  =  1,
    a21  = -1, a22  = -2,
    x1_0 =  5, x2_0 =  4
  )
) -> lf


## ----amat,include=FALSE-------------------------------------------------------
Amatrix <- array(
  coef(lf,c("a11","a21","a12","a22")),
  dim=c(2,2)
)
ev <- eigen(Amatrix)
ic <- coef(lf,c("x1_0","x2_0"))


## ----lf1_plot-----------------------------------------------------------------
lf |>
  as.data.frame() |>
  ggplot(aes(x=x1,y=x2))+
  geom_path()+
  labs(x=expression(x[1]),y=expression(x[2]))


## ----lf2----------------------------------------------------------------------
pars <- parmat(coef(lf),20)
pars["x1_0",1:10] <- seq(-5,5,length=5)
pars["x1_0",11:15] <- 5
pars["x1_0",16:20] <- -5
pars["x2_0",1:5] <- 5
pars["x2_0",6:10] <- -5
pars["x2_0",11:20] <- seq(-5,5,length=5)

lf |>
  trajectory(
    params=pars
  ) -> dat

dat |>
  as.data.frame() |>
  ggplot(aes(x=x1,y=x2,group=.id))+
  geom_path()+
  labs(x=expression(x[1]),y=expression(x[2]))

dat |>
  as.data.frame() |>
  pivot_longer(c(x1,x2)) |>
  ggplot(aes(x=time,y=value,group=.id,color=factor(.id)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name)


## ----sir1---------------------------------------------------------------------
trajectory(
  t0=0, times=seq(0,1,by=0.001),
  params = c(beta=60,gamma=20,N=1000,Sfrac=0.9,Ifrac=0.001),
  skeleton = vectorfield(
    function (S, I, R, beta, gamma, ...) {
      N <- S+I+R
      c(
        S = -beta*S*I/N,
        I = beta*S*I/N - gamma*I,
        R = gamma*I
      )
    }
  ),
  rinit = function (N, Sfrac, Ifrac, ...) {
    c(
      S = Sfrac*N,
      I = Ifrac*N,
      R = (1-Sfrac-Ifrac)*N
    )
  }
) -> sir


## ----sir1_plot1,fig.dim=c(4,3)------------------------------------------------
sir |> plot()


## ----sir1_plot2,fig.dim=c(4,3)------------------------------------------------
sir |>
  as.data.frame() |>
  ggplot(aes(x=S,y=I))+
  geom_path()


## ----rmca1--------------------------------------------------------------------
trajectory(
  t0=0,
  times=seq(0,200,by=0.1),
  rinit=function(n0, p0, ...) {
    c(n=n0,p=p0)
  },
  skeleton=vectorfield(
    function (n, p, alpha, beta, gamma, ...) {
      c(
        n=n*(1-n/gamma)-n*p/(1+n),
        p=beta*p*(n/(1+n)-alpha)
      )
    }
  ),
  params=c(n0=1,p0=0.1,alpha=0.8,beta=1,gamma=2)
) -> RmcA


## ----rmca1_plot---------------------------------------------------------------
RmcA |>
  as.data.frame() |>
  ggplot(aes(x=n,y=p))+
  geom_path()


## ----rmca2,out.height="0.5\\textheight"---------------------------------------
pars <- parmat(coef(RmcA),20)
pars["n0",1:10] <- seq(0,3,length=10)
pars["p0",1:10] <- 5
pars["n0",11:20] <- 3
pars["p0",11:20] <- seq(5,0,length=10)

RmcA |>
  trajectory(
    params=pars,
    format="data.frame"
  ) |>
  ggplot(aes(x=n,y=p,group=.id))+
  geom_path()+
  coord_fixed(ratio=0.5)


## ----rmca3,out.height="0.5\\textheight"---------------------------------------
pars["alpha",] <- 0.65
RmcA |>
  trajectory(
    params=pars,
    format="data.frame"
  ) |>
  ggplot(aes(x=n,y=p,group=.id))+
  geom_path()+
  coord_fixed(ratio=0.5)


## ----rmca4,out.height="0.5\\textheight"---------------------------------------
pars["alpha",] <- 0.5
RmcA |>
  trajectory(
    params=pars,
    format="data.frame"
  ) |>
  ggplot(aes(x=n,y=p,group=.id))+
  geom_path()+
  coord_fixed(ratio=0.5)


## ----rmca5,out.height="0.5\\textheight"---------------------------------------
pars["alpha",] <- 0.33
RmcA |>
  trajectory(
    params=pars,
    times=seq(0,1000,by=0.1),
    format="data.frame"
  ) |>
  ggplot(aes(x=n,y=p,group=.id))+
  geom_path()+
  coord_fixed(ratio=0.5)


## ----rmca6,out.height="0.5\\textheight"---------------------------------------
pars["alpha",] <- 0.1
RmcA |>
  trajectory(
    params=pars,
    format="data.frame"
  ) |>
  ggplot(aes(x=n,y=p,group=.id))+
  geom_path()+
  coord_fixed(ratio=0.5)


## ----rmca7,echo=FALSE---------------------------------------------------------
pars <- parmat(pars,6)
pars["alpha",] <- rep(c(0.8,0.65,0.5,0.33,0.2,0.1),each=20)
colnames(pars) <- sprintf(
  "%s_%02d",
  rep(letters[1:6],each=20),
  rep(1:20,times=6)
)

RmcA |>
  trajectory(
    params=pars,
    format="data.frame"
  ) |> 
  separate(.id,into=c("param","ic")) |>
  left_join(
    bind_cols(
      param=letters[1:6],
      alpha=c(0.8,0.65,0.5,0.33,0.2,0.1)
    ),
    by="param"
  ) |>
  ggplot(aes(x=n,y=p,group=interaction(param,ic)))+
  geom_path()+
  facet_wrap(
    ~alpha,
    labeller=label_bquote(alpha==.(alpha))
  ) -> pl

pl

pl+
  scale_x_log10(limits=c(1e-8,NA))+
  scale_y_log10(limits=c(1e-2,NA))


## ----hp1----------------------------------------------------------------------
trajectory(
  t0=0,
  times=seq(0,200,by=0.1),
  rinit=function(x0, y0, z0, ...) {
    c(x=x0,y=y0,z=z0)
  },
  skeleton=vectorfield(
    function (x, y, z, a1, a2, b1, b2, d1, d2, ...) {
      f1y <- a1*x*y/(1+b1*x)
      f2z <- a2*y*z/(1+b2*y)
      c(
        x=x*(1-x)-f1y,
        y=f1y-f2z-d1*y,
        z=f2z-d2*z
      )
    }
  ),
  params=c(
    a1=5,a2=0.1,
    b1=3,b2=2,
    d1=0.4,d2=0.01,
    x0=0.7,y0=0.15,z0=9.8
  )
) -> HP


## ----hp2----------------------------------------------------------------------
HP |>
  trajectory(
    times=seq(0,6500,by=0.1),
    format="data.frame"
  ) |>
  filter(time>=4000) -> dat


## ----hp2a,fig.dim=c(6,6),echo=FALSE,dev="CairoPNG"----------------------------
library(cowplot)
plot_grid(
  dat |>
    filter(time>=5000) |>
    gather(var,val,x,y,z) |>
    ggplot(aes(x=time,y=val))+
    geom_line()+
    labs(y="")+
    facet_grid(var~.,scales="free_y"),
  dat |>
    ggplot(aes(x=x,y=y))+
    geom_path(),
  dat |>
    ggplot(aes(x=y,y=z))+
    geom_path(),
  dat |>
    ggplot(aes(x=x,y=z))+
    geom_path(),
  ncol=2
)


## ----hp2b,echo=FALSE,dev="CairoPNG"-------------------------------------------
library(plot3D)
with(
  dat,
  lines3D(x,y,z,colvar=NULL,theta=-40,phi=20,r=10)
)

