## ----packages,cache=FALSE,message=FALSE,include=FALSE-------------------------
library(tidyverse)
library(pomp)
library(subplex)
stopifnot(getRversion()>="4.1")
stopifnot(packageVersion("pomp")>="4.0")
library(cowplot)
theme_set(theme_bw(base_size=11,base_family="Helvetica"))
options(scipen=2)
set.seed(228088786)


## ----reset_seed,include=FALSE-------------------------------------------------
set.seed(228088786)


## ----kilpisjarvi1,echo=FALSE,fig.dim=c(6,2)-----------------------------------
library(tidyverse)

read_csv("data/kilpisjarvi_fall_c_rufocanus.csv",comment="#") -> dat

dat |>
  ggplot(aes(x=year,y=ti))+
  geom_line()+
  labs(y="trap index")


## ----kilp1_fit,eval=FALSE,purl=TRUE-------------------------------------------
library(tidyverse)
library(pomp)

read_csv("data/kilpisjarvi_fall_c_rufocanus.csv",comment="#") -> dat

dat |>
  filter(year<=1957) |>
  traj_objfun(
    times="year", t0=1949,
    rinit=Csnippet("
    x = log(N0);
    y = log(P0);
    "),
    skeleton=vectorfield(
      Csnippet("
    double fr = c/(1+a*exp(x));
    Dx = r*(1-exp(x)/K) - fr*exp(y);
    Dy = b*fr*exp(x) - m;
    ")
    ),
    dmeasure=Csnippet("
    lik = -(ti-exp(x))*(ti-exp(x));
    if (!give_log) lik = exp(lik);
    "),
    partrans=parameter_trans(
      log=c("b","r","m","P0")
    ),
    statenames=c("x","y"),
    paramnames=c("a","b","c","r","K","m","N0","P0"),
    params=c(K=25,a=0,r=2,c=1,m=1,b=1,N0=20,P0=1),
    est=c("b","r","m","P0")
  ) -> kilp1

library(subplex)
fit <- subplex(
  fn=kilp1,
  par=coef(kilp1,c("b","r","m","P0"),transform=TRUE)
) # changing the parameter of pomp object in place
kilp1(fit$par); coef(kilp1)



## ----kilp1_plot,echo=FALSE,fig.dim=c(6,4)-------------------------------------
kilp1 |> trajectory(format="d") -> traj1
traj1 |>
  ggplot(aes(x=year,y=exp(x)))+
  geom_line()+
  geom_point(data=dat |> filter(year<=1957),aes(y=ti))+
  labs(y=expression(N==exp(x))) -> pl1
traj1 |>
  ggplot(aes(x=year,y=exp(y)))+
  geom_line()+
  labs(y=expression(P==exp(y))) -> pl2

library(cowplot)
plot_grid(pl1,pl2,labels="auto",align="hv",axis="tlbr",ncol=1)


## ----kilp2_fit,eval=FALSE-----------------------------------------------------
dat |>
  traj_objfun(
    times="year", t0=1949,
    rinit=Csnippet("
      x = log(N0);
      y = log(P0);"),
    skeleton=vectorfield(
      Csnippet("
      double fr = c/(1+a*exp(x));
      Dx = r*(1-exp(x)/K) - fr*exp(y);
      Dy = b*fr*exp(x) - m;")
    ),
    dmeasure=Csnippet("
      lik = -(ti-exp(x))*(ti-exp(x));
      if (!give_log) lik = exp(lik);"),
    partrans=parameter_trans(
      log=c("b","r","m","P0")
    ),
    statenames=c("x","y"),
    paramnames=c("a","b","c","r","K","m","N0","P0"),
    params=coef(kilp1),
    est=c("b","r","m","P0")
  ) -> kilp2

fit <- subplex(
  fn=kilp2,
  par=coef(kilp2,c("b","r","m","P0"),transform=TRUE)
)
kilp2(fit$par); coef(kilp2)



## ----kilp2_plot,echo=FALSE,fig.dim=c(6,4)-------------------------------------
kilp2 |> trajectory(format="d") -> traj2
traj2 |>
  ggplot(aes(x=year,y=exp(x)))+
  geom_line()+
  geom_point(data=dat,aes(y=ti))+
  labs(y=expression(N==exp(x))) -> pl1
traj2 |>
  ggplot(aes(x=year,y=exp(y)))+
  geom_line()+
  labs(y=expression(P==exp(y))) -> pl2
plot_grid(pl1,pl2,labels="auto",align="hv",axis="tlbr",ncol=1)


## ----kilp3_fit,eval=FALSE-----------------------------------------------------
theta <- coef(kilp2)
theta["a"] <- 0.2
kilp2 |>
  traj_objfun(
    partrans=parameter_trans(
      log=c("a","b","r","m","P0")
    ),
    statenames=c("x","y"),
    paramnames=c("a","b","c","r","K","m","N0","P0"),
    est=c("a","b","r","m","P0"),
    params=theta
  ) -> kilp3

fit <- subplex(
  fn=kilp3,
  par=coef(kilp3,c("a","b","r","m","P0"),transform=TRUE)
)
kilp3(fit$par); coef(kilp3)



## ----kilp3_plot,echo=FALSE,fig.dim=c(6,6)-------------------------------------
kilp3 |> trajectory(format="d") -> traj3
traj3 |>
  ggplot(aes(x=year,y=exp(x)))+
  geom_line()+
  geom_point(data=dat,aes(y=ti))+
  labs(y=expression(N==exp(x))) -> pl1
traj3 |>
  ggplot(aes(x=year,y=exp(y)))+
  geom_line()+
  labs(y=expression(P==exp(y))) -> pl2
kilp3 |>
  trajectory(times=seq(1949,1996,by=1/365),format="d") |>
  mutate(N=exp(x),P=exp(y)) |>
  select(year,N,P) |>
  pivot_longer(c(N,P)) |>
  ggplot(aes(x=year,y=value,color=name))+
  geom_line()+
  labs(y="density") -> pl3
plot_grid(pl1,pl2,pl3,labels="auto",align="v",axis="tblr",ncol=1)


## ----kilpis_all_data,echo=FALSE,fig.dim=c(6,2)--------------------------------
read_csv("data/kilpisjarvi_c_rufocanus.csv",comment="#") -> dat1
dat1 |>
  pivot_longer(c(spring,fall)) |>
  ggplot(aes(x=year,y=value,group=name,color=name))+
  geom_line()+
  ##    scale_y_log10()+
  scale_color_manual(values=c(spring="darkgreen",fall="brown"))+
  labs(y="trap index",color="season")


## ----kilp_scape1,echo=FALSE,purl=TRUE,eval=FALSE------------------------------
## ## The following ensures that evaluating 'f' does not change
## ## the parameters stored in 'kilp3'.
## kilp3 |>
##   traj_objfun(partrans=NULL,est=c("r","a")) -> f
## 
## expand_grid(
##   r=seq(0,3,by=0.02),
##   a=seq(0.0,0.25,by=0.002)
## ) |>
##   rowwise() |>
##   mutate(
##     d=f(c(r,a))
##   ) |>
##   ggplot(aes(x=r,y=a,z=d,fill=d))+
##   geom_tile(color=NA)+
##   geom_contour(color="black",alpha=0.1)+
##   geom_point(
##     x=coef(kilp3,"r"),y=coef(kilp3,"a"),
##     shape=3,color='red',size=1
##   )



## ----kilp4_fit,eval=FALSE,purl=TRUE-------------------------------------------
## theta <- coef(kilp3)
## theta[c("r0","r1","phi")] <- c(theta["r"],0.1,0)
## theta <- theta[names(theta)!="r"]
## 
## dat |>
##   traj_objfun(
##     times="year", t0=1949,
##     rinit=Csnippet("
##     x = log(N0);
##     y = log(P0);
##   "),
##   dmeasure=Csnippet("
##     lik = -(ti-exp(x))*(ti-exp(x));
##     if (!give_log) lik = exp(lik);
##   "),
##   skeleton=vectorfield(
##     Csnippet("
##     double fr = c/(1+a*exp(x));
##     double omega = 2*M_PI;
##     Dx = r0*(1+r1*cos(omega*t+phi))*(1-exp(x)/K) - fr*exp(y);
##     Dy = b*fr*exp(x) - m;
##     ")
##   ),
##   partrans=parameter_trans(
##     log=c("a","r0","r1","b","m","P0","N0","K")
##   ),
##   paramnames=c("c","a","r0","r1","phi","K","b","m","P0","N0"),
##   statenames=c("x","y"),
##   params=theta,
##   est=c("a","r0","r1","phi","b","m","P0","N0","K")
##   ) -> kilp4
## 
## fit <- subplex(
##   fn=kilp4,
##   par=coef(kilp4,c("a","r0","r1","phi","b","m","P0","N0","K"),transform=TRUE)
## )
## kilp4(fit$par); coef(kilp4)



## ----kilp4_plot,fig.dim=c(6,6),echo=FALSE-------------------------------------
kilp4 |> trajectory(format="d") -> traj4
traj4 |>
  ggplot(aes(x=year,y=exp(x)))+
  geom_line()+
  geom_point(data=dat,aes(y=ti))+
  labs(y=expression(N==exp(x))) -> pl1
traj4 |>
  ggplot(aes(x=year,y=exp(y)))+
  geom_line()+
  labs(y=expression(P==exp(y))) -> pl2
kilp4 |>
  trajectory(times=seq(1949,1996,by=1/365),format="d") |>
  mutate(N=exp(x),P=exp(y)) |>
  select(year,N,P) |>
  gather(variable,val,N,P) |>
  ggplot(aes(x=year,y=val,color=variable))+
  geom_line()+
  labs(y="density") -> pl3
plot_grid(pl1,pl2,pl3,labels="auto",align="v",axis="tblr",ncol=1)


## ----kilp_scape2,echo=FALSE,purl=TRUE,eval=FALSE------------------------------
## kilp4 |>
##   traj_objfun(partrans=NULL,est="P0") -> f
## expand_grid(
##   P0=exp(seq(log(0.002),log(0.02),length=5000))
## ) |>
##   rowwise() |>
##   mutate(
##     d=f(P0)
##   ) |>
##   ggplot(aes(x=P0,y=d))+
##   geom_line()+
##   scale_x_log10()+
##   labs(x=expression(P[0]),y=expression(d))+
##   geom_vline(
##     xintercept=coef(kilp4,"P0"),
##     color='red'
##   )





## ----kilp4_global,eval=FALSE,purl=TRUE----------------------------------------
## estpars <- kilp4@est
## theta <- coef(kilp4,estpars,transform=TRUE)
## 
## sobol_design(
##   lower=log(0.9)+theta,
##   upper=log(1.1)+theta,
##   nseq=500
## ) -> starts
## 
## library(foreach)
## library(doParallel)
## registerDoParallel()
## 
## foreach (
##   start=iter(starts,"row"),
##   .errorhandling="pass",
##   .inorder=FALSE
## ) %dopar% {
##   f <- kilp4
##   fit <- subplex(
##     fn=f, par=start
##   )
##   f(fit$par)
##   c(coef(f),disc=f(fit$par),conv=fit$convergence)
## } -> fits




## ----kilp4_global2,fig.dim=c(6,6),out.width="90%",echo=FALSE,results="hide"----
parnames <- c(
  K="K",a="a",c="c",
  m="m",b="b",N0="N[0]",
  P0="P[0]",r0="r[0]",r1="r[1]",
  phi="varphi"
)

fits[!sapply(fits,inherits,"error")] |>
  bind_rows() |>
  filter(is.finite(disc)) |>
  arrange(disc) |>
  mutate(
    phi=atan2(sin(phi),cos(phi))
  ) -> fits

fits |> count(conv)

fits |>
  filter(conv>=0) |>
  select(-disc,-conv,-c) |>
  {\(dat)
    GGally::ggpairs(
      dat, aes(),
      columnLabels=parnames[names(dat)],
      labeller=label_parsed,
      upper=list(
        continuous="points",
        combo="facethist",
        discrete="facetbar",
        na="na"
      ),
      lower=list(
        continuous="points",
        combo="facethist",
        discrete="facetbar",
        na="na"
      )
    )+
      theme(
        axis.text.x=element_text(angle=-90),
        strip.background=element_rect(
          fill=NA,
          color=NA
        )
      )
  }()


## ----kilp4_global3,fig.dim=c(6,6),out.width="90%",echo=FALSE,results="hide"----
fits |>
  filter(conv>=0) |>
  filter(rank(disc)<=100) |>
  select(-disc,-conv,-c) |>
  {\(dat)
    GGally::ggpairs(
      dat, aes(),
      columnLabels=parnames[names(dat)],
      labeller=label_parsed,
      upper=list(
        continuous="points",
        combo="facethist",
        discrete="facetbar",
        na="na"
      ),
      lower=list(
        continuous="points",
        combo="facethist",
        discrete="facetbar",
        na="na"
      )
    )+
      theme(
        axis.text.x=element_text(angle=-90),
        strip.background=element_rect(
          fill=NA,
          color=NA
        )
      )
  }()


## ----kilp4_global_best,fig.dim=c(6,4),echo=FALSE------------------------------
fits |>
  filter(conv>=0) |>
  group_by(mode=if_else(r0>5,"high","low")) |>
  filter(disc==min(disc)) |>
  ungroup() |>
  select(-conv) |>
  pivot_longer(-mode) |>
  pivot_wider(names_from=mode) |>
  column_to_rownames("name") |>
  as.matrix() -> pars
kilp4 |>
  trajectory(params=pars,format="d") |>
  mutate(
    N=exp(x),
    P=exp(y)
  ) |>
  select(-x,-y) |>
  pivot_longer(-c(year,.id)) |>
  left_join(dat,by="year") |>
  mutate(ti=if_else(name=="N",ti,NA_real_)) |>
  ggplot(aes(x=year,y=value))+
  geom_line()+
  geom_point(aes(y=ti))+
  labs(y="")+
  facet_grid(
    name~.id,
    scales="free_y",
    labeller=label_bquote(
      rows=.(name),
      cols=.(as.character(.id))~r[0]
    )
  )+
  theme(
    strip.background=element_rect(
      fill=NA,
      color=NA
    )
  )


## ----norm_kilp,results="hide"-------------------------------------------------
fits |> filter(disc==min(disc)) |>
  mutate(rho=1,sigma=sqrt(disc/nrow(dat))) |>
  select(-disc,-conv) -> theta1
est1 <- c("a","b","r0","m","P0","N0","sigma")
kilp4 |>
  traj_objfun(
    ## Data are 'ti'.
    ## 'dnorm' is called with mean = rho*N = rho*exp(x)
    ## and s.d. = sigma = const.
    dmeasure=Csnippet("
      lik = dnorm(ti,rho*exp(x),sigma,give_log);"),
    partrans=parameter_trans(
      log=c("a","r0","b","m","P0","N0","K","rho","sigma")
    ),
    paramnames=c("c","a","r0","K","b","m","P0","N0","rho","sigma"),
    statenames=c("x","y"),
    params=theta1, est=est1
  ) -> norm_kilp

library(subplex)
fit <- subplex(
  fn=norm_kilp,
  par=coef(norm_kilp,est1,transform=TRUE)
)
norm_kilp(fit$par); coef(norm_kilp)


## ----alt_kilp,results="hide"--------------------------------------------------
fits |> filter(disc==min(disc)) |>
  mutate(rho=1,k=mean(dat$ti)/(disc/nrow(dat))) |>
  select(-disc,-conv) -> theta2
est2 <- c("a","b","r0","phi","m","P0","N0","k")
kilp4 |>
  traj_objfun(
    ## Data are 'ti'.
    ## 'dnorm' is called with mean = rho*N = rho*exp(x)
    ## and s.d. = sqrt(variance) = sqrt(rho*N/k)
    dmeasure=Csnippet("
      lik = dnorm(ti,rho*exp(x),sqrt(rho*exp(x)/k),give_log);"),
    partrans=parameter_trans(
      log=c("a","b","r0","m","P0","N0","k")
    ),
    statenames=c("x","y"),
    paramnames=c("a","b","r0","m","N0","P0","k","rho"),
    params=theta2, est=est2
  ) -> alt_kilp

library(subplex)
fit <- subplex(
  fn=alt_kilp,
  par=coef(alt_kilp,est2,transform=TRUE),
  control=list(reltol=1e-12)
)
alt_kilp(fit$par); coef(alt_kilp)

