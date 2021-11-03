## ----packages,cache=FALSE,message=FALSE,include=FALSE-------------------------
library(tidyverse)
library(pomp)
stopifnot(getRversion()>="4.1")
stopifnot(packageVersion("pomp")>="4.0")
library(cowplot)
theme_set(theme_bw(base_size=11,base_family="Helvetica"))
set.seed(228088786)


## ----vandy_data,fig.dim=c(6,2),echo=FALSE-------------------------------------
library(tidyverse)

read_csv("data/vandy_data.csv") |>
  filter(experiment=="Pa") |>
  select(day,rep,count) -> dat

dat |>
  ggplot(aes(x=day,y=count,color=factor(rep)))+
  geom_line()+
  labs(color="replicate")


## ----load_dat,fig.dim=c(6,2)--------------------------------------------------
library(tidyverse)
read_csv("data/vandy_data.csv") |>
  filter(experiment=="Pa") |>
  select(day,rep,count) -> dat
head(dat)


## ----build_pomp---------------------------------------------------------------
library(pomp)
dat |>
  filter(rep==1) |>
  select(-rep) |>
  pomp(
    times="day", t0=0,
    rinit=function (N0, ...) {
      c(N=N0)
    },
    skeleton=vectorfield (
      function (N, r, K, ...) {
        c(N=r*N*(1-N/K))
      }
    )
  ) -> vandy


## ----traj1,fig.dim=c(6,2)-----------------------------------------------------
vandy |>
  trajectory(
    params=c(r=1,K=700,N0=10),
    format="data.frame"
  ) |>
  ggplot(aes(x=day,y=N))+
  geom_line()+
  geom_point()+
  expand_limits(y=800)


## ----dmeasure1----------------------------------------------------------------
vandy |>
  pomp(
    dmeasure=function (count, N, ..., log) {
      d <- (count-N)^2
      if (log) -d else exp(-d)
    }
  ) -> vandy


## ----traj_ofun----------------------------------------------------------------
vandy |>
  traj_objfun(
    params=c(r=1,K=700,N0=10),
    est=c("r")
  ) -> vto


## ----rslice,echo=FALSE--------------------------------------------------------
rvals <- c(0.3,0.65,1.4)
expand_grid(
  r=seq(0.2,2,by=0.1)
) |>
  group_by(r) |>
  mutate(
    d=vto(r)
  ) |>
  ggplot(aes(x=r,y=d/1e6))+
  geom_line()+
  geom_point()+
  geom_vline(xintercept=rvals,color='red')+
  labs(
    x=expression(r),
    y=expression(d(r)/10^6)
  ) -> pl1

## ----rslice2,echo=FALSE-------------------------------------------------------
vandy |>
  trajectory(
    params={
      p <- parmat(coef(vto),3)
      p["r",] <- rvals
      p
    },
    format="d"
  ) |>
  mutate(
    r=rvals[.id]
  ) |>
  left_join(
    dat |> filter(rep==1),
    by="day"
  ) |>
  ggplot(aes(x=day,y=N))+
  geom_line()+
  geom_point(aes(y=count))+
  facet_grid(r~.,labeller=label_bquote(r==.(r))) -> pl2

## ----rsliceplot,cache.dep=c("rslice1","rslice2"),fig.dim=c(6,2.5),echo=FALSE----
library(cowplot)
plot_grid(pl1,pl2,labels="auto",nrow=1,axis="tblr",rel_heights=c(1,2))


## ----traj_match1--------------------------------------------------------------
vto |>
  traj_objfun(
    est=c("r","K","N0"),
    partrans=parameter_trans(log=c("r","K","N0")),
    paramnames=c("r","K","N0")
  ) -> vto

fit <- optim(fn=vto,par=coef(vto,transform=TRUE))


## ----traj_match1a-------------------------------------------------------------
vto(fit$par); coef(vto)


## ----traj_match1b,fig.dim=c(6,2.5)--------------------------------------------
vto |> trajectory(format="d") |>
  left_join(dat |> filter(rep==1), by="day") |>
  ggplot(aes(x=day,y=N))+
  geom_line()+geom_point(aes(y=count))


## ----landscapes,echo=FALSE,fig.dim=c(6,3),dev="CairoPNG"----------------------
thetahat <- coef(vto)
expand_grid(
  r=seq(0.5,1.8,by=0.02),
  K=seq(500,900,by=10),
  N_0=thetahat["N0"]
) |>
  rowwise() |>
  mutate(
    d=vto(log(c(r,K,N_0)))
  ) |>
  ggplot(aes(x=r,y=K,z=d,fill=d))+
  geom_tile(color=NA)+
  geom_contour(color="black")+
  geom_point(x=thetahat["r"],y=thetahat["K"],shape=3,color='red',size=2)+
  expand_limits(fill=2.5e6)+
  theme(legend.position="none")+
  labs(
    x=expression(italic(r)),
    y=expression(italic(K))
  ) -> pl1

expand_grid(
  r=seq(0.7,1.8,by=0.02),
  K=thetahat["K"],
  N_0=seq(0.1,2,by=0.04)
) |>
  rowwise() |>
  mutate(
    d=vto(log(c(r,K,N_0)))
  ) |>
  ggplot(aes(x=r,y=N_0,z=d,fill=d))+
  geom_tile(color=NA)+
  geom_contour(color="black")+
  geom_point(x=thetahat["r"],y=thetahat["N0"],shape=3,color='red',size=2)+
  expand_limits(fill=2.5e6)+
  theme(legend.position="none")+
  labs(
    x=expression(italic(r)),
    y=expression(italic(N[0]))
  ) -> pl2

plot_grid(
  plot_grid(pl1,pl2,labels="auto",nrow=1,align="hv",axis="tblr"),
  get_legend(
    pl1+theme(
          legend.position="right",
          legend.box.margin=margin(l=12)
        )
  ),
  rel_widths=c(1,0.2)
)

