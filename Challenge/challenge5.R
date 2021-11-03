# challenge 5

library(tidyverse)
library(pomp)
library(subplex)


# read_csv("data/kilpisjarvi_fall_c_rufocanus.csv",comment="#") -> dat
read_csv("data/kilpisjarvi_c_rufocanus.csv",comment='#') -> dat_full
dat_full <- dat_full |> na.omit()
spring_data <- dat_full |> select(year, spring) 
fall_data <- dat_full |> select(year, fall) 
start_spring <- spring_data$year[1]
start_fall <- fall_data$year[1]
initpop_spring <- spring_data$spring[1]
initpop_fall <- fall_data$fall[1]

spring_data |>
  traj_objfun(
    times="year", t0=start_spring,
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
      lik = -(spring-exp(x))*(spring-exp(x));
      if (!give_log) lik = exp(lik);
    "),
    partrans=parameter_trans(
      log=c("a", "b","r","m","P0")
    ),
    statenames=c("x","y"),
    paramnames=c("a","b","c","r","K","m","N0","P0"),
    params=c(K=25,a=0.2,r=2,c=1,m=1,b=1,N0=initpop_spring,P0=1),
    est=c("a","b","r","m","P0")
  ) -> RM_spring

spring_data |>
  traj_objfun(
    times="year", t0=start_spring,
    rinit=Csnippet("
      x = log(N0);
      y = log(P0);
      "),
    skeleton=vectorfield(
      Csnippet("
      double fr = c/(1+a*exp(x));
      Dx = r*(1-exp(x)/K) - fr*exp(y);
      Dy = s*(1-q*exp(y)/exp(x));
      ")
    ),
    dmeasure=Csnippet("
      lik = -(spring-exp(x))*(spring-exp(x));
      if (!give_log) lik = exp(lik);
    "),
    partrans=parameter_trans(
      log=c("a", "r", "s", "q", "P0")
    ),
    statenames=c("x","y"),
    paramnames=c("a", "c","r", "K", "s", "q", "N0","P0"),
    params=c(K=25, a=0.2, r=2, c=1, s=1, q=1.5, N0=initpop_spring,P0=1),
    est=c("a","s","r","q","P0")
  ) -> HS_spring


fall_data |>
  traj_objfun(
    times="year", t0=start_fall,
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
      lik = -(fall-exp(x))*(fall-exp(x));
      if (!give_log) lik = exp(lik);
    "),
    partrans=parameter_trans(
      log=c("a", "b","r","m","P0")
    ),
    statenames=c("x","y"),
    paramnames=c("a","b","c","r","K","m","N0","P0"),
    params=c(K=25,a=0.2,r=2,c=1,m=1,b=1,N0=initpop_fall,P0=1),
    est=c("a","b","r","m","P0")
  ) -> RM_fall

fall_data |>
  traj_objfun(
    times="year", t0=start_fall,
    rinit=Csnippet("
      x = log(N0);
      y = log(P0);
      "),
    skeleton=vectorfield(
      Csnippet("
      double fr = c/(1+a*exp(x));
      Dx = r*(1-exp(x)/K) - fr*exp(y);
      Dy = s*(1-q*exp(y)/exp(x));
      ")
    ),
    dmeasure=Csnippet("
      lik = -(fall-exp(x))*(fall-exp(x));
      if (!give_log) lik = exp(lik);
    "),
    partrans=parameter_trans(
      log=c("a", "s","r","q","P0")
    ),
    statenames=c("x","y"),
    paramnames=c("a","s","c","r","K","q","N0","P0"),
    params=c(K=25,a=0.2,r=2,c=1,s=1,q=1.5,N0=initpop_fall,P0=1),
    est=c("a","s","r","q","P0")
  ) -> HS_fall



library(subplex)
fit_spring_RM <- subplex(
  fn=RM_spring,
  par=coef(RM_spring,c("a", "b","r","m","P0"),transform=TRUE)
) # changing the parameter of pomp object in place

RM_spring(fit_spring_RM$par)

fit_spring_HS <- subplex(
  fn=HS_spring,
  par = coef(HS_spring, c("a","s","r","q","P0"), transform=TRUE)
)

HS_spring(fit_spring_HS$par)

fit_fall_RM <- subplex(
  fn=RM_fall,
  par=coef(RM_fall,c("a", "b","r","m","P0"),transform=TRUE)
) # changing the parameter of pomp object in place

RM_fall(fit_fall_RM$par)

fit_fall_HS <- subplex(
  fn=HS_fall,
  par=coef(HS_fall, c("a","s","r","q","P0"),transform=TRUE)
)

HS_fall(fit_fall_HS$par)

coef(RM_spring)
coef(RM_fall)
coef(HS_spring)
coef(HS_fall)

RM_spring |>
  trajectory(times=seq(1953,1996.5,by=0.5),format="d") |>
  mutate(N=exp(x),P=exp(y)) -> spring_prediction

RM_fall |>
  trajectory(times=seq(1953,1996,by=0.5),format="d") |>
  mutate(N=exp(x),P=exp(y)) -> fall_prediction

filter <- rep(c(FALSE, TRUE), times=1996-1953+1)
spring_fall_pred <- spring_prediction[filter, c('year', 'N', 'P')]
spring_spring_pred <- spring_prediction[!filter, c('year', 'N', 'P')]
filter <- filter[-1]
fall_fall_pred <- fall_prediction[filter, c('year', 'N', 'P')]
fall_spring_pred <- fall_prediction[!filter, c('year', 'N', 'P')]






