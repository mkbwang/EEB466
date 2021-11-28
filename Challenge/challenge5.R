# challenge 5

library(tidyverse)
library(pomp)
library(subplex)
library(cowplot)


# read_csv("data/kilpisjarvi_fall_c_rufocanus.csv",comment="#") -> dat
read_csv("data/kilpisjarvi_c_rufocanus.csv",comment='#') -> dat_full
spring_data <- dat_full |> select(year, spring) |> na.omit()
fall_data <- dat_full |> select(year, fall) |> na.omit()
start_spring <- spring_data$year[1]
start_fall <- fall_data$year[1]
initpop_spring <- spring_data$spring[1]
initpop_fall <- fall_data$fall[1]

# load the RM model for fall data fitted in class
prefitted_fall_RM <- readRDS('RM_fall.rds')
prefitted_params <- readRDS('RM_fall_fit_params.rds')
ref_coefs <- coef(prefitted_fall_RM)

prefitted_fall_RM(prefitted_params$par)

# set up RM model for spring data
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
    params=c(K=unname(ref_coefs['K']),a=unname(ref_coefs['a']),r=unname(ref_coefs['r']),c=1,m=unname(ref_coefs['m']),b=unname(ref_coefs['b']),
               N0=initpop_spring,P0=unname(ref_coefs['P0'])/unname(ref_coefs['N0'])*initpop_spring),
    est=c("a","b","r","m","P0")
  ) -> RM_spring

spring_coefs <- coef(RM_spring)

RM_spring |> traj_objfun(
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
  params=c(K=unname(spring_coefs['K']), a=unname(spring_coefs['a']), r=unname(spring_coefs['r']), c=1, s=1, q=20,
           N0=initpop_spring, P0=unname(spring_coefs['P0'])),
  est=c("a","s","r","q","P0")
) -> newmodel

# set up HS model for spring data
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
    params=c(K=unname(spring_coefs['K']), a=unname(spring_coefs['a']), r=unname(spring_coefs['r']), c=1, s=1, q=20,
             N0=initpop_spring, P0=unname(spring_coefs['P0'])),
    est=c("a","s","r","q","P0")
  ) -> HS_spring


# fall_data |>
#   traj_objfun(
#     times="year", t0=start_fall,
#     rinit=Csnippet("
#       x = log(N0);
#       y = log(P0);
#       "),
#     skeleton=vectorfield(
#       Csnippet("
#       double fr = c/(1+a*exp(x));
#       Dx = r*(1-exp(x)/K) - fr*exp(y);
#       Dy = b*fr*exp(x) - m;
#       ")
#     ),
#     dmeasure=Csnippet("
#       lik = -(fall-exp(x))*(fall-exp(x));
#       if (!give_log) lik = exp(lik);
#     "),
#     partrans=parameter_trans(
#       log=c("a", "b","r","m","P0")
#     ),
#     statenames=c("x","y"),
#     paramnames=c("a","b","c","r","K","m","N0","P0"),
#     params=c(K=25,a=0.2,r=2,c=1,m=1,b=1,N0=initpop_fall,P0=1),
#     est=c("a","b","r","m","P0")
#   ) -> RM_fall

# set up HS model for fall data
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
    params=c(K=unname(ref_coefs['K']),a=unname(ref_coefs['a']),r=unname(ref_coefs['r']),c=1,s=1,q=20,
             N0=initpop_fall,P0=unname(ref_coefs['P0'])),
    est=c("a","s","r","q","P0")
  ) -> HS_fall


# fitting models
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

# fit_fall_RM <- subplex(
#   fn=RM_fall,
#   par=coef(RM_fall,c("a", "b","r","m","P0"),transform=TRUE)
# ) # changing the parameter of pomp object in place
# 
# RM_fall(fit_fall_RM$par)

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
  trajectory(times=spring_data$year,format="d") |>
  mutate(N=exp(x),P=exp(y)) -> spring_prediction_RM

prefitted_fall_RM |>
  trajectory(times=fall_data$year,format="d") |>
  mutate(N=exp(x),P=exp(y)) -> fall_prediction_RM

HS_fall |>
  trajectory(times=fall_data$year,format="d") |>
  mutate(N=exp(x),P=exp(y)) -> fall_prediction_HS


HS_spring |>
  trajectory(times=spring_data$year,format="d") |>
  mutate(N=exp(x),P=exp(y)) -> spring_prediction_HS

fall_predictions <- data.frame(year = rep(fall_data$year, times=2),
                               model = rep(c('RM', 'HS'), each=length(fall_data$year)),
                               pred = c(fall_prediction_RM$N, fall_prediction_HS$N))

spring_predictions <- data.frame(year = rep(spring_data$year, times=2),
                                 model = rep(c('RM', 'HS'), each=length(spring_data$year)),
                                 pred = c(spring_prediction_RM$N, spring_prediction_HS$N))

fall_predictions |>
  ggplot(aes(x=year,y=pred))+
  geom_path(aes(linetype=model))+
  scale_linetype_manual(values=c("twodash", "dotted"))+
  geom_point(data=fall_data, aes(y=fall)) + ylab('Fall Prey Number') -> pl1

spring_predictions |>
  ggplot(aes(x=year,y=pred))+
  geom_path(aes(linetype=model))+
  scale_linetype_manual(values=c("twodash", "dotted"))+
  geom_point(data=spring_data, aes(y=spring)) + ylab('Spring Prey Number') -> pl2

plot_grid(pl1,pl2,labels="auto",align="hv",axis="tlbr",ncol=1)


