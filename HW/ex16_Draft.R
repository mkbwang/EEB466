# load packages

library(pomp)
library(tidyverse)
library(cowplot)
library(latex2exp)

## Wiener process code

wiener_mat <- function (times, sigma = 1, n = 1) {
  dt <- diff(c(0,times)) # take in a sequence of time points, calculate the time gaps
  dW <- rnorm(n=n*length(times),mean=0,sd=sigma*sqrt(dt)) # normal distribution vectors
  dW |>
    matrix(nrow=length(times))  -> W
  W # each column is one run of wiener process
}

wiener_df <- function (times, sigma = 1, n = 1) {
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

## exercise 1: geometric brownian process
# generate 100 paths with total time of 20
timestamps <- seq(0,15,by=0.01)
runs <- wiener_mat(times=timestamps,sigma=1,n=10000) |>
  apply(2,cumsum)

repeat_time <- matrix(timestamps, nrow=length(timestamps), ncol=10000, byrow=FALSE)

gbms <- exp(runs - 0.5 * repeat_time)
expectations <- rowMeans(gbms)
probs <- rowMeans(gbms > 0.5)
result_1 <- data.frame(time = timestamps,
                       expectation = expectations,
                       proportion = probs)

plot_expectation <- ggplot(result_1, aes(x=time, y=expectation)) + geom_line() + geom_hline(yintercept=1, linetype="dashed") +
  xlab("Time") + ylab(expression(X[t])) + ggtitle("Average Value across 10000 Runs") + theme_bw()+
  theme(axis.title = element_text(size = 12))

plot_convergence <- ggplot(result_1, aes(x=time, y=proportion)) + geom_line() +
  xlab("Time") + ylab(TeX("Pr$(X_t > 0.5)$")) + ggtitle("Average Value across 10000 Runs") + theme_bw()+
  theme(axis.title = element_text(size = 12))

plot_grid(plot_expectation, plot_convergence, align="v", labels = c('A', 'B'), label_size = 12)


## Q2 plot dX=WdW
set.seed(2021)
wiener_df(times=seq(0,10,by=0.01),sigma=1,n=20) |>
  mutate(X=0.5*(W^2 - t)) |>
  ggplot(aes(x=t,y=X,group=.id,color=factor(.id)))+
  geom_line()+
  guides(color="none")+ theme_bw()+ xlab("Time") + theme(axis.title = element_text(size = 12))+
  labs(y=expression(X[t])) + ggtitle(TeX("Twenty Different Paths for $X_t = X_0 + 0.5(W_t^2 - t)$"))


## Q3 Ou process and Em algorithm

### OU process accurate version
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
  data.frame(
    t=times,
    X=as.numeric(X),
    .id=rep(seq_len(n),each=length(times))
  )
}


### OU process EM approximation
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

set.seed(2021)
ou_analytic <- ou(times=seq(0,36,by=0.6), alpha=1, x0=100, n=100)

ou_em_001 <- ou_em(delta.t=0.01, times=seq(0,36,by=0.6), alpha=1, x0=100, n=100)
ou_em_001$difference <- ou_em_001$X - ou_analytic$X
ou_em_001$Step <- "Step size=0.01"

ou_em_005 <- ou_em(delta.t=0.05, times=seq(0,36,by=0.6), alpha=1, x0=100, n=100)
ou_em_005$difference <- ou_em_005$X - ou_analytic$X
ou_em_005$Step <- "Step size=0.05"

ou_em_01 <- ou_em(delta.t=0.1, times=seq(0,36,by=0.6), alpha=1, x0=100, n=100)
ou_em_01$difference <- ou_em_01$X - ou_analytic$X
ou_em_01$Step <- "Step size=0.1"

ou_em_015 <- ou_em(delta.t=0.15, times=seq(0,36,by=0.6), alpha=1, x0=100, n=100)
ou_em_015$difference <- ou_em_015$X - ou_analytic$X
ou_em_015$Step <- "Step size=0.15"

ou_em_02 <- ou_em(delta.t=0.2, times=seq(0,36,by=0.6), alpha=1, x0=100, n=100)
ou_em_02$difference <- ou_em_02$X - ou_analytic$X
ou_em_02$Step <- "Step size=0.2"


ou_em_03 <- ou_em(delta.t=0.3, times=seq(0,36,by=0.6), alpha=1, x0=100, n=100)
ou_em_03$difference <- ou_em_03$X - ou_analytic$X
ou_em_03$Step <- "Step size=0.3"

ou_em_06 <- ou_em(delta.t=0.6, times=seq(0,36,by=0.6), alpha=1, x0=100, n=100)
ou_em_06$difference <- ou_em_06$X - ou_analytic$X
ou_em_06$Step <- "Step size=0.6"

ou_ems <- rbind(ou_em_001, ou_em_005, 
                ou_em_01, ou_em_02, ou_em_03, ou_em_06)


ou_analytic_plot <- ggplot(ou_analytic, aes(x=t,y=X,group=.id,color=factor(.id)))+
  geom_line()+
  guides(color="none")+
  theme_bw()+ xlab("Time") + theme(axis.title = element_text(size = 12), plot.title=element_text(size = 12))+
  xlab("Time") + ylab(TeX("$X_t$")) + ggtitle("Analytic Solution of Ornstein-Uhlenbeck process")

ou_ems_diff_plot <- ggplot(ou_ems, aes(x=t,y=difference,group=.id,color=factor(.id))) + facet_wrap(vars(Step), nrow = 3)+
  geom_line()+
  guides(color="none") + theme_bw()+ xlab("Time") + theme(axis.title = element_text(size = 12), plot.title=element_text(size = 12))+
  labs(y="Trajectory Difference") + ggtitle("Difference between E-M Simulation and Analytic Solution of O-U Process")

plot_grid(ou_analytic_plot, ou_ems_diff_plot, align="v", labels = c('A', 'B'), label_size = 12)


## use pomp to run Rosensweig-MacArthur Model

### generate the model
simulate(
  t0 = 0, times=seq(0,40,by=0.02),
  rinit=Csnippet("x = log(n0); y=log(p0);"),
  rprocess=euler(
    Csnippet("x += (1 - exp(x)/gamma - exp(y)/(1+exp(x)) - sigma1 * sigma1 / 2)*dt + sigma1 * rnorm(0,sqrt(dt));
             y += (exp(x) / (1+exp(x)) - alpha - sigma2 * sigma2 / 2)*dt + sigma2 * rnorm(0, sqrt(dt));"),
    delta.t=0.01
  ),
  paramnames=c("n0", "p0", "alpha", "gamma", "sigma1", "sigma2"), statenames=c("x", "y"),
  params=c(n0=4, p0=3, alpha=3, gamma=2.5, sigma1=0.02, sigma2=0.02)
) -> RM


### set parameters
starting_points_trajectory <- expand_grid(x=seq(1, 5), y=seq(1, 5)) |> as.matrix() |> t()

coef_mat = parmat(coef(RM), nrep=25)
coef_mat[1:2,] <- starting_points_trajectory

coef_mat['alpha', ] <- 0.75
regime1 <- RM |> simulate(params=coef_mat) |> as.data.frame()
regime1$Regime <- "Regime 1"

coef_mat['alpha', ] <- 0.5
regime2 <- RM |> simulate(params=coef_mat) |> as.data.frame()
regime2$Regime <- "Regime 2"

coef_mat['alpha', ] <- 0.25
regime3 <- RM |> simulate(params=coef_mat) |> as.data.frame()
regime3$Regime <- "Regime 3"

allpaths <- rbind(regime1, regime2, regime3)

ggplot(allpaths, aes(x=exp(x),y=exp(y), group=.id)) + facet_wrap(vars(Regime), nrow = 1)+
  geom_path() + guides(color="none")+
  theme_bw() + labs(y="Predator Population",x="Prey Population") +
  theme(axis.title = element_text(size = 12), plot.title=element_text(size = 12))+
  ggtitle("Sample Paths in Three Stability Regimes")

