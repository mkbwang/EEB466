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
