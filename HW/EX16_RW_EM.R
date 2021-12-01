# load packages

library(pomp)
library(tidyverse)
library(cowplot)
library(latex2exp)

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
  params=c(n0=4, p0=3, alpha=3, gamma=2.5, sigma1=0.1, sigma2=0.1)
) -> RM



### sample paths
timepoints <- seq(0,40,by=0.02)
starting_points_trajectory <- expand_grid(x=seq(1, 5), y=seq(1, 5)) |> as.matrix() |> t()

coef_mat = parmat(coef(RM), nrep=25)
coef_mat[1:2,] <- starting_points_trajectory

coef_mat['alpha', ] <- 0.75
regime1 <- RM |> simulate(params=coef_mat, times=timepoints) |> as.data.frame()
regime1$Regime <- "Regime 1"

coef_mat['alpha', ] <- 0.5
regime2 <- RM |> simulate(params=coef_mat, times=timepoints) |> as.data.frame()
regime2$Regime <- "Regime 2"

coef_mat['alpha', ] <- 0.25
regime3 <- RM |> simulate(params=coef_mat, times=timepoints) |> as.data.frame()
regime3$Regime <- "Regime 3"

allpaths <- rbind(regime1, regime2, regime3)

ggplot(allpaths, aes(x=exp(x),y=exp(y), group=.id)) + facet_wrap(vars(Regime), nrow = 1)+
  geom_path() + guides(color="none")+
  theme_bw() + labs(y="Predator Population",x="Prey Population") +
  theme(axis.title = element_text(size = 12), plot.title=element_text(size = 12), strip.text = element_text(size=10))+
  ggtitle("Sample Paths in Three Stability Regimes")

## set density plots

timepoints <- c(0, 2, 4, 8, 16)
starting_points_density <- expand_grid(x=seq(0.05, 5, 0.05), y=seq(0.05, 5, 0.05)) |> as.matrix() |> t()

coef_mat = parmat(coef(RM), nrep=10000)
coef_mat[1:2, ] <- starting_points_density

coef_mat['alpha', ] <- 0.75
regime1 <- RM |> simulate(params=coef_mat, times=timepoints) |> as.data.frame()
regime1$Regime <- "Regime 1"

coef_mat['alpha', ] <- 0.5
regime2 <- RM |> simulate(params=coef_mat, times=timepoints) |> as.data.frame()
regime2$Regime <- "Regime 2"

coef_mat['alpha', ] <- 0.25
regime3 <- RM |> simulate(params=coef_mat, times=timepoints) |> as.data.frame()
regime3$Regime <- "Regime 3"

allpoints <- rbind(regime1, regime2, regime3)
allpoints$time <- sprintf('t=%d', allpoints$time)
allpoints$time <- factor(allpoints$time, levels=c("t=0", "t=2", "t=4", "t=8", "t=16"))

ggplot(allpoints, aes(x=exp(x), y=exp(y))) + facet_grid(rows=vars(time), cols=vars(Regime))+
  geom_bin2d(bins = 100) + scale_fill_continuous(type="viridis", limits=c(0, 50), oob = scales::squish)+
  labs(y="Predator Population",x="Prey Population") + scale_y_continuous(limits=c(0, 5.5), breaks=seq(0, 5))+
  theme(axis.title = element_text(size = 12), plot.title=element_text(size = 12), strip.text = element_text(size=10))+
  ggtitle("Density Dynamics under Three Stability Regimes")

  
