library(pomp)
library(tidyverse)
library(cowplot)

## euler multinomial model
simulate(
  t0=0, times=seq(0, 1000),
  params=c(birth=0.001, death=0.001, Beta=4, mu_EI=7, mu_IR=1, N0=4000000, I0=2000),
  rprocess = euler(Csnippet("
    double rate[8], dN[8];
    rate[0] = birth; // birth rate
    rate[1] = rate[3] = rate[5] = rate[7] = death; // death rate
    rate[2] = Beta * I / N; // mu_SE
    rate[4] = mu_EI;
    rate[6] = mu_IR;
    // Rprintf(\"%lg %lg %lg %lg %lg\\n\",S,E,I,R,N); 
    dN[0] = rpois(rate[0]*N*dt);
    reulermultinom(2,S,&rate[1],dt,&dN[1]);
    reulermultinom(2,E,&rate[3],dt,&dN[3]);
    reulermultinom(2,I,&rate[5],dt,&dN[5]);
    reulermultinom(1,R,&rate[7],dt,&dN[7]);
    S += dN[0]-dN[1]-dN[2];
    E += dN[2]-dN[3]-dN[4];
    I += dN[4]-dN[5]-dN[6];
    R += dN[6]-dN[7];
    N += dN[0]-dN[1]-dN[3]-dN[5]-dN[7];"), delta.t=0.005),
  rinit = Csnippet("S=nearbyint((N0-I0)*0.3); E=0; I=I0; R=nearbyint((N0-I0)*0.7); N=N0;"),
  paramnames=c("birth", "death", "Beta", "mu_EI", "mu_IR", "N0", "I0"),
  statenames=c("S", "E", "I", "R", "N"),
  nsim=20
) -> em

em |>
  as.data.frame() ->df

# lineplot
df |>
  ggplot(aes(x=time,y=I,group=.id,color=factor(.id)))+
  geom_line()+
  guides(color="none")+
  expand_limits(y=0)+
  scale_y_sqrt()+theme_bw()+
  theme(axis.title = element_text(size = 12), plot.title=element_text(size = 12), strip.text = element_text(size=10))+
  labs(y="Size",x="Time(week)") + ggtitle("Infected Population")


em[[1]] |> simulate(times=c(100, 250, 400, 550, 700, 850),
               nsim=1000) -> em_hist

df <- em_hist |> as.data.frame()

df$time <- sprintf('t=%d', df$time)
df$time <- factor(df$time, levels=c("t=100", "t=250", "t=400", "t=550", "t=700", "t=850"))
df_nozeros <- df |> filter(I != 0 & I < 10000) |> filter(time != "t=100")

# histogram
ggplot(df_nozeros, aes(x=I)) + geom_histogram(color="black", fill="white", binwidth = 10)+facet_grid(rows=vars(time))+
  labs(x="Infected Population Size", y="Count") +theme(axis.title = element_text(size = 12), plot.title=element_text(size = 12), strip.text = element_text(size=10))
  
