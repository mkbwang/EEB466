# challenge 4 code

library(tidyverse)
library(pomp)
library(latex2exp)
library(patchwork)

# linear flow in 2D
## initialize the model
trajectory(
  t0=0, times=seq(0,0.6,by=0.001),
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
    a11  = 1, a12  =  1,
    a21  = 1, a22  = 1,
    x1_0 = 1, x2_0 =  1
  )
) -> lf

## define a function that generates a plot given a coefficient matrix
l2traj <- function(trajobj, a11, a12, a21, a22){
  
  ### redeclare parameters
  params <- parmat(coef(trajobj), 25)
  
  params['a11', ] <- a11
  params['a12', ] <- a12
  params['a21', ] <- a21
  params['a22', ] <- a22
  params['x1_0', ] = rep(c(-2, -1, 0, 1, 2), 5)
  params['x2_0', ] = rep(c(-2, -1, 0, 1, 2), each=5)
  
  ### get all the plots
  trajobj |> trajectory(params=params) |> as.data.frame() -> dat
  
  ### get eigenvalues
  eigen_vals <- eigen(matrix(c(a11, a12, a21, a22), nrow=2))$values
  title_label <- NULL
  
  if(Im(eigen_vals[1])==0){### plot title
    template <- r'($\lambda_{1}=%.1f$ $\lambda_{2}=%.1f$)'
    title_label <- sprintf(template, eigen_vals[1], eigen_vals[2])
  } else{
    template <- r'( $\lambda_{1}=%.1f+%.1fi$ $\lambda_{2}=%.1f+%.1fi$)'
    title_label <- sprintf(template, Re(eigen_vals[1]), Im(eigen_vals[1]),
                           Re(eigen_vals[2]), Im(eigen_vals[2]))
  }
  
  start = data.frame(x1 = rep(c(-2, -1,  0,1, 2), 5),
                     x2 = rep(c(-2, -1, 0,1, 2), each=5))
  
  viz <- ggplot(NULL, aes(x=x1, y=x2)) +geom_path(data=dat, aes(group=.id), show.legend = F)+
    geom_point(data=start, size=1.2, show.legend = F) +xlim(-5, 5) + ylim(-5, 5)+ggtitle(TeX(title_label))+
    theme_bw() + theme(plot.title=element_text(hjust=0.5, size=8),
                       axis.text.x=element_text(size=7),
                       axis.text.y=element_text(size=7),
                       axis.title.x=element_blank(),
                       axis.title.y=element_blank())
  
  viz
}

plots <- vector(mode='list', length=11)
plots[[1]] <- l2traj(lf, 4,3,-2,-1)
plots[[2]] <- l2traj(lf, 1,1,-1,3)
plots[[3]] <- l2traj(lf, -3, 1, -1, -1)
plots[[4]] <- l2traj(lf, -4, 3, -2, 1)
plots[[5]] <- l2traj(lf, 1, 1, 1, -1)
plots[[6]] <- l2traj(lf, 2, 2, -1, -1)
plots[[7]] <- l2traj(lf, -2, 2, -1, 1)
plots[[8]] <- l2traj(lf, 1, -1, 1, -1)
plots[[9]] <- l2traj(lf, 1, 1, -1, 1)
plots[[10]] <- l2traj(lf, -1, 1, -1, -1)
plots[[11]] <- l2traj(lf, 0, 1, -1, 0)

aggregates <- Reduce(`+`, plots) + plot_layout(ncol=4)

aggregates
