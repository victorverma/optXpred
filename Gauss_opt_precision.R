#
#
#  These scripts compute the optimal precision for a Gaussian model
#  as a function of "rho" using numerical integration.
#
get_2dnormal = function(n,rho){
  X = rnorm(n)
  Y = rho*X + sqrt(1-rho^2)*rnorm(n);
  return(cbind(X,Y))
}

f_integrate <- function(u,alpha,rho){
  z_alpha = qnorm(1-alpha);
  return((1-pnorm((z_alpha-rho*qnorm(u))/sqrt(1-rho^2)))/alpha)
}

get_nomal_precision <- function(alpha,rho){
  prec = c()
  for (a in alpha){
    prec = c(prec, integrate(f_integrate,lower=1-a,upper=1,alpha = a, rho = rho)$value)
  }
  return(prec)
}

example_normal_precision <- function(n,rho_s = seq(from=0.8, to=0.999, length.out = 10)){
  for (rho in rho_s){
    Z =get_2dnormal(n,rho)
    p =  seq(from=0.5,to=0.999,length.out = 100)
    out = get_lambda(Z[,1],Z[,2],p=p,main=paste("rho=",rho));
    z_alpha = qnorm(p);
    z_alpha_2 = qnorm(1-(1-p)/2);
    #Implementing Simson's rule after a change of variables
    fa = 1-pnorm(z_alpha*(1 - rho)/sqrt(1-rho^2));
    fab = 1-pnorm((z_alpha - rho* z_alpha_2)/sqrt(1-rho^2));
    fb = 1           
    lines(p,(fa+4*fab+fb)/6,col="red");
    lines(p,get_nomal_precision(alpha = 1-p,rho=rho),col="blue")
    readline();
  }
}

figure_gauss_precision <-function(p=seq(from =0.5, to =0.999, length.out = 200),
                                  rho = seq(from =0, to=1, length.out = 200),monochrome=F){

  library(ggplot2)
  library(metR)
  param_grid <- expand.grid(alpha = 1-p, rho = rho)
  param_grid$precision <- mapply(get_nomal_precision, param_grid$alpha, param_grid$rho)
  param_grid$p = p

  contour_levels <- pretty(range(param_grid$precision), n = 6)
  if (monochrome){
    high_color = "white";
    low_color= "black"
  } else {
    low_color = "red";
    high_color= "blue"
  }
  plt<-ggplot(param_grid, aes(x = rho, y = p, fill = precision)) +
    geom_tile() +
    scale_fill_gradient(
      low = low_color, high = high_color, name = "Precision") +
    geom_contour(
      aes(z = precision),
      color = "white",
      linetype = "dashed",
      breaks = contour_levels
    ) +
    metR::geom_text_contour(
      aes(z = precision),
      breaks = contour_levels,
      stroke = 0.2,
      size = 3,
      check_overlap = FALSE,  # Label all contours
      skip = 0                # Label every contour line (no skipping)
    ) +
    labs(
      title = "Optimal precisions for Gaussian copula",
      x = expression(rho),
      y = expression(p=1-alpha)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  return(plt)
}

plt = figure_gauss_precision(monochrome = F);
ggsave("Gauss_opt_prec.pdf",plot=plt)

plt = figure_gauss_precision(monochrome = T);
ggsave("Gauss_opt_prec_MONOCHROME.pdf",plot=plt)