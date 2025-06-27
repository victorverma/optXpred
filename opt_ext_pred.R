#
# This is a proof-of-concept implementation of the optimal extreme event 
# prediction methodology via homogeneous predictors.
#
# It uses the quantile random forest based estimation of the conditional quantile
# functions implemented in the "ranger" package.
#

#############################################################################
#
#  The following functions implement the methodology for: 
#
#   O P T I M A L  H O M O G E N E O U S  P R E D I C T I O N
#
#     Via Quantile Random Forest 
#
#
#############################################################################
#

library(ranger)

xy_to_utheta <- function(y,x,prob=0.9){
  #
  # input:
  #   y <- n x 1 vector
  #   x <- n x p array of covariates
  #   prop <- threshold proportion (default = 0.9)
  # 
  # output:
  #  u <- m x 1
  #  theta <- n x p 
  # 
  d = dim(x)[2];
  u = abs(y);
  r = u + rowSums(abs(x));
  idx = which(r>= quantile(r,prob))
  u = u[idx]/r[idx];
  keep = which(u<1) #Dropping values of u=1;
  u = u[keep]
  theta = x[idx[keep],]/(matrix((1-u)*r[idx[keep]],length(u),1)%*%matrix(1,1,d));
  return(list("u"=u,"theta"=theta))
}

get_c <- function(u,theta){
  #
  #  Returns the calibration level "c = E[U]"
  #
 return(mean(u))  
}

get_qreg_forest <- function(u, theta, n_trees=2000) {
#
# Builds a random forest for the purpose of computing (weighted) quantile
# regression estimates of the conditional density ~ (1-u)f_{U|\Theta}(u|theta)
# 
# input: 
#    u<- a vector (nx1)
#     theta <- an array (nxd)
#     n_trees <- the number of trees paramter to pass to the function "ranger"
#
# output:
#   rf <- an S3 object containing the data structure storing the random forest 
#         and the methods for computing conditional quantiles
#
#
 theta <- as.data.frame(theta)
#
# IMPORTANT: Performs a *weighted* quantile regression.
#
 weights <- 1 - u
 data <- data.frame(u = u, theta)
 rf <- ranger(u ~ ., data = data, case.weights = weights, quantreg = TRUE, num.trees =n_trees)
return(rf)
}

get_q_alpha <- function(rf, theta_new, alpha,verbose=FALSE) {
#
# input:
#   rf <- an S3 object returned by "ranger"
#   theta_new <- (m x d) array.  IMPORTANT: The column names (if any) must match the 
#             matching the ones used to build the trees in "rf". The columns are treated
#             as covariates.
#
#  alpha <- a number in (0,1)
#
# output:
#  q_alpha <- (m x 1) vector with the estimated conditional alpha-level quantiles 
#              evaluated at "theta_new"
#
  if (verbose){cat("Predicting using random forest quantile regression...")}
  theta_new <- as.data.frame(theta_new)
  pred <- predict(rf, data = theta_new, type = "quantiles", quantiles = alpha)
  if (verbose){cat("Done.")}
  return(pred$predictions)
}

calibrate_alpha <-function(rf,c,u,theta,d,N_theta = d*1e4,a0 = 1e-5,a1 = 1-a0,iter=12,verbose=FALSE){
  if (verbose){
    cat('\n Calibrating using a binary search: ')
    cat('\n Constraint c =',c);
  }
  g0 = get_q_alpha(rf,theta,a0);
  g0 = g0/(1-g0);
  
  g1 = get_q_alpha(rf,theta,a1);
  g1 = g1/(1-g1);
  
  c0 = mean((1-u)*g0);
  c1 = mean((1-u)*g1);
  
  if (verbose){
    cat('\n (c0,c1)=(',c0,",",c1,"); (a0,a1)=(",a0,",",a1,")");
  }
  if ((c<c0)||(c>c1)){
    simpleError("The calibration level is not sandwiched. Try smaller/larger values for alpha0/alpha1.")
    break; 
  }
  
  for (i in c(1:iter)){
   a_new = (a0+a1)/2;
   q_new = get_q_alpha(rf,theta,a_new);
   g_new = q_new/(1-q_new);
   c_new = mean((1-u)*g_new);
   if (c<=c_new){
      a1 = a_new
      c1 = c_new
   } else {
      a0 = a_new;
      c0 = c_new;
   }
   if (verbose){
     cat('\n (c0,c1)=(',c0,",",c1,"); (a0,a1)=(",a0,",",a1,")")
   }
  }
  if(verbose){cat("\nDone calibrating. alpha=",a_new,"\n")}
  return(a_new)
}

get_h_opt <- function(y,x,x_test,prob=0.8,verbose=TRUE, n_trees=2000){
 d=dim(x)[2]
 #
 # Step 1: Conditioning on extreme values.
 # 
 out = xy_to_utheta(y,x,prob = prob);
 u = out$u;
 theta = out$theta;
 #
 # Step 2: Computing the constraint.
 #
 c = get_c(u,theta);
 #
 # Step 3: Computing a random forest for the quantiles
 #
 rf <- get_qreg_forest(u,theta,n_trees = n_trees)
 #
 # Step 4: Calibrating the quantile level alpha to the constraint using 
 #        binary search and a naive Monte Carlo method.
 #
 alpha = calibrate_alpha(rf = rf,c = c, u=u,theta=theta, d= d, verbose=verbose)
 #
 # Step 5: Computing h_opt(x) = |x|*g_opt(x/|x|) 
 #
 norm_x_test = rowSums(x_test);
 n_test = length(norm_x_test);
 theta_test = x_test/(matrix(norm_x_test,n_test,1)%*%matrix(1,1,d))
 return(norm_x_test*get_q_alpha(rf,theta_test,alpha,verbose=verbose))
}
#
#
# G E N E R A L  P U R P O S E  W R A P P E R  F U N C T I O N
#
#
standardize_train_data <- function(y,x){
  #
  # To standard 1-Pareto margins using rank-transform.
  #
  n = length(y);
  if (length(dim(x))==0){
    xr = 1/(1-rank(x)/(n+1));
  } else{
    xr = 1/(1-apply(x,2,rank)/(n+1));
  }
  yr = 1/(1-rank(y)/(n+1));
  return(list("y"=yr,"x"=xr,"y_original"=y,"x_original"=x))
}

opt_pred <- function(y,x,x_test,prob=0.8,verbose=TRUE, extrapolate_tail_CDF = TRUE, n_trees=2000) {
  #
  # This is a very general wrapper function that: 
  #   (1) rank-transforms the training data to standard 1-Pareto
  #   (2) uses a quantile random forest to obtain the optimal homogeneous predictor
  #   (3) computes the optimal predictor over a test set of covariates, 
  #       which is quantile-transformed relative to the training sample.
  #.  (4) returns the rank-transformed training data, test data and the 
  #   predictor (on the 1-Pareto scale)
  #
  x_test = emp_cdf_transform(x_train = x, x_test = x_test, extrapolate_right_tail_via_GPD=extrapolate_tail_CDF);
  out = standardize_train_data(y,x);
  y_train = out$y;
  x_train = out$x;
  colnames(x_test) = colnames(x_train)
  y_pred = get_h_opt(y=y_train, x = x_train, x_test = x_test,prob = prob, verbose = verbose,n_trees=n_trees)
  return(list("y_train"=y_train,"x_train"=x_train,"y_pred"=y_pred,"x_test"=x_test))
}


#
#
#  E N D   O F   T H E F U N C T I O N S   I M P L E M E N T I N G 
#
#   O P T I M A L  H O M O G E N E O U S  P R E D I C T I O N
#
#.    Via Quantile Random Forest
#
##############################################################################

###############################################################################
#
# T E S T S  A N D  E X A M P L E S
#
###############################################################################

##########
#
# A U X I L I A R Y  F U N C T I O N S 
#
#  Simulation of:
#
#    * multivariate logistic (uses the "rstab" function for simulating stable 
#                             subordinators)
#    * Pareto Dirichlet
#
# For the simulation of stable subordinators:
#
source("~/Dropbox/doc/software/grand function respository/rstab.R")
#
# For tail-dependence testing:
#
source("~/Dropbox/doc/software/grand function respository/tail_dependence.R")
#
#
get_frechet <- function(alpha,n){
  return(1/abs(log(runif(n)))^(1/alpha))
}

sim_mvlogistic <- function(beta,s=c(),d,n){
  #
  # Simulates the standard multivariate logistic with CDF
  #
  # F(x) = exp{- (1/x1^beta + ... 1/xd^beta)^(1/beta)}
  #
  if (length(s)==0){
    s = matrix(1,n,d)} 
  else{
    s = matrix(1,n,1)%*%matrix(s,1,d);
  }
  xi = matrix(get_frechet(beta,n*d),n,d)*s;
  Z = matrix(rsub(1/beta,n),n,1)^(1/beta) %*%matrix(1,1,d); 
  return(Z*xi)
}

emp_cdf <- function(sample,x_new){
  #
  # Computes the empirical CDF based on the "sample"
  # applied to a vector of observations "x_new".
  #
  n = length(sample);
  p = length(x_new);
  return(colMeans(matrix(sample,n,1)%*%matrix(1,1,p) <= matrix(1,n,1)%*%matrix(x_new,1,p)));
}

#
#
# TO DO:
# implement PoT-based tail interpolation to compute the CDF for
# x_new basd on sample.
#
#. Implement your own GPD inference method.

gpd.n.loglik <- function(y,par){
  
}

fit.gpd <- function(x,prop=1-1/sqrt(length(x))){
    
}

cdf <- function(sample, x_new, prop=1-1/sqrt(length(unlist(sample)))){
  #
  # Computes the CDF by augmenting the right tail with a GPD fit.  Assumption: the sample is in 
  # the domain of attraction of GPD(xi) for xi>=0.
  #
  source("~/Dropbox/doc/software/grand function respository/tail_inference.R")
  tau = quantile(sample, prop);
  out = tGPD(sample,neg_xi = F, pu=prop);
  
  idx0 = which(x_new<= tau);
  idx1 = which(x_new> tau);
  
  p0 = length(idx0)
  n = length(sample);
  cdf = rep(0,length(x_new))
  
  cdf0 = colMeans(matrix(sample,n,1)%*%matrix(1,1,p0) <= matrix(1,n,1)%*%matrix(x_new[idx0],1,p0))
  
  if (length(idx1)>0){
    cdf1 = prop*gpd.cdf(x_new[idx1]-tau,out$xi,out$sig)
    cdf[idx1]=cdf1;
  } 
  cdf[idx0]=cdf0;
  return(cdf)
}

emp_cdf_transform <- function(x_train,x_test,Pareto_margins = TRUE, extrapolate_right_tail_via_GPD=TRUE){
  n = dim(x_train)[1];
  d = dim(x_train)[2];
  u = c();
  for (i in c(1:d)){
    if (extrapolate_right_tail_via_GPD){
      u = cbind(u, cdf(sample = x_train[,i], x_new = x_test[,i]))
    } else{
      u = cbind(u, emp_cdf(sample = x_train[,i], x_new = x_test[,i]))
    }
  }
  if (Pareto_margins){
    if (extrapolate_right_tail_via_GPD){
      return (1/(1-u))
    } else{
      return(1/(1-n*u/(n+1)))
    }
  } else{
    return(u)
  }
}

oracle_pred_mv_logistic <- function(beta,x_test,x_train=c()){
  #
  # Produces an empirical approximation of a 1-Frechet variable const*g_opt(X).  
  # When non-empty training data is provided it uses the empirical distribution 
  # of g_opt(X_train) to calibrate. Otherwise, it calibrates based on the 
  # rank-transformation of g_opt(X_test).
  #
  d = dim(x_test)[2];
  g_opt_test = 1/rowSums(1/abs(x_test)^beta)^(1/beta);
  # Calibrate relative to the training data, otherwise self-calibrate
  if (length(x_train)>0){ 
    g_opt_train = 1/rowSums(1/abs(x_test)^beta)^(1/beta);
    n = length(g_opt_train);
    u_test = emp_cdf(sample = g_opt_train, x_new = g_opt_train);
    return(1/abs(log((1+u_test*n)/(n+1))))
  } else { ## No training data provided.  Self-calibrate
    u_test = rank(g_opt_test)/(length(g_opt_test)+1);
    return(1/abs(log(u_test)))
  }
}

lambda_opt_mv_logistic <-function(beta,d,iter = 1e6){
  #
  # Computes the optimal extremal precision for the multivariate 
  # logistic model.
  #
  if (beta==1){return(0)}
  else if (beta<1) {return(NA)}
  else{
    G1 = rgamma(iter,shape = 1);
    Gd = rgamma(iter,shape = d);
    c1 = gamma(1-1/beta);
    cd = c1*prod((c(1:(d-1))-1/beta)/c(1:(d-1)));
    return( mean( pmin(G1^(-1/beta)/c1 , Gd^(-1/beta)/cd)))
  }
}


plot_opt_prec_mv_logistic <- function(filename = c(), monochrome=FALSE){
  library(ggplot2)
  library(reshape2) 
  library(metR)
  betas <- seq(1, 5, length.out =200)
  ds <- 2:20
  mat <- outer(betas, ds, Vectorize(lambda_opt_mv_logistic))
  
  colnames(mat) <- ds
  rownames(mat) <- betas
  mat_df <- melt(mat, varnames = c("beta", "d"))
  mat_df$beta <- as.numeric(as.character(mat_df$beta))
  mat_df$d <- as.numeric(as.character(mat_df$d))
  library(ggplot2)
  library(metR)
  contour_levels <- pretty(range(mat_df$value), n = 4)
  if (monochrome){
    low_color = "black";
    high_color= "white"
  } else {
    low_color = "red";
    high_color= "blue"
  }
  
  plt <- ggplot(mat_df, aes(x = beta, y = d)) +
    geom_tile(aes(fill = value)) +
    geom_contour(
      aes(z = value),
      color = "white",
      linetype = "dashed",
      breaks = contour_levels
    ) +
    metR::geom_text_contour(
      aes(z = value),
      breaks = contour_levels,
      stroke = 0.2,
      size = 3,
      check_overlap = FALSE,  # Label all contours
      skip = 0                # Label every contour line (no skipping)
    ) +
    scale_fill_gradient(
      low = low_color, high = high_color,
      name = expression(lambda)
    ) +
    labs(
      x = expression(beta),
      y = "d",
      title = expression("Heatmap of " * lambda^{(opt)} * " for the Gumbel copula")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )
  
  if (length(filename)>0){
    ggsave(filename, plot = plt, width = 5.03, height = 3.98)
  }
  return(list("plot"=plt,"mat"=mat,"mat_df"=mat_df))
}

test_mv_logistic_example <-function(beta=2.1,s=c(),d=10,n_train=1e3,n_test=1e3,prob=0.95){
  x = sim_mvlogistic(beta,s=s,d+1,n =n_train);
  y = x[,1];
  x = x[,-1];
  x_test = sim_mvlogistic(beta,s=s,d+1,n =n_test);
  y_test = x_test[,1];
  x_test = x_test[,-1];
  y_pred <- get_h_opt(y,x,x_test,prob=prob);
  return(list("y_test"=y_test,"y_pred"=y_pred,"x_test"=x_test,"x_train"=x))
}

example_plot_mv_logistic <-function(beta=2.1,s=c(), d=10,n_train=1e4,n_test=1e4,prob=0.95){

  out = test_mv_logistic_example(beta,s=c(),d,n_train=n_train,n_test=n_test,prob=prob);
  oracle = oracle_pred_mv_logistic(beta,x_test = out$x_test,  x_train = c())#out$x_train)
  #p = seq(0.8,0.9999,by=0.01);
  p = seq(0.9,0.999,by=0.001);
  emp_prec = tdep(out$y_test,out$y_pred,p);
  oracle_prec = tdep(out$y_test,oracle,p);
  lambda_opt = lambda_opt_mv_logistic(beta,d)
  
  df <- data.frame(p = p, emp_prec = emp_prec, oracle=oracle_prec, lambda_opt = lambda_opt)
  library(tidyr)
  df_long <- pivot_longer(df, cols = c(emp_prec, oracle, lambda_opt),
                          names_to = "curve", values_to = "value")
  
  plt <- ggplot(df_long, aes(x = p, y = value, linetype = curve, color = curve)) +
    geom_line(size = 1) +
    scale_color_manual(
      name = NULL,
      values = c(
        "emp_prec" = "blue",
        "oracle" = "black",
        "lambda_opt" = "red"
      ),
      labels = c(
        "emp_prec" =  expression(hat(lambda)),
        "oracle" = "Oracle",
        "lambda_opt" = expression(lambda^{(opt)})
      )
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c(
        "emp_prec" = "dashed",
        "oracle" = "solid",
        "lambda_opt" = "dotted"
      ),
      labels = c(
        "emp_prec" = expression(hat(lambda)),
        "oracle" = "Oracle",
        "lambda_opt" = expression(lambda^{opt})
      )
    ) +
    guides(
      color = guide_legend(override.aes = list(linetype = c("solid", "dotdash", "dashed"))),
      linetype = "none"  # suppress the second (redundant) legend
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = bquote("Gumbel copula: " ~ beta == .(beta) ~ ", " ~ d == .(d)),
      x = "p",
      y = "Precision"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5)
    )
  return(plt)
}

###############################################################################
#
#  S P E C T R A L L Y  D I R I C H L E T 
#
#

sim_spec_Dirichlet <- function(n,a,s=c()){
  a_star = sum(a);
  d = length(a);
  G = c();
  if (length(s)==0){s=rep(1,d)}
  for (i in c(1:d)){
    G= cbind(G, s[i]*rgamma(n,shape = a[i]));
  }
  R = matrix(1/runif(n),n,1);
  G_star = matrix(rowSums(G),n,1); 
  e_d = matrix(1,1,d);
  return( (R%*% e_d)*( G/(G_star%*%e_d) ) )
}

test_spec_Dirichlet_example <-function(a = c(1:15)/10,s= c(), n_train=1e4,n_test=1e4,prob=0.95){
  x = sim_spec_Dirichlet(n=n_train,a=a,s=s);
  y = x[,1];
  x = x[,-1];
  x_test = sim_spec_Dirichlet(n=n_test, a=a, s=s);
  y_test = x_test[,1];
  x_test = x_test[,-1];
  y_pred <- get_h_opt(y,x,x_test,prob=prob);
  return(list("y_test"=y_test,"y_pred"=y_pred,"x_test"=x_test,"x_train"=x))
}

lambda_opt_spec_Dirichlet <- function(a0,a1){
  #
  # Computes the extremal precision in the Dirichlet(a) model
  # a0 = a[1]; a1 = a[2]+...+a[d]
  mu= a0/(a0+a1)
  x = (1/mu)*(a0/(a0+a1))*pbeta(mu,a0+1,a1);
  y = (1/(1-mu))*(a1/(a0+a1))*(1-pbeta(mu,a0,a1+1))
  return(x+y)
}

oracle_pred_spec_Dirichlet <- function(x_test,x_train=c()){
  #
  # Produces an empirical approximation of a 1-Pareto variable const*h_opt(X).  
  # When non-empty training data is provided it uses the empirical distribution 
  # of h_opt(X_train) to calibrate. Otherwise, it calibrates based on the 
  # rank-transformation of h_opt(X_test).
  #
  d = dim(x_test)[2];
  h_opt_test =rowSums(abs(x_test));
  # Calibrate relative to the training data, otherwise self-calibrate
  if (length(x_train)>0){ 
    h_opt_train = abs(x_test);
    n = length(g_opt_train);
    u_test = emp_cdf(sample = h_opt_train, x_new = h_opt_train);
    return(1/(1-u_test*n/(n+1)))
  } else { ## No training data provided.  Self-calibrate
    u_test = rank(h_opt_test)/(length(h_opt_test)+1);
    return(1/(1-u_test))
  }
}

example_plot_spec_Dirichlet <-function(a = c(1:25)/10,s = c(), n_train=1e4,n_test=1e4,prob=0.95){
  out = test_spec_Dirichlet_example(a=a,s=s,n_train=n_train,n_test=n_test,prob=prob);
  oracle = oracle_pred_spec_Dirichlet(x_test = out$x_test,  x_train = c())
  p = seq(0.9,0.999,by=0.001);
  emp_prec = tdep(out$y_test,out$y_pred,p);
  lambda_opt = lambda_opt_spec_Dirichlet(a0=a[1],a1=sum(a[-1]));
  oracle_prec = tdep(out$y_test,oracle,p);
 
  df <- data.frame(p = p, emp_prec = emp_prec, oracle=oracle_prec, lambda_opt = lambda_opt)
  library(tidyr)
  df_long <- pivot_longer(df, cols = c(emp_prec, oracle, lambda_opt),
                          names_to = "curve", values_to = "value")
  
  plt <- ggplot(df_long, aes(x = p, y = value, linetype = curve, color = curve)) +
    geom_line(size = c(1)) +
    scale_color_manual(
      name = NULL,
      values = c(
        "emp_prec" = "blue",
        "oracle" = "black",
        "lambda_opt" = "red"
      ),
      labels = c(
        "emp_prec" =  expression(hat(lambda)),
        "oracle" = "Oracle",
        "lambda_opt" = expression(lambda^{(opt)})
      )
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c(
        "emp_prec" = "dashed",
        "oracle" = "solid",
        "lambda_opt" = "dotted"
      ),
      labels = c(
        "emp_prec" = expression(hat(lambda)),
        "oracle" = "Oracle",
        "lambda_opt" = expression(lambda^{opt})
      )
    ) +
    guides(
      color = guide_legend(override.aes = list(linetype = c("solid", "dotdash", "dashed"))),
      linetype = "none"  # suppress the second (redundant) legend
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = paste0("Pareto-Dirichlet model: d=",length(a)-1),
      x = "p",
      y = "Precision"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5)
    )
  return(plt)
}

plot_opt_prec_spec_Dirichlet <- function(filename = c(), monochrome=FALSE){
  library(ggplot2)
  library(reshape2) 
  library(metR)
  a0 <- seq(0.25, 15, length.out = 100);#by = 0.05)
  a1 <- seq(0.25, 15, length.out = 100)
  mat <- outer(a0, a1, Vectorize(lambda_opt_spec_Dirichlet))
  
  colnames(mat) <- a0
  rownames(mat) <- a1
  mat_df <- melt(mat, varnames = c("a0", "a1"))
  mat_df$a0 <- as.numeric(as.character(mat_df$a0))
  mat_df$a1 <- as.numeric(as.character(mat_df$a1))
  library(ggplot2)
  library(metR)
  contour_levels <- pretty(range(mat_df$value), n = 4)
  if (monochrome){
    low_color = "black";
    high_color= "white"
  } else {
    low_color = "red";
    high_color= "blue"
  }
  
  plt <- ggplot(mat_df, aes(x = a0, y = a1)) +
    geom_tile(aes(fill = value)) +
    geom_contour(
      aes(z = value),
      color = "white",
      linetype ="dashed",
      breaks = contour_levels
    ) +
    metR::geom_text_contour(
      aes(z = value),
      breaks = contour_levels,
      stroke = 0.2,
      size = 3,
     
      check_overlap = FALSE,  # Label all contours
      skip = 0                # Label every contour line (no skipping)
    ) +
    scale_fill_gradient(
      low = low_color, high = high_color,
      name = expression(lambda)
    ) +
    labs(
      x = expression(alpha[0]),
      y = expression("||" * alpha* "||"[1]),
      title = expression("Heatmap of " * lambda[G]^{(opt)} * " for the Pareto-Dirichlet")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )
  
  if (length(filename)>0){
    ggsave(filename, plot = plt, width = 5.03, height = 3.98)
  }
  return(list("plot"=plt,"mat"=mat,"mat_df"=mat_df))
}

gen_figure_Dirichlet_pred <- function(dir="~/Dropbox/doc/software/grand function respository/optXpred/figures/"){
  
  plt1 = example_plot_spec_Dirichlet(a = c(1:5),n_train = 1e4,n_test = 1e4);
  plt2 = example_plot_spec_Dirichlet(a = c(1:41)/10,n_train = 1e4,n_test = 1e4)

ggsave(paste0(dir,"Spec_Dirichlet_d4.pdf"), plot=plt1,  width = 5.03, height = 3.98)
ggsave(paste0(dir,"Spec_Dirichlet_d40.pdf"), plot=plt2,  width = 5.03, height = 3.98)

}

gen_figure_mv_logistic_pred <- function(dir="~/Dropbox/doc/software/grand function respository/optXpred/figures/"){
  
  plt1 = example_plot_mv_logistic(beta=1.2,d=10)
  plt2 = example_plot_mv_logistic(beta=3,d=100)
  
  ggsave(paste0(dir,"MV_Logistic_d10.pdf"), plot=plt1,  width = 5.03, height = 3.98)
  ggsave(paste0(dir,"MV_Logistic_d100.pdf"), plot=plt2,  width = 5.03, height = 3.98)
  
}

