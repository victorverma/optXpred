library(arrow)
library(lubridate)
library(dplyr)
library(fields)

source("~/Dropbox/doc/software/grand function respository/tail_dependence.R")
source("~/Dropbox/doc/software/grand function respository/optXpred/opt_ext_pred.R")
setwd("~/Dropbox/doc/software/grand function respository/optXpred/")
get_flux_and_sharp_data <- function(block_max_window_in_hours = 6,lag_hours = 24, scramble=F){
#data <- read_parquet("~/Dropbox-UniversityofMichigan/Stilian Stoev/harp_xray_flux_data/1m_data.parquet")
data <- read_parquet("~/Dropbox/doc/data/harp_xray_flux_data/1m_data.parquet");
flux = bmax(data$flux[-1],12,sliding_blocks = F)
df_flux = data.frame(time = data$time[12*c(1:length(flux))+1],flux)
df_flux$time = as_datetime(df_flux$time)

#df_all <-  read_parquet("~/Dropbox-UniversityofMichigan/Stilian Stoev/harp_xray_flux_data/aggregated_high-qual_near-center-70.parquet")
df_all <- read_parquet("~/Dropbox/doc/data/harp_xray_flux_data/aggregated_high-qual_near-center-70.parquet")
df_all$T_REC <- as_datetime(df_all$T_REC)
df_all <- left_join(x=df_all,y=df_flux,by=c("T_REC"="time"))
df <- df_all[complete.cases(df_all), ] # eliminating rows with missing values

n = nrow(df);
lag = 5*lag_hours;
#lagged_flux = c(rep(NA,1,lag),df$flux[1:(n-lag)])
df$lagged_flux = df$flux;
df$flux[1:(n-lag)] = df$flux[(lag+1):n]
#=  <- df[]
#df <- df[complete.cases(df), ] # eliminating rows with missing values
df <- df[1:(n-lag),]                

N <- nrow(df)
k <- ncol(df)
m <- floor(5 * block_max_window_in_hours)
n <- floor(N / m)

df_old <- df
df <- df[1:n, ]

# Compute T_REC for each block
df$T_REC <- df_old$T_REC[c(1:n)*m]

# Apply bmax to each of the remaining columns (excluding the first one)
for (i in 2:k) {
  # bmax must return at least n values
  df[[i]] <- bmax(df_old[[i]], m, sliding_blocks = FALSE)[1:n]
}
colnames(df)=names(df_old)



if (scramble){
  n = length(df$T_REC);
  k = length(df)
  df[,2:k] = df[order(runif(n)),2:k];
}

return(df)
}

#make_24_hour_max_flux = FALSE
#extrapolate_tail_CDF = FALSE
#if (make_24_hour_max_flux){
#  m = 10;
#  df$flux <- c(bmax(df$flux,m=m,sliding_blocks = TRUE), rep(NA,m-1))
#}
##df$flux = df$flux24
#df <- df[complete.cases(df), ] # eliminating rows with missing values

df = get_flux_and_sharp_data(block_max_window_in_hours = 24,lag_hours = 24,scramble = F)
n_cols = dim(df)[2]
L = c(1);
for (i in c(2:(n_cols))){
  L = c(L,tdep(df[[i]],df$flux,p=0.95))
}
D = matrix(0,n_cols-1,n_cols-1);
for (i in c(2:n_cols)){
  for (j in c(2:n_cols)){
    D[i-1,j-1] = tdep(df[[i]],df[[j]],p=0.95)
  }
}
keep = which(L>=0.25);

pred_flux <- function(df,train_interval,
                         test_interval, block_max_window_in_hours =1,
                         lag_hours = 24, prob=0.95, p = seq(0.5,0.999,by=0.01), to_plot=F,
                         verbose = F, use_flux = T,keep_threshold=0.9, extrapolate_tail_CDF=TRUE,selected_variables=c()){
  
library(lubridate)
#lag = 5*lag_hours; #offsetting the predictors by "lag_hours" hours in the past.
lag = floor(lag_hours/block_max_window_in_hours);

nsec_in_year = 365.24*60*60*24;
i_train = interval(start = as_datetime((train_interval[1]-1970)*nsec_in_year), as_datetime((train_interval[2]-1970)*nsec_in_year))
i_test = interval(start = as_datetime((test_interval[1]-1970)*nsec_in_year), as_datetime((test_interval[2]-1970)*nsec_in_year))

train <- which(df$T_REC %within% i_train);

n_train = length(train)
test <- which( df$T_REC %within% i_test);
n_test = length(test)

idx_lagged_flux = which(names(df)=="lagged_flux")
idx_flux = which(names(df)=="flux")
idx_time =1
variables = colnames(df);

if (length(selected_variables)>1){
  keep_variables = setdiff(c(variables[idx_time],selected_variables),"flux")
} else{
if (use_flux==FALSE){
keep = tdep_partial_select(y = df$flux[train], 
                           x = df[train,variables[-c(idx_time,idx_flux,idx_lagged_flux)]],keep_threshold)
} else {
  keep = tdep_partial_select(y = df$flux[train], 
                             x = df[train,variables[-c(idx_time,idx_flux)]],keep_threshold)
}
keep_variables = c(variables[idx_time], variables[keep]);
}

if (length(use_flux)>0){ #NOTE: if use_flux =c(), then we don't care if lagged_flux is used or not.
if (use_flux==TRUE){
  keep_variables = union(keep_variables,"lagged_flux");
} else {
  keep_variables = setdiff(keep_variables,"lagged_flux");
}
}
if (verbose){
 cat("\n Predicting using variables: ", setdiff(keep_variables,"T_REC"))
}
y_train <- df$flux[train]
x_train <- as.matrix(df[train,keep_variables[-1]])
#y_train = y_train[-c(1:lag)];
x_train <- x_train[1:length(y_train),]

y_test <- df$flux[test];
#y_test = y_test[-c(1:lag)];
x_test <- as.matrix(df[test,keep_variables[-1]])
x_test <- x_test[1:length(y_test),]

out <- opt_pred(y=y_train,x=x_train,x_test = x_test,prob = prob,verbose = verbose, extrapolate_tail_CDF = extrapolate_tail_CDF);
x_test = emp_cdf_transform(x_train = x_train, x_test = x_test);
#simple_pred <- unlist(apply(x_test,1,max))
emp_prec = tdep(y_test,out$y_pred,p);

prop_train_M = mean(y_train < 1e-5);
prop_train_X = mean(y_train < 1e-4);

#out_M = compute_metrics(obs = y_test >1e-5,pred = (1-1/out$y_pred)>0.4)#prop_train_M);
range_metrics_M = PAM_metrics(obs = y_test >1e-5,pred = (1-1/out$y_pred));

out_M = compute_metrics(obs = y_test >1e-5,pred = (1-1/out$y_pred)>prop_train_M);
PAM_C = PAM_metrics(obs= y_test >1e-6,pred = (1-1/out$y_pred))
PAM_M = PAM_metrics(obs= y_test >1e-5,pred = (1-1/out$y_pred))

TSS_M = out_M$TSS
prec_M = out_M$prec
F1_M = out_M$F1;
missed_M = out_M$missed;
alarm_M = out_M$alarm;
#out_X = compute_metrics(obs = y_test >1e-4,pred = (1-1/out$y_pred)>0.75)#prop_train_X);
out_X = compute_metrics(obs = y_test >1e-4,pred = (1-1/out$y_pred)>prop_train_X);
PAM_X = PAM_metrics(obs= y_test >1e-4,pred = (1-1/out$y_pred))
TSS_X = out_X$TSS
prec_X = out_X$prec
F1_X = out_X$F1;
missed_X = out_X$missed;
alarm_X = out_X$alarm
#emp_prec_simple_pred =  tdep(y_test,simple_pred,p)
#emp_prec =  tdep(y_test,rowMeans(x_test),p);

TSS_pred = TSS(emp_prec,p);
p_M_class = mean(df$flux[train]<1e-5);
p_X_class = mean(df$flux[train]<1e-4);

if (to_plot){

plot(p,emp_prec,type="l",ylim=c(0,1), cex=0.25, ylab="Precision & TSS", col="black",
     main=paste0("Train interval: ",i_train,
                 "\nTest interval:",i_test))
lines(p,TSS_pred,col="black",lty=2)
abline(v=c(p_M_class,p_X_class),col=c("blue","red"),lty=2)
legend("topright",legend = c("prec","TSS","prop_M","prop_X"),col=c("black","black","blue","red"),lty=c(1,2,2,2),lwd=c(1,1,1,1))

}


if (max(p)< p_M_class){
  TSS_M = NA
  prec_M_tdep = NA
} else {
  ip = min(which(p_M_class<=p));
  prec_M_tdep = emp_prec[ip];
}

return(list("TSS_M"=TSS_M,"prec_M"=prec_M,"F1_M"=F1_M, "TSS_X"=TSS_X,"prec_X"=prec_X,"F1_X"=F1_X,
            "missed_M"=missed_M,"missed_X"=missed_X,"alarm_M"=alarm_M,"alarm_X"=alarm_X,
            "emp_prec"=emp_prec,"variables" = keep_variables[-1],"TSS_pred"=TSS_pred,"p_M_class"=p_M_class,
            "p_X_class"=p_X_class,
            "PAM_M"=PAM_M,"PAM_X"=PAM_X,"PAM_C"=PAM_C,"n_train"=n_train,"n_test"=n_test))
}

scrambled_examples <- function(lag_hours_block_maxima = 24, lag_hours_prediction = 24, iter = 100, use_flux =TRUE, 
                                extrapolate_tail_CDF =FALSE, to_plot = F, 
                                selected_variables = c("TOTUSJH","AREA_ACR","MEANALP","ABSNJZH",
                                                       "MEANPOT","lagged_flux","NPIX","SAVNCPP","SHRGT45")){
  
  prec_M = c();
  prec_X = c();
  emp_tss_M = c();
  emp_tss_X = c();
  F1_M = c();
  F1_X = c();
  variables = list();
  for (i in c(1:iter)){
    if (to_plot==TRUE){
     cat('.')
    }
  df = get_flux_and_sharp_data(block_max_window_in_hours = lag_hours_block_maxima,lag_hours = lag_hours_prediction,scramble = T)
  out = pred_flux(df,c(2011,2020),c(2021,2023),to_plot = to_plot,verbose = F,use_flux = T, extrapolate_tail_CDF = extrapolate_tail_CDF,
  keep_threshold = 0.9,prob = 0.9, selected_variables = selected_variables); 
  prec_M = c(prec_M, out$prec_M);
  prec_X = c(prec_X, out$prec_X);
  emp_tss_M = c(emp_tss_M,out$TSS_M)
  emp_tss_X = c(emp_tss_X,out$TSS_X)
  F1_M = c(F1_M, out$F1_M)
  F1_X = c(F1_X, out$F1_X)
  variables = c(variables,list(out$variables))
  }
  x = cbind(prec_M,prec_X, emp_tss_M, emp_tss_X,F1_M,F1_X)
  colnames(x)= c("prec(M)","prec(X)","TSS(M)","TSS(X)","F1(M)","F1(X)");
  if (to_plot){
   boxplot(x,ylim=c(0,1))
   abline(h = 0.5,lty=2)
  }
  return(list("res"=x,"variables"=variables))
}

#
#
# Good sets of variables:
#
#


good_for_X = c("AREA_ACR","MEANJZH","ABSNJZH","SHRGT45",
               "SAVNCPP","MEANGBZ","MEANSHR","NACR",
               "TOTUSJZ","SIZE","lagged_flux");

good_for_M_and_X = c("TOTUSJH", "AREA_ACR","ABSNJZH",
               "SAVNCPP", "TOTPOT");

best_2_M <- c("ABSNJZH","AREA_ACR","lagged_flux")
best_2_X <- c("ABSNJZH","AREA","lagged_flux")


best_3_M <- c("TOTUSJZ","TOTUSJH","ABSNJZH","lagged_flux")
best_3_X <- c("MEANPOT","TOTPOT","MEANSHR","lagged_flux" )

figure_good_for_M_and_X <- function(){
  
   
   
   out_fixed_24_best_for_X = scrambled_examples(lag_hours_block_maxima = 24, lag_hours_prediction = 24, 
                                                selected_variables = good_for_X,
                                                use_flux = T,iter = 20,extrapolate_tail_CDF = F)
  
   out_fixed_24_best_for_M = scrambled_examples(lag_hours_block_maxima = 24, lag_hours_prediction = 24,
                                                selected_variables = good_for_M, to_plot = T,
                                                use_flux = T,iter = 20,extrapolate_tail_CDF = F)
   
   out_6 = scrambled_examples(lag_hours_block_maxima = 6, lag_hours_prediction = 6,
                                                selected_variables = c("TOTUSJH", "AREA_ACR","ABSNJZH",
                                                                       "SAVNCPP", "TOTPOT"),
                                                use_flux = T,iter = 20,extrapolate_tail_CDF = F)
   out_12 = scrambled_examples(lag_hours_block_maxima = 12, lag_hours_prediction = 12,
                              selected_variables = c("TOTUSJH", "AREA_ACR","ABSNJZH",
                                                     "SAVNCPP", "TOTPOT"),
                              use_flux = T,iter = 20,extrapolate_tail_CDF = F)
   
   out_24 = scrambled_examples(lag_hours_block_maxima = 24, lag_hours_prediction = 24,
                              selected_variables = c("TOTUSJH", "AREA_ACR","ABSNJZH",
                                                     "SAVNCPP", "TOTPOT"),
                              use_flux = T,iter = 20,extrapolate_tail_CDF = F)
   out_48 = scrambled_examples(lag_hours_block_maxima = 48, lag_hours_prediction = 48,
                            selected_variables = c("TOTUSJH", "AREA_ACR","ABSNJZH",
                                                     "SAVNCPP", "TOTPOT"),
                              use_flux = T,iter = 20,extrapolate_tail_CDF = F)
   
   res_M = cbind(out_6$res[,3],out_12$res[,3],out_24$res[,3],out_48$res[,3])
   boxplot(res_M,ylim=c(0,1))
   
   out_X6 = scrambled_examples(lag_hours_block_maxima = 6, lag_hours_prediction = 6, 
                                                selected_variables = c("AREA_ACR","MEANJZH","ABSNJZH","SHRGT45",
                                                                       "SAVNCPP","MEANGBZ","MEANSHR","NACR",
                                                                       "TOTUSJZ","SIZE","lagged_flux"),
                                                use_flux = T,iter = 20,extrapolate_tail_CDF = F)
   out_X12 = scrambled_examples(lag_hours_block_maxima = 12, lag_hours_prediction = 12, 
                               selected_variables = c("AREA_ACR","MEANJZH","ABSNJZH","SHRGT45",
                                                      "SAVNCPP","MEANGBZ","MEANSHR","NACR",
                                                      "TOTUSJZ","SIZE","lagged_flux"),
                               use_flux = T,iter = 20,extrapolate_tail_CDF = F)
   out_X24 = scrambled_examples(lag_hours_block_maxima = 24, lag_hours_prediction = 24, 
                               selected_variables = c("AREA_ACR","MEANJZH","ABSNJZH","SHRGT45",
                                                      "SAVNCPP","MEANGBZ","MEANSHR","NACR",
                                                      "TOTUSJZ","SIZE","lagged_flux"),
                               use_flux = T,iter = 20,extrapolate_tail_CDF = F)
   out_X48 = scrambled_examples(lag_hours_block_maxima = 48, lag_hours_prediction = 48, 
                               selected_variables = c("AREA_ACR","MEANJZH","ABSNJZH","SHRGT45",
                                                      "SAVNCPP","MEANGBZ","MEANSHR","NACR",
                                                      "TOTUSJZ","SIZE","lagged_flux"),
                               use_flux = T,iter = 20,extrapolate_tail_CDF = F)
   
}


next_k_tuple <- function(a,n){
  #
  # a = (a_1,..,a_k) is a vector of non-negative spacings such that sum(a)+k <=n 
  # the function returns the "next" vector of spacings obtained lexicographically 
  # or "NA" if none exists.
  #
  k = length(a);
  if (a[1]==n-k){
    return(c())}
  else{
     if (sum(a)==n-k){
       i0 = max(which(a>0));
       a[i0]=0
       a[i0-1] = a[i0-1]+1;
     } else {
       a[k] = a[k]+1;
     }
    return(a)
  }
}

test_k_tuple <- function(k = 3, n=10){
  as = rep(0,k);
  a = as;
  while(length(a)>0){
    a = next_k_tuple(a,n);
    as = rbind(as,a)
  }
  return(as)
}

find_best_variables_for_M_and_X <- function(df, k = 3){
  variable_names = names(df)
  idx_flux = which(variable_names=="flux");
  variable_names = variable_names[-c(1,idx_flux)];
  n_vars = length(variable_names);
  a = rep(0,k);
  res = c();
  while(length(a)>0){
    out = scrambled_examples(lag_hours_block_maxima = 24, 
                             lag_hours_prediction = 24,
                             selected_variables = variable_names[cumsum(a+1)],
                             use_flux = c(),iter = 20,to_plot = F, extrapolate_tail_CDF = F);
    TSS.M =  mean(out$res[,3]);
    TSS.X = mean(out$res[,4]);
    
    vars =  variable_names[cumsum(a+1)];
    res = rbind(res,data.frame("TSS.M" =TSS.M , "TSS.X" = TSS.X,variables = I(list(out$variables[[1]])),
                "tss.M.sample" =I(list(out$res[,3])), "tss.X.sample" =I(list(out$res[,4]))))
    cat("\n", vars," : TSS.M ",TSS.M," TSS.X ",TSS.X);
    a = next_k_tuple(a,n_vars)
    save(file=paste0("best_flux_pred_vars_k",k,".RData"),list = c("res"))
  }
  return(res)
}

#
#
#
#

some_examples_BROKEN <- function(train_period = c(2010+5/12,2010+7/12),keep_threshold=0.95,to_plot=T,p = seq(0.4,0.99,by=0.01)){
  #keep = tdep_partial_select(y = df$flux[year(df$T_REC)==year], 
  #                          x = df[year(df$T_REC)==year,-c(1,24)],keep_threshold)
  #names(df)[keep]
  #[1]
  #keep = c(1, 8, 11, 12, 13, 15, 24);
  #keep =c(1,1+keep,24)
  #keep_variables = names(df)[keep];
  #cat(keep_variables)
#  keep_variables
#[1] "T_REC"   "TOTUSJZ" "TOTUSJH" "ABSNJZH" "SAVNCPP" "TOTPOT"  "flux"   
  out_with_flux = pred_flux(df,train_period,test_period,lag_hours=24,to_plot=T,prob = 0.99, 
                  p = seq(0.4,0.99,by=0.01),use_flux = T);
  out_no_flux = pred_flux(df,year,month,year,month,year,month+1,lag_hours=24,to_plot=T,prob = 0.99, 
                            p = seq(0.4,0.99,by=0.01),use_flux = F);
 if(to_plot){
  plot(p,out_with_flux$emp_prec,type="l",cex=0.25,lwd=2,ylim=c(-.5,3),ylab="Precision & TSS",
       main=year);
  lines(p,out_with_flux$emp_prec_simple_pred,col="blue")
  lines(p,out_no_flux$emp_prec,lty=2,lwd=2)
  #cat(out_no_flux$emp_prec_simple_pred)
  cat(TSS(out_with_flux$emp_prec,p=p))
  lines(p,out_no_flux$emp_prec_simple_pred,col="blue",lty=2)
  #lines(p,TSS(out_with_flux$emp_prec,p=p),col="red",lty=1,lwd=2);
  #lines(p,TSS(out_no_flux$emp_prec,p=p),col="red",lty=2,lwd=2);
  legend("topleft",legend=c("emp prec","simple", "emp prec (NO flux)","simple (NO FLUX)"),
                            #,"TSS", "TSS NO flux"),
         lwd=c(1,1,1,1),
         col=c("black","blue","black","blue"),lty=c(1,1,2,2))
  abline(h=0,lty=3,lwd=0.5) 
  abline(h =0.5, lty=3,lwd=0.5) 
  abline(v=out_with_flux$M_class,lty=3,lwd=2)
}}

compute_over_time_BROKEN <- function(df){
ym <- df %>% mutate(year=year(df$T_REC),month=month(df$T_REC))
ym <- unique(cbind(ym$year,ym$month))

n_ym = length(ym[,1]);
lag_hours = 12;
k=3;
res = c();
pb <- txtProgressBar(min = 1, max = n_ym, style = 3)
for (i in c(1:(n_ym-k))){
  res = cbind(res,pred_flux(df,train_year0=ym[i,1],train_month0=ym[i,2],
                               train_year1=ym[i+(k-1),1],train_month1=ym[i+(k-1),2],
                               test_year=ym[i+k,1],test_month = ym[i+k,2],
                            lag_hours = lag_hours))
  setTxtProgressBar(pb, i)
}
close(pb)
return(res)
}

TSS <- function(emp_prec,p){
  return((emp_prec-(1-p))/p)
}

compute_tss <- function(obs, pred) { # From ChatGPT
  if (length(obs) != length(pred)) {
    stop("Observed and predicted vectors must be of the same length.")
  }
  
  # Convert to logical for clarity
  obs <- as.logical(obs)
  pred <- as.logical(pred)
  
  # Confusion matrix components
  TP <- sum(obs & pred)
  TN <- sum(!obs & !pred)
  FP <- sum(!obs & pred)
  FN <- sum(obs & !pred)
  
  sensitivity <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  specificity <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  
  TSS <- sensitivity + specificity - 1
  return(list("TSS"=TSS,"prec"=TP/(TP+FP),"FNR"=(FN)/(FN+TN)))
}


compute_metrics <- function(obs, pred) {
  if (length(obs) != length(pred)) {
    stop("Observed and predicted vectors must be of the same length.")
  }
  
  obs <- as.integer(obs)
  pred <- as.integer(pred)
  
  # Confusion matrix components
  TP <- sum(obs == 1 & pred == 1)
  TN <- sum(obs == 0 & pred == 0)
  FP <- sum(obs == 0 & pred == 1)
  FN <- sum(obs == 1 & pred == 0)
  
  # TSS components
  sensitivity <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  specificity <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  TSS <- sensitivity + specificity - 1
  # F1 components
  precision <- if ((TP + FP) == 0) NA else TP / (TP + FP)
  recall <- sensitivity
  F1 <- if (is.na(precision) || is.na(recall) || (precision + recall == 0)) NA else {
    2 * precision * recall / (precision + recall)
  }
  
  return(list(
    TSS = TSS,
    prec = precision,
    recall = recall,
    missed = 1-TP/(TP+FN),
    alarm = mean(pred),
    F1 = F1,
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN
  ))
}

PAM_metrics <- function(obs, pred, p_level =seq(from=0.01,to=0.99,length.out=100)){
#
# computes the precision, alarm, and missed metrics based on the prob threshold "p_level"
#
  
precision = c();
alarm = c();
missed = c();
TSS = c()
for (p in p_level){
  out = compute_metrics(obs,pred>p);
  precision = c(precision, out$prec);
  alarm = c(alarm, out$alarm)
  missed = c(missed, out$missed)
  TSS = c(TSS,out$TSS)
}
return(list("prec"=precision,"alarm"=alarm,"missed"=missed,"TSS"=TSS,"p_level"=p_level))
}

PAM_plot <-function(PAM,title_text=""){
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  # Convert data to long format
  df_long <- as.data.frame(PAM) %>%
    select(p_level, prec, alarm, TSS, missed) %>%
    pivot_longer(cols = -p_level, names_to = "Metric", values_to = "Value")
  
  # Define custom line types and colors
  metric_colors <- c(
    prec = "black",
    alarm = "red",
    TSS = "black",
    missed = "blue"
  )
  
  metric_linetypes <- c(
    prec = "solid",
    alarm = "solid",
    TSS = "dashed",
    missed = "solid"
  )
  
  metric_sizes <- c(
    prec = 0.5,
    alarm = 1,
    TSS = 1.5,
    missed = 1
  )
  
  # Plot
  plt = ggplot(df_long, aes(x = p_level, y = Value, color = Metric, linetype = Metric, size = Metric)) +
    geom_line() +
    scale_color_manual(values = metric_colors) +
    scale_linetype_manual(values = metric_linetypes) +
    scale_size_manual(values = metric_sizes) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = title_text,
      x = "Prob threshold",
      y = "Metric Value",
      color = "Metric",
      linetype = "Metric"
    ) +
    theme_minimal(base_size = 12)
  return(plt)
}

PAM_plot_by_alarm <- function(PAM, title_text = "") {
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  # Convert data to long format with "alarm" as x-axis
  df_long <- as.data.frame(PAM) %>%
    select(alarm, prec, TSS, missed) %>%
    pivot_longer(cols = -alarm, names_to = "Metric", values_to = "Value")
  
  # Define custom line types and colors
  metric_colors <- c(
    prec = "blue",
    TSS = "black",
    missed = "red"
  )
  
  metric_linetypes <- c(
    prec = "solid",
    TSS = "dashed",
    missed = "solid"
  )
  
  metric_sizes <- c(
    prec = 0.5,
    TSS = 1.5,
    missed = 1
  )
  
  # Plot
  plt <- ggplot(df_long, aes(x = alarm, y = Value, color = Metric, linetype = Metric, size = Metric)) +
    geom_line() +
    scale_color_manual(values = metric_colors) +
    scale_linetype_manual(values = metric_linetypes) +
    scale_size_manual(values = metric_sizes) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = title_text,
      x = "Alarm Rate",
      y = "Metric Value",
      color = "Metric",
      linetype = "Metric"
    ) +
    theme_minimal(base_size = 12)
  
  return(plt)
}


 gen_figures_xray_flux_pred <- function(dir="~/Dropbox/doc/software/grand function respository/optXpred/figures/",via_alarm=TRUE){
   df = get_flux_and_sharp_data(24,24,F)
   good_for_M_and_X = c("TOTUSJH", "AREA_ACR","ABSNJZH","SAVNCPP", "TOTPOT");
   #good_for_X = c("USFLUX","TOTUSJH", "AREA_ACR","ABSNJZH","SAVNCPP", "TOTPOT");
   #good_for_X = c("USFLUX","AREA_ACR","ABSNJZH","SAVNCPP", "TOTPOT");
   
   set.seed(2)
   out = pred_flux(df,c(2010,2017),test_interval = c(2017.01,2025), to_plot = T, selected_variables = good_for_M_and_X,prob = 0.97,extrapolate_tail_CDF = F); 
  ## out = pred_flux(df,c(2010,2017),test_interval = c(2017.01,2025), to_plot = T, selected_variables = good_for_X,prob = 0.97,extrapolate_tail_CDF = F); 
  ## out = pred_flux(df,c(2010,2017),test_interval = c(2017.01,2025), to_plot = T, selected_variables = c(), prob = 0.97,extrapolate_tail_CDF = F); 
   
   if(via_alarm){
     plt_C = PAM_plot_by_alarm(out$PAM_C,"Train (2010-2017) Test (2017-2025): C+ class")
     plt_M = PAM_plot_by_alarm(out$PAM_M,"Train (2010-2017) Test (2017-2025): M+ class")
     plt_X = PAM_plot_by_alarm(out$PAM_X,"Train (2010-2017) Test (2017-2025): X+ class")
     
   }else{
   plt_C = PAM_plot(out$PAM_C,"Train (2010-2017) Test (2017-2025): C+ class")
   plt_M = PAM_plot(out$PAM_M,"Train (2010-2017) Test (2017-2025): M+ class")
   plt_X = PAM_plot(out$PAM_X,"Train (2010-2017) Test (2017-2025): X+ class")
   }
    if (length(dir)>0){
     string = "";
     if (via_alarm){string="VIA_ALARM"}
     ggsave(paste0(dir,"Xray_flux_flares_C_class_",string,".pdf"), plot = plt_C, width = 5.03, height = 3.98)
     ggsave(paste0(dir,"Xray_flux_flares_M_class_",string,".pdf"), plot = plt_M, width = 5.03, height = 3.98)
     ggsave(paste0(dir,"Xray_flux_flares_X_class_",string,".pdf"), plot = plt_X, width = 5.03, height = 3.98)
   }
 }

 gen_figure_GOES_flux <- function(dir = "~/Dropbox/doc/software/grand function respository/optXpred/figures/"){
   df = get_flux_and_sharp_data(24,24,F);
   library(ggplot2)
   
   plt = ggplot(df, aes(x = T_REC, y = flux)) +
     geom_line(linewidth = 0.5, color = "black") +
     geom_hline(yintercept = 1e-4, linetype = "dashed", color = "red") +
     #scale_y_log10() +  # Optional: if your flux spans several orders of magnitude
     labs(
       title = "GOES X-ray flux time series: SC24 and SC25",
       y = "W/m²",
       x = NULL
     ) +
     theme_minimal(base_size = 12)
  
   ggsave(paste0(dir,"Xray_flux_time_series.pdf"), plot = plt, width = 5.03, height = 3.98)
   
 }

