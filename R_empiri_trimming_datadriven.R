setwd("C:\\Users\\yanni\\OneDrive\\Documenten\\Universiteit\\EOR\\Year 3\\Thesis") 

library(stats)
library(Matrix)
library(pracma)

#source('global.R', echo=FALSE)
rm(list=ls())
data <- read.csv("SPF_Final.csv")
test_set <- 16

add_data <- function(df, date_range) {
  for (i in 1:nrow(df)) {
    time <- date_range[df[i, 'TIME_PERIOD']]
    df[i, 'Quarter'] <- time
  }
  return(df)
}

date_range <- seq(as.Date("1999-10-01"), as.Date("2018-09-30"), by = "quarter")
df <- add_data(data, date_range)
### Now we have proper df with quarterly info

proper_df <- function(df, name, horizon) {
  df_name <- df[df$FCT_TOPIC == name,]
  df <- df_name[df_name$FCT_HORIZON == horizon,]
  return(df)
}

make_pred_error_matrix <- function(df){
  Merror <- matrix(0, nrow=length(unique(df$FCT_SOURCE)), ncol=max(df$TIME_PERIOD))
  forecasters <- unique(df$FCT_SOURCE)
  for (i in 1:length(forecasters)) {
    df_1 <- df[df$FCT_SOURCE == forecasters[i],]
    for (j in 1:nrow(df_1)) {
      t <- df_1[j,5]  
      Merror[i,t] <- df_1[j,8] ## ERR_VALUE 
      if (df_1[j,8] == 0) {
        Merror[i,t] <- 0.0001 ## ERR_VALUE  so we can still see if forecaster participated
      }
    }
  }
  return(t(Merror))
}

trim <- function(weight, c, spec){
  if (spec == 1){
    weight[weight < c] <- c
    alpha <- 1/sum(weight)
    return(alpha*weight)
  }
  if (spec == 2){
    idw <- which(weight > c)
    count <- length(weight[weight < c])
    weight[weight < c] <- c
    alpha <- (1 - count * c)/sum(weight[idw])
    weight[idw] <- weight[idw] * alpha
    return(weight)
  }
  if (spec == 3){
    minw <- min(weight)
    idw <- which(weight > c)
    temp <- 0
    if (length(idw)<length(weight)){
      temp <- sum(weight[weight < c] * (c/minw))  #afterwards more difficult to find these values
      weight[weight < c] <- weight[weight < c] * (c/minw)
    }
    alpha <- (1 - temp)/sum(weight[idw])
    weight[idw] <- weight[idw] * alpha
    return(weight)
  }
}

fminvar_optim <- function(w, Mcov){
  return(t(w)%*% Mcov %*% w)
}    # combination variance 

trim4 <- function(Mcov,c){
  fminvar <- function(w) { fminvar_optim(w, Mcov) }
  n <- nrow(Mcov)
  Aeq <- matrix(1,1,n)                            # equality constraint : weights sum up to 1
  beq <- 1
  w0 <- runif(n,-5,5)                             # Initial random weights
  lb <- c 
  opt <- fmincon(w0,fminvar,Aeq=Aeq,beq=beq,lb=lb)  # optimal solution with norm1 constraint
  w_opt <- matrix(opt$par,nrow = 1,ncol=n)
  return(w_opt)
}

trim5 <- function(Mcov,c){
  fminvar <- function(w) { fminvar_optim(w, Mcov) }
  n <- nrow(Mcov)
  Aeq <- matrix(1,1,n)                            # equality constraint : weights sum up to 1
  beq <- 1
  w0 <- runif(n,-5,5)                           # Initial random weights                                  # TR5
  c1 <-  1 + 2*c                                # inequality constraint
  nonlcon <- function(x) {
    c <- norm(x,type="O") - c1        
    return(list(c=c,ceq=NULL))}     
  w_opt <- solnl(w0,fminvar,Aeq=Aeq,Beq=beq, confun = nonlcon)
  return(matrix(w_opt$par,nrow = 1,ncol=n))
}

est_cov <- function(Merror, t) {
  Nforc <- ncol(Merror)
  T <- nrow(Merror)
  Mcov <- matrix(0, nrow = Nforc, ncol = Nforc)
  for (i in 1:Nforc) {
    for (j in 1:Nforc) {
      v1 <- Merror[1:(T-t), i]
      v2 <- Merror[1:(T-t), j]
      
      nonzero_indices_i <- which(v1 != 0)
      nonzero_indices_j <- which(v2 != 0)
      common_forecasts <- intersect(nonzero_indices_i, nonzero_indices_j)
      if (length(common_forecasts) != 0) {
        selected_elements_i <- v1[common_forecasts]
        selected_elements_j <- v2[common_forecasts]
        Mcov[i, j] <- sum(selected_elements_i * selected_elements_j) / length(common_forecasts)
      }
    }
  }
  Mcov2 <- nearPD(Mcov)  
  Mcov2 <- as.matrix(Mcov2$mat)
  nonzero_forecaster <- which(Merror[T-t+1,] != 0)
  Mcov3 <- Mcov2[nonzero_forecaster, nonzero_forecaster]
  return(list(Mcov3,nonzero_forecaster))
}

calc_trimmed_weights <- function(weights,c,trimvalue,Mcov){
  if (trimvalue %in% c(1, 2, 3)) {
    trimmed_w <- trim(weights,c,trimvalue) 
  } else if(trimvalue == 4) {
    trimmed_w <- trim4(Mcov,c)
  }else if(trimvalue == 5) {
    success <- FALSE # catch any errors
    n_tries <- 0
    while (!success && n_tries < 10) {
      n_tries <- n_tries + 1
      tryCatch({
        trimmed_w <- trim5(Mcov,-c)  # sometimes not a pd matrix
        success <- TRUE
      }, error = function(e) {
      })
    }
    if (!success) {
      trimmed_w <- NaN
      print('failure weight 5')
    }
  }
  return(trimmed_w)
}

calc_ratio <- function(df, name, horizon, test_set, trimvalue,values_c){
  df1 <- proper_df(df, name, horizon)
  Merror <- make_pred_error_matrix(df1)
  T <- nrow(Merror)
  
  dymMSE <- numeric(test_set)
  eqqMSE <- 0
  trimfcterror  <- array(0,c(2,test_set))
  cvector  <- numeric(test_set)
  wvector <- numeric(0)
  for (t in test_set:1){
    Results <- est_cov(Merror,t)
    Mcov <- Results[[1]]
    nonzero_forecaster <- Results[[2]]
    iota <- matrix(1,nrow(Mcov),1)
    weights1 <- t(iota)%*%solve(Mcov) * as.numeric((1/(t(iota)%*%solve(Mcov)%*%iota)))
    eqweights1 <- (1/length(weights1))*matrix(1,nrow(weights1),ncol(weights1))
    
    eqqMSE <- eqqMSE + (eqweights1%*%Merror[(T-t+1),nonzero_forecaster])^2
    trimfcterror[2,t] <- eqweights1%*%Merror[(T-t+1),nonzero_forecaster]
    MSE <- numeric(length(values_c))
    for (i in 1:length(values_c))
    {
      tr_w <- calc_trimmed_weights(weights1,values_c[i], trimvalue,Mcov) 
      for (tim in 1:(T-t-1))
      { 
        MSE[i] = MSE[i] + (sum(tr_w*t(Merror[tim,nonzero_forecaster]), na.rm=TRUE))^2 
      }
    } 
    min_mse <- which(MSE==min(MSE))
    tr_w <- calc_trimmed_weights(weights1,values_c[min_mse[1]], trimvalue,Mcov) 
    dymMSE[t] <- (tr_w%*%Merror[(T-t+1),nonzero_forecaster])^2
    trimfcterror[1,t] <- tr_w%*%Merror[(T-t+1),nonzero_forecaster]
    cvector[t] <- values_c[min_mse[1]]
    wvector <- c(wvector,tr_w)
  }
  ratio <- sum(dymMSE)/eqqMSE
  return(list(ratio,trimfcterror,cvector,wvector))
}

calc_fterror <- function(trimfcterror){
  e1 = trimfcterror[1,]
  e2 = trimfcterror[2,]
  dm = dm.test(e1,e2,alternative = c("two.sided"))
  fterror = dm$p.value
  return(fterror)
}

calc_all_ratio <- function(df, name, horizon, test_set,values_c,trimvalue){
  tot_ratio <- numeric(length(trimvalue))
  fterror_pvalue <- numeric(0)
  wvector <- numeric(0)
  cvector <- numeric(0)
  errorvector <- numeric(0)
  for (i in 1:length(trimvalue)) {
    results <- calc_ratio(df, name, horizon, test_set, trimvalue[i],values_c)
    tot_ratio[i] <- results[[1]]
    fterror <-calc_fterror(results[[2]])
    cvector <- cbind(cvector,results[[3]])
    wvector <- cbind(wvector,results[[4]])
    Verror <- results[[2]][1, ]
    errorvector <- rbind(errorvector,Verror)
  }
  Verror_eq <- results[[2]][2, ]
  errorvector <- rbind(errorvector,Verror_eq)
  return(list(tot_ratio,fterror_pvalue,cvector,wvector,errorvector))
}

test_set <- 16
values_c <- seq(-2,0,0.1)
trimvalue <- seq(from = 1, to = 5, by = 1)

results1 <- calc_all_ratio(df, 'HICP', 1, test_set,values_c,trimvalue)
dynamic_ratio_HIPC_1 <- results1[[1]]
p_val_dr_HIPC_1 <- results1[[2]]
vError1 <- results1[[5]]
rownames(vError1) <- c("Trim 1","Trim 2","Trim 3","Trim 4", "Trim 5", "Trim eq")
write.csv(vError1, file <- "dynamic_vError_HICP1.csv")



results2 <- calc_all_ratio(df, 'HICP', 2, test_set,values_c,trimvalue)
dynamic_ratio_HIPC_2 <- results2[[1]]
p_val_dr_HIPC_2 <- results2[[2]]
vError2 <- results2[[5]]
rownames(vError2) <- c("Trim 1","Trim 2","Trim 3","Trim 4", "Trim 5", "Trim eq")
write.csv(vError2, file <- "dynamic_vError_HICP2.csv")

results3 <- calc_all_ratio(df, 'RGDP', 1, test_set,values_c,trimvalue)
dynamic_ratio_RGDP_1 <- results3[[1]]
p_val_dr_RGDP_1 <- results3[[2]]
vError3 <- results3[[5]]
rownames(vError3) <- c("Trim 1","Trim 2","Trim 3","Trim 4", "Trim 5", "Trim eq")
write.csv(vError3, file <- "dynamic_vError_RGDP1.csv")

results4 <- calc_all_ratio(df, 'RGDP', 2, test_set,values_c,trimvalue)
dynamic_ratio_RGDP_2 <- results4[[1]]
p_val_dr_RGDP_2 <- results4[[2]]
vError4 <- results4[[5]]
rownames(vError4) <- c("Trim 1","Trim 2","Trim 3","Trim 4", "Trim 5", "Trim eq")
write.csv(vError4, file <- "dynamic_vError_RGDP2.csv")

results5 <- calc_all_ratio(df, 'UNEM', 1, test_set,values_c,trimvalue)
dynamic_ratio_UNEM_1 <- results5[[1]]
p_val_dr_UNEM_1 <- results5[[2]]
vError5 <- results5[[5]]
rownames(vError5) <- c("Trim 1","Trim 2","Trim 3","Trim 4", "Trim 5", "Trim eq")
write.csv(vError5, file <- "dynamic_vError_UNEM1.csv")

results6 <- calc_all_ratio(df, 'UNEM', 2, test_set,values_c,trimvalue)
dynamic_ratio_UNEM_2 <- results6[[1]]
p_val_dr_UNEM_2 <- results6[[2]]
vError6 <- results6[[5]]
rownames(vError6) <- c("Trim 1","Trim 2","Trim 3","Trim 4", "Trim 5", "Trim eq")
write.csv(vError6, file <- "dynamic_vError_UNEM6.csv")

write.csv(dynamic_ratio_HIPC_1, file <- "dynamic_ratio_0.2_HIPC_1.csv")
write.csv(dynamic_ratio_HIPC_2, file <- "dynamic_ratio_0.2_HIPC_2.csv")
write.csv(dynamic_ratio_RGDP_1, file <- "dynamic_ratio_0.2_RGDP_1.csv")
write.csv(dynamic_ratio_RGDP_2, file <- "dynamic_ratio_0.2_RGDP_2.csv")
write.csv(dynamic_ratio_UNEM_1, file <- "dynamic_ratio_0.2_UNEM_1.csv")
write.csv(dynamic_ratio_UNEM_2, file <- "dynamic_ratio_0.2_UNEM_2.csv")

write.csv(p_val_dr_HIPC_1, file <- "p_val_dr_0.2_HIPC_1.csv")
write.csv(p_val_dr_HIPC_2, file <- "p_val_dr_0.2_HIPC_2.csv")
write.csv(p_val_dr_RGDP_1, file <- "p_val_dr_0.2_RGDP_1.csv")
write.csv(p_val_dr_RGDP_2, file <- "p_val_dr_0.2_RGDP_2.csv")
write.csv(p_val_dr_UNEM_1, file <- "p_val_dr_0.2_UNEM_1.csv")
write.csv(p_val_dr_UNEM_2, file <- "p_val_dr_0.2_UNEM_2.csv")

msfe_HIPC_1_prev <- read.csv("Results_dynamic//dynamic_ratio_HIPC_1.csv", header = TRUE)
msfe_HIPC_2_prev <- read.csv("Results_dynamic//dynamic_ratio_HIPC_2.csv", header = TRUE)
msfe_RGDP_1_prev <- read.csv("Results_dynamic//dynamic_ratio_RGDP_1.csv", header = TRUE)
msfe_RGDP_2_prev <- read.csv("Results_dynamic//dynamic_ratio_RGDP_2.csv", header = TRUE)
msfe_UNEM_1_prev <- read.csv("Results_dynamic//dynamic_ratio_UNEM_1.csv", header = TRUE)
msfe_UNEM_2_prev <- read.csv("Results_dynamic..dynamic_ratio_UNEM_2.csv", header = TRUE)



