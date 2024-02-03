
AJS <- function(X, M, Y, Z, Delta, Model){
library("survival")
n <- dim(M)[1]
p <- dim(M)[2]
weight <- 1
cutoff <- weight*sqrt(n)/log(n)
if (Model=="Linear"){
  if (Z[1]=="null"){
    alpha_est <- matrix(0,1,p)
    alpha_SE <- matrix(0,1,p)
    MX <- cbind(M, X)
    fit_Y <- lsfit(MX,Y,intercept = TRUE) # M-Y
    beta_est <- matrix(coef(fit_Y))[2:(p+1)]
    beta_SE <- ls.diag(fit_Y)$std.err[2:(p+1)]
    for (k in 1:p){
      fit_M <- lsfit(X,M[,k],intercept = TRUE) # x-M
      alpha_est[k] <- matrix(coef(fit_M))[2]
      alpha_SE[k] <- ls.diag(fit_M)$std.err[2]
    }
    P_alpha <- 2*(1-pnorm(abs(alpha_est)/alpha_SE,0,1))
    P_beta  <- 2*(1-pnorm(abs(beta_est)/beta_SE,0,1))
    P_ab <- matrix(0,2,p)
    P_ab[1,] <- P_alpha
    P_ab[2,] <- P_beta
    P_max <- apply(P_ab,2,max)
    T_alpha <- abs(alpha_est/alpha_SE)
    T_beta <- abs(beta_est/beta_SE)
    P_AJS <- matrix(0,1,p)
    for (j in 1:p){
      if (T_alpha[j] < cutoff&&T_beta[j]< cutoff){
        P_AJS[j] <- (P_max[j])^2
      }else{
        P_AJS[j] <- P_max[j]
      }
    }
  }else{
    alpha_est <- matrix(0,1,p)
    alpha_SE <- matrix(0,1,p)
    MXZ <- cbind(M,X,Z)
    fit_Y <- lsfit(MXZ,Y,intercept = TRUE) # M-Y
    beta_est <- matrix(coef(fit_Y))[2:(p+1)]
    beta_SE <- ls.diag(fit_Y)$std.err[2:(p+1)]
    XZ <- cbind(X,Z)
    for (k in 1:p){
      fit_M <- lsfit(XZ,M[,k],intercept = TRUE) # x-M
      alpha_est[k] <- matrix(coef(fit_M))[2]
      alpha_SE[k] <- ls.diag(fit_M)$std.err[2]
    }
    P_alpha <- 2*(1-pnorm(abs(alpha_est)/alpha_SE,0,1))
    P_beta  <- 2*(1-pnorm(abs(beta_est)/beta_SE,0,1))
    P_ab <- matrix(0,2,p)
    P_ab[1,] <- P_alpha
    P_ab[2,] <- P_beta
    P_max <- apply(P_ab,2,max)
    T_alpha <- abs(alpha_est/alpha_SE)
    T_beta <- abs(beta_est/beta_SE)
    P_AJS <- matrix(0,1,p)
    for (j in 1:p){
      if (T_alpha[j] < cutoff&&T_beta[j]< cutoff){
        P_AJS[j] <- (P_max[j])^2
      }else{
        P_AJS[j] <- P_max[j]
      }
    }
  }
} else if (Model=="Logistic"){
  if (Z[1]=="null"){
    alpha_est <- matrix(0,p,1)
    alpha_SE <- matrix(0,p,1)
    MX <- cbind(M, X)
    random_data <- data.frame(Y,MX)
    fit_Y <- glm(formula = Y~MX,data = random_data, family = binomial(link = "logit"))
    beta_est <- matrix(summary(fit_Y)$coefficients[2:(1+p),1])
    beta_SE <-  matrix(summary(fit_Y)$coefficients[2:(1+p),2])

    for (k in 1:p){
      fit_M <- lsfit(X,M[,k],intercept = TRUE) # x-M
      alpha_est[k] <- matrix(matrix(coef(fit_M))[2])
      alpha_SE[k] <-  matrix(ls.diag(fit_M)$std.err[2])
    }
    P_alpha <- 2*(1-pnorm(abs(alpha_est)/alpha_SE,0,1))
    P_beta  <- 2*(1-pnorm(abs(beta_est)/beta_SE,0,1))
    P_ab <- matrix(0,2,p)
    P_ab[1,] <- P_alpha
    P_ab[2,] <- P_beta
    P_max <- apply(P_ab,2,max)
    T_alpha <- abs(alpha_est/alpha_SE)
    T_beta <- abs(beta_est/beta_SE)
    P_AJS <- matrix(0,1,p)
    for (j in 1:p){
      if (T_alpha[j] < cutoff&&T_beta[j]< cutoff){
        P_AJS[j] <- (P_max[j])^2
      }else{
        P_AJS[j] <- P_max[j]
      }
    }
  }else{
    alpha_est <- matrix(0,p,1)
    alpha_SE <- matrix(0,p,1)
    MXZ <- cbind(M, X, Z)
    random_data <- data.frame(Y,MXZ)
    fit_Y <- glm(formula = Y~MXZ,data = random_data, family = binomial(link = "logit"))
    beta_est <- matrix(summary(fit_Y)$coefficients[2:(1+p),1])
    beta_SE <-  matrix(summary(fit_Y)$coefficients[2:(1+p),2])
    XZ <- cbind(X, Z)
    for (k in 1:p){
      fit_M <- lsfit(XZ,M[,k],intercept = TRUE) # x-M
      alpha_est[k] <- matrix(matrix(coef(fit_M))[2])
      alpha_SE[k] <-  matrix(ls.diag(fit_M)$std.err[2])
    }
    P_alpha <- 2*(1-pnorm(abs(alpha_est)/alpha_SE,0,1))
    P_beta  <- 2*(1-pnorm(abs(beta_est)/beta_SE,0,1))
    P_ab <- matrix(0,2,p)
    P_ab[1,] <- P_alpha
    P_ab[2,] <- P_beta
    P_max <- apply(P_ab,2,max)
    T_alpha <- abs(alpha_est/alpha_SE)
    T_beta <- abs(beta_est/beta_SE)
    P_AJS <- matrix(0,1,p)
    for (j in 1:p){
      if (T_alpha[j] < cutoff&&T_beta[j]< cutoff){
        P_AJS[j] <- (P_max[j])^2
      }else{
        P_AJS[j] <- P_max[j]
      }
    }
  }
}else if (Model=="Cox"){
  if (Z[1]=="null"){
    OT= Y
    status <- Delta > 0
    y = Surv(OT, status)
    alpha_est <- matrix(0,p,1)
    alpha_SE <- matrix(0,p,1)
    MX <- cbind(M, X)
    fit_Y <-  coxph(y ~ MX)
    beta_est <- matrix(summary(fit_Y)$coefficients[1:p,1])
    beta_SE <-  matrix(summary(fit_Y)$coefficients[1:p,3])
    for (k in 1:p){
      fit_M <- lsfit(X,M[,k],intercept = TRUE) # x-M
      alpha_est[k] <- matrix(matrix(coef(fit_M))[2])
      alpha_SE[k] <-  matrix(ls.diag(fit_M)$std.err[2])
    }
    P_alpha <- 2*(1-pnorm(abs(alpha_est)/alpha_SE,0,1))
    P_beta  <- 2*(1-pnorm(abs(beta_est)/beta_SE,0,1))
    P_ab <- matrix(0,2,p)
    P_ab[1,] <- P_alpha
    P_ab[2,] <- P_beta
    P_max <- apply(P_ab,2,max)
    T_alpha <- abs(alpha_est/alpha_SE)
    T_beta <- abs(beta_est/beta_SE)
    P_AJS <- matrix(0,1,p)
    for (j in 1:p){
      if (T_alpha[j] < cutoff&&T_beta[j]< cutoff){
        P_AJS[j] <- (P_max[j])^2
      }else{
        P_AJS[j] <- P_max[j]
      }
    }
  }else{
    OT= Y
    status <- Delta > 0
    y = Surv(OT, status)
    alpha_est <- matrix(0,p,1)
    alpha_SE <- matrix(0,p,1)
    MXZ <- cbind(M, X,Z)
    fit_Y <-  coxph(y ~ MXZ)
    beta_est <- matrix(summary(fit_Y)$coefficients[1:p,1])
    beta_SE <-  matrix(summary(fit_Y)$coefficients[1:p,3])
    XZ <- cbind(X,Z)
    for (k in 1:p){
      fit_M <- lsfit(XZ,M[,k],intercept = TRUE) # x-M
      alpha_est[k] <- matrix(matrix(coef(fit_M))[2])
      alpha_SE[k] <-  matrix(ls.diag(fit_M)$std.err[2])
    }
    P_alpha <- 2*(1-pnorm(abs(alpha_est)/alpha_SE,0,1))
    P_beta  <- 2*(1-pnorm(abs(beta_est)/beta_SE,0,1))
    P_ab <- matrix(0,2,p)
    P_ab[1,] <- P_alpha
    P_ab[2,] <- P_beta
    P_max <- apply(P_ab,2,max)
    T_alpha <- abs(alpha_est/alpha_SE)
    T_beta <- abs(beta_est/beta_SE)
    P_AJS <- matrix(0,1,p)
    for (j in 1:p){
      if (T_alpha[j] < cutoff&&T_beta[j]< cutoff){
        P_AJS[j] <- (P_max[j])^2
      }else{
        P_AJS[j] <- P_max[j]
      }
    }
  }
}# 

OUT <- list(alpha_est=as.numeric(alpha_est), alpha_SE=as.numeric(alpha_SE), beta_est=as.numeric(beta_est), 
            beta_SE=as.numeric(beta_SE),P_AJS=as.numeric(P_AJS))
return(OUT)

} 



