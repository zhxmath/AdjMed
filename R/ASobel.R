
ASobel <- function(X, M, Y, Z, Delta, Model, tau){
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
    T_alpha <- abs(alpha_est/alpha_SE)
    T_beta <- abs(beta_est/beta_SE)
    T_sobel_1 <- alpha_est*beta_est
    T_sobel_2 <- alpha_est^2*(beta_SE^2) + beta_est^2*(alpha_SE^2)
    q_sobel <- abs(qnorm(tau/(2),0,1))
    
    CI_sobel <- matrix(0,p,2)
    CI_sobel[,1] <- alpha_est*beta_est - q_sobel*sqrt(T_sobel_2)
    CI_sobel[,2] <- alpha_est*beta_est + q_sobel*sqrt(T_sobel_2)
    
    q_Asobel <- abs(qnorm(tau/(2),0,1/2))
    CI_Asobel <- matrix(0,p,2)
    CI_Asobel[,1] <- alpha_est*beta_est - q_Asobel*sqrt(T_sobel_2)
    CI_Asobel[,2] <- alpha_est*beta_est + q_Asobel*sqrt(T_sobel_2)
    
    for (j in 1:p){
      if (max(T_alpha[j],T_beta[j]) > cutoff){
        
        CI_Asobel[j,1] <- CI_sobel[j,1]
        CI_Asobel[j,2] <- CI_sobel[j,2]
        
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
    T_alpha <- abs(alpha_est/alpha_SE)
    T_beta <- abs(beta_est/beta_SE)
    T_sobel_1 <- alpha_est*beta_est
    T_sobel_2 <- alpha_est^2*(beta_SE^2) + beta_est^2*(alpha_SE^2)
    q_sobel <- abs(qnorm(tau/(2),0,1))
    
    CI_sobel <- matrix(0,p,2)
    CI_sobel[,1] <- alpha_est*beta_est - q_sobel*sqrt(T_sobel_2)
    CI_sobel[,2] <- alpha_est*beta_est + q_sobel*sqrt(T_sobel_2)
    
    q_Asobel <- abs(qnorm(tau/(2),0,1/2))
    CI_Asobel <- matrix(0,p,2)
    CI_Asobel[,1] <- alpha_est*beta_est - q_Asobel*sqrt(T_sobel_2)
    CI_Asobel[,2] <- alpha_est*beta_est + q_Asobel*sqrt(T_sobel_2)
    
    for (j in 1:p){
      if (max(T_alpha[j],T_beta[j]) > cutoff){
        
        CI_Asobel[j,1] <- CI_sobel[j,1]
        CI_Asobel[j,2] <- CI_sobel[j,2]
        
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
    T_alpha <- abs(alpha_est/alpha_SE)
    T_beta <- abs(beta_est/beta_SE)
    T_sobel_1 <- alpha_est*beta_est
    T_sobel_2 <- alpha_est^2*(beta_SE^2) + beta_est^2*(alpha_SE^2)
    q_sobel <- abs(qnorm(tau/(2),0,1))
    
    CI_sobel <- matrix(0,p,2)
    CI_sobel[,1] <- alpha_est*beta_est - q_sobel*sqrt(T_sobel_2)
    CI_sobel[,2] <- alpha_est*beta_est + q_sobel*sqrt(T_sobel_2)
    
    q_Asobel <- abs(qnorm(tau/(2),0,1/2))
    CI_Asobel <- matrix(0,p,2)
    CI_Asobel[,1] <- alpha_est*beta_est - q_Asobel*sqrt(T_sobel_2)
    CI_Asobel[,2] <- alpha_est*beta_est + q_Asobel*sqrt(T_sobel_2)
    
    for (j in 1:p){
      if (max(T_alpha[j],T_beta[j]) > cutoff){
        
        CI_Asobel[j,1] <- CI_sobel[j,1]
        CI_Asobel[j,2] <- CI_sobel[j,2]
        
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
    T_alpha <- abs(alpha_est/alpha_SE)
    T_beta <- abs(beta_est/beta_SE)
    T_sobel_1 <- alpha_est*beta_est
    T_sobel_2 <- alpha_est^2*(beta_SE^2) + beta_est^2*(alpha_SE^2)
    q_sobel <- abs(qnorm(tau/(2),0,1))
    
    CI_sobel <- matrix(0,p,2)
    CI_sobel[,1] <- alpha_est*beta_est - q_sobel*sqrt(T_sobel_2)
    CI_sobel[,2] <- alpha_est*beta_est + q_sobel*sqrt(T_sobel_2)
    
    q_Asobel <- abs(qnorm(tau/(2),0,1/2))
    CI_Asobel <- matrix(0,p,2)
    CI_Asobel[,1] <- alpha_est*beta_est - q_Asobel*sqrt(T_sobel_2)
    CI_Asobel[,2] <- alpha_est*beta_est + q_Asobel*sqrt(T_sobel_2)
    
    for (j in 1:p){
      if (max(T_alpha[j],T_beta[j]) > cutoff){
        
        CI_Asobel[j,1] <- CI_sobel[j,1]
        CI_Asobel[j,2] <- CI_sobel[j,2]
        
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
    T_alpha <- abs(alpha_est/alpha_SE)
    T_beta <- abs(beta_est/beta_SE)
    T_sobel_1 <- alpha_est*beta_est
    T_sobel_2 <- alpha_est^2*(beta_SE^2) + beta_est^2*(alpha_SE^2)
    q_sobel <- abs(qnorm(tau/(2),0,1))
    
    CI_sobel <- matrix(0,p,2)
    CI_sobel[,1] <- alpha_est*beta_est - q_sobel*sqrt(T_sobel_2)
    CI_sobel[,2] <- alpha_est*beta_est + q_sobel*sqrt(T_sobel_2)
    
    q_Asobel <- abs(qnorm(tau/(2),0,1/2))
    CI_Asobel <- matrix(0,p,2)
    CI_Asobel[,1] <- alpha_est*beta_est - q_Asobel*sqrt(T_sobel_2)
    CI_Asobel[,2] <- alpha_est*beta_est + q_Asobel*sqrt(T_sobel_2)
    
    for (j in 1:p){
      if (max(T_alpha[j],T_beta[j]) > cutoff){
        
        CI_Asobel[j,1] <- CI_sobel[j,1]
        CI_Asobel[j,2] <- CI_sobel[j,2]
        
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
    T_alpha <- abs(alpha_est/alpha_SE)
    T_beta <- abs(beta_est/beta_SE)
    T_sobel_1 <- alpha_est*beta_est
    T_sobel_2 <- alpha_est^2*(beta_SE^2) + beta_est^2*(alpha_SE^2)
    q_sobel <- abs(qnorm(tau/(2),0,1))
    
    CI_sobel <- matrix(0,p,2)
    CI_sobel[,1] <- alpha_est*beta_est - q_sobel*sqrt(T_sobel_2)
    CI_sobel[,2] <- alpha_est*beta_est + q_sobel*sqrt(T_sobel_2)
    
    q_Asobel <- abs(qnorm(tau/(2),0,1/2))
    CI_Asobel <- matrix(0,p,2)
    CI_Asobel[,1] <- alpha_est*beta_est - q_Asobel*sqrt(T_sobel_2)
    CI_Asobel[,2] <- alpha_est*beta_est + q_Asobel*sqrt(T_sobel_2)
    
    for (j in 1:p){
      if (max(T_alpha[j],T_beta[j]) > cutoff){
        
        CI_Asobel[j,1] <- CI_sobel[j,1]
        CI_Asobel[j,2] <- CI_sobel[j,2]
        
      }
    }
}
}# 

OUT <- list(alpha_est=as.numeric(alpha_est), alpha_SE=as.numeric(alpha_SE), beta_est=as.numeric(beta_est), 
            beta_SE=as.numeric(beta_SE),CI_Asobel=CI_Asobel)
return(OUT)

} 



