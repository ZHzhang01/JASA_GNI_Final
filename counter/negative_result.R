library(igraph)
library(tidyverse)
library(stringdist)
if (!requireNamespace("Matrix", quietly = TRUE)) {
  install.packages("Matrix")
}
library(Matrix)

library(quadprog)


args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
  script_dir <- dirname(script_path)
  repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
} else {
  repo_root <- normalizePath(".", mustWork = TRUE)
}

n <- 1000

pscore1 <-  runif(n, 0.1, 0.9)
pscore0 <- 1 - pscore1

C <- 0.5 * rnorm(n) 
Y_1 <- pscore1/pscore0 * C 
Y_0 <- -C

errors <- runif(n, 1, 5)
errors <- rnorm(n) 

X <- C 
X <- scale(X,center = T, scale = F)

data1 <- read.table(file.path(repo_root, "counter", "linear_parameter_X_counter.txt"), header = FALSE)
data2 <- read.table(file.path(repo_root, "counter", "linear_parameter_Y0_counter.txt"), header = FALSE)
data3 <- read.table(file.path(repo_root, "counter", "linear_parameter_Y1_counter.txt"), header = FALSE)
data4 <- read.table(file.path(repo_root, "counter", "linear_parameter_noise_counter.txt"), header = FALSE)
data5 <- read.table(file.path(repo_root, "counter", "linear_parameter_ps1_counter.txt"), header = FALSE)
data6 <- read.table(file.path(repo_root, "counter", "linear_parameter_ps0_counter.txt"), header = FALSE)

X <- as.matrix(data1[, -1])  
Y_0 <- as.matrix(data2[, -1])  
Y_1 <- as.matrix(data3[, -1])  
errors <- as.vector(data4[, -1])   
pscore1 <-  as.vector(data5[, -1])   
pscore0 <-  as.vector(data6[, -1])   

# This counterexample constructs observed outcomes from the two potential outcomes and additive noise.
get_Y <- function(Z){
  Y <- Z * Y_1 + (1-Z) * Y_0 + 0.5 * errors
  return(Y)
}

tau <- map_dbl(1:100000, ~{
  Z <- rbinom(n, size = 1, prob = pscore1); Y <- get_Y(Z)
  T_vec <- Z; D <- Y*T_vec/pscore1-Y*(1-T_vec)/pscore0
  return(mean(D))
}) %>% mean()

tau

tau_haj <- map_dbl(1:100000, ~{
  Z <- rbinom(n, size = 1, prob = pscore1); Y <- get_Y(Z)
  T_vec <- Z; D <- Y*T_vec/(pscore1*mean(T_vec/pscore1))-Y*(1-T_vec)/(pscore0*mean((1-T_vec)/pscore0))
  return(mean(D))
}) %>% mean()

tau_haj

sum <- c()

sim_res<- map_dfr(1:100000, ~{
  # A is identity here, so the design is treated as independent across units.
  A <- diag(n)
  sum <- c(sum,1); print("sum:"); print(length(sum))
  Z <- rbinom(n, size = 1, prob = pscore1); Y <- get_Y(Z); X_aug <- X
  T_vec <- (Z)
  w <- T_vec/pscore1-(1-T_vec)/pscore0
  w_haj <- T_vec/(pscore1*mean(T_vec/pscore1))-(1-T_vec)/(pscore0*mean((1-T_vec)/pscore0))
  w_1 <- T_vec/pscore1
  w_0 <- (1-T_vec)/pscore0
  w_haj_1 <- T_vec/(pscore1*mean(T_vec/pscore1))
  w_haj_0 <- (1-T_vec)/(pscore0*mean((1-T_vec)/pscore0))
  D <- Y*w
  D_haj <- Y*w_haj
  # D is the inverse-probability score; the quadratic forms below deliver the variance analogs under A = I.
  Leung <- mean(D); 
  var_Leung <- t(D-Leung)%*%A%*%(D-Leung)/n^2 %>% as.vector();
  var_Leung_naive <- t(D-Leung)%*%diag(n)%*%(D-Leung)/n^2 %>% as.vector()
  coverage_Leung <- abs(Leung-tau)<=qnorm(0.975)*sqrt(var_Leung)
  coverage_Leung_naive <- abs(Leung-tau)<=qnorm(0.975)*sqrt(var_Leung_naive)
  X_db <- X
  D_2 <- (X_db * w)
  # Solve the A-weighted normal equations for the variance-reducing adjustment coefficient.
  hbeta_2_haj <- solve(t(D_2)%*%A%*%(D_2),   t(D_2)%*%A%*%(D- (mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) ))
  Ours_X_haj <- mean((Y-(X_db) %*% hbeta_2_haj)*w_haj)
  var_Ours_X_haj <- t(D- D_2%*%hbeta_2_haj  -(mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) ) %*%A %*% (D- D_2%*%hbeta_2_haj  -(mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) )/n^2 %>% as.vector()
  var_Ours_X_haj_naive <- t(D- D_2%*%hbeta_2_haj  -(mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) ) %*%diag(n) %*% (D- D_2%*%hbeta_2_haj  -(mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) )/n^2 %>% as.vector()
  coverage_Ours_X_haj <- abs(Ours_X_haj - tau)<=qnorm(0.975)*sqrt(var_Ours_X_haj)
  coverage_Ours_X_haj_naive <- abs(Ours_X_haj - tau)<=qnorm(0.975)*sqrt(var_Ours_X_haj_naive)
  X_db <- X
  X_db <- cbind(X_db * T_vec, X_db * (1-T_vec))
  D_2 <- (X_db * w)
  V <- (D- (mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0)        )
  hbeta_2_haj <- solve(t(D_2)%*%A%*%(D_2),   t(D_2)%*%A%*%V)
  Ours_X_haj_lin <- mean((Y-X_db%*%hbeta_2_haj)*w_haj)
  var_Ours_X_haj_lin <- t(V - D_2%*%hbeta_2_haj  ) %*%A %*% (V - D_2%*%hbeta_2_haj )/n^2 %>% as.vector()
  var_Ours_X_haj_lin_naive <- t(V - D_2%*%hbeta_2_haj  ) %*%diag(n) %*% (V - D_2%*%hbeta_2_haj )/n^2 %>% as.vector()
  coverage_Ours_X_haj_lin <- abs(Ours_X_haj_lin - tau)<=qnorm(0.975)*sqrt(var_Ours_X_haj_lin)
  coverage_Ours_X_haj_lin_naive <- abs(Ours_X_haj_lin - tau)<=qnorm(0.975)*sqrt(var_Ours_X_haj_lin_naive)
  X_db <- X
  X_db <- cbind(X_db * T_vec, X_db * (1-T_vec))
  X_db <- cbind(X * pscore1, X * pscore0)
  D_2 <- (X_db * w)
  V <- (D- (mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0)        )
  hbeta_2_haj <- solve(t(D_2)%*%A%*%(D_2),   t(D_2)%*%A%*%V)
  Ours_X_haj_lin_phi <- mean((Y-X_db%*%hbeta_2_haj)*w_haj)
  var_Ours_X_haj_lin_phi <- t(V - D_2%*%hbeta_2_haj  ) %*%A %*% (V - D_2%*%hbeta_2_haj )/n^2 %>% as.vector()
  var_Ours_X_haj_lin_naive_phi <- t(V - D_2%*%hbeta_2_haj  ) %*%diag(n) %*% (V - D_2%*%hbeta_2_haj )/n^2 %>% as.vector()
  coverage_Ours_X_haj_lin_phi <- abs(Ours_X_haj_lin_phi - tau)<=qnorm(0.975)*sqrt(var_Ours_X_haj_lin_phi)
  coverage_Ours_X_haj_lin_naive_phi <- abs(Ours_X_haj_lin_phi - tau)<=qnorm(0.975)*sqrt(var_Ours_X_haj_lin_naive_phi)
  lm_haj <- lm(Y~1+T_vec+X+T_vec:X, w = T_vec/pscore1+(1-T_vec)/pscore0)
  e_haj <- lm_haj %>% resid()
  Gao_L <- lm_haj %>% coef() %>% .[2];
  X_db <- X
  X_db <- cbind(X_db * (1-T_vec), X_db * (T_vec))
  hbeta_2_haj <- c(lm_haj %>% coef() %>% .[3],   lm_haj %>% coef() %>% .[3] + lm_haj %>% coef() %>% .[4]   )
  w = T_vec/pscore1 - (1-T_vec)/pscore0
  D_2 <- X_db * w
  var_Gao_L <- t(D- D_2%*%hbeta_2_haj  - (mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) ) %*%A %*% (D- D_2%*%hbeta_2_haj  -(mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) )/n^2 %>% as.vector()
  var_Gao_L_naive <- t(D- D_2%*%hbeta_2_haj  - (mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) ) %*%diag(n) %*% (D- D_2%*%hbeta_2_haj  -(mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) )/n^2 %>% as.vector()
  coverage_Gao_L <- abs(Gao_L-tau)<=qnorm(0.975)*sqrt(var_Gao_L)
  coverage_Gao_L_naive <- abs(Gao_L-tau)<=qnorm(0.975)*sqrt(var_Gao_L_naive)
  lm_haj <- lm(Y~1+T_vec+X, w = T_vec/pscore1+(1-T_vec)/pscore0)
  e_haj <- lm_haj %>% resid()
  Gao_F <- lm_haj %>% coef() %>% .[2]
  X_db <- X
  w <- T_vec/pscore1 - (1-T_vec)/pscore0
  D_2 <- X_db * w
  hbeta_2_haj <- (lm_haj %>% coef() %>% .[3]  )
  var_Gao_F <- t(D- D_2 * hbeta_2_haj  - (mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) ) %*%A %*% (D- D_2 * hbeta_2_haj  -(mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) )/n^2 %>% as.vector()
  var_Gao_F_naive <- t(D- D_2 * hbeta_2_haj  - (mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) ) %*%diag(n) %*% (D- D_2 * hbeta_2_haj  -(mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) )/n^2 %>% as.vector()
  coverage_Gao_F <- abs(Gao_F-tau)<=qnorm(0.975)*sqrt(var_Gao_F)
  coverage_Gao_F_naive <- abs(Gao_F-tau)<=qnorm(0.975)*sqrt(var_Gao_F_naive)
  lm_haj <- lm(Y~1+T_vec, w = T_vec/pscore1+(1-T_vec)/pscore0)
  w <- T_vec/pscore1+(1-T_vec)/pscore0
  e_haj <- lm_haj %>% resid(); Gao <- lm_haj %>% coef() %>% .[2]
  C_haj <- cbind(1,T_vec)
  Gao_plus <- Gao
  var_Gao <- t(D  - (mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) ) %*%A %*% (D  -(mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) )/n^2 %>% as.vector()
  var_Gao_naive <- t(D - (mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) ) %*%diag(n) %*% (D  -(mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) )/n^2 %>% as.vector()
  coverage_Gao <- abs(Gao-tau)<=qnorm(0.975)*sqrt(var_Gao)
  coverage_Gao_naive <- abs(Gao-tau)<=qnorm(0.975)*sqrt(var_Gao_naive)
  return(tibble(Leung, Gao, Gao_F, Gao_L,  Ours_X_haj,   Ours_X_haj_lin,               Ours_X_haj_lin_phi,
                var_Leung, var_Gao, var_Gao_F, var_Gao_L,     var_Ours_X_haj,    var_Ours_X_haj_lin,  var_Ours_X_haj_lin_phi,  
                coverage_Leung, coverage_Gao, coverage_Gao_F, coverage_Gao_L,  coverage_Ours_X_haj,     coverage_Ours_X_haj_lin, coverage_Ours_X_haj_lin_phi, 
                var_Leung_naive, var_Gao_naive, var_Gao_F_naive, var_Gao_L_naive,    var_Ours_X_haj_naive,   var_Ours_X_haj_lin_naive, var_Ours_X_haj_lin_naive_phi,
                coverage_Leung_naive, coverage_Gao_naive, coverage_Gao_F_naive, coverage_Gao_L_naive,  coverage_Ours_X_haj_naive,     coverage_Ours_X_haj_lin_naive, coverage_Ours_X_haj_lin_naive_phi))
})

index <- 7

print("tau:"); print(tau)

print("Estimation:"); print((sim_res %>% summarise_all(mean)  %>% as.data.frame() )[1:index]) 
print("oracle sd:"); print((sim_res %>% summarise_all(sd)  %>% as.data.frame() )[1:index]) 
print("practical sd:"); print(sqrt(sim_res %>% summarise_all(mean)  %>% as.data.frame() )[(index+1):(index*2)]) 
print("naive sd:"); print(sqrt(sim_res %>% summarise_all(mean)  %>% as.data.frame() )[(3*index+1):(4*index)]) 
print("practical coverage:"); print((sim_res %>% summarise_all(mean)  %>% as.data.frame() )[(2*index+1):(3*index)]) 
print("naive coverage:"); print((sim_res %>% summarise_all(mean)  %>% as.data.frame() )[(4*index+1):(5*index)])

oracle <- c()
oracle_plus <- c()

sim_res_oracle_cover <- map_dfr(1:1, ~{
  variables <- c('Leung', 'Gao', 'Gao_F', 'Gao_L',   'Ours_X_haj',    'Ours_X_haj_lin', 'Ours_X_haj_lin_phi')
  for (var in variables) {
    o_coverage <- mean(abs(sim_res[[var]] - tau) <= qnorm(0.975) * sd(sim_res[[var]]) )
    cat("Oracle Coverage for", var, ":", o_coverage, "\n")
    oracle <- c(oracle, o_coverage)
  }
  return(tibble(oracle, oracle_plus))
})

file_path <- file.path(repo_root, "counter", "result_synthetic.txt")
write.table(tau, file = file_path, col.names = FALSE)
write.table((sim_res %>% summarise_all(mean)  %>% as.data.frame() )[1:index], file = file_path, col.names = FALSE, append = TRUE)
write.table((sim_res %>% summarise_all(sd)  %>% as.data.frame() )[1:index], file = file_path, col.names = FALSE, append = TRUE)
write.table(sqrt(sim_res %>% summarise_all(mean)  %>% as.data.frame() )[(index+1):(index*2)], file = file_path, col.names = FALSE, append = TRUE)
write.table(sqrt(sim_res %>% summarise_all(mean)  %>% as.data.frame() )[(3*index+1):(4*index)], file = file_path, col.names = FALSE, append = TRUE)
write.table((sim_res_oracle_cover$oracle)  , file = file_path, col.names = FALSE, append = TRUE)
write.table((sim_res %>% summarise_all(mean)  %>% as.data.frame() )[(2*index+1):(3*index)], file = file_path, col.names = FALSE, append = TRUE)
write.table((sim_res %>% summarise_all(mean)  %>% as.data.frame() )[(4*index+1):(5*index)], file = file_path, col.names = FALSE, append = TRUE)

cat("Data has been successfully written to the text file (bn = 3):", file_path, "\n")

tempt <- mean(pscore0 * pscore1 * (Y_1/pscore1 - mean(Y_1)/pscore1 + Y_0/pscore0 - mean(Y_0)/pscore0 )^2)
sigma_haj <- sqrt(tempt / n)
tempt <- mean(pscore0 * pscore1 * (Y_1/pscore1 - 0/pscore1 + Y_0/pscore0 - 0/pscore0 )^2)
sigma_ht <- sqrt(tempt / n)
beta_L_1 <- sum(X * X)^{-1} * sum(X * Y_1)
beta_L_0 <- sum(X * X)^{-1} * sum(X * Y_0)
beta_F <- (beta_L_0+beta_L_1)/2
tempt <- mean(pscore0 * pscore1 * (Y_1/pscore1 - mean(Y_1)/pscore1 + Y_0/pscore0 - mean(Y_0)/pscore0 - X * beta_L_1/pscore1 - X * beta_L_0/pscore0 )^2)
sigma_L <- sqrt(tempt / n)
tempt <- mean(pscore0 * pscore1 * (Y_1/pscore1 - mean(Y_1)/pscore1 + Y_0/pscore0 - mean(Y_0)/pscore0 - X * beta_F/(pscore1*pscore0)  )^2)
sigma_F <- sqrt(tempt / n)
tempt1 <- sum(pscore0/pscore1 * X * X)
tempt2 <- sum(X * X)
tempt3 <- tempt1 <- sum(pscore1/pscore0 * X * X)
tempt_left <- matrix(c(tempt1, tempt2, tempt2, tempt3), ncol = 2)
tempt1 <- sum(pscore0/pscore1 * X * Y_1 + X * Y_0)
tempt2 <- sum(pscore1/pscore0 * X * Y_0 + X * Y_1)
tilde_beta_1 <- (solve(tempt_left) %*% (matrix(cbind(tempt1, tempt2))))[1,] 
tilde_beta_0 <- (solve(tempt_left) %*% (matrix(cbind(tempt1, tempt2))))[2,] 
delta <- mean(pscore1 * pscore0 * (X * tilde_beta_1/pscore1 + X * tilde_beta_0/pscore0)^2)
delta_1 <- mean(pscore1 * pscore0 * (X * (tilde_beta_1-beta_L_1)/pscore1 + X * (tilde_beta_0- beta_L_0)/pscore0)^2)
sighaj_sigL <- delta - delta_1
sighaj_sigL / (sigma_haj^2 * n)
sqrt(abs(sighaj_sigL)/n)/sigma_haj 
delta_2 <- mean(pscore1 * pscore0 * (X * (tilde_beta_1-beta_F)/pscore1 + X * (tilde_beta_0- beta_F)/pscore0)^2)
sighaj_sigF <- delta - delta_2
sighaj_sigF / (sigma_haj^2 * n)
beta_F <- (beta_L_0+beta_L_1)/2


tempt <- mean(pscore0 * pscore1 * (Y_1/pscore1 - mean(Y_1)/pscore1 + Y_0/pscore0 - mean(Y_0)/pscore0 - X * beta_L_1/pscore1 - X * beta_L_0/pscore0 )^2)


sigma_L <- sqrt(tempt / n)


tempt <- mean(pscore0 * pscore1 * (Y_1/pscore1 - mean(Y_1)/pscore1 + Y_0/pscore0 - mean(Y_0)/pscore0 - X * beta_F/(pscore1*pscore0)  )^2)


sigma_F <- sqrt(tempt / n)


tempt1 <- sum(pscore0/pscore1 * X * X)
tempt2 <- sum(X * X)
tempt3 <- tempt1 <- sum(pscore1/pscore0 * X * X)
tempt_left <- matrix(c(tempt1, tempt2, tempt2, tempt3), ncol = 2)


tempt1 <- sum(pscore0/pscore1 * X * Y_1 + X * Y_0)
tempt2 <- sum(pscore1/pscore0 * X * Y_0 + X * Y_1)




tilde_beta_1 <- (solve(tempt_left) %*% (matrix(cbind(tempt1, tempt2))))[1,] 
tilde_beta_0 <- (solve(tempt_left) %*% (matrix(cbind(tempt1, tempt2))))[2,] 


delta <- mean(pscore1 * pscore0 * (X * tilde_beta_1/pscore1 + X * tilde_beta_0/pscore0)^2)


delta_1 <- mean(pscore1 * pscore0 * (X * (tilde_beta_1-beta_L_1)/pscore1 + X * (tilde_beta_0- beta_L_0)/pscore0)^2)


sighaj_sigL <- delta - delta_1


sighaj_sigL / (sigma_haj^2 * n)
sqrt(abs(sighaj_sigL)/n)/sigma_haj 


delta_2 <- mean(pscore1 * pscore0 * (X * (tilde_beta_1-beta_F)/pscore1 + X * (tilde_beta_0- beta_F)/pscore0)^2)


sighaj_sigF <- delta - delta_2


sighaj_sigF / (sigma_haj^2 * n)









