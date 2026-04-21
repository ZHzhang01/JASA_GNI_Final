library(igraph)
library(tidyverse)
library(stringdist)


args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
  script_dir <- dirname(script_path)
  repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
} else {
  repo_root <- normalizePath(".", mustWork = TRUE)
}

n <- 3000
r1 <- 0.5 

ball_vol <- function(d,r){
  volume <- pi^(d/2) * r^d / gamma(d/2+1)
  return(volume)
}

positions <- matrix(runif(2 * n), nrow = n, ncol = 2)
dist_matrix <- as.matrix(dist(positions, method = "euclidean"))
r <- (1.5 / (ball_vol(2, 1) * n))^{1/2}
E <- (dist_matrix <= r)
E <- (E) * 1
for (i in 1:nrow(E)){
  E[i,i] <- 0
}

g <- graph_from_adjacency_matrix(as.matrix(E), mode = "undirected")
avg_path_length <- mean_distance(g)








connected_edges <- which(E != 0, arr.ind = TRUE)

nodes <- unique(c(connected_edges[, 1], connected_edges[, 2]))

E_new <- matrix(0, nrow = length(nodes), ncol = length(nodes))

for (i in 1:nrow(connected_edges)) {
  row <- connected_edges[i, 1]
  col <- connected_edges[i, 2]
  E_new[match(row, nodes), match(col, nodes)] <- 1
  E_new[match(col, nodes), match(row, nodes)] <- 1
}



print("Symmetric adjacency matrix E_new after removing edges with no connected nodes:")
print(E_new)

E <- E_new





positions <- matrix(runif(2 * ncol(E)), nrow = ncol(E), ncol = 2)
num_nb <- rowSums(E)
G <- E/rowSums(E)
n <- ncol(G)


errors <- rnorm(n) + (positions[,1] - 0.5)
X <-rnorm(n) %>% scale(scale = FALSE); epsilon <- errors







data1 <- read.table(file.path(repo_root, "syn", "linear_parameter_X.txt"), header = FALSE)
data2 <- read.table(file.path(repo_root, "syn", "linear_parameter_G.txt"), header = FALSE)
data3 <- read.table(file.path(repo_root, "syn", "linear_parameter_noise.txt"), header = FALSE)

X <- as.matrix(data1[, -1])
G <- as.matrix(data2[, -1])
E <- apply(G, c(1, 2), function(x) ifelse(x != 0, 1, 0))
num_nb <- rowSums(E)
errors <- as.vector(data3[, -1])







n <- length(X)
r1 <- 0.5








G <- as.matrix(G)

# Iterate the network fixed-point map until outcomes stabilize under assignment Z.
get_Y <- function(Z){
  
  alpha1 <- 0.1
  alpha2 <- 1
  alpha3 <- 1
  alpha4 <- 1
  alpha5 <- 1
  
  Y_0 <- (-1 + alpha2 * G%*%Z +  alpha3 * Z + alpha4 * X + alpha5 * errors); Y_0 <- as.numeric(Y_0)
  
  Y_next <- (-1 +  alpha1 * (G%*% (Y_0))  +  alpha2 * G%*%Z +  alpha3 * Z + alpha4 * X + alpha5 * errors); Y_next <- as.numeric(Y_next)
  
  while( max( abs(as.numeric(Y_0) - as.numeric(Y_next)) ) > 0.0001){
    Y_0 <- Y_next; Y_0 <- as.numeric(Y_0)
    Y_next <- (-1 +  alpha1 * (G%*%Y_0)  +  alpha2 * G%*%Z +  alpha3 * Z + alpha4 * X + alpha5 * errors); Y_next <- as.numeric(Y_next)
  }
  return(Y_0)
}

get_Y_heter <- function(Z){
  
  term1 <- solve(diag(length(X)) - mat1 %*% G) 
  term2 <-  (-diag(length(X)) %*% rep(1, length(X)) + (mat2 %*% G + mat3) %*% Z + mat4 %*% X + errors)
  Y <- term1 %*% term2 
  
  return(Y)
}







get_Y_nonlinear <- function(Z){
  Y_0 <- (-1 + 1* G%*%Z + Z + X + errors) > 0
  Y_next <- Y_0 + 1
  while( max( abs(Y_0 - Y_next) ) > 0.0001){
    Y_0 <- Y_next
    Y_next <- (-1 + 1.5* G%*%Y_0  + 1* G%*%Z + Z + X + errors) >0
  }
  return(Y_0)
}

threshold <- 0.5 * num_nb

# Map the individual assignment Z to the exposure indicator induced by neighbors.
get_T <- function(Z){
  return(drop(  floor(E %*% Z) > threshold ) *1)
}

pscore0 <- pbinom(threshold ,size = num_nb, prob = r1); pscore1 <- 1-pscore0

g <- graph_from_adjacency_matrix(E, mode = "undirected")
shortest_paths_mat <- shortest.paths(g, mode = "all")


if (avg_path_length < 2*log(n)/log(sum(E)/n)){
  bn <- avg_path_length/2
}
if (avg_path_length >= 2*log(n)/log(sum(E)/n)){
  bn <- avg_path_length^{1/3}
}

K<-1
bn<- round(max(bn, 2*K)) 
# Scripts in this family are nearly identical except for bn, so users can modify bn here instead of maintaining separate copies.
bn <- 4
# A records which pairs are within graph distance bn and is the dependence matrix used in variance calculations.
A <- (shortest_paths_mat <= bn); temp <- eigen(A);  A_p <- (temp$vectors)%*%diag((temp$values)*(temp$values >= 0))%*%solve(temp$vectors) 

# Build the augmented regressors (Z, GZ, X, GX) used for orthogonalization.
get_X <- function(X,Z,G){
  return(matrix(c(Z,drop(G%*%Z),X,drop(G%*%X)), nrow=n, ncol = 4))
}


pscore0 <- as.numeric(pscore0)
pscore1 <- as.numeric(pscore1)

# Monte Carlo moments estimate projection coefficients that residualize augmentation terms.
mom_mat <- matrix(0, nrow = n, ncol = 5)
mom_mat_haj <- matrix(0, nrow = n, ncol = 5)
mom_mat_lin <- matrix(0, nrow = n, ncol = 5+4)
mom_mat_haj_lin <- matrix(0, nrow = n, ncol = 5+4)
for(i in 1:10000){
  Z <- rbinom(n, size = 1, prob = r1);  
  X_aug <- get_X(X,Z,G)
  T_vec <- get_T(Z);
  X_aug_lin <- cbind(X_aug * T_vec, X_aug * (1-T_vec))
  
  w <- T_vec/pscore1-(1-T_vec)/pscore0
  w_haj <- T_vec/(pscore1*mean(T_vec/pscore1))-(1-T_vec)/(pscore0*mean((1-T_vec)/pscore0))
  mom_mat <- mom_mat + c(w^2, X_aug*w)
  mom_mat_haj <- mom_mat_haj + c(w_haj^2, X_aug*w_haj)
  mom_mat_lin <- mom_mat_lin + c(w^2, X_aug_lin*w) 
  mom_mat_haj_lin <- mom_mat_haj_lin + c(w_haj^2, X_aug_lin*w_haj) 
}
orth_coef <- mom_mat[, 2:5] / mom_mat[, 1]
orth_coef_haj <- mom_mat_haj[, 2:5] / mom_mat_haj[, 1]
orth_coef_lin <- mom_mat_lin[, 2:(5+4)] / mom_mat_lin[, 1]
orth_coef_haj_lin <- mom_mat_haj_lin[, 2:(5+4)] / mom_mat_haj_lin[, 1]










paste("orth_coef:", orth_coef)

tau <- map_dbl(1:10000, ~{
  Z <- rbinom(n, size = 1, prob = r1); Y <- get_Y(Z)
  T_vec <- get_T(Z); D <- Y*T_vec/pscore1-Y*(1-T_vec)/pscore0
  return(mean(D))
}) %>% mean()

tau_haj <- map_dbl(1:10000, ~{
  Z <- rbinom(n, size = 1, prob = r1); Y <- get_Y(Z)
  T_vec <- get_T(Z); D <- Y*T_vec/(pscore1*mean(T_vec/pscore1))-Y*(1-T_vec)/(pscore0*mean((1-T_vec)/pscore0))
  return(mean(D))
}) %>% mean()

sum <- c()

sim_res<- map_dfr(1:10000, ~{
  
  sum <- c(sum,1); print("sum:"); print(length(sum))
  Z <- rbinom(n, size = 1, prob = r1); Y <- get_Y(Z); X_aug <- get_X(X,Z,G)
  T_vec <- get_T(Z)
  w <- T_vec/pscore1-(1-T_vec)/pscore0
  w_haj <- T_vec/(pscore1*mean(T_vec/pscore1))-(1-T_vec)/(pscore0*mean((1-T_vec)/pscore0))
  w_1 <- T_vec/pscore1
  w_0 <- (1-T_vec)/pscore0
  w_haj_1 <- T_vec/(pscore1*mean(T_vec/pscore1))
  w_haj_0 <- (1-T_vec)/(pscore0*mean((1-T_vec)/pscore0))
  
  # D is the inverse-probability score; Leung uses it directly with network dependence matrix A.
  D <- Y*w
  D_haj <- Y*w_haj
  
  Leung <- mean(D); 
  var_Leung <- t(D-Leung)%*%A%*%(D-Leung)/n^2 %>% as.vector();
  var_Leung_naive <- t(D-Leung)%*%diag(n)%*%(D-Leung)/n^2 %>% as.vector()
  coverage_Leung <- abs(Leung-tau)<=qnorm(0.975)*sqrt(var_Leung)
  coverage_Leung_naive <- abs(Leung-tau)<=qnorm(0.975)*sqrt(var_Leung_naive)
  
  
  
  
  
  
  X_db <- X_aug - w_haj * (orth_coef_haj)
  D_2 <- (X_db * w)
  
  D_2 <- D_2 - (as.matrix(w_1) %*%  t(as.matrix(colMeans(X_db*w_haj_1)) ) - as.matrix(w_0) %*%  t(as.matrix(colMeans(X_db*w_haj_0)) ))
  
  # Solve the A-weighted normal equations for the variance-reducing adjustment coefficient.
  V <- (D- (mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0)        )
  hbeta_2_haj <- solve(t(D_2)%*%A%*%(D_2),   t(D_2)%*%A%*%V)
  
  Ours_G_haj <- mean((Y-X_db%*%hbeta_2_haj)*w_haj)
  var_Ours_G_haj <- t(V - D_2%*%hbeta_2_haj) %*%A %*% (V- D_2%*%hbeta_2_haj)/n^2 %>% as.vector()
  var_Ours_G_haj_naive <- t(V - D_2%*%hbeta_2_haj) %*% diag(n) %*% (V- D_2%*%hbeta_2_haj)/n^2 %>% as.vector()
  coverage_Ours_G_haj <- abs(Ours_G_haj - tau)<=qnorm(0.975)*sqrt(var_Ours_G_haj)
  coverage_Ours_G_haj_naive <- abs(Ours_G_haj - tau)<=qnorm(0.975)*sqrt(var_Ours_G_haj_naive)
  
  
  
  
  X_db <- X
  D_2 <- (X_db * w)
  
  
  
  hbeta_2_haj <- solve(t(D_2)%*%A%*%(D_2),   t(D_2)%*%A%*%(D- (mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) ))
  Ours_X_haj <- mean((Y-X_db%*%hbeta_2_haj)*w_haj)
  var_Ours_X_haj <- t(D- D_2%*%hbeta_2_haj  -(mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) ) %*%A %*% (D- D_2%*%hbeta_2_haj  -(mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) )/n^2 %>% as.vector()
  var_Ours_X_haj_naive <- t(D- D_2%*%hbeta_2_haj  -(mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) ) %*%diag(n) %*% (D- D_2%*%hbeta_2_haj  -(mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) )/n^2 %>% as.vector()
  coverage_Ours_X_haj <- abs(Ours_X_haj - tau)<=qnorm(0.975)*sqrt(var_Ours_X_haj)
  coverage_Ours_X_haj_naive <- abs(Ours_X_haj - tau)<=qnorm(0.975)*sqrt(var_Ours_X_haj_naive)
  
  
  
  
  
  
  X_db <- cbind(X_aug * T_vec, X_aug * (1-T_vec))
  X_db <- X_db - w_haj * (orth_coef_haj_lin)
  D_2 <- (X_db * w)
  
  
  D_2 <- D_2 - (as.matrix(w_1) %*%  t(as.matrix(colMeans(X_db*w_haj_1)) ) - as.matrix(w_0) %*%  t(as.matrix(colMeans(X_db*w_haj_0)) ))
  
  
  
  D <- Y * w
  
  V <- (D- (mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0)        )
  hbeta_2_haj <- solve(t(D_2)%*%A%*%(D_2),   t(D_2)%*%A%*%V)
  
  Ours_G_haj_lin <- mean((Y-X_db%*%hbeta_2_haj)*w_haj)
  var_Ours_G_haj_lin <- t(V- D_2%*%hbeta_2_haj) %*%A %*% (V- D_2%*%hbeta_2_haj )/n^2 %>% as.vector()
  var_Ours_G_haj_lin_naive <- t(V- D_2%*%hbeta_2_haj) %*% diag(n) %*% (V- D_2%*%hbeta_2_haj )/n^2 %>% as.vector()
  coverage_Ours_G_haj_lin <- abs(Ours_G_haj_lin - tau)<=qnorm(0.975)*sqrt(var_Ours_G_haj_lin)
  coverage_Ours_G_haj_lin_naive <- abs(Ours_G_haj_lin - tau)<=qnorm(0.975)*sqrt(var_Ours_G_haj_lin_naive)
  
  
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  X_db <- cbind(X_aug * T_vec, X_aug * (1-T_vec))
  X_db <- X_db - w_haj * (orth_coef_haj_lin)
  lm_haj <- lm(Y~1+T_vec+X_db, w = T_vec/pscore1+(1-T_vec)/pscore0)
  e_haj <- lm_haj %>% resid()
  Gao_L_abla <- lm_haj %>% coef() %>% .[2];
  
  hbeta_2_haj <- c(lm_haj %>% coef() %>% .[3],  lm_haj %>% coef() %>% .[4], lm_haj %>% coef() %>% .[5], lm_haj %>% coef() %>% .[6], lm_haj %>% coef() %>% .[7], lm_haj %>% coef() %>% .[8], lm_haj %>% coef() %>% .[9], lm_haj %>% coef() %>% .[10])
  w = T_vec/pscore1 - (1-T_vec)/pscore0
  D_2 <- X_db * w
  D_2 <- D_2 - (as.matrix(w_1) %*%  t(as.matrix(colMeans(X_db*w_haj_1)) ) - as.matrix(w_0) %*%  t(as.matrix(colMeans(X_db*w_haj_0)) ))
  var_Gao_L_abla <- t(D- D_2%*%hbeta_2_haj  - (mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) ) %*%A %*% (D- D_2%*%hbeta_2_haj  -(mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) )/n^2 %>% as.vector()
  var_Gao_L_naive_abla <- t(D- D_2%*%hbeta_2_haj  - (mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) ) %*%diag(n) %*% (D- D_2%*%hbeta_2_haj  -(mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) )/n^2 %>% as.vector()
  coverage_Gao_L_abla <- abs(Gao_L_abla-tau)<=qnorm(0.975)*sqrt(var_Gao_L_abla)
  coverage_Gao_L_naive_abla <- abs(Gao_L_abla-tau)<=qnorm(0.975)*sqrt(var_Gao_L_naive_abla)
  
  
  
  X_db <- X_aug - w_haj * (orth_coef_haj)
  D_2 <- (X_db * w)
  D_2 <- D_2 - (as.matrix(w_1) %*%  t(as.matrix(colMeans(X_db*w_haj_1)) ) - as.matrix(w_0) %*%  t(as.matrix(colMeans(X_db*w_haj_0)) ))
  
  lm_haj <- lm(Y~1+T_vec+X_db, w = T_vec/pscore1+(1-T_vec)/pscore0)
  e_haj <- lm_haj %>% resid()
  Gao_F_abla <- lm_haj %>% coef() %>% .[2]
  w <- T_vec/pscore1 - (1-T_vec)/pscore0
  hbeta_2_haj <- c(lm_haj %>% coef() %>% .[3], lm_haj %>% coef() %>% .[4], lm_haj %>% coef() %>% .[5],lm_haj %>% coef() %>% .[6]  )
  var_Gao_F_abla <- t(D- D_2 %*% hbeta_2_haj  - (mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) ) %*%A %*% (D- D_2 %*% hbeta_2_haj  -(mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) )/n^2 %>% as.vector()
  var_Gao_F_naive_abla <- t(D- D_2 %*% hbeta_2_haj  - (mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) ) %*%diag(n) %*% (D- D_2 %*% hbeta_2_haj  -(mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0) )/n^2 %>% as.vector()
  coverage_Gao_F_abla <- abs(Gao_F_abla-tau)<=qnorm(0.975)*sqrt(var_Gao_F_abla)
  coverage_Gao_F_naive_abla <- abs(Gao_F_abla-tau)<=qnorm(0.975)*sqrt(var_Gao_F_naive_abla)
  
  
  
  X_db <- X_aug
  D_2 <- (X_db * w)
  
  D_2 <- D_2 - (as.matrix(w_1) %*%  t(as.matrix(colMeans(X_db*w_haj_1)) ) - as.matrix(w_0) %*%  t(as.matrix(colMeans(X_db*w_haj_0)) ))
  
  V <- (D- (mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0)        )
  hbeta_2_haj <- solve(t(D_2)%*%A%*%(D_2),   t(D_2)%*%A%*%V)
  Ours_G_haj_abla <- mean((Y-X_db%*%hbeta_2_haj)*w_haj)
  var_Ours_G_haj_abla <- t(V - D_2%*%hbeta_2_haj) %*%A %*% (V- D_2%*%hbeta_2_haj)/n^2 %>% as.vector()
  var_Ours_G_haj_naive_abla <- t(V - D_2%*%hbeta_2_haj) %*% diag(n) %*% (V- D_2%*%hbeta_2_haj)/n^2 %>% as.vector()
  coverage_Ours_G_haj_abla <- abs(Ours_G_haj_abla - tau)<=qnorm(0.975)*sqrt(var_Ours_G_haj_abla)
  coverage_Ours_G_haj_naive_abla <- abs(Ours_G_haj_abla - tau)<=qnorm(0.975)*sqrt(var_Ours_G_haj_naive_abla)
  
  
  X_db <- cbind(X_aug * T_vec, X_aug * (1-T_vec))
  D_2 <- (X_db * w)
  D_2 <- D_2 - (as.matrix(w_1) %*%  t(as.matrix(colMeans(X_db*w_haj_1)) ) - as.matrix(w_0) %*%  t(as.matrix(colMeans(X_db*w_haj_0)) ))
  D <- Y * w
  
  V <- (D- (mean(Y*w_haj_1)*w_1-mean(Y*w_haj_0)*w_0)        )
  hbeta_2_haj <- solve(t(D_2)%*%A%*%(D_2),   t(D_2)%*%A%*%V)
  
  Ours_G_haj_lin_abla <- mean((Y-X_db%*%hbeta_2_haj)*w_haj)
  var_Ours_G_haj_lin_abla <- t(V- D_2%*%hbeta_2_haj) %*%A %*% (V- D_2%*%hbeta_2_haj )/n^2 %>% as.vector()
  var_Ours_G_haj_lin_naive_abla <- t(V- D_2%*%hbeta_2_haj) %*% diag(n) %*% (V- D_2%*%hbeta_2_haj )/n^2 %>% as.vector()
  coverage_Ours_G_haj_lin_abla <- abs(Ours_G_haj_lin_abla - tau)<=qnorm(0.975)*sqrt(var_Ours_G_haj_lin_abla)
  coverage_Ours_G_haj_lin_naive_abla <- abs(Ours_G_haj_lin_abla - tau)<=qnorm(0.975)*sqrt(var_Ours_G_haj_lin_naive_abla)
  
  
  
  return(tibble(Leung, Gao, Gao_F, Gao_L, Gao_F_abla, Gao_L_abla,     Ours_X_haj, Ours_G_haj, Ours_G_haj_abla,   Ours_X_haj_lin,Ours_G_haj_lin,Ours_G_haj_lin_abla,  
                var_Leung, var_Gao, var_Gao_F, var_Gao_L, var_Gao_F_abla, var_Gao_L_abla,      var_Ours_X_haj, var_Ours_G_haj, var_Ours_G_haj_abla,    var_Ours_X_haj_lin,  var_Ours_G_haj_lin, var_Ours_G_haj_lin_abla,
                coverage_Leung, coverage_Gao, coverage_Gao_F, coverage_Gao_L, coverage_Gao_F_abla, coverage_Gao_L_abla,       coverage_Ours_X_haj, coverage_Ours_G_haj, coverage_Ours_G_haj_abla,     coverage_Ours_X_haj_lin, coverage_Ours_G_haj_lin,coverage_Ours_G_haj_lin_abla,
                var_Leung_naive, var_Gao_naive, var_Gao_F_naive, var_Gao_L_naive, var_Gao_F_naive_abla, var_Gao_L_naive_abla,       var_Ours_X_haj_naive, var_Ours_G_haj_naive,var_Ours_G_haj_naive_abla,   var_Ours_X_haj_lin_naive,  var_Ours_G_haj_lin_naive, var_Ours_G_haj_lin_naive_abla,
                coverage_Leung_naive, coverage_Gao_naive, coverage_Gao_F_naive, coverage_Gao_L_naive,coverage_Gao_F_naive_abla, coverage_Gao_L_naive_abla,      coverage_Ours_X_haj_naive, coverage_Ours_G_haj_naive, coverage_Ours_G_haj_naive_abla,     coverage_Ours_X_haj_lin_naive, coverage_Ours_G_haj_lin_naive, coverage_Ours_G_haj_lin_naive_abla,
                
                
                
                
  ))
  
})



index <- 12


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
  variables <- c('Leung', 'Gao', 'Gao_F', 'Gao_L', 'Gao_F_abla', 'Gao_L_abla',  'Ours_G_haj',  'Ours_X_haj','Ours_G_haj_abla',  'Ours_G_haj_lin',  'Ours_X_haj_lin', 'Ours_G_haj_lin_abla')
  
  for (var in variables) {
    o_coverage <- mean(abs(sim_res[[var]] - tau) <= qnorm(0.975) * sd(sim_res[[var]]) )
    cat("Oracle Coverage for", var, ":", o_coverage, "\n")
    oracle <- c(oracle, o_coverage)
  }
  
  
  
  return(tibble(oracle, oracle_plus))
  
})




file_path <- file.path(repo_root, "modify", "bn_4_result_synthetic_linear.txt")

write.table(tau, file = file_path, col.names = FALSE)

write.table((sim_res %>% summarise_all(mean)  %>% as.data.frame() )[1:index], file = file_path, col.names = FALSE, append = TRUE)
write.table((sim_res %>% summarise_all(sd)  %>% as.data.frame() )[1:index], file = file_path, col.names = FALSE, append = TRUE)
write.table(sqrt(sim_res %>% summarise_all(mean)  %>% as.data.frame() )[(index+1):(index*2)], file = file_path, col.names = FALSE, append = TRUE)
write.table(sqrt(sim_res %>% summarise_all(mean)  %>% as.data.frame() )[(3*index+1):(4*index)], file = file_path, col.names = FALSE, append = TRUE)
write.table((sim_res_oracle_cover$oracle)  , file = file_path, col.names = FALSE, append = TRUE)
write.table((sim_res %>% summarise_all(mean)  %>% as.data.frame() )[(2*index+1):(3*index)], file = file_path, col.names = FALSE, append = TRUE)
write.table((sim_res %>% summarise_all(mean)  %>% as.data.frame() )[(4*index+1):(5*index)], file = file_path, col.names = FALSE, append = TRUE)





cat("Data has been successfully written to the text file (bn = 3):", file_path, "\n")
