LD_eas <- read.csv("/Users/jingruimu/Desktop/COMP/Final Project/viprs_x_files/eas/LD.csv", header = FALSE, sep = ",")
X_test_eas <- read.csv("/Users/jingruimu/Desktop/COMP/Final Project/viprs_x_files/eas/X_test.csv", header = FALSE, sep = ",")
y_test_eas <- read.csv("/Users/jingruimu/Desktop/COMP/Final Project/viprs_x_files/eas/y_test.csv", header = FALSE, sep = ",")

X_train_eas <- read.csv("/Users/jingruimu/Desktop/COMP/Final Project/viprs_x_files/eas/X_train.csv", header = FALSE, sep = ",")
y_train_eas <- read.csv("/Users/jingruimu/Desktop/COMP/Final Project/viprs_x_files/eas/y_train.csv", header = FALSE, sep = ",")

beta_marginal_eas <- read.csv("/Users/jingruimu/Desktop/COMP/Final Project/viprs_x_files/eas/beta_marginal.csv", header = FALSE, sep = ",")
beta_marginal_eas <- as.vector(beta_marginal_eas)
beta_marginal_eas <- beta_marginal_eas$V1

beta_marginal_eur <- read.csv("/Users/jingruimu/Desktop/COMP/Final Project/viprs_x_files/eur/beta_marginal.csv", header = FALSE, sep = ",")
beta_marginal_eur <- as.vector(beta_marginal_eur)
beta_marginal_eur <- beta_marginal_eur$V1


LD_eur <- read.csv("/Users/jingruimu/Desktop/COMP/Final Project/viprs_x_files/eur/LD.csv", header = FALSE, sep = ",")
X_test_eur <- read.csv("/Users/jingruimu/Desktop/COMP/Final Project/viprs_x_files/eur/X_test.csv", header = FALSE, sep = ",")
y_test_eur <- read.csv("/Users/jingruimu/Desktop/COMP/Final Project/viprs_x_files/eur/y_test.csv", header = FALSE, sep = ",")

X_train_eur <- read.csv("/Users/jingruimu/Desktop/COMP/Final Project/viprs_x_files/eur/X_train.csv", header = FALSE, sep = ",")
y_train_eur <- read.csv("/Users/jingruimu/Desktop/COMP/Final Project/viprs_x_files/eur/y_train.csv", header = FALSE, sep = ",")


E_step <- function(X, Y, M, R, N, K, P, sigma_epsilon_current, pii_current,  
                   eta_current, mu_jk_current, Sigma_k_current, gamma_jk_current, beta_marginal){
  Sigma_jk_update <- vector("list", M)
  for(j in 1:M){
    Sigma_jk_update[[j]] <- array(data = NA, dim = c(P,P,K))
  }
  Sigma_jk_update_inv <- vector("list",M)
  for(j in 1:M){
    Sigma_jk_update_inv[[j]] <- array(data = NA, dim = c(P,P,K))
  }
  # eta_current <- eta_update
  mu_jk_update <- mu_jk_current
  # u_jk_update <- matrix(data = NA, nrow = M, ncol = K)
  gamma_jk_update <- gamma_jk_current
  eta_update <- eta_current
  temp_mu_jk <- rep(0,P)
  ujk_update <- matrix(data = NA, nrow = K, ncol = M)
  temp_gamma_jk_update <- matrix(data = NA, nrow = K, ncol = M)
  temp_Sigma_jk <- rep(0,P)
  temp_mu_jk <- matrix(data= NA, nrow = P, ncol = M)
  for(j in 1:M){
    for(k in 1:K){
      for(p in 1:P){
        temp_Sigma_jk[p] <- (t(X[[p]][,j])%*%X[[p]][,j])/sigma_epsilon_current[p]
      }
      Sigma_jk_update_inv[[j]][,,k] <- diag(temp_Sigma_jk,nrow = P, ncol = P) + ginv(Sigma_k_current[,,k])
      Sigma_jk_update[[j]][,,k] <- ginv(Sigma_jk_update_inv[[j]][,,k])
    }
  }
  for(j in 1:M){
    for(k in 1:K){
      for(p in 1:P){
        temp_mu_jk[p,j] <- (N[p] * (beta_marginal[j,p] -  t(R[-j,j,p])%*% eta_update[-j,p]))/sigma_epsilon_current[p]
      }
      mu_jk_update[,k,j] <- Sigma_jk_update[[j]][,,k] %*% temp_mu_jk[,j]
    }
  }
  for(j in 1:M){
    for(k in 1:K){
      ujk_update[k,j] <- log(pii_current[k]) + 0.5 * log(abs(det(Sigma_jk_update[[j]][,,k]))/abs(det(Sigma_k_current[,,k]))) + 0.5 * tr(Sigma_jk_update_inv[[j]][,,k] %*% mu_jk_update[,k,j] %*%t(mu_jk_update[,k,j]))
      # if(ujk_update[k,j] < log(0.01) | ujk_update[k,j] > log(0.99)){
      #   ujk_update[k,j] <- log(pii_current[k])
      # }
    }
  }
  for(j in 1:M){
    for(k in 1:K){
      temp_gamma_jk_update[k,j] <- exp(ujk_update[k,j]) / sum(exp(ujk_update[,j]))
      # gamma_jk_update[k,j] <- temp_gamma_jk_update[k,j]
      if(temp_gamma_jk_update[k,j] < 0.01){
        gamma_jk_update[k,j] <- 0.01
      }else if(temp_gamma_jk_update[k,j] > 0.99){
        gamma_jk_update[k,j] <- 0.99
      }else{
        gamma_jk_update[k,j] <- temp_gamma_jk_update[k,j]
      }
    }
  }
  for(j in 1:M){
    for(p in 1:P){
      eta_update[j,p] <- gamma_jk_update[,j] %*% mu_jk_update[p,,j]
    }
  }
  return(list(Sigma_jk = Sigma_jk_update, mu_jk = mu_jk_update, gamma_jk = gamma_jk_update, eta = eta_update, ujk = ujk_update))
}

M_step <- function(gamma_jk, mu_jk, Sigma_jk, eta, X, Y, M, N, K, P){
  temp22 <- array(data = NA, dim = c(P,P,M))
  temp33 <- array(data = NA, dim = c(P,K,M))
  Var_betajp <- matrix(data= NA, nrow = M, ncol = P)
  zeta <- matrix(data = NA, nrow = M, ncol = P)
  sigma_epsilon <- rep(0, P)
  pii <- rep(0, K)
  Sigma_k_update <- array(data = NA, dim = c(P,P,K))
  for(k in 1:K){
    for(j in 1:M){
      temp22[,,j] <- gamma_jk[k,j] * (mu_jk[,k,j] %*% t(mu_jk[,k,j]) + Sigma_jk[[j]][,,k])
    }
    Sigma_k_update[,,k] <- rowSums(temp22, dims = 2) / sum(gamma_jk[k,])
    pii[k] <- (1/M) * sum(gamma_jk[k,])
  }
  for(p in 1:P){
    for(j in 1:M){
      for(k in 1:K){
        temp33[p,k,j] <- mu_jk[p,k,j]^2 + Sigma_jk[[j]][p,p,k]^2
      }
      Var_betajp[j,p] <- t(gamma_jk[,j]) %*% temp33[p,,j]  - eta[j,p]^2
      zeta[j,p] <- t(X[[p]][,j]) %*% X[[p]][,j] %*% Var_betajp[j,p]
    }
    # zeta[j,p] <- t(X[,j,p]) %*% X[,j,p] %*% Var_betajp[j,p]
    sigma_epsilon[p] <- 1/N[p] * (t(Y[[p]] - X[[p]]%*%eta[,p])%*%(Y[[p]] - X[[p]]%*%eta[,p]) + sum(zeta[,p]))
  }
  return(list(Sigma_k = Sigma_k_update, sigma_epsilon = sigma_epsilon, pii = pii, zeta = zeta))
}

ELBO <- function(N, sigma_epsilon, Y, X, zeta, gamma_jk, pii, Sigma_jk, Sigma_k, mu_jk, eta, M, K){
  part1 <- rep(0,P)
  part2 <- rep(0,M)
  part3 <- rep(0,M)
  temp_part2 <- matrix(data = NA, nrow = K, ncol = M)
  temp_part3 <- matrix(data = NA, nrow = K, ncol = M)
  for(p in 1:P){
    part1[p] <- N[p] * log(2*pi*sigma_epsilon[p]) + 1/sigma_epsilon[p]*(t(Y[[p]] - X[[p]]%*%eta[,p])%*%(Y[[p]] - X[[p]]%*%eta[,p]) + sum(zeta[,p]))
  }
  for(j in 1:M){
    for(k in 1:K){
      temp_part2[k,j] <- log(gamma_jk[k,j] / pii[k])
    }
    part2[j] <- t(gamma_jk[,j]) %*% temp_part2[,j]
  }
  for(j in 1:M){
    for(k in 1:K){
      temp_part3[k,j] <- 0.5*gamma_jk[k,j]*(P + log(det(Sigma_jk[[j]][,,k])/det(Sigma_k[,,k])) - tr(ginv(Sigma_k[,,k]) %*% (mu_jk[,k,j] %*% t(mu_jk[,k,j]) + Sigma_jk[[j]][,,k])))
    }
    part3[j] <- sum(temp_part3[,j])
  }
  ELBO_value <- -0.5 * sum(part1) - sum(part2) + sum(part3)
  return(ELBO_value)
}

VIPRS_x <- function(maxiter, X, Y, R, M, N, P, K, Sigma_k_initial, Sigma_jk_initial, sigma_epsilon_initial, eta_initial, pii_initial, mu_jk_initial, u_jk_initial, gamma_jk_initial, zeta_initial){
  Sigma_jk_path <- vector("list", maxiter)
  for(iter in 1:maxiter){
    Sigma_jk_path[[iter]] <- vector("list", M)
    for(j in 1:M){
      Sigma_jk_path[[iter]][[j]] <- array(data = NA, dim = c(P,P,K-1))
    }
  }
  Sigma_k_path <- vector("list", maxiter)
  for(iter in 1:maxiter){
    Sigma_k_path[[iter]] <- array(data = NA, dim = c(P,P,K-1))
  }
  mu_jk_path <- vector("list", maxiter)
  for(iter in 1:maxiter){
    mu_jk_path[[iter]] <- array(data = NA, dim = c(P,K-1,M))
  }
  # u_jk_path <- array(data = NA, dim = c(K,M,maxiter))
  gamma_jk_path <- array(data = NA, dim = c(K-1, M, maxiter))
  eta_path <- array(data = NA, dim = c(M,P,maxiter))
  sigma_epsilon_path <- matrix(data = NA, nrow = P, ncol = maxiter)
  pii_path <- matrix(data = NA, nrow = K-1, ncol = maxiter)
  zeta_path <- array(data = NA, dim = c(M, P, maxiter))
  ELBO_path <- rep(0, maxiter)
  
  # Sigma_jk_path[[1]] <- Sigma_jk_initial
  Sigma_k_path[[1]] <- Sigma_k_initial
  mu_jk_path[[1]] <- mu_jk_initial
  # u_jk_path[,,1] <- u_jk_initial
  gamma_jk_path[,,1] <- gamma_jk_initial
  eta_path[,,1] <- eta_initial
  sigma_epsilon_path[,1] <- sigma_epsilon_initial
  pii_path[,1] <- pii_initial
  zeta_path[,,1] <- zeta_initial
  
  for(iter in 2:maxiter){
    E_list <- E_step(X = X, Y = Y, M = M, R = R, N = N, K = K, P = P, sigma_epsilon_current = sigma_epsilon_path[,iter-1], pii_current = pii_path[,iter-1],  
                     eta_current = eta_path[,,iter-1], mu_jk_current = mu_jk_path[[iter-1]], Sigma_k_current = Simga_k_path[[iter-1]], gamma_jk_current = gamma_jk_path[,,iter-1])
    Sigma_jk_path[[iter]] = E_list$Sigma_jk
    mu_jk_path[[iter]] = E_list$mu_jk
    gamma_jk_path[,,iter] = E_list$gamma_jk
    eta_path[,,iter] = E_list$eta
    M_list <- M_step(gamma_jk = gamma_jk_path[,,iter], mu_jk = mu_jk_path[[iter]], Sigma_jk = Sigma_jk_path[[iter]], 
                     eta = eta_path[,,iter], X = X, Y = Y, M = M, N = N, K = K, P = P)
    
    Sigma_k_path[[iter]] <- M_list$Sigma_k
    sigma_epsilon_path[,iter] <- M_list$sigma_epsilon
    pii_path[,iter] <- M_list$pii
    zeta_path[,,iter] <- M_list$zeta
    
    ELBO_path[iter] <- ELBO(N = N, sigma_epsilon = sigma_epsilon_path[,iter], Y = Y, X = X, zeta = zeta_path[,,iter], gamma_jk = gamma_jk_path[,,iter], pii = pii_path[,iter], Sigma_jk = Sigma_jk_path[[iter]], Sigma_k = Sigma_k[[iter]], mu_jk = mu_jk_path[[iter]], eta = eta_path[,,iter], M = M, K = K)
    
    if(abs(ELBO_path[iter] - ELBO_path[iter-1]) < 1e-4){
      numiter <- iter
      break
    }
  }
  
  return(list(ELBO_out = ELBO_path[2:numiter],
              eta_out = eta_path[,,1:numiter],
              Sigma_jk_out = Sigma_jk_path[[1:numiter]],
              mu_jk_out = mu_jk_path[[1:numiter]],
              gamma_jk_out = gamma_jk_path[,,1:numiter],
              sigma_epsilon_out = sigma_epsilon_path[,1:numiter],
              pii_out <- pii_path[,1:numiter],
              zeta_out <- zeta_path[,,1:numiter]
  ))
}

X_data <- vector("list",2)
scaled.X_train_eas <- scale(X_train_eas)
scaled.y_train_eas <- scale(y_train_eas)

scaled.X_train_eur <- scale(X_train_eur)
scaled.y_train_eur <- scale(y_train_eur)
X_data[[1]] <- scaled.X_train_eas
X_data[[2]] <- scaled.X_train_eur
Y_data <- vector("list", 2)
Y_data[[1]] <- scaled.y_train_eas
Y_data[[2]] <- scaled.y_train_eur
R_data <- array(data = NA, dim = c(100,100,2))
LD_eas <- as.matrix(LD_eas)
LD_eur <- as.matrix(LD_eur)
R_data[,,1] <- LD_eas
R_data[,,2] <- LD_eur
P <- 2
beta_marginal <- matrix(data = NA, nrow = M, ncol = P)
# beta_marginal[,1] <- beta_marginal_eas
# beta_marginal[,2] <- beta_marginal_eur
for(p in 1:P){
  for(j in 1:M){
    beta_marginal[j,p] <- t(Y_data[[p]][,1]) %*% X_data[[p]][,j] / N[p]
  }
}

Sigma_k_initial <- array(data = NA, dim = c(2,2,2))
Sigma_k_initial[,,1] <- matrix(data = c(0.08^2, 0.8, 0.8, 0.08^2), nrow = 2, ncol = 2)
Sigma_k_initial[,,2] <- matrix(data = c(1e-5^2, 0.8, 0.8, 1e-5^2), nrow = 2, ncol = 2)
sigma_epsilon_initial <- rep(0.99,2)
pii_initial <- c(0.2, 0.8)
eta_initial <- matrix(data = rep(0,200), nrow = 100, ncol = 2)
mu_jk_initial <- array(data = rep(0,400), dim = c(2,2,100))
gamma_jk_initial <- matrix(data = c(rep(0.2,100), rep(0.8,100)), nrow = 2, ncol = 100, byrow = TRUE)

lol <- E_step(X = X_data, Y = Y_data, M = 100, R = R_data, N = c(nrow(X_train_eas), nrow(X_train_eur)), K = 2, P = 2, sigma_epsilon_current = sigma_epsilon_initial, pii_current = pii_initial,  
                  eta_current = eta_initial, mu_jk_current = mu_jk_initial, Sigma_k_current = Sigma_k_initial, gamma_jk_current = gamma_jk_initial, beta_marginal = beta_marginal)


ujk_initial <- lol$ujk
Sigma_jk_initial <- lol$Sigma_jk
zeta_initial <- matrix(data = rep(0,200), nrow = 100, ncol = 2)

P <- 2; M <- 100; K <-2
N <- c(352, 273)
X <- X_data
Y <- Y_data
R <- R_data
maxiter <- 100
Sigma_jk_path <- vector("list", maxiter)
for(iter in 1:maxiter){
  Sigma_jk_path[[iter]] <- vector("list", M)
  for(j in 1:M){
    Sigma_jk_path[[iter]][[j]] <- array(data = NA, dim = c(P,P,K))
  }
}
Sigma_k_path <- vector("list", maxiter)
for(iter in 1:maxiter){
  Sigma_k_path[[iter]] <- array(data = NA, dim = c(P,P,K))
}
mu_jk_path <- vector("list", maxiter)
for(iter in 1:maxiter){
  mu_jk_path[[iter]] <- array(data = NA, dim = c(P,K,M))
}
u_jk_path <- array(data = NA, dim = c(K,M,maxiter))
gamma_jk_path <- array(data = NA, dim = c(K, M, maxiter))
eta_path <- array(data = NA, dim = c(M,P,maxiter))
sigma_epsilon_path <- matrix(data = NA, nrow = P, ncol = maxiter)
pii_path <- matrix(data = NA, nrow = K, ncol = maxiter)
zeta_path <- array(data = NA, dim = c(M, P, maxiter))
ELBO_path <- rep(0, maxiter)

Sigma_jk_path[[1]] <- Sigma_jk_initial
Sigma_k_path[[1]] <- Sigma_k_initial

mu_jk_path[[1]] <- mu_jk_initial
u_jk_path[,,1] <- ujk_initial
gamma_jk_path[,,1] <- gamma_jk_initial
eta_path[,,1] <- eta_initial
sigma_epsilon_path[,1] <- sigma_epsilon_initial
pii_path[,1] <- pii_initial
zeta_path[,,1] <- zeta_initial

for(iter in 2:50){
  E_list <- E_step(X = X, Y = Y, M = M, R = R, N = N, K = K, P = P, sigma_epsilon_current = sigma_epsilon_path[,iter-1], pii_current = pii_path[,iter-1],  
                   eta_current = eta_path[,,iter-1], mu_jk_current = mu_jk_path[[iter-1]], Sigma_k_current = Sigma_k_path[[iter-1]], gamma_jk_current = gamma_jk_path[,,iter-1], beta_marginal = beta_marginal)
  Sigma_jk_path[[iter]] = E_list$Sigma_jk
  mu_jk_path[[iter]] = E_list$mu_jk
  gamma_jk_path[,,iter] = E_list$gamma_jk
  eta_path[,,iter] = E_list$eta
  u_jk_path[,,iter] = E_list$ujk
  M_list <- M_step(gamma_jk = gamma_jk_path[,,iter], mu_jk = mu_jk_path[[iter]], Sigma_jk = Sigma_jk_path[[iter]],
                   eta = eta_path[,,iter], X = X, Y = Y, M = M, N = N, K = K, P = P)

  Sigma_k_path[[iter]] <- M_list$Sigma_k
  sigma_epsilon_path[,iter] <- M_list$sigma_epsilon
  pii_path[,iter] <- M_list$pii
  zeta_path[,,iter] <- M_list$zeta

  ELBO_path[iter] <- ELBO(N = N, sigma_epsilon = sigma_epsilon_path[,iter], Y = Y, X = X, zeta = zeta_path[,,iter], gamma_jk = gamma_jk_path[,,iter], pii = pii_path[,iter], Sigma_jk = Sigma_jk_path[[iter]], Sigma_k = Sigma_k_path[[iter]], mu_jk = mu_jk_path[[iter]], eta = eta_path[,,iter], M = M, K = K)
  # 
  # if(abs(ELBO_path[iter] - ELBO_path[iter-1]) < 1e-4){
  #   numiter <- iter
  #   break
  # }
}
y_train_predicted_eas <- X_data[[1]] %*% eta_path[,1,50]

data1 <- data.frame(Predicted.phenotype = y_train_predicted_eas, True.phenotype = Y_data[[1]][,1])
p1 <- ggplot(data1, aes(x = Predicted.phenotype, y = True.phenotype)) +
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "red")+
  xlab("Predicted phenotype")+
  ylab("True phenotype")
p1

y_train_predicted_eur <- X_data[[2]] %*% eta_path[,2,50]
data2 <- data.frame(Predicted.phenotype = y_train_predicted_eur, True.phenotype = Y_data[[2]][,1])
p2 <- ggplot(data2, aes(x = Predicted.phenotype, y = True.phenotype)) +
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "red")+
  xlab("Predicted phenotype")+
  ylab("True phenotype")
p2

library(ggplot2)
library(gridExtra)

# Create function to generate the plot for each dataset
create_plot <- function(data, title, corr) {
  ggplot(data, aes(x = Predicted.phenotype, y = True.phenotype)) +
    geom_point(size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(
      title = paste0(title, " (Pearson Correlation Coef. = ", corr, ")"),
      x = "Predicted phenotype",
      y = "True phenotype"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_line(color = "gray90"),
      panel.grid.major = element_line(color = "gray90"),
      plot.title = element_text(hjust = 0.5)
    ) +
    coord_cartesian(xlim = c(-0.001, 0.001), ylim = c(-3, 3))
}
train_PPC_eas <- cor(data1$Predicted.phenotype, data1$True.phenotype, method = "pearson")
train_PPC_eur <- cor(data2$Predicted.phenotype, data2$True.phenotype, method = "pearson")

# Create individual plots
p1 <- create_plot(data1, "East Asian", round(train_PPC_eas,4))
p2 <- create_plot(data2, "Europe", round(train_PPC_eur,4))

# Combine plots
combined_plot <- grid.arrange(p1, p2, ncol = 2)

# # Add overall title
# title <- grid.text("Figure 2: PRS prediction on training and testing set.",
#                    y = 0.02, 
#                    gp = gpar(fontsize = 12))
X_data_test <- vector("list",2)
scaled.X_test_eas <- scale(X_test_eas)
scaled.y_test_eas <- scale(y_test_eas)

scaled.X_test_eur <- scale(X_test_eur)
scaled.y_test_eur <- scale(y_test_eur)
X_data_test[[1]] <- scaled.X_test_eas
X_data_test[[2]] <- scaled.X_test_eur
Y_data_test <- vector("list", 2)
Y_data_test[[1]] <- scaled.y_test_eas
Y_data_test[[2]] <- scaled.y_test_eur

y_test_predicted_eas <- X_data_test[[1]] %*% eta_path[,1,50]
data1_test <- data.frame(Predicted.phenotype = y_test_predicted_eas, True.phenotype = Y_data_test[[1]][,1])
test_PPC_eas <- cor(data1_test$Predicted.phenotype, data1_test$True.phenotype, method = "pearson")

y_test_predicted_eur <- X_data_test[[2]] %*% eta_path[,2,50]
data2_test <- data.frame(Predicted.phenotype = y_test_predicted_eur, True.phenotype = Y_data_test[[2]][,1])
test_PPC_eur <- cor(data2_test$Predicted.phenotype, data2_test$True.phenotype, method = "pearson")

# Create individual plots
p1 <- create_plot(data1_test, "East Asian", round(test_PPC_eas,4))
p2 <- create_plot(data2_test, "Europe", round(test_PPC_eur,4))

# Combine plots
combined_plot <- grid.arrange(p1, p2, ncol = 2)

###################
X = X; Y = Y; M = M; R = R; N = N; K = K; P = P; 
sigma_epsilon_current = rep(0.95,2); pii_current = pii_path[,3];  
eta_current = eta_path[,,3]; mu_jk_current = mu_jk_path[[3]]; 
Sigma_k_current = Sigma_k_path[[3]]; gamma_jk_current = gamma_jk_path[,,3]


Sigma_jk_update <- vector("list", M)
for(j in 1:M){
  Sigma_jk_update[[j]] <- array(data = NA, dim = c(P,P,K))
}
Sigma_jk_update_inv <- vector("list",M)
for(j in 1:M){
  Sigma_jk_update_inv[[j]] <- array(data = NA, dim = c(P,P,K))
}
# eta_current <- eta_update
mu_jk_update <- mu_jk_current
# u_jk_update <- matrix(data = NA, nrow = M, ncol = K)
gamma_jk_update <- gamma_jk_current
eta_update <- eta_current
temp_mu_jk <- rep(0,P)
ujk_update <- matrix(data = NA, nrow = K, ncol = M)
temp_gamma_jk_update <- matrix(data = NA, nrow = K, ncol = M)
temp_Sigma_jk <- matrix(data= NA, nrow = P, ncol = M)
temp_mu_jk <- matrix(data= NA, nrow = P, ncol = M)
for(j in 1:M){
  for(k in 1:K){
    for(p in 1:P){
      temp_Sigma_jk[p,j] <- (t(X[[p]][,j])%*%X[[p]][,j])/sigma_epsilon_current[p]
    }
    Sigma_jk_update_inv[[j]][,,k] <- diag(temp_Sigma_jk[,j],nrow = P, ncol = P) + ginv(Sigma_k_current[,,k])
    Sigma_jk_update[[j]][,,k] <- ginv(Sigma_jk_update_inv[[j]][,,k])
  }
}
for(j in 1:M){
  for(p in 1:P){
    temp_mu_jk[p,j] <- (N[p] * (beta_marginal[j,p] -  t(R[-j,j,p])%*% eta_update[-j,p]))/sigma_epsilon_current[p]
  }
  for(k in 1:K){
    mu_jk_update[,k,j] <- Sigma_jk_update[[j]][,,k] %*% temp_mu_jk[,j]
  }
}
for(j in 1:M){
  for(k in 1:K){
    ujk_update[k,j] <- log(pii_current[k]) + 0.5 * log(abs(det(Sigma_jk_update[[j]][,,k]))/abs(det(Sigma_k_current[,,k]))) + 0.5 * tr(Sigma_jk_update_inv[[j]][,,k] %*% mu_jk_update[,k,j] %*%t(mu_jk_update[,k,j]))
    if(ujk_update[k,j] < log(0.01) | ujk_update[k,j] > log(0.99)){
      ujk_update[k,j] <- log(pii_current[k])
    }
  }
}
for(j in 1:M){
  for(k in 1:K){
    temp_gamma_jk_update[k,j] <- exp(-ujk_update[k,j]) / sum(exp(-ujk_update[,j]))
    # gamma_jk_update[k,j] <- temp_gamma_jk_update[k,j]
    if(temp_gamma_jk_update[k,j] < 0.01){
      gamma_jk_update[k,j] <- 0.01
    }else if(temp_gamma_jk_update[k,j] > 0.99){
      gamma_jk_update[k,j] <- 0.99
    }else{
      gamma_jk_update[k,j] <- temp_gamma_jk_update[k,j]
    }
  }
}
for(j in 1:M){
  for(p in 1:P){
    eta_update[j,p] <- gamma_jk_update[1,j] * mu_jk_update[p,1,j]
  }
}



