
sigma_epsilon_current <- sigma_epsilon_initial
mu_jk_current <- mu_jk_initial
gamma_jk_current <- gamma_jk_initial
eta_current <- eta_initial
# Sigma_jk_update <- vector("list", 100)
# for(i in 1:100){
#   Sigma_jk_update[[i]] <- array(data = NA, dim = c(2,2,2))
# }
# for(i in 1:100){
#   Sigma_jk_update[[i]][,,1] <- matrix(data = c(0.08, 0.8, 0.8, 0.08), nrow = 2, ncol = 2)
#   Sigma_jk_update[[i]][,,2] <- matrix(data = c(1e-5, 0.8, 0.8, 1e-5), nrow = 2, ncol = 2)
# }
X <- X_data
Y <- Y_data
P <- 2; M <- 100; K <- 2
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
      temp_mu_jk[p] <- (t(X[[p]][,j])%*%Y[[p]] - t(eta_update[-j,p]) %*% (R[-j,j,p] * N[p]))/sigma_epsilon_current[p]
    }
    mu_jk_update[,k,j] <- Sigma_jk_update[[j]][,,k] %*% temp_mu_jk
  }
}
for(j in 1:M){
  for(k in 1:K){
    ujk_update[k,j] <- log(pii_current[k]) + 0.5 * log(abs(det(Sigma_jk_update[[j]][,,k]))/abs(det(Sigma_k_current[,,k]))) + 0.5 * tr(Sigma_jk_update_inv[[j]][,,k] %*% mu_jk_update[,k,j] %*%t(mu_jk_update[,k,j]))
  }
}
for(j in 1:M){
  for(k in 1:K){
    temp_gamma_jk_update[k,j] <- exp(ujk_update[k,j]) / sum(exp(ujk_update[,j]))
    # gamma_jk_update[k,j] <- temp_gamma_jk_update[k,j]
    if(temp_gamma_jk_update[k,j] < 0.001){
      gamma_jk_update[k,j] <- 0.001
    }else if(temp_gamma_jk_update[k,j] > 0.999){
      gamma_jk_update[k,j] <- 0.999
    }else{
      gamma_jk_update[k,j] <- temp_gamma_jk_update[k,j]
    }
  }
}
for(j in 1:M){
  for(p in 1:P){
    eta_update[j,p] <- t(gamma_jk_update[,j]) %*% mu_jk_update[p,,j]
  }
}
