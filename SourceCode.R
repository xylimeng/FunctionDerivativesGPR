library(RandomFieldsUtils)
library(MASS)
library(Rcpp)
sourceCpp("dmvrm.cpp")
par(mai = c(0.5, 0.5, 0.1, 0.1))

# Second-order Sobolev kernel
sob = function(x1, x2){
  n1 = length(x1)
  n2 = length(x2)
  XX = t(matrix(rep(x1, n2), ncol = n2))
  xx = matrix(rep(x2, n1), ncol = n1)
  return(1 + XX * xx + pmin(XX, xx)^2 * (3 * pmax(XX, xx) - pmin(XX, xx)) / 6)
}

sob_prime = function(x1, x2){
  n1 = length(x1)
  n2 = length(x2)
  XX = t(matrix(rep(x1, n2), ncol = n2))
  xx = matrix(rep(x2, n1), ncol = n1)
  return(XX + pmin(XX, xx) * XX - pmin(XX, xx)^2 / 2)
}

sob_2prime = function(x1, x2){
  n1 = length(x1)
  n2 = length(x2)
  XX = t(matrix(rep(x1, n2), ncol = n2))
  xx = matrix(rep(x2, n1), ncol = n1)
  return(1 + pmin(XX, xx))
}

# Squared exponential kernel
SE = function(x1, x2){
  n1 = length(x1)
  n2 = length(x2)
  XX = matrix(rep(x2, n1), ncol = n1) - t(matrix(rep(x1, n2), ncol = n2))
  return(matrix(exp(-XX^2), ncol = n1))
}

SE_prime = function(x1, x2){
  n1 = length(x1)
  n2 = length(x2)
  XX = matrix(rep(x2, n1), ncol = n1) - t(matrix(rep(x1, n2), ncol = n2))
  return(-2 * XX * matrix(exp(-XX^2), ncol = n1))
}

SE_2prime = function(x1, x2){
  n1 = length(x1)
  n2 = length(x2)
  XX = matrix(rep(x2, n1), ncol = n1) - t(matrix(rep(x1, n2), ncol = n2))
  return(2 * matrix(exp(-XX^2), ncol = n1) - 4 * XX^2 * matrix(exp(-XX^2), ncol = n1))
}

# Matern kernel
mat = function(x1, x2, nu){
  n1 = length(x1)
  n2 = length(x2)
  XX = matrix(rep(x2, n1), ncol = n1) - t(matrix(rep(x1, n2), ncol = n2))
  return(matrix(matern(abs(XX), nu, derivative = 0), ncol = n1))
}

mat_prime = function(x1, x2, nu){
  n1 = length(x1)
  n2 = length(x2)
  XX = matrix(rep(x2, n1), ncol = n1) - t(matrix(rep(x1, n2), ncol = n2))
  return(matrix(matern(abs(XX), nu, derivative = 1), ncol = n1) * sign(XX))
}

mat_2prime = function(x1, x2, nu){
  n1 = length(x1)
  n2 = length(x2)
  XX = matrix(rep(x2, n1), ncol = n1) - t(matrix(rep(x1, n2), ncol = n2))
  return(-matrix(matern(abs(XX), nu, derivative = 2), ncol = n1))
}

# Choose the best smoothness parameter for Matern kernel via leave-one-out cross validation 
get_loo = function(x, y, nu, lambda){
  MSE = rep(0, n)
  for(i in 1:n){
    test_y = y[i]
    test_x = x[i]
    train_y = y[-i]
    train_x = x[-i]
    
    K_XX = mat(train_x, train_x, nu)
    #XX = abs(matrix(rep(trainx, n-1), ncol = n-1) - t(matrix(rep(trainx, n-1), ncol = n-1)))
    #K_XX = matrix(matern(abs(XX), nu, derivative = 0), ncol = n-1)
    K_middle = t(chol(K_XX + (n-1) * lambda * diag(n-1)))
    K_Xx = mat(train_x, test_x, nu)
    #Xx = matrix(rep(testx, n-1), ncol = n-1) - t(matrix(rep(trainx, 1), ncol = 1))
    #K_Xx = matrix(matern(abs(Xx), nu, derivative = 0), ncol = n-1)
    f_hat = as.vector(crossprod(forwardsolve(K_middle, t(K_Xx)), forwardsolve(K_middle, train_y)))
    MSE[i] = mean((f_hat - test_y)^2)
  }
  return(mean(MSE))
}

get_func_est = function(x, x_new, y, lambda, kernel, nu = NULL){
  # Computer kernel matrix and derivatives
  if(kernel == "Sobolev"){
    K_XX = sob(x, x)
    K_Xx = sob(x, x_new)
    K_xx = sob(x_new, x_new)
    K10_Xx = sob_prime(x, x_new)
    K11_xx = sob_2prime(x_new, x_new)
  }
  if(kernel == "SE"){
    K_XX = SE(x, x)
    K_Xx = SE(x, x_new)
    K_xx = SE(x_new, x_new)
    K10_Xx = SE_prime(x, x_new)
    K11_xx = SE_2prime(x_new, x_new)
  }
  if(kernel == "Matern"){
    K_XX = mat(x, x, nu)
    K_Xx = mat(x, x_new, nu)
    K_xx = mat(x_new, x_new, nu)
    K10_Xx = mat_prime(x, x_new, nu)
    K11_xx = mat_2prime(x_new, x_new, nu)
  }
  
  # Compute posterior mean and its derivative
  K_middle = t(chol(K_XX + n * lambda * diag(n)))
  sigma_hat = sqrt(lambda * as.vector(crossprod(forwardsolve(K_middle, y))))
  f_hat = as.vector(crossprod(forwardsolve(K_middle, t(K_Xx)), forwardsolve(K_middle, y)))
  f_prime_hat = as.vector(crossprod(forwardsolve(K_middle, t(K10_Xx)), forwardsolve(K_middle, y)))
  
  # Compute posterior variance
  cov = sigma_hat^2 / (n * lambda) * (K_xx - crossprod(forwardsolve(K_middle, t(K_Xx))))
  
  # Compute the L-infinity simultaneous credible band for posterior mean
  pf_hat = mvrnorm(n = 1000, mu = rep(0, n_new), Sigma = cov)
  pf_max = apply(abs(pf_hat), 1, max)
  radius = quantile(pf_max, 1 - 0.05)
  f_min = f_hat - radius
  f_max = f_hat + radius
  
  # Plot for estimating regression function
  plot(x_new, f0_new, type = "l", lwd = 3, xlab = "", ylab = "",
       ylim = c(min(c(f_min, f0_new)), max(c(f_max, f0_new))))
  lines(x_new, f_hat, lty = 3, lwd = 3)
  lines(x_new, f_min, lty = 5, lwd = 3)
  lines(x_new, f_max, lty = 5, lwd = 3)
  
  # Compute the L-infinity simultaneous credible band for derivative of posterior mean
  cov_11 = sigma_hat^2 / (n * lambda) * (K11_xx - crossprod(forwardsolve(K_middle, t(K10_Xx))))
  pf_prime_hat = mvrnorm(n = 1000, mu = rep(0, n_new), Sigma = cov_11)
  pr_prime_max = apply(abs(pf_prime_hat), 1, max)
  radius_prime = quantile(pr_prime_max, 1 - 0.05)
  f_prime_min = f_prime_hat - radius_prime
  f_prime_max = f_prime_hat + radius_prime
  
  # Plot for estimating derivative of regression function
  plot(x_new, f0_new_prime, type = "l", lwd = 3, xlab = "", ylab = "",
       ylim = c(min(c(f_prime_min, f0_new_prime)), max(c(f_prime_max, f0_new_prime))))
  lines(x_new, f_prime_hat, lty = 3, lwd = 3)
  lines(x_new, f_prime_min, lty = 5, lwd = 3)
  lines(x_new, f_prime_max, lty = 5, lwd = 3)
  
  # Report RMSE
  return(c(sqrt(sum((f_hat - f0_new)^2)), sqrt(sum((f_prime_hat - f0_new_prime)^2))))
}

# Log marginal posterior likelihood
log_post = function(x, y, lambda, kernel, nu = NULL){
  if(kernel == "Sobolev") K_XX = sob(x, x)
  if(kernel == "SE") K_XX = SE(x, x)
  if(kernel == "Matern") K_XX = mat(x, x, nu)
  K_middle = t(chol(K_XX + n * lambda * diag(n)))
  mar_var = 1/n * as.vector(crossprod(forwardsolve(K_middle, y))) * (K_XX + n * lambda * diag(n))
  log_p = dmvnrm_arma_fast(t(y), rep(0, n), mar_var, TRUE)
  return(log_p)
}

# Main function
get_GPR = function(kernel){
  # Matern kernel
  if(kernel == "Matern"){
    # Candidates of smoothness parameter
    nu_min = 2; nu_max = 3
    seq_nu = seq(nu_min, nu_max, by = 0.5)
    n_seq = length(seq_nu)
    seq_lambda = mseBF = rep(0, n_seq)
    
    # Choose the best smoothness parameter via leave-one-out cross validation
    for (i in 1:n_seq){
      # Compute the MMLE of lambda
      opt = optim(par = 1/n,
                  fn = log_post,
                  control = list(fnscale = -1),
                  method = c("L-BFGS-B"), lower = 1e-6, upper = 1,
                  x = x, y = y, kernel = "Matern", nu = seq_nu[i])
      seq_lambda[i] = opt$par
      mseBF[i] = get_loo(x = x, y = y, nu = seq_nu[i], lambda = seq_lambda[i])
    }
    opt_nu = seq_nu[which(mseBF==min(mseBF))]
    lambda = seq_lambda[which(mseBF==min(mseBF))]
    
    return(get_func_est(x, x_new, y, lambda, kernel, nu = opt_nu))
    
    # Sobolev kernel of SE kernel
  } else{
    # Compute the MMLE of lambda
    opt = optim(par = 1/n,
                fn = log_post,
                control = list(fnscale = -1),
                method = c("L-BFGS-B"), lower = 1e-6, upper = 1,
                x = x, y = y, kernel = "SE")
    lambda = opt$par
    
    return(get_func_est(x, x_new, y, lambda, kernel))
  }
}