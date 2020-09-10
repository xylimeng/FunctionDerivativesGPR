# Function Derivatives via GPR 
Nonparametric estimation of function and function derivatives using Gaussian Process Regression 

## Example
The following code is available in the file "GPR.R".

Define the true regression function and its derivative.
```R
# True regression function
f0 = function(x){
  sum(sqrt(2) * I ^ (- 4) * sin(I) * cos((I - 0.5) * pi * x))
}

# Derivative of true regression function
f0_prime = function(x){
  sum(sqrt(2) * I ^ (- 4) * sin(I) * (0.5 - I) * pi * sin((I - 0.5) * pi * x))
}

set.seed(1000)
n = 1000 # number of data points
x = sort(runif(n, 0, 1)) # random design
I = 1:9999
y0 = sapply(x, f0) # true value at samples
sigma0 = sqrt(0.1) # true regression standard deviation
y = y0 + sigma0 * rnorm(n)  #generate some data

n_new = 100 # number of grid points
x_new = seq(0, 1,, n_new) # grid points
f0_new = sapply(x_new, f0) # true value at grid points
f0_new_prime = sapply(x_new, f0_prime) # true derivative at grid points
```
