#' ---
#' title: "Functional Quantile Regression"
#' author: "Meng Li and Kehui Wang"
#' date: "5/27/2021"
#' output: "html_document"
#' ---
#' 

#+ message = FALSE
library(quantreg) # R package for multi-covariate quantile regression 
library(refund) # R package for FPCA 

### sim.fQR simulates data to test proposed methods -----
# n.sub - number of subjects
# gamma - a parameter used in the model to control the deviation from the null hypothesis
# n.tp - number of time points (equal-spaced in [0, 1]) 
# eigen.basis: basis function to gerate functional covariates
# eigen.value: decreasing and non-zero 
# noise.x: standard deviation of the white noise in the functional covariates 
# missing.rate: missing rate in the functional covariates 

sim.fQR<- function(n.sub, gamma = 0, n.tp = 100, 
                   eigen.basis = "DFT", eigen.value = c(1, 1/2, 1/4), noise.x = 1, missing.rate = 0){
  
  if (eigen.basis == "DFT"){
    eigen.basis.fun <- function(t, degree = 0){
      length.t= length(t)
      if(degree==0){temp.lp=rep(1, length.t)}
      if ((degree > 0) && (degree %% 2 == 1)){
        temp.lp = sqrt(2) * cos(2 * pi * (degree + 1)/2 * t)
      }
      if ((degree > 0) && (degree %% 2 == 0)){
        temp.lp = sqrt(2) * sin(2 * pi * degree/2 * t)
      }
      return(temp.lp)
    }
  }
  
  if (eigen.basis == "legendre.polynomials"){
    eigen.basis.fun =function(t, degree = 0){
      length.t= length(t)
      if(degree==0){temp.lp = rep(1, length.t)}
      if(degree==1){temp.lp =sqrt(3)*(2*t-1)}
      if(degree==2){temp.lp =sqrt(5)*(6*t^2  - 6*t +1) }
      if(degree==3){temp.lp =sqrt(7)*(20*t^3 - 30*t^2  +12 *t -1) } 
      return(temp.lp)
    }
  }
  
  eigen.value = eigen.value[eigen.value > 0]
  eigen.value = sort(eigen.value, decreasing = TRUE)
  cat(sprintf("Functional covariates: \n eigen.value = %s \n eigen.function = %s",
              paste(eigen.value, collapse = " "), eigen.basis))
  
  tp <- seq(from = 0, to = 1, length = n.tp)  # dense design   
  K.lambda = length(eigen.value) 
  phi.basys.system = sapply(c(1:K.lambda), function(k) eigen.basis.fun(tp,degree = k))
  # start from 1, so we don't use the intercept to generate data
  
  cov.true =  phi.basys.system %*% diag(eigen.value)  %*% t(phi.basys.system) + diag(noise.x^2, n.tp)      # white noise
  scores.matrix = matrix(rnorm(n.sub*K.lambda), ncol=n.sub)
  var.pointwise = diag(phi.basys.system %*% diag(eigen.value)  %*% t(phi.basys.system) )
  noise = noise.x * matrix(rnorm(n.sub*length(tp)), nrow=n.sub)    
  X = t(phi.basys.system %*% diag(sqrt(eigen.value)) %*% scores.matrix ) + noise    
  
  ## next we generate the response according to quantile model 
  b <- tp
  btx <- apply(X, 1, function(k) mean(k* b))    
  y <- btx +  rnorm(n.sub) * (1 + apply(X, 1, function(k) mean(k* gamma *tp^2)) )  
  
  ## possible missing in the functional covariates 
  if (missing.rate > 0){
    for (i in 1:n.sub){
      idx2missing = sample(n.tp, size = floor(n.tp * missing.rate), replace = FALSE)
      X[i, idx2missing] <- NA
    }
  }   
  
  ret = list(Y = y, X = X)
  return(ret)
}



### fQR returns functional quantile regression estimation -----
# Y: vector of responses
# X: matrix of covariates; NA's allowed for missing values
# taus: quantile levels as a vector
# PVE: percentage variance explained used in FPCA (default value is 0.95)
fQR <- function(Y, X, taus, pve = 0.95){
  x.fpca <- fpca.sc(X, pve = pve, var=T)
  x.new <- x.fpca$scores
  
  # wald test function 
  wald.func = function(y,x,taus) {
    # x is the design matrix WITHOUT intercept 
    n= length(y)
    X= cbind(rep(1,n),x) # design matrix 
    K= length(taus)
    p= ncol(X)
    f_K = matrix(0,n,K)
    Omega = outer(taus, taus, pmin) - outer(taus, taus)
    J = crossprod(X)/n
    Hinv = list(length = K)
    for (k in 1:K){
      eps = .Machine$double.eps^(2/3)
      h <- bandwidth.rq(taus[k], n, hs = TRUE)
      if (taus[k] + h > 1) {h = 1 - taus[k]; warning("tau + h > 1")}
      if (taus[k] - h < 0) {h = taus[k]; warning("tau - h < 0")}
      bhi <- rq.fit.fnb(X, y, tau = taus[k] + h)$coef
      blo <- rq.fit.fnb(X, y, tau = taus[k] - h)$coef
      dyhat <- X %*% (bhi - blo)         # X is design matrix with column of 1's
      if (any(dyhat <= 0)) 
        warning(paste(sum(dyhat <= 0), "non-positive fis"))
      f_K[,k] <- pmax(0, (2 * h)/(dyhat - eps))
      fxxinv <- diag(p)
      fxxinv <- backsolve(qr(sqrt(f_K[,k]) * X)$qr[1:p, 1:p, drop = FALSE], fxxinv)
      Hinv[[k]] <- n*fxxinv %*% t(fxxinv) # inv(H) in Page 74
    }
    V = matrix(0,K*p,K*p)
    for (i in 1:K){
      for (j in 1:K){
        V[((i-1)*p+1):(i*p), ((j-1)*p+1):(j*p)] = Omega[i,j]*Hinv[[i]]%*%J%*%Hinv[[j]]
      }
    }
    
    coef = as.vector(rq(y~x,taus)$coefficients)
    D <- kronecker(diff(diag(K)), cbind(0, diag(p - 1)))  # LHS matrix in null hypothesis
    ndf <- (p - 1) * (K - 1)
    Tn <- as.numeric(n*t(D %*% coef) %*% solve(D %*% V %*% t(D), D %*% coef))/ndf
    ddf <- n * K - (p - 1) * (K - 1)
    pvalue.f <- 1-pf(Tn,ndf,ddf) # F test 
    # pvalue.c <- 1-pchisq(ndf*Tn,ndf) # Chi-squared test 
    
    coef = matrix(coef, nrow = ncol(x) + 1)
    beta.hat = x.fpca$efunctions %*% coef[-1, ]
    
    out=list(pvalue = pvalue.f, Tn=Tn, coef=coef, p = p, K = K, x.fpca = x.fpca)
    return(out)
  }
  
  # 
  ret = wald.func(y = Y, x = x.new, taus)
  coef = matrix(ret$coef, nrow = ncol(x.new) + 1)
  beta.hat = x.fpca$efunctions %*% coef[-1, ]
  
  return(list(beta.hat = beta.hat, coef = coef, pvalue = ret$pvalue, p = ret$p, 
              K = ret$K, x.fpca = x.fpca))
}




### Composite Quantile regression ---- 


# Method 2: CRQ 
# Reference: Wang, K. and Wang, H. J. (2015). Optimally combined estimation for tail quantile regression. Statistica Sinica, To appear (doi:10.5705/ss.2014.051)

rq.fit.fnb2 <- function (x, y, taus, beta = 0.99995, eps = 1e-06)
{
  # modified rq.fit.fnb function
  # taus is a vector, taus=(tau_1,...,tau_n)
  # minimize \sum \rho_{\tau_i}(y_i-x_i^T\theta)
  
  n <- length(y)
  p <- ncol(x)
  if (n != nrow(x))
    stop("x and y don't match n")
  #if (tau < eps || tau > 1 - eps)
  #    stop("No parametric Frisch-Newton method.  Set tau in (0,1)")
  rhs <- t(x)%*%(1-taus)
  d <- rep(1, n)
  u <- rep(1, n)
  wn <- rep(0, 10 * n)
  wn[1:n] <- (1 - taus)
  z <- .Fortran("rqfnb", as.integer(n), as.integer(p), a = as.double(t(as.matrix(x))),
                c = as.double(-y), rhs = as.double(rhs), d = as.double(d),
                as.double(u), beta = as.double(beta), eps = as.double(eps),
                wn = as.double(wn), wp = double((p + 3) * p), it.count = integer(3),
                info = integer(1), PACKAGE = "quantreg")
  if (z$info != 0)
    stop(paste("Error info = ", z$info, "in stepy: singular design"))
  coefficients <- -z$wp[1:p]
  names(coefficients) <- dimnames(x)[[2]]
  residuals <- y - x %*% coefficients
  list(coefficients = coefficients, taus = taus, residuals = residuals)
}



crq.fit.func = function(y, x, taus)
{
  # obtain composite estimate assuming common slopes at taus, equal weights at all quantiles
  # x does not include the intercept
  n = length(y)
  ntau = length(taus)
  Taus = rep(taus, each=n)
  y2 = rep(y, ntau)
  x2 = cbind(diag(ntau)%x%rep(1,n), rep(1,ntau)%x%x)
  tmp = rq.fit.fnb2(x2, y2, taus = Taus, beta = 0.99995, eps = 1e-06)
  Coef = tmp$coef
  intersp = Coef[1:ntau]
  slope = matrix(rep(Coef[-(1:ntau)], each=ntau), byrow=T, ncol=ntau)
  Coef = rbind(intersp, slope)
  return(Coef)
}

