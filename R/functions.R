#' @importFrom stats as.formula
#' @importFrom stats lm
#' @importFrom stats pnorm
#' @importFrom stats rnorm
#' @importFrom utils combn
#' @importFrom lmf rmnorm
#' @importFrom lmf inv
#' @import VGAM


#' @export fit_all_subset_linear_models
fit_all_subset_linear_models <- function(y, X, intercept=c(TRUE,FALSE)){

  n <- nrow(X)
  p <- ncol(X)

  if(intercept == TRUE){
    X=cbind(rep(1,n),X)
  }

  best_adj_r_squared <- -Inf
  best_model <- NULL
  k <- vector("list", length = p)
  R_M_k <- vector("list", length = p)
  kappa_M_k <- vector("list", length = p)
  sst <- sum((y-mean(y))^2)
  for (i in 1:p) {
    k[[i]] <- list()
    R_M_k[[i]] <- list()
    kappa_M_k[[i]] <- list()
    subsets <- combn(p, i)
    for (j in 1:ncol(subsets)) {
      selected_cols <- subsets[,j]
      Construct_adj_r <- construct_adj_r_squared(X,selected_cols,y,n,intercept,sst)
      subset_X <- Construct_adj_r$X_M_k
      R_M_k_current <- Construct_adj_r$R_M_k
      kappa_M_k_current <- Construct_adj_r$kappa_M_k
      adj_r_squared <- Construct_adj_r$adj_r_squared

      if (adj_r_squared > best_adj_r_squared) {
        best_adj_r_squared <- adj_r_squared
        if(intercept == FALSE){
          model_formula <- as.formula(paste("y ~ 0 +", paste(paste0("X[,", selected_cols, "]"), collapse = "+")))
          best_model <- lm(model_formula, data = data.frame(y, subset_X))
        }
        if(intercept == TRUE){
          model_formula <- as.formula(paste("y ~ ", paste(paste0("X[,", selected_cols+rep(1,length(selected_cols)), "]"), collapse = "+")))
          best_model <- lm(model_formula, data = data.frame(y, subset_X[,-1]))
        }
        phat <- selected_cols
        X_M_phat <- subset_X
        R_M_phat <- R_M_k_current
        kappa_M_phat <- kappa_M_k_current
      }
      k[[i]][[j]] <- as.vector(selected_cols)
      R_M_k[[i]][[j]] <- as.matrix(R_M_k_current)
      kappa_M_k[[i]][[j]] <- kappa_M_k_current
    }
  }
  return(list(k = k, best_model = best_model, phat = phat, X_M_phat = X_M_phat, best_adj_r_squared = best_adj_r_squared,
              R_M_phat = R_M_phat, kappa_M_phat = kappa_M_phat, R_M_k = R_M_k, kappa_M_k = kappa_M_k))
}


#' @export fit_specified_size_subset_linear_models
fit_specified_size_subset_linear_models <- function(y, X, size, intercept=c(TRUE,FALSE)){

  n <- nrow(X)
  p <- ncol(X)
  subsets <- combn(p, size)
  K <- ncol(subsets)

  if(intercept == TRUE){
    X=cbind(rep(1,n),X)
  }

  best_adj_r_squared <- -Inf
  best_model <- NULL
  k <- vector("list", length = K)
  R_M_k <- vector("list", length = K)
  kappa_M_k <- vector("list", length = K)
  sst <- sum((y-mean(y))^2)
  for (i in 1:K){
    R_M_k[[i]] <- list()
    kappa_M_k[[i]] <- list()
    selected_cols <- subsets[,i]
    Construct_adj_r <- construct_adj_r_squared(X,selected_cols,y,n,intercept,sst)
    subset_X <- Construct_adj_r$X_M_k
    R_M_k_current <- Construct_adj_r$R_M_k
    kappa_M_k_current <- Construct_adj_r$kappa_M_k
    adj_r_squared <- Construct_adj_r$adj_r_squared

    if (adj_r_squared > best_adj_r_squared) {
      best_adj_r_squared <- adj_r_squared
      if(intercept == FALSE){
        model_formula <- as.formula(paste("y ~ 0 +", paste(paste0("X[,", selected_cols, "]"), collapse = "+")))
        best_model <- lm(model_formula, data = data.frame(y, subset_X))
      }
      if(intercept == TRUE){
        model_formula <- as.formula(paste("y ~ ", paste(paste0("X[,", selected_cols+rep(1,length(selected_cols)), "]"), collapse = "+")))
        best_model <- lm(model_formula, data = data.frame(y, subset_X[,-1]))
      }
      phat <- selected_cols
      X_M_phat <- subset_X
      R_M_phat <- R_M_k_current
      kappa_M_phat <- kappa_M_k_current
    }
    k[[i]] <- as.vector(selected_cols)
    R_M_k[[i]] <- as.matrix(R_M_k_current)
    kappa_M_k[[i]] <- kappa_M_k_current
  }
  return(list(k = k, best_model = best_model, phat = phat, X_M_phat = X_M_phat, best_adj_r_squared = best_adj_r_squared,
              R_M_phat = R_M_phat, kappa_M_phat = kappa_M_phat, R_M_k = R_M_k, kappa_M_k = kappa_M_k))
}


#' @export construct_test_statistic
construct_test_statistic <- function(j,X_M_phat,y,phat,Sigma,intercept=c(TRUE,FALSE)){

  n=nrow(y)
  ej=c()

  if(intercept == FALSE){

    for(each_j in 1:length(phat)){
      if(phat[each_j]==j){
        ej[each_j]=1
      }
      else{
        ej[each_j]=0
      }
    }

    inv=inv(t(X_M_phat)%*%X_M_phat)
    X_M_phatinv=X_M_phat%*%inv
    etaj=X_M_phatinv%*%ej
    etajTy=t(etaj)%*%y

    sq_norm <- norm(t(etaj)%*%Sigma%*%etaj,type='F') # Frobenius norm of matrix is L2-norm of vector
    a <- (diag(n)-(Sigma%*%etaj%*%t(etaj))/sq_norm)%*%y
    b <- Sigma%*%etaj/sq_norm
  }

  if(intercept == TRUE){
    phat=c(1,phat+rep(1,length(phat)))
    for(each_j in 1:length(phat)){
      if(phat[each_j]==j){
        ej[each_j]=1
      }
      else{
        ej[each_j]=0
      }
    }
    inv=inv(t(X_M_phat)%*%X_M_phat)
    X_M_phatinv=X_M_phat%*%inv
    etaj=X_M_phatinv%*%ej
    etajTy=t(etaj)%*%y

    sq_norm <- norm(t(etaj)%*%Sigma%*%etaj,type='F') # Frobenius norm of matrix is L2-norm of vector
    a <- (diag(n)-(Sigma%*%etaj%*%t(etaj))/sq_norm)%*%y
    b <- Sigma%*%etaj/sq_norm
  }

  return(list('etaj'= etaj,
              'etajTy'= etajTy,
              'a'= a,
              'b'= b))
}


#' @export construct_adj_r_squared
construct_adj_r_squared <- function(X,k,y,n,intercept=c(TRUE,FALSE),sst){

  if(intercept == FALSE){
    X_M_k <- as.matrix(X[,k])
    kappa_M_k <- (n-length(k))^(-1)
  }
  if(intercept == TRUE){
    X_M_k <- cbind(rep(1,n),as.matrix(X[,k+rep(1,length(k))]))
    kappa_M_k <- (n-length(k)-1)^(-1)
  }
  P_M_k <- as.matrix(X_M_k%*%inv((t(X_M_k)%*%X_M_k))%*%t(X_M_k))
  R_M_k <- as.matrix(diag(n) - P_M_k)

  adj_r_squared <- 1-((t(y)%*%R_M_k%*%y*kappa_M_k)/(sst/(n-1)))

  return(list('X_M_k'=X_M_k,
              'P_M_k'=P_M_k,
              'R_M_k'=R_M_k,
              'kappa_M_k'=kappa_M_k,
              'adj_r_squared'=adj_r_squared))
}


#' @export construct_selection_event
construct_selection_event <- function(a,b,R_M_k,kappa_M_k,R_M_phat,kappa_M_phat){

  c_k <- (t(b)%*%R_M_k%*%b)*kappa_M_k-(t(b)%*%R_M_phat%*%b)*kappa_M_phat
  d_k <- 2*((t(b)%*%R_M_k%*%a)*kappa_M_k-(t(b)%*%R_M_phat%*%a)*kappa_M_phat)
  e_k <- (t(a)%*%R_M_k%*%a)*kappa_M_k-(t(a)%*%R_M_phat%*%a)*kappa_M_phat

  return(list('c_k'=c_k,
              'd_k'=d_k,
              'e_k'=e_k))
}


#' @export solve_selection_event
solve_selection_event <- function(a,b,R_M_k,kappa_M_k,R_M_phat,kappa_M_phat,k){

  i <- 1

  intervals_list <- list()  # Initialize a list to store intervals

  while(i <= length(k)){

    if(class(k[[i]])=="list"){ # If "fit_all_subset_linear_models" was used, we have to loop over two dimensions

      j <- 1

      while(j <= length(k[[i]])){

        Construct_event <- construct_selection_event(a,b,R_M_k[[i]][[j]],kappa_M_k[[i]][[j]],R_M_phat,kappa_M_phat)
        c_k=Construct_event$c_k
        d_k=Construct_event$d_k
        e_k=Construct_event$e_k

        # Calculate the discriminant
        D <- d_k^2 - 4 * c_k * e_k

        # Find the roots of the quadratic equation
        roots <- NULL
        if (D >= 1e-3) {
          root1 <- (-d_k + sqrt(D)) / (2 * c_k)
          root2 <- (-d_k - sqrt(D)) / (2 * c_k)
          roots <- sort(c(root1, root2))
        }

        # Function to check if the inequality is satisfied for a value of Z
        satisfies_inequality <- function(Z) {
          return(c_k * Z^2 + d_k * Z + e_k >= 0)
        }

        # Find intervals where the inequality holds
        intervals <- list()
        if (D < 1e-3) {
          intervals[[1]] <- c(-Inf, Inf)
          intervals[[2]] <- c(-Inf, Inf)
        } else {
          if (roots[1]  > -Inf) {
            intervals[[1]] <- c(-Inf, roots[1])
          }
          for (l in 1:(length(roots) - 1)) {
            if (satisfies_inequality((roots[l] + roots[l + 1]) / 2)) { # Fills in root[1] +  root[2]/2 for Z to see if area between two roots satisfies inequality
              intervals[[l + 1]] <- c(roots[l], roots[l + 1]) # If so, it returns this interval
            }
          }
          if (roots[2] < Inf) {
            intervals[[2]] <- c(roots[length(roots)], Inf)
          }
        }

        # Store intervals for this step of the inner loop
        if (D >= 1e-3) {
          intervals_list <- c(intervals_list, intervals)
        }
        j <- j+1
      }
    }

    else {

      Construct_event <- construct_selection_event(a,b,R_M_k[[i]],kappa_M_k[[i]],R_M_phat,kappa_M_phat)
      c_k=Construct_event$c_k
      d_k=Construct_event$d_k
      e_k=Construct_event$e_k

      # Calculate the discriminant
      D <- d_k^2 - 4 * c_k * e_k

      # Find the roots of the quadratic equation
      roots <- NULL
      if (D >= 1e-3) {
        root1 <- (-d_k + sqrt(D)) / (2 * c_k)
        root2 <- (-d_k - sqrt(D)) / (2 * c_k)
        roots <- sort(c(root1, root2))
      }

      # Function to check if the inequality is satisfied for a value of Z
      satisfies_inequality <- function(Z) {
        return(c_k * Z^2 + d_k * Z + e_k >= 0)
      }

      # Find intervals where the inequality holds
      intervals <- list()
      if (D < 1e-3) {
        intervals[[1]] <- c(-Inf, Inf)
        intervals[[2]] <- c(-Inf, Inf)
      } else {
        if (roots[1]  > -Inf) {
          intervals[[1]] <- c(-Inf, roots[1])
        }
        for (l in 1:(length(roots) - 1)) {
          if (satisfies_inequality((roots[l] + roots[l + 1]) / 2)) { # Fills in root[1] +  root[2]/2 for Z to see if area between two roots satisfies inequality
            intervals[[l + 1]] <- c(roots[l], roots[l + 1]) # If so, it returns this interval
          }
        }
        if (roots[2] < Inf) {
          intervals[[2]] <- c(roots[length(roots)], Inf)
        }
      }

      # if (D >= 1e-3) {
      intervals_list <- c(intervals_list, intervals)
      # }
    }

    j <- 1
    i <- i + 1
  }

  # Pair up the intervals
  paired_intervals <- split(intervals_list, rep(1:(length(intervals_list)/2), each = 2))

  paired_intervals_clean <- list()
  for(i in 1:length(paired_intervals)){
    if(paired_intervals[[i]][[1]][[2]]<0 & paired_intervals[[i]][[2]][[1]]>0){ # Take intersection of intervals
      paired_intervals_clean <- c(paired_intervals_clean, list(paired_intervals[[i]]))
    }
  }

  if(length(paired_intervals_clean)==0){
    intervals[[1]] <- c(-Inf, Inf)
    intervals[[2]] <- c(-Inf, Inf)
    paired_intervals_clean=c(paired_intervals_clean,list(intervals))
  }

  # Initialize vectors for interval bounds
  interval_1_ub <- c()
  interval_2_lb <- c()

  # Populate interval bounds vectors
  for (i in 1:length(paired_intervals_clean)) {
    interval_1_ub <- c(interval_1_ub, paired_intervals_clean[[i]][[1]][[2]])
    interval_2_lb <- c(interval_2_lb, paired_intervals_clean[[i]][[2]][[1]])
  }

  # Calculate intervals
  min_interval_1_ub <- min(interval_1_ub)
  max_interval_2_lb <- max(interval_2_lb)

  # Create z_interval with specified intervals
  z_interval <- list(
    c(-Inf, min_interval_1_ub),
    c(max_interval_2_lb, Inf)
  )

  return(list('intervals'=intervals_list,
              'z_interval'=z_interval))
}


#' @export pivot_with_specified_interval
pivot_with_specified_interval <- function(z_interval, etaj, etajTy, tn_mu, tn_sigma){

  numerator=0
  denominator=0

  for(i in 1:length((z_interval))){

    al=z_interval[[i]][1]
    ar=z_interval[[i]][2]

    denominator=denominator+pnorm((ar-tn_mu)/tn_sigma)-pnorm((al-tn_mu)/tn_sigma)

    if(etajTy>=ar){ #if test statistic larger than upper bound of truncation cdf = 1 and p = 0 will be returned
      numerator=numerator+pnorm((ar-tn_mu)/tn_sigma)-pnorm((al-tn_mu)/tn_sigma)
    }
    else if(etajTy>=al && etajTy<ar){
      numerator=numerator+pnorm((etajTy-tn_mu)/tn_sigma)-pnorm((al-tn_mu)/tn_sigma)
    }
  }

  if(denominator!=0){
    return(numerator/denominator)
  }
  else{
    return(Inf)
  }
}


#' @export postselp_value_specified_interval
postselp_value_specified_interval <- function(z_interval, etaj, etajTy, tn_mu, tn_sigma){

  value=pivot_with_specified_interval(z_interval, etaj, etajTy, tn_mu, tn_sigma)

  return(2*min(1-value,value))
}


# DATA GENERATING #

#' @export datagen.norm
datagen.norm <- function(seed, n, p, rho, beta_vec){
  set.seed(seed)
  X=rmnorm(n=n,mean=rep(0,p), varcov=rho*(matrix(1,p,p)-diag(p))+diag(p))
  true_y=t(beta_vec%*%t(X))
  noise=rnorm(n,mean=0,sd=1)
  y=t(beta_vec%*%t(X))+noise
  return(list('X'=X,'y'=y,'true_y'=true_y))
}


#' @export datagen.norm.intercept
datagen.norm.intercept <- function(seed, n, p, rho, beta_vec){
  set.seed(seed)
  q=length(beta_vec)
  X=cbind(rep(1,n),rmnorm(n=n,mean=rep(0,q-1), varcov=rho*(matrix(1,q-1,q-1)-diag(q-1))+diag(q-1)))
  true_y=t(beta_vec%*%t(X))
  noise=rnorm(n,mean=0,sd=1)
  y=t(beta_vec%*%t(X))+noise
  return(list('X'=X,'y'=y,'true_y'=true_y))
}


# CONFIDENCE INTERVALS #

#' @export f
f <- function(z_interval,etajTy,mu,tn_sigma){

  numerator=0
  denominator=0

  for(i in 1:length(z_interval)){

    al=z_interval[[i]][1]
    ar=z_interval[[i]][2]

    denominator=denominator+pnorm((ar-mu)/tn_sigma)-pnorm((al-mu)/tn_sigma)

    if(etajTy>=ar){ #if test statistic larger than upper bound of truncation, CDF should be 1 of course
      numerator=numerator+pnorm((ar-mu)/tn_sigma)-pnorm((al-mu)/tn_sigma)
    }
    else if(etajTy>=al & etajTy<ar){
      numerator=numerator+pnorm((etajTy-mu)/tn_sigma)-pnorm((al-mu)/tn_sigma)
    }
  }

  if(denominator!=0){
    return(numerator/denominator)
  }
  else{
    return(Inf)
  }
}


#' @export find_root
find_root <- function(z_interval,etajTy,tn_sigma,y,lb,ub,tol=1e-6){

  a=lb
  b=ub
  fa=f(z_interval,etajTy,a,tn_sigma) # truncated normal CDF with mean lb
  fb=f(z_interval,etajTy,b,tn_sigma) # truncated normal CDF with mean ub

  # assume a < b
  if(fa>y & fb>y){
    while(fb>y){
      b=b+(b-a)
      fb=f(z_interval,etajTy,b,tn_sigma)
    }
  }
  else if (fa<y & fb<y){
    while(fa<y){
      a=a-(b-a)
      fa=f(z_interval,etajTy,a,tn_sigma)
    }
  }

  max_iter=ceiling((log(tol) - log(b-a)) / log(0.5))

  for(i in 1:range(max_iter)[1]){
    c=(a+b)/2
    fc=f(z_interval,etajTy,c,tn_sigma)
    if(fc>y){
      a=c
    }
    else if(fc<y){
      b=c
    }
  }
  return(c)
}


#' @export equal_tailed_interval
equal_tailed_interval <- function(z_interval,etajTy,alpha,tn_mu,tn_sigma){

  lb=tn_mu-20*tn_sigma
  ub=tn_mu+20*tn_sigma

  L=find_root(z_interval,etajTy,tn_sigma,1.0-0.5*alpha,lb,ub)
  U=find_root(z_interval,etajTy,tn_sigma,0.5*alpha,lb,ub)

  return(list('L'=L,'U'=U))
}


#' @export compute_ci_with_specified_interval
compute_ci_with_specified_interval <- function(z_interval,etaj,etajTy,Sigma,tn_mu,alpha){

  tn_sigma=sqrt((t(etaj)%*%Sigma)%*%etaj)

  ci=equal_tailed_interval(z_interval,etajTy,alpha,tn_mu,tn_sigma)

  return(as.numeric(ci))
}
