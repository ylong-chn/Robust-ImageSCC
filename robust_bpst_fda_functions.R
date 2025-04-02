#' Function library for robust FDA with BPST 
#' 
#' 0. Preliminary functions 
#' 0.0.1 signal functions
#' @param z coordinates of pixels, must be a npix by 2 matrix
#' @param type 'quadratic', 'exponential', 'cubic', or 'sine'
true_mean <- function(z, type) {
  if (type == "quadratic") 20*((z[,1]-0.5)^2+(z[,2]-0.5)^2)
  else if (type == "exponential") 5*exp(-15*((z[,1]-0.5)^2+(z[,2]-0.5)^2))+0.5
  else if (type == "cubic") 3.2*(-z[,1]^3+z[,2]^3)+2.4
  else if (type == "sine") 3*sin(5*pi*(z[,1]+0.22))+3*sin(5*pi*(z[,2]-0.18))+4.8
  # else if (type == "sine") 3/2*sin(3*pi*(z[,1]+0.22))+3/2*sin(3*pi*(z[,2]-0.18))+4.8
  else if (type == "rdnn") (-8)/(1+exp(cot(z[,1]^2)*cos(2*pi*z[,2])))
  else if (type == "horseshoe") mgcv::fs.test(z[,1], z[,2], exclude = FALSE)
}
#' 0.0.2 eigenfunctions
psi_1 <- function(z) 0.9878812*sin(pi*z[,1])+0.5
psi_2 <- function(z) 2.157205*(cos(pi*z[,2])-0.03903665)
#' 0.0.3 generate grid of pixel location
#' @param n1,n2 desired resolution of generated images
generate_grid <- function(n1, n2) expand.grid(seq(0, 1, length=n1), seq(0, 1, length=n2)) %>% arrange_all() %>% as.matrix()
#' 0.0.4 generate xi 
#' @param n number of images to generate 
#' @param lam1,lam2 eigenvalues to adjust subject-level variation 
generate_xi <- function(n, lam1, lam2) cbind(rnorm(n, mean=0, sd=sqrt(lam1)), rnorm(n, mean=0, sd=sqrt(lam2)))
#' 0.0.5 generate within-image dependence 
#' @param type type of within-image dependence, either use 'bpst' from Wang (2019) or 'rdnn' from Wang and Cao (2022)
#' @param z coordinates of pixels, must be a npix by 2 matrix
#' @param xi vector of 2, random number following normal distribution
generate_wi_dep <- function(type, z, xi) {
  if (type == 'bpst') psi_1(z)*xi[1] + psi_2(z)*xi[2]
  else if (type == 'rdnn') {
    n1 <- length(unique(z[,1]))
    n2 <- length(unique(z[,2]))
    npix <- n1 * n2
    mvrnorm(1, rep(0, npix), cos(2*pi*outer(z[,1], z[,1], "-")) + cos(2*pi*outer(z[,2], z[,2], "-")))
  }
}
#'
#' 1. Create simulation sample, with various types of contamination 
#' @param z coordinates of pixels, must be a npix by 2 matrix
#' @param ind_inside a vector of indices for pixels inside domain
#' @param n number of images to generate 
#' @param type type of true functions, see `true_mean()` for options
#' @param error_type within-image dependence, either use 'bpst' from Wang (2019) or 'rdnn' from Wang and Cao (2022)
#' @param lam1,lam2 eigenvalues to adjust subject-level variation 
#' @param con_rate contamination rate, only applies to 'stripe' and 'square' cont. types
#' @param con_type contamination type, 'stripe', 'square', 'cauchy', 'slash'
#' @param clean_sd standard deviation to generate clean observations
#' @param a_0 parameters specifying location of 'stripe' contamination (can be a vector of values to have multiple stripes)
#' @param a_1,a_2 parameter specifying location of 'square' contamination
#' @param mix_weight weight of the fat-tailed distributions for mix-dist. types of contam.
#' @param seed random seed
generate_sample <- function(z, ind_inside, n, type, 
                            error_type='bpst', lam1, lam2, 
                            con_rate=0, con_type=NULL,
                            clean_sd=0.2,
                            a_0=0.3, 
                            # segment=FALSE,
                            a_1=0.1, a_2=0.5,
                            a_square=matrix(c(0.1, 0.2, 0.8, 0.9), nrow=2, byrow=TRUE),
                            rand_cont_ct=10,
                            cont_dist_min=10,
                            cont_dist_max=20,
                            mix_weight=0.5,
                            seed=NULL) {
  
  if (!is.null(seed)) set.seed(seed=seed)
  
  # location info
  n1 <- length(unique(z[,1]))
  n2 <- length(unique(z[,2]))
  npix <- n1 * n2
  
  # generate true image data 
  beta_true <- true_mean(z=z, type=type)
  ind_outside <- setdiff(1:npix, ind_inside)
  beta_true[ind_outside] <- NA
  
  # generate xi matrix for within image variations 
  xi_mat <- generate_xi(n=n, lam1=lam1, lam2=lam2)
  eps_mat <- matrix(rnorm(n*npix, mean=0, sd=clean_sd), nrow=n)
  
  # index of realizations without outliers 
  ind_conta <- NULL
  if (con_rate > 0) {
    # if the contamination is distributional, all images have the same distributed error
    dist_flag <- con_type == "cauchy" || con_type == "slash"
    if (dist_flag) ind_conta <- seq(n)
    else ind_conta <- sample(n, n*con_rate)
    ind_clean <- c(1:n)[!1:n %in% ind_conta]
  } else ind_clean <- 1:n
  
  # clean observations are needed for clean sample and stripe and square type contaminated samples
  if (length(ind_clean) > 0) {
    clean_sample <- lapply(ind_clean, function(i) {
      # realizations
      beta_true + generate_wi_dep(type=error_type, z=z, xi=xi_mat[i,]) + eps_mat[i,]
      
    })
  } else clean_sample <- c()
  
  # add errors for comparison 
  
  if (con_rate > 0) {
    
    # stripe type of contamination
    if (con_type == "stripe") {
      if (is.null(a_0)) stop('Value of a_0 is missing')
      
      # stripe location
      con_n1 <- floor(a_0*n1) 
      z_cont <- as.vector(sapply(con_n1, function(i) ((i-1)*n2+1):(i*n2)))
      # # if only segments, instead of entire stripe is contaminated
      # if (segment) {
      #   # segment complies with the paper
      #   partition <- split(n2_levels, cut(unique(z[,2]), breaks=seq(0,1,by=1/10), include.lowest = TRUE, right= FALSE))
      #   z_cont <- z_cont[z[z_cont,2] %in% do.call(c, partition[c(TRUE, FALSE)])]
      # }
      # generate sample 
      contaminated_sample <- lapply(ind_conta, function(i) {
        # vector of epsilon
        # overwrites the entries epsilon matrix with the outliers 
        eps_mat[i, z_cont] <- runif(length(z_cont), min=cont_dist_min, max=cont_dist_max)
        # generate images
        beta_true + generate_wi_dep(type=error_type, z=z, xi=xi_mat[i,]) + eps_mat[i,]
      })
    }
    
    
    
    # square type of contamination
    if (con_type == "square") {
      if (is.null(a_1) || is.null(a_2)) stop('Value of a_1 or a_2 is missing')
      
      # square location
      z1_range <- diff(range(z[,1]))
      z2_range <- diff(range(z[,2]))
      z_cont <- c(1:npix)[as.logical(between(z[,1], a_1 * z1_range + min(z[,1]), a_2 * z1_range + min(z[,1])) * 
                                       between(z[,2], a_1 * z2_range + min(z[,2]), a_2 * z2_range + min(z[,2])))]
      # generate sample 
      contaminated_sample <- lapply(ind_conta, function(i) {
        # vector of epsilon
        eps_mat[i, z_cont] <- runif(length(z_cont), min=cont_dist_min, max=cont_dist_max)
        # generate images
        beta_true + generate_wi_dep(type=error_type, z=z, xi=xi_mat[i,]) + eps_mat[i,]
      })
    }
    
    # additional square type of contamination
    if (con_type == "square2") {
      if (is.null(a_square)) stop('Value of a_square is missing')
      
      # square location
      z1_range <- diff(range(z[,1]))
      z2_range <- diff(range(z[,2]))
      
      ind_cont <- rep(FALSE, npix)
      for (i in 1:nrow(a_square)) {
        a_1 <- a_square[i, 1]
        a_2 <- a_square[i, 2]
        ind_cont <- ind_cont +
          between(z[,1], a_1 * z1_range + min(z[,1]), a_2 * z1_range + min(z[,1])) * 
          between(z[,2], a_1 * z2_range + min(z[,2]), a_2 * z2_range + min(z[,2]))
      }
      z_cont <- c(1:npix)[as.logical(ind_cont)]
      # generate sample 
      contaminated_sample <- lapply(ind_conta, function(i) {
        # vector of epsilon
        eps_mat[i, z_cont] <- runif(length(z_cont), min=cont_dist_min, max=cont_dist_max)
        # generate images
        beta_true + generate_wi_dep(type=error_type, z=z, xi=xi_mat[i,]) + eps_mat[i,]
      })
    }
    
    # square type of contamination
    if (con_type == "random") {
      
      z_cont <- sample(ind_inside, rand_cont_ct)
      # generate sample 
      contaminated_sample <- lapply(ind_conta, function(i) {
        # vector of epsilon
        eps_mat[i, z_cont] <- runif(length(z_cont), min=cont_dist_min, max=cont_dist_max)
        # generate images
        beta_true + generate_wi_dep(type=error_type, z=z, xi=xi_mat[i,]) + eps_mat[i,]
      })
    }
    
    
    
    # contamination with Cauchy error
    if (con_type == "cauchy") {
      # generate sample 
      contaminated_sample <- lapply(ind_conta, function(i) {
        # vector of epsilon
        eps_mat[i, ] <- mix_weight*rcauchy(npix, location = 0, scale = 0.5)+(1-mix_weight)*rnorm(npix, mean=0, sd=clean_sd)
        # generate images
        beta_true + generate_wi_dep(type=error_type, z=z, xi=xi_mat[i,]) + eps_mat[i, ]
      })
    }
    
    
    
    # contamination with Slash error
    if (con_type == "slash") {
      # generate sample 
      contaminated_sample <- lapply(ind_conta, function(i) {
        # vector of epsilon
        eps_mat[i, ] <- mix_weight*rslash(npix, mu = 0, sigma = 0.5)+(1-mix_weight)*rnorm(npix, mean=0, sd=clean_sd)
        # generate images
        beta_true + generate_wi_dep(type=error_type, z=z, xi=xi_mat[i,]) + eps_mat[i, ]
      })
    }
  }
  
  if (con_rate > 0) {
    if (!dist_flag) {
      y <- do.call(rbind, c(clean_sample, contaminated_sample))
      y <- y[match(1:n,c(ind_clean, ind_conta)),]
    }
    else y <- do.call(rbind, contaminated_sample)
  }
  else y <- do.call(rbind, clean_sample)
  
  # output true function and generated sample
  list(true_func = beta_true,
       ind_cont = ind_conta,
       ind_clean = ind_clean,
       xi_mat = xi_mat,
       eps_mat = eps_mat,
       y = y)
}
#'
#' 2. M-estimator for mean function using B-spline smoothing
#' @param B matrix storing values of basis functions
#' @param ind_inside a vector of indices for pixels inside domain
#' @param Q2 Q2 matrix obtained from QR decomposition
#' @param K energy function a.k.a. penalty matrix 
#' @param Y vector of values of response variable
#' @param lambda vector of penalty parameters 
#' @param huber_const Huber tuning constant
#' @param maxit maximum number of iterations for IRWLS algorithm 
#' @param tol tolerance level for convergence 
#' @returns a list containing estimates of coefficient and mean function, and number of iterations
mean_func_m_est <- function(
    B, ind_inside, Q2, K, Y, lambda=10^(-5:2), 
    huber_const=1,
    maxit=100, tol=1e-6
) {
  
  # some preliminary matrix computation
  B_til <- B%*%Q2
  D <- t(Q2)%*%K%*%Q2
  J <- ncol(Q2)
  WW <- crossprod(B_til)
  
  Y_inside <- Y[,ind_inside]
  
  lambda <- lambda[order(lambda, decreasing = TRUE)]
  
  result <- MEstLambdaPath(
    as.matrix(B_til),
    t(Y_inside),
    as.matrix(D),
    lambda,
    huber_const,
    maxit,
    tol
  )
  
  W_final <- result$WeightMatFinal
  
  # output
  list(
    coefficients = result$Coeff,
    est_func_robust = drop(B_til%*%result$Coeff),
    iteration = result$Iteration,
    gcv_all = result$GCVall,
    gcv_final = result$GCV,
    lambda_gcv = result$LambdaCV,
    weight_mat_final = W_final
  )
  
}

#' weighted m estimation for weighted bootstrap 
#' @param sample_weight vector weights for the images (observations)
mean_func_m_est_weighted <- function(
    B_til, D,
    ind_inside, Y, 
    sample_weight,
    lambda=10^(-6:3), 
    huber_const=1,
    maxit=100, tol=1e-6
) {
  
  # some preliminary matrix computation
  Y_inside <- Y[,ind_inside]
  
  lambda <- lambda[order(lambda, decreasing = TRUE)]
  
  result <- MEstLambdaPathWeighted(
    as.matrix(B_til),
    t(Y_inside),
    as.matrix(D),
    sample_weight,
    lambda,
    huber_const,
    maxit,
    tol
  )
  
  W_final <- result$WeightMatFinal
  
  # output
  list(
    coefficients = result$Coeff,
    est_func_robust = drop(B_til%*%result$Coeff),
    iteration = result$Iteration,
    gcv_all = result$GCVall,
    gcv_final = result$GCV,
    lambda_gcv = result$LambdaCV,
    weight_mat_final = W_final
  )
  
}

#' Implementation of weighted bootstrap to construct RSCC
#' @param nboot number of bootstrap replications
#' @param alpha level of significance 
robust_imagescc <- function(
    B.est, ind.inside, Q2.est, K.est, Y, 
    lambda=lambda, huber_const=huber_const,
    nboot=1000, alpha=0.05,
    maxit=100, tol=1e-6
  ) {
  
  # mean function estimation 
  mfit0 <- mean_func_m_est(
    B.est, ind.inside, Q2.est, K.est, Y, lambda=lambda, 
    huber_const=huber_const,
    maxit = maxit, tol = tol
  )
  mu.hat <- mfit0$est_func_robust
  
  # compute pointwise variance 
  resid_mat <- t(Y_inside) - mu.hat
  delta <- rowMeans(abs(resid_mat) <= huber_const)
  
  lambda_cv <- mfit0$lambda_gcv
  N <- length(ind.inside)
  L_n <- n/N * t(B_til) %*% diag(delta) %*% B_til + lambda_cv/N * D
  L_n_inv <- solve(L_n)

  psi_mat_est <- robustbase::Mpsi(resid_mat, cc = huber_const, psi = 'huber', deriv = 0)
  G <- 1/N^2 * crossprod(t(psi_mat_est))
  LGL <- L_n_inv %*% (t(B_til)%*%G%*%B_til) %*% L_n_inv
  Sigma <- diag(B_til %*% LGL %*% t(B_til))
  
  # weighted bootstrap implementation 
  bootstrap_func <- c()
  for (ib in 1:n_b) {
    if (ib == 1 | ib %% 100 == 0) cat('Bootstrapping sample No.', ib, '... \n')
    
    # use exponential with parameter 1
    sample_wts <- rexp(n)
    
    mfit_weight_bs <- mean_func_m_est_weighted(
      B_til, D,
      ind.inside, Y, 
      sample_weight=sample_wts,
      lambda=lambda, 
      huber_const=huber_const,
      maxit = 1000, tol = 1e-6
    )
    
    bootstrap_func <- rbind(
      bootstrap_func,
      mfit_weight_bs$est_func_robust
    )
  }
  
  # compute the (1-alpha) quantile of the absolute maximal deviations 
  all_dev <- abs((t(bootstrap_func)-mu.hat)/sqrt(Sigma))
  max_dev <- apply(all_dev, 2, max)
  c_b <- quantile(max_dev, 1-alpha)
  
  # construct the upper and lower RSCC
  cc.l <- mu.hat - c_b * sqrt(Sigma)
  cc.u <- mu.hat + c_b * sqrt(Sigma)
  
  list(cc.l=cc.l, cc.u=cc.u, bw=mean(c_b * sqrt(Sigma)*2))
}

#' Modified `fit.image()` function from the GitHub repository `ImageSCC` (https://github.com/FIRST-Data-Lab/ImageSCC) package to take a given set of basis results 
fit.image <- function(Y,Z,V.est,Tr.est,d.est=5,r=1,lambda=10^(-6:3),proj.matrix=FALSE,Bfull=NULL, Y_inside=NULL){
  
  n <- nrow(Y)
  # Basis for estimating mu: Tr.est,V.est,d.est,r
  if(is.null(Bfull)) Bfull.est <- basis(V.est,Tr.est,d.est,r,Z)
  else Bfull.est <- Bfull
  B <- Bfull.est$B
  ind.inside <- Bfull.est$Ind.inside
  Q2 <- Bfull.est$Q2
  K <- Bfull.est$K
  if(is.null(Y_inside)) Y <- matrix(Y[,ind.inside],nrow=n)
  else {
    Y_inside_mean <- matrix(apply(Y_inside,2,mean),nrow=1)
    Y <- matrix(Y_inside_mean, nrow=n)
  }
  lambda <- as.matrix(lambda)
  
  this.call <- match.call()
  n <- nrow(Y)
  npix <- ncol(Y)
  J <- ncol(Q2)
  
  W <- as.matrix(B%*%Q2)
  WW <- crossprod(W,W)
  rhs <- crossprod(W,t(Y))
  D <- crossprod(t(crossprod(Q2,as.matrix(K))),Q2)
  D <- as.matrix(D)
  
  flag <- (rankMatrix(WW)<J)
  if(!flag){
    Ainv <- chol(WW,pivot=TRUE)
    A <- solve(t(Ainv))
    ADA <- A%*%D%*%t(A)
    eigs <- eigen(ADA)
    Cval <- eigs$values
  }
  
  nl <- length(lambda)
  gcv_all <- sapply(lambda,FUN=function(Lam){  # eventually a n*nl matrix
    Dlam <- Lam*D
    lhs <- WW+Dlam
    lhs.inv <- chol2inv(chol(lhs));
    theta <- crossprod(t(lhs.inv),rhs)
    gamma <- crossprod(t(Q2),theta)
    Yhat <- crossprod(t(W),theta)
    res <- t(Y)-Yhat
    sse <- apply(res^2,2,sum)
    if(!flag){
      df <- sum(1/(1+Cval*Lam))
    }
    if(flag){
      Hmtx <- crossprod(t(crossprod(t(W),lhs.inv)),t(W))
      df <- sum(diag(Hmtx))
    }
    gcv <- npix*sse/(npix-df)^2
  })
  gcv_all <- matrix(gcv_all,nrow=n)
  lam.ind <- apply(gcv_all,1,which.min)
  lambdac <- lambda[lam.ind]
  theta <- c()
  gamma <- c()
  Yhat <- c()
  df <- c()
  for (i in 1:n){
    lamc.tmp <- lambdac[i]
    Dlam <- lamc.tmp*D
    lhs <- WW+Dlam
    lhs.inv <- chol2inv(chol(lhs));
    rhs.tmp <- as.matrix(rhs[,i],ncol=1)
    theta.tmp <- crossprod(t(lhs.inv),rhs.tmp)
    theta <- cbind(theta,theta.tmp)
    gamma.tmp <- crossprod(t(Q2),theta.tmp)
    gamma <- cbind(gamma,gamma.tmp)
    Yhat.tmp <- crossprod(t(W),theta.tmp)
    Yhat <- cbind(Yhat,Yhat.tmp)
    if(!flag){
      df.tmp <- sum(1/(1+Cval*lamc.tmp))
    }
    if(flag){
      Hmtx <- crossprod(t(crossprod(t(W),lhs.inv)),t(W))
      df.tmp <- sum(diag(Hmtx))
    }
    df <- c(df,df.tmp)
  }
  if(proj.matrix==TRUE){   # projection matrix
    if(n==1){
      Dlam <- lambdac*D
      lhs <- WW+Dlam
      lhs.inv <- chol2inv(chol(lhs))
      Slam <- crossprod(t(W),lhs.inv%*%t(W))
    }
    if(n>1){
      Slam <- lapply(lambdac,FUN=function(lamc){
        Dlam <- lamc*D
        lhs <- WW+Dlam
        lhs.inv <- chol2inv(chol(lhs))
        Slam <- crossprod(t(W),lhs.inv%*%t(W))
        return(Slam)
      })
    }
  }else{Slam <- NA}
  Yhat <- t(Yhat)
  res <- Y-Yhat
  sse <- apply(res^2,1,sum)
  gcv <- npix*sse/(npix-df)^2
  bic <- log(sse/npix)+df*log(npix)/npix
  mfit <- list()
  mfit$theta <- theta; mfit$gamma <- gamma;
  mfit$Yhat <- Yhat; mfit$res <- res; mfit$sse <- sse;
  mfit$gcv <- gcv; mfit$bic <- bic;
  mfit$lambdac <- lambdac; mfit$Slam <- Slam;
  mfit$scc <- NULL;
  mfit$Y <- Y; mfit$Z <- Z; mfit$V.est.a <- V.est; mfit$Tr.est.a <- Tr.est;
  mfit$V.est.b <- NULL; mfit$Tr.est.b <- NULL;
  mfit$d.est <- d.est; mfit$r <- r; mfit$ind.inside <- ind.inside;
  mfit$call <- this.call;
  class(mfit) <- "image"
  return(mfit)
}

#' The implementing function of the RDNN method, modified from the GitHub repository `FDADNN` (https://github.com/FDASTATAUBURN/FDADNN) to only take the points inside the domain as the input
rFDADNN=function(Data=NULL, Y_ins=NULL, z_ins=NULL, d, Grid=NULL, N, n, L, p, s, epoch, batch, huber_const=NULL, quantile=NULL){
  
  if (is.null(z_ins) && is.null(Y_ins)) {
    #Process grid points
    if(d==1){
      x_train=Grid[[1]]
      x_train.all=rep(x_train, n)
    }else if(d==2){
      x1=rep(rep(Grid[[1]],N[2]))
      x2=rep(Grid[[2]],each=N[1])
      x_train=cbind(x1,x2)
      x1.all=rep(x1, n)
      x2.all=rep(x2, n)
      x_train.all=cbind(x1.all, x2.all)
    }else if(d==3){
      x1=rep(Grid[[1]],N[2]*N[3])
      x2=rep(rep(Grid[[2]],each=N[1]),N[3])
      x3=rep(Grid[[3]],each=N[1]*N[2])
      x_train=cbind(x1,x2,x3)
      x1.all=rep(x1, n)
      x2.all=rep(x2, n)
      x3.all=rep(x3, n)
      x_train.all=cbind(x1.all, x2.all, x3.all)
    }else if(d==4){
      x1=rep(rep(Grid[[1]],N[2]*N[3]*N[4]))
      x2=rep(rep(Grid[[2]],each=N[1]),N[3]*N[4])
      x3=rep(rep(Grid[[3]],each=N[1]*N[2]), N[4])
      x4=rep(Grid[[4]],each=N[1]*N[2]*N[3])
      x_train=cbind(x1,x2,x3, x4)
      x1.all=rep(x1, n)
      x2.all=rep(x2, n)
      x3.all=rep(x3, n)
      x4.all=rep(x4, n)
      x_train.all=cbind(x1.all, x2.all, x3.all, x4.all)
    }
    y_train.raw=matrix(NA, n, base::prod(N))
    for(i in 1:n){
      y_train.raw[i, ]=as.vector(Data[[i]])
    }
    y_train=as.vector(t(y_train.raw))
  }
  else {
    x_train <- z_ins
    x_train.all <- do.call(rbind, replicate(n, x_train, simplify = FALSE))
    y_train <- as.vector(t(Y_ins))
  }
  
  if (is.null(huber_const)) huber_const = 1
  
  model <- keras_model_sequential() 
  
  model %>% 
    layer_dense(units=p, activation = "relu", input_shape = c(d), kernel_initializer = "normal")%>%
    layer_dropout(rate = s) 
  for(xx in 1:L){
    model %>% 
      layer_dense(units=p, activation = "relu", kernel_initializer = "normal")%>%
      layer_dropout(rate = s) 
  }
  model %>% layer_dense(units=1)
  
  model %>% compile(
    loss = loss_huber(delta=huber_const),
    optimizer = optimizer_adam(),
    metrics = list("mse")
  )
  
  
  
  history <- model %>% fit(
    x_train.all, y_train, 
    epochs = epoch, batch_size = batch,
    verbose=0
  )
  y.reg=model %>% stats::predict(x_train)
  y.reg.all=model %>% stats::predict(x_train.all)
  # if(d==1){
  #   estimation=array(y.reg, N)
  # }else if(d==2){
  #   estimation=array(y.reg, c(N[1], N[2]))
  # }else if(d==3){
  #   estimation=array(y.reg, c(N[1], N[2], N[3]))
  # }else if(d==4){
  #   estimation=array(y.reg, c(N[1], N[2], N[3], N[4]))
  # }
  pse=mean((y.reg.all-y_train)^2)
  # if(d==1){
  #   plot.dnn=plot(x_train, y.reg, main="DNN estimation of 1d curve", type="l", xlab="Grids", ylab="Estimation")
  # }else if(d==2){
  #   plot.dnn=plotly::plot_ly(
  #     x = Grid[[1]], 
  #     y = Grid[[2]], 
  #     z = matrix(y.reg, N[1], N[2]), 
  #     type = "contour",
  #     colorbar=list(tickfont=list(size=18)))
  # }else if(d==3){
  #   plot.dnn=plot3D::scatter3D(x_train[,1],x_train[,2],x_train[,3],colvar = y.reg, pch=16, phi=0, theta=30, colkey=FALSE, xlab="x1", ylab="x2", zlab="x3")
  # }
  # list(plot.dnn=plot.dnn, estimation=estimation, pse=pse)
  list(estimation=y.reg, pse=pse)
}
