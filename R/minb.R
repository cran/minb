minb <- function(X,y,pars_init=NULL,lambda1_set=NULL,lambda2_set=NULL,
                 ntune=10,maxiter=200,tol=1e-03,vrbs=FALSE){

  n<-dim(X)[1]
  X<-cbind(rep(1,n),X)
  p<-dim(X)[2]

  if(!is.null(pars_init)){
    kappa.init<-pars_init$kappa
    omega.init<-pars_init$omega
    beta.init<-pars_init$beta
    phi.init<-pars_init$phi
    rho1<-1.0/omega.init
    rho2<-1.0/abs(beta.ini)
    rho2[1]<-0
    if(length(kappa.init)!=length(omega.init))
      stop("the number of 'kappa' should be equal with the number of 'omega'")
    if(length(beta)!=p)
      stop("error with the length of starting 'beta'")
    if(is.null(phi))
      stop("value of 'phi' is not supplied ")
  }else{
    kappa.ini<-sort(unique(y))
    omega.ini<-sapply(kappa.ini, function(x) sum(y == x))/n
    omega.ini<-omega.ini - min(omega.ini) / 2
    omega.ini[which(omega.ini>0.5)]<-0.5
    rho1<-1.0/abs(omega.ini)

    nb.fit<-NULL
    tryCatch({ nb.fit <- glm.nb(y~., data = data.frame(X[,-1],y),
                               control=glm.control(maxit=500))
    nb <- list(beta = nb.fit$coefficients,
              phi = 1.0/nb.fit$theta)
    beta.ini<-nb$beta
    phi.ini<-nb$phi
    if(phi.ini<1e-5) phi.ini<-1e-5
    rho2 <- 1.0/abs(beta.ini)
    rho2[1] <- 0
    },
    error=function(e){print('fitting NB distribution for initialization failed, using ZINB distribution instead')},
    finally={})

    if(is.null(nb.fit)){
      zinb <- pscl::zeroinfl(y ~ . | 1, data = data.frame(X[,-1],y),
                             dist = "negbin",link = "probit",
                             control = zeroinfl.control(method = "CG"))
      beta.ini <- zinb$coefficients$count
      phi.ini <- 1.0/zinb$theta
      rho2 <- 1.0/abs(beta.ini)
      rho2[1] <- 0 } }

  names(omega.ini) <- names(rho1) <- paste0("omega", kappa.ini)
  names(beta.ini) <- names(rho2) <- paste0("beta", 0:(p-1))

  init_set <- list(beta = beta.ini,
                  phi = phi.ini,
                  omega = omega.ini,
                  kappa = kappa.ini)

  tune_n <- 0.3*n/1000
  if(is.null(lambda1_set)){
    lambda1_set <- c(seq(from=1,to=3,length.out=ntune))
  }
  if(is.null(lambda2_set)){
    lambda2_set <- c(seq(from=1+tune_n,to=3+tune_n,length.out=ntune))
  }

  lll<-1
  l1<-length(lambda1_set)
  l2<-length(lambda2_set)
  para<-matrix(0,l1*l2,3)
  for (m1 in 1:l1){
    for (m2 in 1:l2){
      para[lll,1]<-lll
      para[lll,2]<-lambda1_set[m1]
      para[lll,3]<-lambda2_set[m2]
      lll<-lll+1
    }
  }

  bic<-matrix(1e+20,l1*l2,1)
  result_temp<-list()
  lll<-1
  for (m1 in 1:l1){
    if(vrbs) print(paste0('The ',m1,'th lambda1'))
    lambda1<-lambda1_set[m1]
    for (m2 in 1:l2){
      if(vrbs) print(paste0('The ',m2,'th lambda2'))
      lambda2<-lambda2_set[m2]
      tryCatch({ temp<-minb_onetune(X,y,init_set,lambda1,lambda2,rho1,rho2,maxiter=maxiter)
      beta<-temp$beta
      phi<-temp$phi
      omega<-temp$omega
      kappa<-temp$kappa
      gamma<-gammas(X,y,beta,phi,omega,kappa)
      result_temp[[lll]]<-temp
      lik<-L(X,y,beta,phi,omega,kappa,gamma)
      bic[lll] <- -2*lik+(sum(beta!=0)+sum(omega!=0))*log(n)
      },
      error=function(e){print('MINB one_tune failed')},
      finally={})
      lll<-lll+1
    }
  }
  id<-which.min(bic)
  output<-result_temp[[id]]

  beta<-output$beta
  indicator<-which(beta!=0)
  if(length(indicator)>1){
    X1<-X[,indicator]
    beta<-beta[indicator]
    phi<-output$phi
    omega<-output$omega
    kappa<-output$kappa
    pars<-list(beta=beta,phi=phi,omega=omega,kappa=kappa)
    output<-minb_refit(X=X1,y=y,pars=pars,maxiter=maxiter)
  }


  return(output)
}



minb_onetune <- function(X,y,init,lambda1,lambda2,rho1,rho2,maxiter=200,tol=1e-03){

  beta <- init$beta
  phi <- init$phi
  omega <- init$omega
  kappa <- init$kappa

  ui <- rbind(c(1),c(-1))
  ci <- c(1e-5,-100)
  diff_a <- i <- 1

  while(diff_a > tol && i< maxiter){
    tmpbeta <- beta
    tmpphi <- phi
    tmpomega <- omega
    gamma <- gammas(X,y,beta,phi,omega,kappa)
    beta <- beta_upd(X,y,beta,phi,lambda2,rho2,gamma,maxiter)
    try(phi <- constrOptim(phi,f,g,ui,ci,method='BFGS',
                           X=X,y=y,beta=beta,omega=omega,gamma=gamma)$par,silent = TRUE)
    omega <- omega_upd(X,y,omega,kappa,lambda1,rho1,gamma,maxiter)
    if(sum(omega)>=1) omega<-omega-min(omega)/2

    diff_a <-max(max(abs(beta-tmpbeta)),max(abs(phi-tmpphi)),max(abs(omega-tmpomega)))
    b <- which(omega>0)
    kappa <- as.numeric(gsub('omega','',names(b)))
    omega <-  omega[b]
    rho1 <- rho1[b]

    i <- i + 1
  }

  omega[abs(omega)<1e-3]<-0
  b <- which(omega>0)
  kappa <- as.numeric(gsub('omega','',names(b)))
  omega <-  omega[b]
  beta[abs(beta)<1e-3]<-0

  output<-list(beta=beta,phi=phi,omega=omega,kappa=kappa)
}

minb_refit <- function(X,y,pars,maxiter=200,tol=1e-03){

  beta <- pars$beta
  phi <- pars$phi
  omega <- pars$omega
  kappa <- pars$kappa

  ui <- rbind(c(1),c(-1))
  ci <- c(1e-05,-100)
  diff_r <- i <- 1

  while(diff_r > tol && i< maxiter){
    tmpbeta <- beta
    tmpphi <- phi
    tmpomega <- omega
    gamma <- gammas(X,y,beta,phi,omega,kappa)
    beta <- beta_upd(X,y,beta,phi,lambda2=NULL,rho2=NULL,gamma,penalize=FALSE,maxiter)
    try(phi <- constrOptim(phi,f,g,ui,ci,method='BFGS',
                           X=X,y=y,beta=beta,omega=omega,gamma=gamma)$par,silent = TRUE)
    if(length(omega)){
      omega <- omega_upd(X,y,omega,kappa,lambda1=NULL,rho1=NULL,gamma,penalize=FALSE,maxiter)
      if(sum(omega)>=1) omega<-omega-min(omega)/2

      diff_r <-max(max(abs(beta-tmpbeta)),max(abs(phi-tmpphi)),max(abs(omega-tmpomega)))
    }else{ diff_r <-max(max(abs(beta-tmpbeta)),max(abs(phi-tmpphi))) }
    i <- i + 1
  }

  output<-list(beta=beta,phi=phi,omega=omega,kappa=kappa)
}
L<-function(X,y,beta,phi,omega,kappa,gamma){

  n <- dim(X)[1]
  miu <- exp(drop(X%*%beta))
  J <- length(omega)
  a <- lgamma(y+1.0/phi)-lgamma(y+1)-lgamma(1.0/phi)+y*log(phi*miu)-(y+1.0/phi)*log(1 + phi*miu)
  a <- gamma[,J+1]* (log(1-sum(omega)) + a)

  for(i in 1:n){
    temp <- which(y[i]==kappa)
    if(length(temp)) a[i] <- a[i] + gamma[i,temp]*log(omega[temp])
  }

  sum(a)
}

deltas<-function(y,kappa){
  J <- length(kappa)
  n <- length(y)
  delta <- matrix(0,n,length(kappa)+1)
  for(i in 1:n){
    delta[i, ifelse(sum(y[i]==kappa), which(y[i]==kappa), J+1)]  <- 1
  }
  delta
}

gammas<-function(X,y,beta,phi,omega,kappa){
  n <- length(y)
  J <- length(kappa)
  miu <- exp(drop(X%*%beta))
  delta <- deltas(y,kappa)
  gamma <- matrix(0,nrow=n,ncol=J+1)

  temp  <- lgamma(y+1.0/phi)-lgamma(y+1)-lgamma(1.0/phi)+y*log(phi*miu)-(y+1.0/phi)*log(1 + phi*miu)
  temp <- (1-sum(omega))*exp(temp)
  temp[which(temp==0)] <-1e-20
  gamma[,J+1] <- temp

  if(J>0){
    for( i in 1:n){
      gamma[i,1:J] <- delta[i,1:J]*omega
      gamma[i,] <- gamma[i,]/sum(gamma[i,])
    }}
  return(gamma)

}


g <- function(X,y,beta,phi,omega,gamma){
  J <- length(omega)
  miu <- exp(drop(X%*%beta))
  g <- digamma(y+1.0/phi)*(-1.0/phi^2) - digamma(1.0/phi)*(-1.0/phi^2) + y/phi +
    1.0/phi^2*log(1+miu*phi)-(y+1.0/phi)*(miu/(1+miu*phi))
  g <- -1.0*gamma[,J+1]*g
  return(sum(g))
}

f <-function(X,y,beta,phi,omega,gamma){

  J <- length(omega)
  miu <- exp(drop(X%*%beta))
  a <-   lgamma(y+1.0/phi)-lgamma(y+1)-lgamma(1.0/phi)+y*log(phi*miu)-(y+1.0/phi)*log(1 + phi*miu)
  a <- a + log(1-sum(omega))
  obj <- -1.0*sum(gamma[,J+1]*a)
  return(sum(obj))
}

beta_upd <- function(X,y,beta,phi,lambda2,rho2,gamma,penalize=TRUE,maxiter=200){
  n<-dim(X)[1]
  p<-dim(X)[2]
  z_vec <- gamma[,ncol(gamma)]
  beta_old <- beta
  iter_beta <- 0
  repeat{
    abs_sum <- 0
    for( j  in 1:p){
      mu_vec <- exp(drop(X%*%beta_old) )

      X[,j][X[,j]==0] <- 1e-04#0.5*(1 + 2/pi*atan(X[,j][X[,j]==0]^2/10^-4))
      A_vec <- (y - mu_vec)/(mu_vec *(1 + mu_vec * phi )) * mu_vec * X[,j]
      B_vec <-   (y - mu_vec)/(mu_vec * (1 + mu_vec * phi )) * mu_vec * X[,j]^2 + (-2 * phi * y * mu_vec + mu_vec^2 * phi - y)/
        (mu_vec^2 * (1 + mu_vec*phi)^2) *mu_vec^2 * X[,j]^2

      omega_vec <- B_vec
      tau_vec <- beta_old[j] - A_vec/B_vec

      beta_tilde <- sum( z_vec * omega_vec * tau_vec)/sum(z_vec * omega_vec)

      if(penalize){betanew <- (sum( z_vec * omega_vec * tau_vec) + lambda2 * rho2[j] *sign(beta_tilde))/sum(z_vec * omega_vec)
      betanew <- ifelse( (sign(beta_tilde) * (betanew)) >0,  betanew, 0 )

      abs_sum <- abs_sum + abs(beta_old[j]-betanew)
      beta_old[j] <- betanew
      }else{
        abs_sum <- abs_sum + abs(beta_old[j]-beta_tilde)
        beta_old[j]<-beta_tilde
      }
    }
    #print(beta_old)
    iter_beta <- iter_beta + 1
    if(iter_beta > maxiter){break}
    if( abs_sum < 10^(-3)){break}
  }
  return(beta_old)
}

omega_upd <- function(X,y,omega,kappa,lambda1,rho1,gamma,penalize=TRUE,maxiter=100){

  J <- length(kappa)
  diff_o<-i<-1
  gamma_hat <- apply(gamma,2,sum)
  if(penalize){
    while(diff_o>=1e-04 && i<maxiter)
    {
      omega_old <- omega
      omega <- c(lambda1*rho1*omega,0)-gamma_hat
      omega[which(omega*omega[J+1]<0)] <- 0
      omega <- omega/sum(omega)
      omega <- omega[-length(omega)]
      diff_o <- max(abs(omega_old-omega))
      i<-i+1
    }}else{
      name <- names(omega)
      omega <- gamma_hat/sum(gamma_hat)
      omega <- omega[-length(omega)]
      names(omega)<-name
    }

  return(omega)

}



