#' gibbs_cox
#'
#' Function to perform factorized regression with the Cox Proportional Hazards Model
#' @import magrittr

#' @param y vector of responses
#' @param X matrix of predictors
#' @param nu vector of censoring indicators
#' @param X_test test set matrix (if NULL, this is X)
#' @param nrun number of iterations of gibbs sampler
#' @param burn number of initial samples to discard
#' @param thin thinning amount
#' @param delta_eta eta initial acceptance ratio
#' @param delta_omega omega initial acceptance ratio
#' @param k desired number of factors
#' @param verbose description (not sure)
#'
#' @return Matrix (Lambda) of factor loadings
#' @export
#'

gibbs_cox <- function(y, X, nu, X_test = NULL, nrun = 1000, burn = 500, thin = 1,
                     delta_eta = 0.2, delta_omega = 0.2, k = NULL, verbose = T){

  n = nrow(X)

  if(is.null(k)) k = floor(log(ncol(X))*3)
  if(is.null(X_test)) X_test <- X
  if(length(y) != n) stop("Mismatching input lengths")
  at = ceiling(nrun/100)
  if(verbose) {
    pb = utils::txtProgressBar(style = 3)
  }


  n_test <- nrow(X_test)
  p <- ncol(X)                    # collect data attributes
  y_sort <- sort(y, index.return=TRUE)
  y_ord <- y_sort$x; y_ind <- y_sort$ix



  sp <- floor((nrun - burn)/thin)

  # hyperparameters
  as <- 1; bs <- 1

  # initial values
  eta <- matrix(0, n, k)
  eta_test <- matrix(0, n_test, k)

  ps = stats::rgamma(p,as,bs)
  Sigma = diag(1/ps)
  Lambda = matrix(0,p,k)
  Plam = matrix(1,p,k)
  omega = numeric(k)
  acp_eta_st <- list()
  acp_omega_st <- list()
  eta_omega <- eta%*%omega
  # alpha = 0         # intercept

  # storage
  omega_st = matrix(0,nrow = sp, ncol = k)
  beta_st = matrix(0,nrow = sp, ncol = p)   # induced main effects
  Lambda_st = list()
  eta_st = list()
  ps_st = matrix(0, nrow = sp, ncol = p)
  acp_eta = numeric(n)
  acp_omega = numeric(k)

  count = 1
  count_100 = 1
  for(i in 1:nrun){


    # --- Update eta --- #
    eta_omega <- eta%*%omega %>% c()
    for (h in 1:n){                # Metropolis hastings step

      # proposal distribution
      eta_star <- eta
      eta_star[h,] <- bayesSurv::rMVNorm(1, eta[h,], diag(k)*delta_eta)

      cox_term <- 0
      etastar_omega <- eta_star%*%omega %>% c()
      for(l in 1:n){
        # https://en.wikipedia.org/wiki/Proportional_hazards_model
        # this can be made much more efficient
        # min(l+1, n) is to avoid the case (n+1):n
        if(nu[y_ind[l]] != 0){
          cox_term <- cox_term + t(omega)%*%eta_star[y_ind[l],] - t(omega)%*%eta[y_ind[l],] -
            log( sum( exp( etastar_omega[y_ind[min(l+1,n):n]] ) ) ) + log( sum( exp( eta_omega[y_ind[min(l+1,n):n]] ) ) )
        }
      }

      # MH ratio: likelihood + prior (both logged)
      logr <- cox_term +
        mvtnorm::dmvnorm((X[h,] - Lambda%*%eta_star[h,]) %>% as.vector(), sigma = Sigma, log = T) -
        mvtnorm::dmvnorm((X[h,] - Lambda%*%eta[h,]) %>% as.vector(), sigma = Sigma, log = T) +
        mvtnorm::dmvnorm(eta_star[h,], log = T) -
        mvtnorm::dmvnorm(eta[h,], log = T)

      logu = stats::runif(1) %>% log()

      if (logr > logu){

        eta = eta_star
        acp_eta[h] = acp_eta[h] + 1
        eta_omega <- etastar_omega

      }
    }

    # --- Update omega --- #
    # modify here as well, first two lines
    for (j in 1:k){                # Metropolis hastings step

      omega_star <- omega
      omega_star_j <- stats::rnorm(1, omega[j], delta_omega)
      omega_star[j] <- omega_star_j

      cox_term <- 0
      eta_omegastar <- eta%*%omega_star %>% c()
      for(l in 1:n){
        # https://en.wikipedia.org/wiki/Proportional_hazards_model
        # this can be made much more efficient
        # min(l+1, n) is to avoid the case (n+1):n
        if(nu[y_ind[l]] != 0){
          cox_term <- cox_term + t(omega_star)%*%eta[y_ind[l],] - t(omega)%*%eta[y_ind[l],] -
            log( sum( exp( eta_omegastar[y_ind[min(l+1,n):n]] ) ) ) + log( sum( exp( eta_omega[y_ind[min(l+1,n):n]] ) ) )
        }
      }

      logr = cox_term -
        mvtnorm::dmvnorm(omega_star, log = T) -
        mvtnorm::dmvnorm(omega, log = T)

      logu = stats::runif(1) %>% log()

      if (logr > logu){

        eta_omega <- eta_omegastar
        omega[j] = omega_star_j
        acp_omega[j] = acp_omega[j] + 1

      }
    }


    # --- Update Lambda --- #
    eta2 = t(eta)%*%eta
    zlams = stats::rnorm(k*p)       # generate normal draws all at once
    Lambda.T = t(Lambda)

    for(j in 1:p) {
      Llamt = chol(diag(Plam[j,]) + ps[j]*eta2)
      Lambda[j,] = t(solve(Llamt,
                           zlams[1:k + (j-1)*k]) +
                       solve(Llamt,
                             solve(t(Llamt),
                                   ps[j] * t(eta) %*% X[,j])))
    }

    # --- Update Sigma --- #
    Xtil = X - eta%*%t(Lambda)
    ps = stats::rgamma(p,as+0.5*n,1)
    ps = (1/(bs + 0.5*apply(Xtil^2,2,sum)))*ps
    Sigma = diag(1/ps)


    #store posterior samples
    if ((i %% thin == 0) & (i > burn)){

      Sigma_inv <- diag(ps)
      V_n = solve(t(Lambda)%*%Sigma_inv%*%Lambda + diag(rep(1,ncol(Lambda))))
      A_n = V_n%*%t(Lambda)%*%Sigma_inv
      beta_st[count,] = omega%*%A_n
      # Lambda_st[count,,] = dsVX %*% Lambda
      Lambda_st[[count]] = Lambda
      eta_st[[count]] = eta
      ps_st[count,] = ps

      count = count + 1

    }

    if(verbose & (i %% at == 0)) utils::setTxtProgressBar(pb, i / nrun)

    if (i%%100==0){

      # print(i)

      # adapt deltas
      acp_eta_mean = mean(acp_eta)/100
      acp_omega_mean = mean(acp_omega)/100

      if(acp_eta_mean > 0.3){
        delta_eta = delta_eta*2
      }else if(acp_eta_mean < 0.2){
        delta_eta = delta_eta*2/3
      }
      acp_eta_st[[count_100]] <- acp_eta_mean
      acp_eta = numeric(n)

      if(acp_omega_mean > 0.3){
        delta_omega = delta_omega*2
      }else if(acp_omega_mean < 0.2){
        delta_omega = delta_omega*2/3
      }
      acp_omega_st[[count_100]] <- acp_omega_mean
      acp_omega = numeric(k)

      count_100 <- count_100 + 1

    }
  }

  return(list(Lambda = Lambda_st,
              ps = ps_st,
              acp_eta = acp_eta_st,
              acp_omega = acp_omega_st,
              beta = beta_st,
              eta = eta_st))
}
