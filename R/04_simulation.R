# Tthree functions in this file: ---------------------

## - Simulate Dirichlet-Multinomial Regression Model
## - Simulate Dirichlet-Multinomial Tree Regression Model
## - Simulate DMLMbvs model Doubly Multivariate Linear Model

## 4.1 simulate_DM {{{---------------

#' Title: Simulate Dirichlet-Multinomial Regression Model
#'
#' @description # This function can be used to simulate DM data.
#' Taken from Wadsworth (2017)
#' An integrative Bayesian Dirichlet-Multinomial regression model for
#' the analysis of taxonomic abundances in microbiome data
#'
#' @param n_obs `integer` Number of observations
#' @param n_vars `integer` Number of variables
#' @param n_taxa `integer` Number of taxa
#' @param n_relevant_vars `integer` Number of variables that are important
#' @param n_relevant_taxa `integer` Number of important taxa
#' @param beta_min `numeric` The desired min beta value
#' @param beta_max `numeric` The desired max beta value
#' @param signoise `numeric` Default is 1.
#' @param n_reads_min `numeric` Default is 5000. Used to obtain a random subject from n_reads_min
#'                              to n_reads_max from the Dirichlet Multinomial Distribution using dirmult::simPop()
#' @param n_reads_max `numeric` Default is 10000. Used to obtain a random subject from n_reads_min
#'                              to n_reads_max from the Dirichlet Multinomial Distribution using dirmult::simPop()
#' @param theta0
#' @param rho `numeric` Used for covariate correlation between 0 & 1. Default is NULL.
#'                      Need either rho or Sigma to be input.
#' @param Sigma `matrix` A symmetric positive definite matrix used for covariate
#'                       correlation structure columns must match number of covariates.
#'                       Default is NULL. Need either rho or Sigma to be input.
#'
#' @return A list of outcomes including:
#'   - X = XX
#'   - Y = YY
#'   - alphas = intercept,
#'   - betas = betas,
#'   - n_reads_min = n_reads_min,
#'   - n_reads_max = n_reads_max,
#'   - theta0 = theta0,
#'   - phi = phi,
#'   - rho = rho,
#'   - signoise = signoise,
#'   - Sigma = Sigma
#'
#' @export

simulate_DM <- function(n_obs = 100,
                        n_vars = 30,
                        n_taxa = 75,
                        n_relevant_vars = 4,
                        n_relevant_taxa = 4,
                        beta_min = 1,
                        beta_max = 1.25,
                        signoise = 1.0,
                        n_reads_min = 5000,
                        n_reads_max = 10000,
                        theta0 = 0.01,
                        rho = NULL,
                        Sigma = NULL) {

  # check for required packages
  rlang::check_installed("dirmult", reason = "to use `simulate_DM()`")
  rlang::check_installed("MASS", reason = "to use `simulate_DM()`")
  rlang::check_installed("matrixcalc", reason = "to use `simulate_DM()`")

  # Defense
  if (!is.null(rho) & !is.null(Sigma)) {
    stop("Please just provide either rho or Sigma for covariate correlation structure.")
  }
  if (!is.null(rho)) {
    if (rho > 1 | rho < 0) {
      stop("Please provide rho between 0 and 1.")
    }
  }
  if (!is.null(Sigma)) {
    if (!matrixcalc::is.positive.definite(Sigma)) {
      stop("Covaraince matrix must be a symetric positive definite matrix.")
    }
  }
  if (!is.null(Sigma)) {
    if (ncol(Sigma) != n_vars) {
      stop("Please provide covariance matrix to match the number of covariates")
    }
  }

  # covariance matrix for predictors
  if (!is.null(rho)) {
    Sigma <- matrix(0, n_vars, n_vars)
    Sigma <- rho^abs(row(Sigma) - col(Sigma))
  }


  # include the intercept
  XX <- cbind(rep(1, n_obs),
              scale(MASS::mvrnorm(n = n_obs,
                                  mu = rep(0, n_vars),
                                  Sigma = Sigma)))
  # Mon May  2 11:56:48 2022 ------------------------------
  # using empty matrix
  YY <- matrix(0, n_obs, n_taxa)
  betas <- matrix(0, n_taxa, n_vars)
  phi <- matrix(0, n_obs, n_taxa)

  # parameters with signs alternating
  st <- 0
  low_side <- beta_min
  high_side <- beta_max

  if (n_relevant_taxa != 1) {
    # warning if the lengths don't match
    coef <- suppressWarnings(seq(low_side,
                                 high_side,
                                 len = n_relevant_taxa) * c(1, -1))
    } else {
      coef <- (low_side + high_side) / 2
    }


  coef_g <- rep(1.0, len = n_relevant_vars)


  for (ii in 1:n_relevant_vars) {
    # overlap species
    betas[(st:(st + n_relevant_taxa - 1)) %% n_taxa + 1, 3 * ii - 2] <-
      coef_g[ii] * sample(coef)[((ii - 1):(ii + n_relevant_taxa - 2)) %% n_relevant_taxa + 1]
    st <- st + 1
  }

  # -2.3 and 2.3 so that the intercept varies over three orders of magnitude
  intercept <- runif(n_taxa, -2.3, 2.3)
  Beta <- cbind(intercept, signoise * betas)

  # row totals
  ct0 <- sample(n_reads_min:n_reads_max, n_obs, rep = T)

  for (ii in 1:n_obs) {
    thisrow <- as.vector(exp(Beta %*% XX[ii, ]))
    phi[ii, ] <- thisrow / sum(thisrow)
    YY[ii, ] <- dirmult::simPop(J = 1,
                                n = ct0[ii],
                                pi = phi[ii, ],
                                theta = theta0)$data[1, ]
  }

  return(list(
    X = XX,
    Y = YY,
    alphas = intercept,
    betas = Beta,
    n_reads_min = n_reads_min,
    n_reads_max = n_reads_max,
    theta0 = theta0,
    phi = phi,
    rho = rho,
    signoise = signoise,
    Sigma = Sigma
  ))
}
## }}}-----------------


## 4.2 Xsim {{{---------------
#' Title: Simulate Design matrix and Indicators
#'
#' @param subject_sim `integer` Number of subjects to simulate
#' @param tree `phylo` The function is set to simulate the tree structure for the
#'                     simulated data. Default is NULL.
#' @param num_leaf`integer` Number of nodes in tree to simulate
#' @param covariates_sim `integer` Number of covariates to simulate
#' @param rho `numeric` Used for covariate correlation between 0 & 1. Default is NULL.
#'                      Need either rho or Sigma to be input.
#' @param Sigma `matrix` A symmetric positive definite matrix used for covariate
#'                       correlation structure columns must match number of covariates.
#'                       Default is NULL. Need either rho or Sigma to be input.
#' @param num_branch `integer` Number of tree branches
#' @param num_cov `integer` Number of associated covariates to simulate
#' @param seed `integer` The desired seed for reproducibility. Default is 555.
#'
#' @return A list of outcomes including:
#'   - X = The design matrix
#'   - Zeta = The matrix indicating which covariates are significant
#'   - Sigma = The simulated covariance matrix
#'   - tree = The simulated tree, `phylo` object generated by the `ape` package
#' @export
#'
Xsim<- function(subject_sim = 100,
                tree = NULL,
                num_leaf = 10,
                covariates_sim = 50,
                rho = NULL,
                Sigma = NULL,
                num_branch = 3,
                num_cov = 5,
                seed = 555) {

  set.seed(seed)

  if (!is.null(rho)) {
    Sigma <- matrix(0, covariates_sim, covariates_sim)
    Sigma <- rho^abs(row(Sigma) - col(Sigma))
  }

  tree.ex <- if (is.null(tree)) {
    ape::rtree(n = num_leaf)
  } else {tree}

  # Get dimensions from tree
  # Set number of parent nodes = #subtrees = #parentheses sets
  V <- tree.ex$Nnode

  # Set number of child nodes for each parent node
  Cv <- table(tree.ex$edge[, 1])

  # Set number of leaves (tips) in the tree
  K <- length(tree.ex$tip.label)

  # Set parameters
  B_sim <- sum(Cv)

  # Simulate covariates for selection
  # Set true inclusion indicators
  zeta_sim <- matrix(0, B_sim, covariates_sim)
  for (i in sample(seq(1, B_sim), num_branch)) {
    select <- sample(1:covariates_sim, num_cov)
    zeta_sim[i, select] <- 1
  }

  X <- scale(mvtnorm::rmvnorm(subject_sim,
                            rep(0, covariates_sim),
                            Sigma))
  zeta_sim
  return(list(X = X,
              zeta = zeta_sim,
              Sigma = Sigma,
              tree = tree.ex))
}
## }}}-----------------------



## 4.3 simulate_DTM {{{---------------
#
#' Title: Simulate Dirichlet-Multinomial-Tree Regression Model
#' @description Wrapper function for the Rcpp code to simulate DTM data
#' Taken from Wadsworth (2017)
#'
#' @param subject_sim `integer` Number of subjects to simulate, default = 100
#' @param tree `phylo` Tree structure for simulated data, default is NULL.
#' @param num_leaf `integer` The number of leaves in tree to simulate
#' @param covariates_sim `integer` Number of covariates to simulate, default = 50
#' @param rho `numeric` Correlation structure between 0 & 1, default = 0.20
#' @param Sigma `matrix` Variance-covariance matrix must be symmetric positive semi-definite
#' @param num_branch `integer` Number of branches the associated covariates are in, default = 3
#' @param num_cov `integer` Number of associated covariates to simulate, default = 5
#' @param phi_min
#' @param phi_max
#' @param seed `integer` Seed for reproducibility, default = 555
#'
#' @return A list of results:
#' - Y: the compositional counts
#' - X: the covariates design matrix
#' - phi_sim: true regression coefficients
#' - zeta_sim: true inclusion indicators
#' - tree: the random tree, `phylo` object generated by the `ape` package
#' - Sigma: the simulated covariance matrix
#'
#' @export
#'
#'
simulate_DTM <- function(subject_sim = 100,
                         tree = NULL,
                         num_leaf = 10,
                         covariates_sim = 50,
                         rho = NULL,
                         Sigma = NULL,
                         X = X,
                         zeta_sim = zeta,
                         rep = rep,
                         num_branch = 3,
                         num_cov = 5,
                         phi_min = 0.9,
                         phi_max = 1.2,
                         seed = 555) {

  # Set seed for replication
  set.seed(seed)

  # use_package("mvtnorm")
  # use_package("MCMCpack")
  # use_package("ape")
  # use_package("GGMselect")
  # use_package("rlang")

  rlang::check_installed("dirmult", reason = "to use `simulate_DM()`")
  rlang::check_installed("mvtnorm", reason = "to use `simulate_DM()`")
  rlang::check_installed("MCMCpack", reason = "to use `simulate_DM()`")
  rlang::check_installed("ape", reason = "to use `simulate_DM()`")
  rlang::check_installed("GGMselect", reason = "to use `simulate_DM()`")

  #### Defense {{{{-----------------
  # Simulate data if there is no tree, X, or Y
  if (subject_sim %% 1 != 0 | subject_sim < 0) {
    stop("Bad input: number of subjects should be a positive integer")
  }
  if (num_leaf %% 1 != 0 | num_leaf < 0) {
    stop("Bad input: number of leaves should be a positive integer")
  }
  if (covariates_sim %% 1 != 0 | covariates_sim < 0) {
    stop("Bad input: number of covariates should be a positive integer")
  }
  if (num_branch %% 1 != 0 | num_branch < 0) {
    stop("Bad input: number of branches for associated covariates should be a positive integer")
  }
  if (num_cov %% 1 != 0 | num_cov < 0) {
    stop("Bad input: number of associated covariates should be a positive integer")
  }
  if (num_cov >= covariates_sim) {
    stop("Bad input: number of associated covariates should be less than total number of covariates")
  }
  if (!is.null(rho) & !is.null(Sigma)) {
    stop("Bad input: Please provide either rho or Sigma for covariate correlation structure.")
  }
  if (is.null(rho) & is.null(Sigma)) {
    stop("Bad input: Please provide either rho or Sigma for covariate correlation structure.")
  }
  if (!is.null(rho)) {
    if (rho > 1 | rho < 0) {
      stop("Bad input: Please provide rho between 0 and 1.")
    }
  }
  if (!is.null(Sigma)) {
    if (!matrixcalc::is.positive.definite(Sigma)) {
      stop("Bad input: Please provide positive definite covariance matrix.")
    }
  }
  if (!is.null(Sigma)) {
    if (ncol(Sigma) != covariates_sim) {
      stop("Bad input: Please provide covariance matrix to match the number of covariates")
    }
  }

  # covariance matrix for predictors
  if (!is.null(rho)) {
    Sigma <- matrix(0, covariates_sim, covariates_sim)
    Sigma <- rho^abs(row(Sigma) - col(Sigma))
  }
  ####}}}}-----------------


  set.seed(seed)

  # Simulate random DTM with phylogenetic tree
  # Relies on 'ape' package
  tree.ex <- if (is.null(tree)) {
    ape::rtree(n = num_leaf)
  } else {tree}

  # Get dimensions from tree
  # Set number of parent nodes = #subtrees = #parentheses sets
  V <- tree.ex$Nnode

  # Set number of child nodes for each parent node
  Cv <- table(tree.ex$edge[, 1])

  # Set number of leaves (tips) in the tree
  K <- length(tree.ex$tip.label)

  # Set parameters
  B_sim <- sum(Cv)

  # Simulate covariates for selection
  # Set true inclusion indicators

  truth <<- which(zeta_sim == 1)

  # Simulate true alpha parameters and regression coefficients phi
  alpha_sim <- matrix(1,
                      nrow = subject_sim,
                      ncol = 1) %*%
    (runif(n = B_sim, -1.3, 1.3))

  true_cov <- which(zeta_sim == 1)
  phi_sim <- matrix(0, B_sim, covariates_sim)
  phi_sim[true_cov] <- runif(sum(zeta_sim), phi_min, phi_max) * # What is phi used for?
    sample(c(-1, 1), sum(zeta_sim), replace = TRUE)

  # Used for count probabilities
  inside_sim <- exp(alpha_sim + X %*% t(phi_sim))

  # Look through the tree and separate to
  # get dirichlet parameters and then simulated

  Y_list <- list()

  for (rm in seq_along(1:rep)) {
      node_counts <- matrix(0,
                            nrow = subject_sim,
                            ncol = (V + K))

      ## Why using those two values? 7500, 10000???-------
      node_counts[, (K + 1)] <- sample(seq(7500, 10000),
                                       subject_sim)

      for (b in (K + 1):(sum(Cv) + 1)) {
        node <- which(tree.ex$edge[, 1] == b)

        # Split inside by each subtree
        inside_branches <- inside_sim[, node]

        # Simulate probabilities for each subtree
        prob_sim <- apply(inside_branches,
                          1,
                          function(x) {dirmult::rdirichlet(1, x)})

        # Simulate count data
        for (i in 1:subject_sim) {
          y <- t(rmultinom(1, node_counts[i, b], t(prob_sim)[i, ]))
          node_counts[i, (tree.ex$edge[node, 2])] <- y
        }
      }
      Y <- node_counts[, 1:K]
      Y_list[[rm]]<- Y
  }

  X <- scale(X)

  return(list(
    Y = Y_list,
    X = X,
    phi_sim = phi_sim,
    zeta_sim = zeta_sim,
    tree = tree.ex,
    Sigma = Sigma
  ))
}
## }}}----------------



## 4.4 simulate_DMLM {{{------------
# Code to simulate data for DMLMbvs model (Doubly Multivariate Linear Model)
#' Title simulate data for DMLMbvs model
#' @param subject_sim `integer` Number of subjects to simulate
#' @param B_sim `integer` Number of betas
#' @param n_vars `integer` Number of covariates to simulate
#' @param active_cov `integer` Number of associated covariates
#' @param rho `numeric` Used for covariate correlation between 0 & 1. Default is NULL.
#'                      Need either rho or Sigma to be input.
#' @param Sigma `matrix` A symmetric positive definite matrix used for covariate
#'                       correlation structure columns must match number of covariates.
#'                       Default is NULL. Need either rho or Sigma to be input.
#' @param seed `integer` Seed for reproducibility, default = 555
#'
#' @return a list of outcomes:
#' - Y = Y,
#' - Z = Z,
#' - X = X,
#' - true_cov = true_cov,
#' - true_coeff = true_coeff,
#' - true_coeff_beta = true_coeff_beta)
#' @export
#'
simulate_DMLM <- function(subject_sim = 50,
                          B_sim = 50,
                          n_vars = 50,
                          active_cov = 10,
                          rho = NULL,
                          Sigma = NULL,
                          seed = 111) {

  # Call libraries
  rlang::check_installed("MCMCpack", reason = "to use `simulate_DM()`")
  rlang::check_installed("mvtnorm", reason = "to use `simulate_DM()`")

  # Defense
  if (!is.null(rho) & !is.null(Sigma)) {
    stop("Bad input: Please provide either rho or Sigma for covariate correlation structure.")
  }
  if (is.null(rho) & is.null(Sigma)) {
    stop("Bad input: Please provide either rho or Sigma for covariate correlation structure.")
  }
  if (!is.null(rho)) {
    if (rho > 1 | rho < 0) {
      stop("Bad input: Please provide rho between 0 and 1.")
    }
  }
  if (!is.null(Sigma)) {
    if (!matrixcalc::is.positive.definite(Sigma)) {
      stop("Bad input: Please provide positive definite covariance matrix.")
    }
  }
  if (!is.null(Sigma)) {
    if (ncol(Sigma) != n_vars) {
      stop("Bad input: Please provide covariance matrix to match the number of covariates")
    }
  }

  # covariance matrix for predictors
  if (!is.null(rho)) {
    Sigma <- matrix(0, n_vars, n_vars)
    Sigma <- rho^abs(row(Sigma) - col(Sigma))
  }

  # Set seed
  set.seed(seed)

  # Set covariance strucuture for covariates
  X <- scale(rmvnorm(subject_sim,
                     rep(0, n_vars),
                     Sigma))

  zeta_sim <- matrix(0, B_sim, n_vars)
  true_cov <- cbind(sample(seq(1, active_cov),
                           active_cov,
                           replace = T),
                    sample(seq(1, n_vars),
                           active_cov))
  zeta_sim[true_cov] <- 1
  alpha_sim <-
    matrix(1, nrow = subject_sim,
           ncol = 1) %*%
    (runif(n = B_sim, -2.3, 2.3))

  true_coeff <- runif(10, 0.75, 1.25)
  phi_sim <- matrix(0, B_sim, n_vars)
  phi_sim[true_cov] <-
    true_coeff * sample(c(1, -1), 10, replace = TRUE)

  inside_sim <- exp(alpha_sim + X %*% t(phi_sim))
  psi_sim <- inside_sim / rowSums(inside_sim)
  psi_sim_overdispersed <- psi_sim * (1 - 0.01) / 0.01

  # Simulate probabilities for each branch
  prob_sim <- apply(psi_sim_overdispersed, 1,
                    function(x) {rdirichlet(1, x)})

  # Simulate count data for each subtree branch
  # (simplified for bifurcating tree in this example)
  Z <-
    t(apply(t(prob_sim), 1,
            function(x) {rmultinom(1, sample(seq(2500, 7500), 1), x)}))

  # adjust for zero counts

  Z_norm <- Z
  Z_norm[Z_norm == 0] <- 0.5
  Z_norm <- Z_norm / rowSums(Z_norm)

  prob_sim[prob_sim < 0.5 / 7500] <- 0.5 / 7500
  prob_sim2 <- t(prob_sim) / colSums(prob_sim)

  # Calculate Balances
  Balances <- ilr_fun_cpp(prob_sim2)

  # Select Balances
  xi_sim <- rep(0, B_sim - 1)
  true_cov_beta <- seq(1:10)

  xi_sim[true_cov_beta] <- 1
  true_coeff_beta <- runif(length(true_cov_beta), 1.25, 1.75)

  beta_sim <- rep(0, B_sim - 1)

  beta_sim[true_cov_beta] <-
    true_coeff_beta * sample(c(1, -1),
                             length(true_cov_beta),
                             replace = TRUE)

  Y <- Balances %*% beta_sim + rnorm(subject_sim)

  return(list(Y = Y,
              Z = Z,
              X = X,
              true_cov = true_cov,
              true_coeff = true_coeff,
              true_coeff_beta = true_coeff_beta))
}


## }}} ----------






## simulate_DTM {{{---------------
#
#' @title Simulate Dirichlet-Multinomial-Tree Longitudinal Regression Model
#' @description This is the adaptation of previous simulation function to the
#' longitudinal setting, where we can have repeated measure for each of the individual
#'
#' @param subject_sim `integer` the number of subjects to simulate, default = 100
#' @param n_measurement `integer` number of measurements for each subject, default = 3
#' @param X design matrix of covariate
#' @param zeta_sim true inclusion indicators (must be provided if X is imported in argument)
#' @param tree `phylo` a tree structure for simulated data
#' @param num_leaf `integer` the number of leaves in tree to simulate, default = 5
#' @param covariates_sim `integer` the number of covariates to simulate, default = 50
#' @param rho `numeric` the correlation structure, default = 0.20
#' @param Sigma `matrix` the variance-covariance matrix must be symmetric positive semi-definite
#' @param global_random `logi` whether all the taxa share the same random effect or not, default = TRUE
#' @param random_time `logi` whether to include random slope term for time or not, default = FALSE
#' @param random_sim_mean a vector of mean used to simulate random effect, length = 1 if random intercept only
#' @param random_sim_cov a covariance matrix used to simulate random effect, length = 1 if random intercept only
#' @param random_sim_mean_local a vector of mean if local random effect is in used, length = number of taxa
#' @param random_sim_cov_local a covariance matrix if local random effect is in used, length = number of taxa
#' @param num_branch `integer` number of branches the associated covariates are in, default = 3
#' @param num_cov `integer` the number of associated covariates to simulate, default = 5
#' @param phi_min
#' @param phi_max
#' @param seed `integer` seed setup, default = 555
#'
#' @return A list of results:
#' - Y: the compositional counts
#' - X: the covariates design matrix
#' - phi_sim: true regression coefficients
#' - zeta_sim: true inclusion indicators
#' - tree: the random tree, `phylo` object generated by the `ape` package
#' - Sigma: the simulated covariance matrix
#'
#' @examples # random intercept only
#' data_random_intercept_global <- simulate_correlated_DTM(rho = 0.2, random_sim_mean = 0, random_sim_cov = 1, rep = 100)
#' # random intercept and slope (uncorrelated)
#' data_random_intercept_slope_global <- simulate_correlated_DTM(rho = 0.2, X = data_random_intercept_global$X, zeta_sim = data_random_intercept_global$zeta_sim, random_sim_mean = c(1,1), random_sim_cov = diag(1,2), rep = 100, random_time = TRUE)
#' # random intercept and random slope (correlated)
#' data_random_intercept_slope_correlated_global <- simulate_correlated_DTM(rho = 0.2, X = data_random_intercept_global$X, zeta_sim = data_random_intercept_global$zeta_sim,random_sim_mean = c(1,1), random_sim_cov =  matrix(c(1,0.2,0.2,1), nrow = 2), rep = 100, random_time = TRUE, seed = 100)
#'
#' @export
#'
#'





simulate_correlated_DTM <- function(subject_sim = 100,
                                    n_measurement = 3,
                                    X = NULL,
                                    zeta_sim = NULL,
                                    tree = NULL,
                                    num_leaf = 10,
                                    covariates_sim = 50,
                                    rho = NULL,
                                    Sigma = NULL,
                                    global_random = TRUE,
                                    random_time = FALSE,
                                    random_sim_mean = NULL,
                                    random_sim_cov = NULL,
                                    random_sim_mean_local = NULL,
                                    random_sim_cov_local = NULL,
                                    rep = rep,
                                    num_branch = 3,
                                    num_cov = 5,
                                    phi_min = 0.9,
                                    phi_max = 1.2,
                                    seed = 555) {

  # Set seed for replication
  set.seed(seed)

  # use_package("mvtnorm")
  # use_package("MCMCpack")
  # use_package("ape")
  # use_package("GGMselect")
  # use_package("rlang")



  rlang::check_installed("dirmult", reason = "to use `simulate_DM()`")
  rlang::check_installed("mvtnorm", reason = "to use `simulate_DM()`")
  rlang::check_installed("MCMCpack", reason = "to use `simulate_DM()`")
  rlang::check_installed("ape", reason = "to use `simulate_DM()`")
  rlang::check_installed("GGMselect", reason = "to use `simulate_DM()`")

  #### Defense {{{{-----------------
  # Simulate data if there is no tree, X, or Y
  if (subject_sim %% 1 != 0 | subject_sim < 0) {
    stop("Bad input: number of subjects should be a positive integer")
  }

  if (num_leaf %% 1 != 0 | num_leaf < 0) {
    stop("Bad input: number of leaves should be a positive integer")
  }

  if (covariates_sim %% 1 != 0 | covariates_sim < 0) {
    stop("Bad input: number of covariates should be a positive integer")
  }

  if (num_branch %% 1 != 0 | num_branch < 0) {
    stop("Bad input: number of branches for associated covariates should be a positive integer")
  }

  if (num_cov %% 1 != 0 | num_cov < 0) {
    stop("Bad input: number of associated covariates should be a positive integer")
  }

  if (num_cov >= covariates_sim) {
    stop("Bad input: number of associated covariates should be less than total number of covariates")
  }

  if (!is.null(rho) & !is.null(Sigma)) {
    stop("Bad input: Please provide either rho or Sigma for covariate correlation structure.")
  }

  if (is.null(rho) & is.null(Sigma)) {
    stop("Bad input: Please provide either rho or Sigma for covariate correlation structure.")
  }

  if (!is.null(rho)) {
    if (rho > 1 | rho < 0) {
      stop("Bad input: Please provide rho between 0 and 1.")
    }
  }

  if (!is.null(Sigma)) {
    if (!matrixcalc::is.positive.definite(Sigma)) {
      stop("Bad input: Please provide positive definite covariance matrix.")
    }
  }

  if (!is.null(Sigma)) {
    if (ncol(Sigma) != covariates_sim) {
      stop("Bad input: Please provide covariance matrix to match the number of covariates")
    }
  }

  # covariance matrix for predictors
  if (!is.null(rho)) {
    Sigma <- matrix(0, covariates_sim, covariates_sim)
    Sigma <- rho^abs(row(Sigma) - col(Sigma))
  }
  ####}}}}-----------------


  set.seed(seed)

  # Simulate random DTM with phylogenetic tree
  # Relies on 'ape' package
  tree.ex <- if (is.null(tree)) {
    ape::rtree(n = num_leaf)
  } else {tree}

  # Get dimensions from tree
  # Set number of parent nodes = #subtrees = #parentheses sets
  V <- tree.ex$Nnode

  # Set number of child nodes for each parent node
  Cv <- table(tree.ex$edge[, 1])

  # Set number of leaves (tips) in the tree
  K <- length(tree.ex$tip.label)

  # Set parameters
  B_sim <- sum(Cv)

  # Simulate covariates for selection
  # Set true inclusion indicators

  # simulate design matrix of X (if X is not provided)
  if (is.null(X))
  {
    X_sim <- Xsim(subject_sim = subject_sim * n_measurement,
                  tree = tree,
                  num_leaf = num_leaf,
                  covariates_sim = covariates_sim,
                  rho = rho,
                  Sigma = Sigma,
                  num_branch = num_branch,
                  num_cov = num_cov,
                  seed = seed)
    X <- X_sim$X
    zeta_sim <- X_sim$zeta
  }


  # create a vector for measurement time
  time <- rep(1:n_measurement, subject_sim)

  # create a vector for subject id
  id <- rep(1:subject_sim, each = n_measurement)

  # create a vector

  truth <<- which(zeta_sim == 1)

  # Simulate true alpha parameters and regression coefficients phi
  alpha_sim <- matrix(1,
                      nrow = subject_sim * n_measurement,
                      ncol = 1) %*%
    (runif(n = B_sim, -1.3, 1.3))

  true_cov <- which(zeta_sim == 1)
  phi_sim <- matrix(0, B_sim, covariates_sim)
  phi_sim[true_cov] <- runif(sum(zeta_sim), phi_min, phi_max) *
    sample(c(-1, 1), sum(zeta_sim), replace = TRUE)



  # adding global or local random effect to the count probability (Z matrix)
  if (global_random == TRUE)
  {

    # if we have random slope for time
    if (random_time == TRUE)
    {
      # randomly simulate coefficient for random terms
      b_coef <- MASS::mvrnorm(n = subject_sim, mu = random_sim_mean, Sigma = random_sim_cov,
                              tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

      b_matrix <- apply(b_coef, 2, function(x)rep(x, each = n_measurement))
      # obtain design matrix
      Z <- cbind(rep(1, n_measurement * subject_sim), rep(1:n_measurement, subject_sim))

      # obtain the total random effect
      random_effect_sum <- apply(Z * b_matrix, 1, sum)
    }else{
      b_matrix <- rep(rnorm(subject_sim), mean = random_sim_mean,
                      sd = random_sim_cov, each = n_measurement)

      # obtain the total random effect
      random_effect_sum <- b_matrix
    }



  }else{


    # if we have random slope for time
    if (random_time == TRUE)
    {
      # randomly simulate coefficient for random terms
      b_coef <- MASS::mvrnorm(n = subject_sim, mu = random_sim_mean_local %x% random_sim_mean, Sigma = random_sim_cov_local %x% random_sim_cov,
                              tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

      # create an empty matrix for random effect
      random_effect_sum <- matrix(0, nrow = subject_sim * n_measurement, ncol = B_sim)

      # loop through each vc
      for (i in seq(from = 1, to = ncol(b_coef) - 1, by = 2))
      {
        b_coef <- cbind(b_coef[,i], b_coef[,i+1])
        b_matrix <- apply(b_coef, 2, function(x)rep(x, each = n_measurement))
        # obtain design matrix
        Z <- cbind(rep(1, n_measurement * subject_sim), rep(1:n_measurement, subject_sim))
        # obtain the total random effect
        random_effect_sum[,i] <- apply(Z * b_matrix, 1, sum)
      }


    }else{
      # randomly simulate coefficient for random terms
      b_coef <- MASS::mvrnorm(n = subject_sim, mu = random_sim_mean_local, Sigma = random_sim_cov_local,
                              tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
      # create an empty matrix for random effect
      random_effect_sum <- matrix(0, nrow = subject_sim * n_measurement, ncol = B_sim)

      # loop through each vc
      for (i in 1:ncol(b_coef))
      {
        b_coef <- b_coef[,i]
        b_matrix <- rep(b_coef, each = n_measurement)

        # obtain the total random effect
        random_effect_sum[,i] <- b_matrix
      }



    }



  }


  # Used for count probabilities
  inside_sim <- exp(alpha_sim + X %*% t(phi_sim) + random_effect_sum)

  # Look through the tree and separate to
  # get dirichlet parameters and then simulated

  Y_list <- list()

  for (rm in seq_along(1:rep)) {
    node_counts <- matrix(0,
                          nrow = subject_sim * n_measurement,
                          ncol = (V + K))

    node_counts[, (K + 1)] <- sample(seq(7500, 10000),
                                     subject_sim * n_measurement)

    for (b in (K + 1):(sum(Cv) + 1)) {
      node <- which(tree.ex$edge[, 1] == b)

      # Split inside by each subtree
      inside_branches <- inside_sim[, node]

      # Simulate probabilities for each subtree
      prob_sim <- apply(inside_branches,
                        1,
                        function(x) {dirmult::rdirichlet(1, x)})

      # Simulate count data
      for (i in 1:subject_sim) {
        y <- t(rmultinom(1, node_counts[i, b], t(prob_sim)[i, ]))
        node_counts[i, (tree.ex$edge[node, 2])] <- y
      }
    }
    Y <- node_counts[, 1:K]
    Y_list[[rm]]<- Y
  }

  X <- scale(X)

  return(list(
    Y = Y_list,
    X = X,
    b = b,
    random_effect = b_matrix,
    total_random_effect = random_effect_sum,
    phi_sim = phi_sim,
    zeta_sim = zeta_sim,
    tree = tree.ex,
    Sigma = Sigma,
    n_measurement = n_measurement,
    id = id,
    time = time
  ))
}
## }}}---------------























