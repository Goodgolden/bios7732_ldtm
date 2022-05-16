## 5.0 Single fold {{{---------------

#' Title
#'
#' @param obs
#' @param k
#'
#' @return
#' @export
#'

singlefold <- function(obs, k) {
  if (k == 1) {
    return(rep(1, obs))
  }
  else {
    i <- obs/k
    if (i < 1) {
      stop("insufficient records:", obs,
           ", with k=", k)
    }
    i <- round(c(0, i * 1:(k - 1), obs))
    times = i[-1] - i[-length(i)]
    group <- c()
    for (j in 1:(length(times))) {
      group <- c(group, rep(j, times = times[j]))
    }
    r <- order(runif(obs))
    return(group[r])
  }
}
## }}}------------------



## 5.1 K-fold dismo {{{-----------------------

#' Title
#' @details  adapted from R package dismo
#' ??dismo # Species Distribution Modeling
# library(dismo)
# ?dismo
#' @param x
#' @param k
#' @param by
#'
#' @return
#' @export
#'
#' @examples
#'
#'

kfold.dismo <- function (x, k = 5, by = NULL) {

  if (is.vector(x)) {
    if (length(x) == 1) {
      if (x > 1) {
        x <- 1:x
      }
    }
    obs <- length(x)
  }
  else if (inherits(x, "Spatial")) {
    if (inherits(x, "SpatialPoints")) {
      obs <- nrow(coordinates(x))
    }
    else {
      obs <- nrow(x@data)
    }
  }
  else {
    obs <- nrow(x)
  }
  if (is.null(by)) {
    return(singlefold(obs, k))
  }
  by = as.vector(as.matrix(by))
  if (length(by) != obs) {
    stop("by should be a vector with the same number of records as x")
  }
  un <- unique(by)
  group <- vector(length = obs)
  for (u in un) {
    i = which(by == u)
    kk = min(length(i), k)
    if (kk < k)
      warning("lowered k for by group: ", u, "  because the number of observations was  ",
              length(i))
    group[i] <- singlefold(length(i), kk)
  }
  return(group)
}


## 5.2 cv_T_group_fun{{{-----------------
#' Title
#'
#' @param Y
#' @param X
#' @param treeinfo
#' @param alpha
#' @param cutoff
#' @param model
#' @param err.conv
#' @param iter.max
#' @param L.init
#' @param lambda.max
#' @param lambda.min
#' @param nfolds
#' @param CV
#' @param cv.ind
#'
#' @return
#' @export
#'
cv_T_group_fun <-
  function(Y, X,
           treeinfo,
           alpha = 0.5,
           cutoff = 0.8,
           model = "dirmult",
           err.conv = 1e-3,
           iter.max = 30,
           L.init = NULL,
           lambda.max = NULL,
           lambda.min = NULL,
           nfolds = 10,
           CV = T,
           cv.ind = NULL) {

    output <- T_group_path(Y, X,
                           treeinfo,
                           alpha,
                           cutoff,
                           model,
                           err.conv,
                           iter.max,
                           L.init,
                           lambda.max,
                           lambda.min)

    # cross validation starts here

    if (CV == F) {
      return(output)
    }

    ## cross validation
    else {
      output.path <- output$output.path
      path.length <- length(output.path)
      lambda.max <- output.path[[1]]$lambda
      lambda.min <- output.path[[path.length]]$lambda

      E <- matrix(NA, nrow = nfolds, ncol = path.length)
      if (is.null(cv.ind)) cv.ind <- kfold.dismo(1 : dim(X)[1], nfolds)

      for (i in 1 : nfolds) {

        cat("Starting CV fold #", i, sep = "", "\n")
        X1 <- X[cv.ind != i, , drop = FALSE]
        Y1 <- Y[cv.ind != i, , drop = FALSE]
        X2 <- X[cv.ind == i, , drop = FALSE]
        Y2 <- Y[cv.ind == i, , drop = FALSE]

        output.i <- T_group_path(Y1, X1, treeinfo,
                                 alpha, cutoff = 1,
                                 model, err.conv, iter.max,
                                 L.init, lambda.max,
                                 lambda.min)

        output.path.i <- output.i$output.path
        path.length.i <- length(output.path.i)
        for (j in 1 : path.length.i) {

          pred.Y2.i.j <- pred_Y(Y2, X2, output.path.i[[j]]$B, treeinfo)
          E[i, j] <- sum((Y2 - pred.Y2.i.j)^2 / pred.Y2.i.j / rowSums(Y2))
        }
      }

      cve <- apply(E, 2, mean, na.rm = T)
      cvse <- apply(E, 2, sd, na.rm = T) / sqrt(nfolds)

      min.i <- which.min(cve)[1]
      min.1se <- which(cve <= cve[min.i] + cvse[min.i])[1]

      output.cv <- list(cve = cve,
                        cvse = cvse,
                        min.i = min.i,
                        min.1se = min.1se,
                        E = E,
                        output = output)
      return(output.cv)
    }
  }
## }}}---------------------------


## 5.3 IC_fun {{{----------------------

#' Title Getting the information criteria
#'
#' @param output
#' @param Y
#' @param X
#' @param treeinfo
#' @param model
IC_fun <- function(output,
                   Y, X,
                   treeinfo,
                   model = "dirmult") {

  print("lasso and group_lasso_only!")

  n <- dim(X)[1]
  Ytree <- Y_tree(Y, treeinfo)

  Aic.list <- Bic.list <-
    Df.list <- list()
  Aic <- Bic <- Df.Aic <-
    Df.Bic <- min.i.Aic <- min.i.Bic <-
    rep(NA, length(output))

  for (j in 1 : length(output)) {

    output.path.j <- output[[j]]$output.path
    alpha.j <- output[[j]]$alpha
    Aic.j <- Bic.j <- Df.j <- NULL

    for (k in 1 : length(output.path.j)) {

      B.j.k <- output.path.j[[k]]$B
      loglik.j.k <- Loglik_Dirichlet_tree(Ytree, X, B.j.k, model = model)
      dev.j.k <- -2 * loglik.j.k * n

      if (alpha.j == 1) {
        Df.j.k <- sum(colSums(B.j.k[, -1] != 0) != 0) * dim(B.j.k[, -1])[1]
      } else {
        Df.j.k <- sum(B.j.k[, -1] != 0)
      }

      Aic.j.k <- dev.j.k + Df.j.k * 2
      Bic.j.k <- dev.j.k + Df.j.k * log(n)

      Df.j <- c(Df.j, Df.j.k)
      Aic.j <- c(Aic.j, Aic.j.k)
      Bic.j <- c(Bic.j, Bic.j.k)
    }

    Df.list[[j]] <- Df.j
    Aic.list[[j]] <- Aic.j
    Bic.list[[j]] <- Bic.j

    min.i.Aic[j] <- which.min(Aic.j)[1]
    Df.Aic[j] <- Df.j[min.i.Aic[j]]
    Aic[j] <- Aic.j[min.i.Aic[j]]

    min.i.Bic[j] <- which.min(Bic.j)[1]
    Df.Bic[j] <- Df.j[min.i.Bic[j]]
    Bic[j] <- Bic.j[min.i.Bic[j]]
  }

  return(list(Aic = Aic,
              Bic = Bic,
              Df.Aic = Df.Aic,
              Df.Bic = Df.Bic,
              min.i.Aic = min.i.Aic,
              min.i.Bic = min.i.Bic,
              Aic.list = Aic.list,
              Bic.list = Bic.list,
              Df.list = Df.list))
}
## }}}-----------------
