# standard group lasso fit
group_lasso_fit <- function(X, y, grp, w, lambda, maxit = 4000, tol = 1e-9) {
  n <- nrow(X)
  gs <- split(seq_len(ncol(X)), grp)
  XtX <- crossprod(X) / n
  Xty <- as.numeric(crossprod(X, y) / n)
  Lip <- max(eigen(XtX, symmetric = TRUE, only.values = TRUE)$values) + 1e-8
  b <- numeric(ncol(X))
  z <- b
  tk <- 1
  obj_old <- Inf
  for (it in seq_len(maxit)) {
    bn <- z - as.numeric(XtX %*% z - Xty) / Lip
    for (k in seq_along(gs)) {
      wk <- w[[names(gs)[k]]]
      if (wk > 0) {
        id <- gs[[k]]
        nb <- sqrt(sum(bn[id]^2))
        th <- lambda * wk / Lip
        bn[id] <- if (nb > th) (1 - th / nb) * bn[id] else 0
      }
    }
    tk2 <- (1 + sqrt(1 + 4 * tk^2)) / 2
    z <- bn + ((tk - 1) / tk2) * (bn - b)
    b <- bn
    tk <- tk2
    if (it %% 20 == 0) {
      gn <- vapply(gs, function(id) sqrt(sum(b[id]^2)), 0)
      obj <- sum((y - X %*% b)^2) / (2 * n) + lambda * sum(unlist(w[names(gs)]) * gn)
      if (abs(obj_old - obj) < tol * max(1e-10, obj)) break
      obj_old <- obj
    }
  }
  b
}
