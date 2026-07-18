# 
fit_destination <- function(dest, ridx, U, Ulag, E0, P, L = 3,
                             eps_w = 0.005, hold_fr = 0.25) {
  nodes <- colnames(U)
  src <- nodes[E0[, dest]]
  if (!length(src)) return(NULL)
  ridx <- ridx[ridx > L]
  y <- U[ridx, dest]
  Xo <- sapply(1:L, function(l) Ulag[[l]][ridx, dest])
  Xs <- do.call(cbind, lapply(src, function(s)
          sapply(1:L, function(l) Ulag[[l]][ridx, s])))
  X <- cbind(Xo, Xs)
  grp <- c(rep("own", L), rep(src, each = L))
  sds <- apply(X, 2, sd)
  sdy <- sd(y)
  keep <- sds > 1e-8
  if (sdy < 1e-8 || !any(keep[grp != "own"])) return(NULL)
  X <- scale(X[, keep, drop = FALSE])
  grp <- grp[keep]
  y <- as.numeric(scale(y))
  n <- length(y)
  sds <- sds[keep]
  src <- unique(grp[grp != "own"])
  wpen <- 1 / sqrt(P[src, dest] + eps_w)
  wpen <- wpen / mean(wpen)
  w <- as.list(c(0, wpen))
  names(w) <- c("own", src)
  own_id <- which(grp == "own")
  r0 <- if (length(own_id)) lm.fit(X[, own_id, drop = FALSE], y)$residuals else y
  lmax <- max(vapply(src, function(s) {
    id <- which(grp == s)
    sqrt(sum((crossprod(X[, id, drop = FALSE], r0) / n)^2)) / w[[s]]
  }, 0))
  grid <- exp(seq(log(lmax * 0.98), log(lmax * 0.02), length.out = 20))
  i0 <- floor(0.6 * n)
  bl <- split((i0 + 1):n, cut((i0 + 1):n, 4, labels = FALSE))
  cv <- sapply(grid, function(lam) mean(vapply(bl, function(b) {
    tr <- 1:(min(b) - 1)
    bb <- group_lasso_fit(X[tr, , drop = FALSE], y[tr], grp, w, lam)
    mean((y[b] - X[b, , drop = FALSE] %*% bb)^2)
  }, 0)))
  lam <- grid[which.min(cv)]
  beta <- group_lasso_fit(X, y, grp, w, lam)
  tr <- 1:floor((1 - hold_fr) * n)
  ho <- (max(tr) + 1):n
  bh <- group_lasso_fit(X[tr, , drop = FALSE], y[tr], grp, w, lam)
  mse_full <- mean((y[ho] - X[ho, , drop = FALSE] %*% bh)^2)
  out <- NULL
  for (s in src) {
    id <- which(grp == s)
    bs <- beta[id]
    if (sqrt(sum(bs^2)) <= 1e-6) next
    bt <- bh
    bt[id] <- 0
    mse_wo <- mean((y[ho] - X[ho, , drop = FALSE] %*% bt)^2)
    out <- rbind(out, data.frame(
      from = s, to = dest, P_ji = round(P[s, dest], 4),
      b1 = round(bs[1], 3), b2 = round(bs[2], 3), b3 = round(bs[3], 3),
      B_sd = round(sum(bs), 3),
      B_elast = round(sum(bs * sdy / sds[id]), 3),
      D_dR2 = round((mse_wo - mse_full) / var(y[ho]), 4)))
  }
  out
}
