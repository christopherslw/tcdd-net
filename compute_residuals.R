compute_residuals <- function(Y, hours, q_own = 3) {
  ylog <- log1p(Y)
  lt <- as.POSIXlt(hours)
  hod <- as.integer(lt$hour)
  dow <- as.integer(((lt$wday + 6) %% 7) + 1)
  S <- model.matrix(~ factor(hod) + factor(dow))
  lagv <- function(x, l) c(rep(NA, l), x[1:(length(x) - l)])
  rows <- (q_own + 1):nrow(Y)
  p <- ncol(Y)
  E <- matrix(NA_real_, length(rows), p, dimnames = list(NULL, colnames(Y)))
  for (i in 1:p) {
    Xi <- cbind(S, sapply(1:q_own, function(l) lagv(ylog[, i], l)))
    f <- lm.fit(Xi[rows, ], ylog[rows, i])
    E[, i] <- f$residuals
  }
  U <- E - (rowSums(E) - E) / (p - 1)
  list(U = U, hours = hours[rows], ylog = ylog)
}
