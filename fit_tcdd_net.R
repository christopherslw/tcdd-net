# final fit on the trajectory constrained dependence network
# Granger causality interpretation under suitable conditions
# general interpretation: an edge i->j means i contains information that significantly
# improves prediction of mobility activity at node j.
fit_tcdd_net <- function(U, E0, P, regimes, L = 3, eps_w = 0.005, hold_fr = 0.25) {
  nodes <- colnames(U)
  lagv <- function(x, l) c(rep(NA, l), x[1:(length(x) - l)])
  Ulag <- lapply(1:L, function(l) apply(U, 2, lagv, l = l))
  edges <- NULL
  for (rg in names(regimes))
    for (dest in nodes) {
      e <- fit_destination(dest, regimes[[rg]], U, Ulag, E0, P, L, eps_w, hold_fr)
      if (!is.null(e)) edges <- rbind(edges, cbind(regime = rg, e))
    }
  edges[order(edges$regime, -abs(edges$B_sd)), ]
}

