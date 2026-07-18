# builds the raw graph based on visitation history
# edge candidates for the final Granger causality graph
build_visitation_graph <- function(events, gap_max = Inf, c0 = 0.01) {
  ev <- events[order(events$user, events$time, method = "radix"), ]
  nodes <- sort(unique(ev$location), method = "radix")
  nu <- ev$user[-1] == ev$user[-nrow(ev)]
  same_ep <- if (!is.null(ev$episode)) ev$episode[-1] == ev$episode[-nrow(ev)] else rep(FALSE, nrow(ev)-1)
  gap <- as.numeric(ev$time[-1]) - as.numeric(ev$time[-nrow(ev)])
  linked <- nu & (same_ep | (gap > 0 & gap <= gap_max))
  from <- ev$location[-nrow(ev)][linked]
  to <- ev$location[-1][linked]
  keep <- from != to
  from <- from[keep]
  to <- to[keep]
  Ntr <- table(factor(from, levels = nodes), factor(to, levels = nodes))
  Nj <- table(factor(ev$location, levels = nodes))
  P <- sweep(matrix(as.numeric(Ntr), length(nodes), length(nodes),
                     dimnames = list(nodes, nodes)),
             1, pmax(as.numeric(Nj), 1), "/")
  E0 <- P >= c0
  diag(E0) <- FALSE
  list(P = P, E0 = E0, nodes = nodes)
}
