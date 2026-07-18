# build the graph
activated_edges <- function(edges, from_regime, to_regime, b_floor = 0.02) {
  key <- function(e) paste(e$from, "->", e$to)
  e0 <- edges[edges$regime == from_regime & abs(edges$B_sd) >= b_floor, ]
  e1 <- edges[edges$regime == to_regime & abs(edges$B_sd) >= b_floor, ]
  list(baseline = e0, disrupted = e1,
       activated = e1[!key(e1) %in% key(e0), ],
       dropped = e0[!key(e0) %in% key(e1), ])
}