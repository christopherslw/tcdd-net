# test the tajectory-constrained dynamic dependence network approach on real mobility data 
# from Kincade Fire, Sonoma County, 2019
# source: Xu et al., 2022
# https://github.com/EvacuationBehavior/Highway-Routing-Data-Processing


source("build_visitation_graph.R")
source("compute_residuals.R")
source("fit_destination.R")
source("fit_tcdd_net.R")
source("group_lasso_fit.R")
source("activated_edges.R")

# config parameters #
tz_ <- "America/Los_Angeles"
L <- 3
q_own <- 3
gap_max <- 8 * 3600
c0 <- 0.01
eps_w <- 0.005
hold_fr <- 0.25
set.seed(1)
ignition <- as.POSIXct("2019-10-23 21:27", tz = tz_) #approximate ignition point according to SFTimes, 2019
evac <- as.POSIXct("2019-10-26 04:00", tz = tz_)  # evacuation orders according to Xu et al, 2022
contained <- as.POSIXct("2019-11-06 19:00", tz = tz_)

d <- read.csv("data/Highway_Vehicle_Routing_Data_R1.csv", stringsAsFactors=F)
ms2t <- function(x) as.POSIXct(x / 1000, origin = "1970-01-01", tz = tz_)
d$t_in <- ms2t(d$Entrance_TIMESTAMP)
d$t_out <- ms2t(d$Exit_TIMESTAMP)
node_of <- function(hwy, lat, lon) {
  ifelse(hwy == "Hwy 101", ifelse(lat >= 38.55, "101 N", ifelse(lat >= 38.35, "101 C", "101 S")),
  ifelse(hwy == "Hwy 12", ifelse(lon < -122.65, "12 W", "12 E"),
  ifelse(hwy == "Hwy 116", ifelse(lon < -122.90, "116 W", "116 E"), hwy)))
}
d$node_in <- node_of(d$Entrance_HWY, d$Entrance_LAT, d$Entrance_LON)
d$node_out <- node_of(d$Exit_HWY, d$Exit_LAT, d$Exit_LON)
d$episode <- seq_len(nrow(d))
events <- rbind(
  data.frame(user = d$ID, location = d$node_in, time = d$t_in, episode = d$episode,
             lat = d$Entrance_LAT, lon = d$Entrance_LON),
  data.frame(user = d$ID, location = d$node_out, time = d$t_out, episode = d$episode,
             lat = d$Exit_LAT, lon = d$Exit_LON))
nodes <- sort(unique(events$location), method = "radix")
p <- length(nodes)
cat("Locations and event counts:\n")
print(table(events$location))
expected <- c("101 C", "101 N", "101 S", "116 E", "116 W", "12 E", "12 W", "Hwy 1", "Hwy 128", "Hwy 37")
xy <- t(sapply(nodes, function(nd)
  c(lon = mean(events$lon[events$location == nd]), lat = mean(events$lat[events$location == nd]))))


## hourly activity panel Y 
hsec <- 3600
h0 <- floor(min(as.numeric(events$time)) / hsec)
hidx <- floor(as.numeric(events$time) / hsec) - h0 + 1
Tn <- max(hidx)
hours <- as.POSIXct((h0 + 0:(Tn - 1)) * hsec, origin = "1970-01-01", tz = tz_)
tab <- table(factor(hidx, levels = 1:Tn), factor(events$location, levels = nodes))
Y <- matrix(as.numeric(tab), Tn, p, dimnames = list(NULL, nodes))

# step 1
g <- build_visitation_graph(events[events$time < ignition, ], gap_max, c0)
cat(sprintf("Opportunity graph: %d of %d ordered pairs eligible (c0=%.2f)\n",
            sum(g$E0), p * (p - 1), c0))

# step 2
res <- compute_residuals(Y, hours, q_own)
regimes <- list(pre = which(res$hours < ignition), fire = which(res$hours >= ignition))
cat(sprintf("Regime sizes (hours): pre=%d, fire+recovery=%d\n",
            length(regimes$pre), length(regimes$fire)))

# step 3
edges <- fit_tcdd_net(res$U, g$E0, g$P, regimes, L, eps_w, hold_fr)

# step 4
r <- activated_edges(edges, "pre", "fire", b_floor = 0.02) 
print(r$baseline, row.names=F)  # pre fire baseline
print(r$disrupted, row.names=F)  # post fire network
cat(sprintf("\nEdges Activated in fire: %d | dropped: %d\n", nrow(r$activated), nrow(r$dropped)))


