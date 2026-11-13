###################Started by Mao####################
###################Oct 18,2026######################

rm(list = ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

library("dplR")
library("dplyr")
library("ggplot2")
library("rstan")
setwd("C:/PhD/Project/PhD_thesis/seeingTrees/analyses")
util <- new.env()
source("mcmc_analysis_tools_rstan.R", local=util)
source("mcmc_visualization_tools.R", local=util)

# Read in the WLDAS climate data
wldas <- readRDS("C:/PhD/Project/PhD_thesis/seeingTrees/data/climate/wldas_climpredictors_updated.rds")

# Read in the climateNA climate data
climateNA <- read.csv("C:/PhD/Project/PhD_thesis/seeingTrees/data/climate/climateNA/climateNAMORAStands.csv")

# Read in the ring width data
AB08THSE <- read.rwl("C:/PhD/Project/PhD_thesis/seeingTrees/data/treeCores/AB08_TSHE_dated.rwl", format = "auto")

treeID <- sub("\\(.*\\)$", "", colnames(AB08THSE))

# Average across cores from the same tree
ring <- sapply(split(seq_along(treeID), treeID), function(i) {
  if (length(i) == 1) {
    AB08THSE[, i]
  } else {
    rowMeans(AB08THSE[, i, drop = FALSE], na.rm = TRUE)
  }
})

rownames(ring) <- rownames(AB08THSE)

# Get data for only one tree
oneTree <- data.frame(
  year  = as.numeric(rownames(ring)),
  width = ring[, "4195"]
)
plot(width~year, oneTree, pch=16, cex=0.25)
lines(width~year, oneTree)

# subset to years with climate data
all_years <- 1980:2023
oneTree_1980 <- oneTree[oneTree$year %in% all_years, ]

x <- oneTree_1980$year-1000
y <- oneTree_1980$width*10


## subset to AB08
wldas_AB08 <- wldas[wldas$plotname == "AB08",]
gdd_wldas <- wldas_AB08$gdd_ings/100


data <- list(
  N_obs = length(x),
  x_obs = x,
  y_obs = y,
  temp_obs = gdd_wldas,
  
  N_pred = length(x),
  x_pred = x)


modelhere <- stan_model('stan/gpTreeGrowth.stan')
fit <- sampling(modelhere, data = data, iter = 4000, warmup = 3000, chains = 4, cores = 4)

diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples, c('rho', 'gamma', 'sigma'))
util$check_all_expectand_diagnostics(base_samples)

names <- paste0('y_pred[',1:data$N_pred,']')
util$plot_conn_pushforward_quantiles(samples, names, data$x_pred)
points(data$x_obs, data$y_obs, pch=16, cex=1, col="white")
points(data$x_obs, data$y_obs, pch=16, cex=0.5, col="black")



# Hard to see if plot the whole prediction, make it 800
idx <- which(data$x_pred >= 800)

names_sub <- paste0('y_pred[', idx, ']')
util$plot_conn_pushforward_quantiles(samples, names_sub, data$x_pred[idx])
points(data$x_obs, data$y_obs, pch=16, cex=1, col="white")
points(data$x_obs, data$y_obs, pch=16, cex=0.5, col="black")

# Only taking the range of real data
idx <- which(data$x_pred >= min(data$x_obs))

names_sub <- paste0('y_pred[', idx, ']')
shift <- 1000
util$plot_conn_pushforward_quantiles(samples, names_sub, data$x_pred[idx] + shift)
points(data$x_obs + shift, data$y_obs, pch=16, cex=1, col="white")
points(data$x_obs + shift, data$y_obs, pch=16, cex=0.5, col="black")


# Posterior predictive check ----
ring <- data.frame(
  year = rep(as.integer(rownames(ring)), times = ncol(ring)),
  tree = rep(colnames(ring), each = nrow(ring)),
  rw   = as.vector(as.matrix(ring))
)
ring <- ring[!is.na(ring$rw), ]

climateNA_AB08 <- climateNA[climateNA$stand_id == "AB08",]
climateNA_AB08$DD5 <- climateNA_AB08$DD5/100

gdd_climateNA <- data.frame(year = as.integer(climateNA_AB08[["year"]]), temp = climateNA_AB08[["DD5"]])

n_before <- nrow(ring)
ring <- merge(ring, gdd_climateNA, by = "year")
ring <- ring[!is.na(ring$temp), ]
ring <- ring[order(ring$tree, ring$year), ]

tree_id <- unique(ring$tree)
n_tree <- as.integer(table(factor(ring$tree, levels = tree_id)))

end <- cumsum(n_tree)
start <- end - n_tree + 1

all_year_obs <- seq_len(max(n_tree))

stan_data <- list(
  N_tree = length(n_tree),
  N_obs = nrow(ring),
  N_all_year_obs = length(all_year_obs),
  temp_obs = ring$temp,
  start = as.integer(start),
  end = as.integer(end),
  N_obs_tree = as.integer(n_tree),
  all_year_obs = as.integer(all_year_obs),
  year_obs = as.numeric(ring$year),
  log_rw_obs = log(as.numeric(ring$rw))
)

modelhere <- stan_model('stan/gpTreeGrowthMulti.stan')
fit <- sampling(modelhere, data = stan_data, iter = 4000, warmup = 3000, chains = 4, cores = 4)

diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples, c('rho', 'gamma', 'sigma'))
util$check_all_expectand_diagnostics(base_samples)

names(samples)

names <- c(grep('alpha', names(samples), value = TRUE),
           grep('beta', names(samples), value = TRUE),
           grep('rho', names(samples), value = TRUE),
           grep('gamma', names(samples), value = TRUE),
           grep('sigma', names(samples), value = TRUE))

print(fit, pars = names)

## Posterior predictive check ----
pdf("figures/climateNAAB08.pdf", width = 15, height = 15)
par(mfrow = c(2, 2), mar = c(4,4,1,1))

for (i in seq_len(stan_data$N_tree)) {
  idx <- stan_data$start[i]:stan_data$end[i]
  names_log <- paste0("log_rw_pred[", idx, "]")
  names_natural <- paste0("rw_pred[", idx, "]")
  
  for (k in seq_along(idx)){
    samples[[names_natural[k]]] <- exp(samples[[names_log[k]]])
  }
  util$plot_conn_pushforward_quantiles(samples, names_natural, stan_data$year_obs[idx],xlab = "Year", ylab = "Ring width", main = tree_id[i])
  
  rw_natural <- exp(stan_data$log_rw_obs[idx])
  points(stan_data$year_obs[idx], rw_natural, pch = 16, cex = 1,   col = "white")
  points(stan_data$year_obs[idx], rw_natural, pch = 16, cex = 0.5, col = "black")
}

dev.off()

## Posterior for the gp parameters ----
pdf("figures/climateNAAB08TsheParameters.pdf", width = 9, height = 3)
par(mfrow = c(1, 2), mar = c(4,4,1,1))

util$plot_expectand_pushforward(samples[["rho"]], 50,
                                display_name = "rho",
                                flim = c(-10, 10))

util$plot_expectand_pushforward(samples[["gamma"]], 25,
                                display_name = "gamma",
                                flim = c(-1, 1))

dev.off()