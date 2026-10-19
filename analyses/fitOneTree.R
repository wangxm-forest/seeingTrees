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
class(ring)

# Get data for only one tree
oneTree <- data.frame(
  year  = as.numeric(rownames(ring)),
  width = ring[, "4195"]
)
plot(width~year, oneTree, pch=16, cex=0.25)
lines(width~year, oneTree)

x <- oneTree$year-1000
y <- oneTree$width*10

data <- list(
  N_obs = length(x),
  x_obs = x,
  y_obs = y,
  
  N_pred = 100,
  x_pred = seq(100, 1000, length.out = 100)
)

modelhere <- stan_model('stan/gpTreeGrowth.stan')
fit <- sampling(modelhere, data = data, chains = 4, cores = 4)

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
util$plot_conn_pushforward_quantiles(samples, names_sub, data$x_pred[idx])
points(data$x_obs, data$y_obs, pch=16, cex=1, col="white")
points(data$x_obs, data$y_obs, pch=16, cex=0.5, col="black")

