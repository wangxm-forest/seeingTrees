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
# Read in the climate data

climate <- readRDS("C:/PhD/Project/PhD_thesis/seeingTrees/data/climate/wldas_climpredictors_updated.rds")

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
climate_AB08 <- climate[climate$plotname == "AB08",]
gdd <- climate_AB08$gdd_ings/100


data <- list(
  N_obs = length(x),
  x_obs = x,
  y_obs = y,
  temp_obs = gdd,
  
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


