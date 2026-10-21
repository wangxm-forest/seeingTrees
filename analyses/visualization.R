###################Started by Mao####################
###################Oct 19,2026######################

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

## subset to AB08
climate_AB08 <- wldas_climpredictors_updated[wldas_climpredictors_updated$plotname == "AB08",]

oneTree$width_std <- as.numeric(scale(oneTree$width))
climate_AB08$gdd_std <- as.numeric(scale(climate_AB08$gdd_ings))
climate_AB08$soil_std <- as.numeric(scale(climate_AB08$soilmoist_ings))

plot(width_std ~ year, data = oneTree,
     pch = 16, cex = 0.25,
     col = "black",
     xlab = "Year", ylab = "Standardized value",
     ylim = range(c(oneTree$width_std, climate_AB08$gdd_std),
                  na.rm = TRUE))

lines(width_std ~ year, data = oneTree, col = "black")

# gdd
points(gdd_std ~ year, data = climate_AB08,
       pch = 16, cex = 0.5, col = "red")

lines(gdd_std ~ year, data = climate_AB08, col = "red")

# soil moisture
points(soil_std ~ year, data = climate_AB08,
       pch =16, cex = 0.5, col = "lightgreen")

lines(soil_std ~ year, data = climate_AB08, col = "lightgreen")
