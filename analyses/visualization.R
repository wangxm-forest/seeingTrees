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


# Plotting with Ailene's data
## Asked AI to help me make this into a loop so I can run for all the datasets
datasets <- c("AE10Abam", "AE10Tsme", "AE10Xano", "AG05Abam", "AG05Psme", "AG05Thpl", "AG05Tshe", "AM16Abam", "AM16Tshe", "AM16Tsme", "AM16Xano", "AR07Abam", "AR07Tsme", "AR07Xano", "AV06Abam", "AV06Psme", "AV06Thpl", "AV06Tshe", "AX15Psme", "AX15Thpl", "AX15Tshe", "PARAAbam", "PARATsme", "PARAXano", "TB13Abam", "TB13Psme", "TB13Thpl", "TB13Tshe", "TO04Abam", "TO04Psme", "TO04Thpl", "TO04Tshe") 

input_dir <- "C:/PhD/Project/PhD_thesis/seeingTrees/data/ailene/rwlFiles"
output_dir <- "figures"

files <- list.files(input_dir, pattern = "\\.rwl$", full.names = TRUE)
file_stands <- substr(basename(files), 1, 4)


stands <- unique(file_stands)

for (stand in stands) {
  
  stand_files <- files[file_stands == stand]
  
  pdf(file.path(output_dir, paste0(stand, ".pdf")),
      width = 8, height = 6)
  par(mfrow = c(2, 3))
  
  for (file_path in stand_files) {
    
    species <- substr(basename(file_path), 6, 9)
    
    rwl_data <- read.rwl(file_path)

    treeID <- sub("[ab]$", "", colnames(rwl_data))
    
    ring <- sapply(split(seq_along(treeID), treeID), function(i) {
      if (length(i) == 1) {
        rwl_data[, i]
      } else {
        rowMeans(rwl_data[, i, drop = FALSE], na.rm = TRUE)
      }
    })
    
    rownames(ring) <- rownames(rwl_data)
    ring <- as.data.frame(ring, check.names = FALSE)
    
    growth_mat <- t(as.matrix(ring))
    years <- as.numeric(colnames(growth_mat))
    
    for (i in 1:nrow(growth_mat)) {
      plot(
        years,
        growth_mat[i, ],
        type = "l",
        xlab = "Year",
        ylab = "Growth (ring width mm)",
        main = paste(rownames(growth_mat)[i], species)
      )
    }
  }
  
  dev.off()
}
