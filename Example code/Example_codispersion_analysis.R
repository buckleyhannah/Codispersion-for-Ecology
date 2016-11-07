## Code for running codispersion analysis
## November 2016
## Two examples are given:
## 1. A point pattern with two marks where we show the codispersion of the two marks
## 2. Two raster patterns where we show the codispersion of the two rasters

## Load required packages and source code
library(ggplot2)
library(tidyr)
library(dplyr)
library(geoR)

source("Codispersion_analysis_functions.R")  # code that runs the codispersion kernal functions

### Example 1: Codispersion analysis of a point pattern with two marks (diameter and ant richness)
## Simulate a tree map with 250 randomly distributed trees of spatially random diameters in a 100 x 100 m plot; each tree has an associated ant species richness measurement (random). The relationship between tree diameter and ant richness is random (codispersion is close to zero)
set.seed(1)
dat <- data.frame(x = round(runif(250, min = 0, max = 100), 1), 
                  y = round(runif(250, min = 0, max = 100), 1), 
                  tree.diameter = rlnorm(n = 250, meanlog = 1.5, sdlog = 0.15), 
                  ant.richness = as.integer(runif(n = 250, min = 0, max = 10)) ) 

# Draw a map of the trees
ggplot(dat, aes(x, y, colour = tree.diameter, size = ant.richness)) + geom_point(alpha = 0.8) + coord_fixed(ratio = 1) + theme_bw() + ggtitle("random ant richness")

# Create geoR objects as input for codispersion analysis
geo.tree <- as.geodata(dat, coords.col = 1:2, data.col = 3)
geo.antr <- as.geodata(dat, coords.col = 1:2, data.col = 4)

# Run the codispersion analysis
codisp.output.ants <- codisp.fn(geo.tree, geo.antr, k = c(4, 4, 4), max.window.size = 100/4, lx = 20, ly = 10)

# Plot the codispersion graph
ggplot(codisp.output.ants, aes(x = X, y = Y, fill = Codispersion)) + scale_fill_gradientn(colours = c('#0000FF', '#FFFFFF', '#FF6666'),limits=c(-1, 1)) + coord_fixed(ratio = 1) + geom_tile() + theme_bw(base_size = 14) + ggtitle('Tree diameter vs. ant species richness') + xlab('Spatial lag in X (m)') + ylab('Spatial lag in Y (m)')

### Example 2: Codispersion analysis of a raster pattern (basal area of tree sp A and tree sp B)
## Simulate two raster patterns representing the abundance of two tree species in 20 x 20 m grid cells in a 300 x 300 m forest plot. The relationship between the two species is anisotropic across the plot.
X <- seq(from = 0, to = 280, by = 20)
Y <- seq(from = 0, to = 280, by = 20)
gridxy <- expand.grid(x = X , y = Y)
gridxy$spA <- 1 + (2 * gridxy$x + 5) / 10
gridxy$spB <- 1 + rev(((gridxy$x + 1)^2 + (gridxy$y + 1)^2) / 3000)

# Draw maps of the rasters
gather(gridxy, key = Species, value = Abundance, -x, -y) %>% ggplot(aes(x, y, fill = Abundance)) + geom_tile(col = 'grey50') + coord_fixed(ratio = 1) + facet_grid(~ Species) + theme_bw(base_size = 14) 

# Create geoR objects as input for codispersion analysis
geo.spA <- as.geodata(gridxy, coords.col = 1:2, data.col = 3)
geo.spB <- as.geodata(gridxy, coords.col = 1:2, data.col = 4)

# Run the codispersion analysis
codisp.output.trees <- codisp.fn(geo.spA, geo.spB, k = c(20, 20, 20), max.window.size = 300/4, lx = 20, ly = 10)

# Plot the codispersion graph
ggplot(codisp.output.trees, aes(x = X, y = Y, fill = Codispersion)) + scale_fill_gradientn(colours = c('#0000FF', '#FFFFFF', '#FF6666'),limits=c(-1, 1)) + coord_fixed(ratio = 1) + geom_tile() + theme_bw(base_size = 14) + ggtitle('Tree abundance: spA vs. spB') + xlab('Spatial lag in X (m)') + ylab('Spatial lag in Y (m)')




