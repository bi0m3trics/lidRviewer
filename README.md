![Github](https://img.shields.io/badge/Github-1.0.0-green.svg) ![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg)

# **My Fork** - lidRviewer - **My Fork**

This is a point cloud viewer for R. The primary goal of this package is to serve as an alternative backend for displaying point clouds in the [lidR](https://github.com/Jean-Romain/lidR) package, replacing `rgl`.

While `rgl` is a powerful package, it has some limitations when it comes to handling large point clouds. The `lidRviewer` package is designed to efficiently display arbitrarily large in-memory point clouds, with tested cases including over 880 million points.

<img width="2028" height="700" alt="image" src="https://github.com/user-attachments/assets/316e6c1d-a86a-4b14-a29d-823e498b255a" />

### Advantages of `lidRviewer`:

* Capable of handling hundred millions of points efficiently if they fits in memory, unlike `rgl`, which may struggle with a few million.
* More memory efficient: it requires less additional memory than `rgl` that may require gigabytes to store a full copy of the point cloud.
* Provides keyboard controls for changing coloring attributes and point size on-the-fly.
* Features eyes-dome lighting for better visualization.

### Drawbacks:

* Limited functionality: only displays point clouds and does not offer additional features. It is not intended to replace `rgl` for other visualization needs and custom rendering.
* Limited to only one windows and thus one rendering.

## Installation

```r
install.packages('lidRviewer', repos = c('https://r-lidar.r-universe.dev'))
```

### MacOS

```
brew install sdl2 mesa mesa-glu
```

## Usage

### Point Clouds

```r
library(lidR)
library(lidRviewer)
LASfile <- system.file("extdata", "Megaplot.laz", package="lidR")
las <- readLAS(LASfile)
view(las)
# Close the view to continue working in your R session
```

### 3D Polygon Visualization (NEW - what I've added)

You can now visualize sf polygon objects in 3D by extracting Z values from a SpatRaster DEM and you can colorize the las by any attribute:

<img width="1296" height="700" alt="image" src="https://github.com/user-attachments/assets/39357898-4f31-4135-8c25-3b55a62ccb7c" />

```r
library(lidRviewer)
library(lidR
library(sf)
library(terra)

LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
las <- readLAS(LASfile, select = "xyz"
las <- classify_ground(las, csf())
dem <- rasterize_terrain(las, tin(), res=1)
las<-normalize_height(las, tin())
chm <- rasterize_canopy(las, 0.5, p2r(0.3))
ttops <- locate_trees(chm, lmf(4, 2))
las <- segment_trees(las, dalponte2016(chm, ttops))
metrics <- crown_metrics(las, .stdtreemetrics, geom = "convex")

# ===================================================
# NEW APPROACH: Separate LAS and SF coloring
# ===================================================
# 1. Color the LAS data separately using colorizeLAS
las_by_height <- colorizeLAS(las, "Z", lidR::height.colors)
las_by_trees <- colorizeLAS(las, "treeID", function(n) random.colors(500))

# 2. Use the new sf_objects + sf_styling approach for multiple SF objects
view(las_by_height,
     sf_objects = list(
       polygons = metrics,
       points = ttops
     ),
     sf_styling = list(
       polygons = list(color_by = "treeID", palette = rainbow(50)),
       points = list(color_by = "Z", palette = "#FFFFFF")
     ),
     dem = dem)
```

The polygons will be displayed as white wireframes in 3D space, with Z coordinates extracted from the provided DEM at polygon vertex locations.

- Rotate with left mouse button
- Zoom with mouse wheel
- Pan with right mouse button
- Keyboard <kbd>r</kbd> or <kbd>g</kbd> or <kbd>b</kbd> to color with RGB
- Keyboard <kbd>z</kbd> to color with Z
- Keyboard <kbd>i</kbd> to color with Intensity
- Keyboard <kbd>c</kbd> to color with Classification
- Keyboard <kbd>+</kbd> or <kbd>-</kbd> to change the point size (on the numberpad)
- Keyboard <kbd>l</kbd> to enable/disable eyes-dome lightning

## Technical details

`lidRviewer` is based on [Markus Schultz thesis](https://www.cg.tuwien.ac.at/research/publications/2016/SCHUETZ-2016-POT/) with some adaptation and variation. One of the main difference is that Potree spatially indexes the point cloud in dedicated and optimized files on-disk file for an out-of-core rendering. `lidRviewer` on its side creates a nested octree on-the-fly on an in-memory `data.frame` without modifying the original data (not sorting, no data layout optimization).
