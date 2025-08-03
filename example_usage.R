# Test specific functions with NEW SEPARATED COLORING approach
# ==============================================================

# Load required libraries
library(lidRviewer)
library(lidR)
library(sf)
library(terra)
library(RColorBrewer)

# Load and process example data
LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
las <- readLAS(LASfile, select = "xyz")

# Create DEM, CHM and tree metrics
dem <- rasterize_terrain(las, tin(), res=1)
chm <- rasterize_canopy(las, 0.5, p2r(0.3))
ttops <- locate_trees(chm, lmf(4, 2))
las <- segment_trees(las, dalponte2016(chm, ttops))
metrics <- crown_metrics(las, .stdtreemetrics, geom = "convex")

# ===================================================
# NEW APPROACH: Separate LAS and SF coloring
# ===================================================

# 1. Color the LAS data separately using colorizeLAS
las_by_height <- colorizeLAS(las, "Z", lidR::height.colors)
las_by_trees <- colorizeLAS(las, "treeID", rainbow)

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

# ===================================================
# LEGACY APPROACH: Still works for backward compatibility
# ===================================================

# Test 1: Basic functionality (works)
cat("1. Testing basic plot_sf_3d...\n")
plot_sf_3d(metrics, dem)
plot_sf_3d(ttops, chm)

# Test 2: colorizeLAS (works)
cat("2. Testing colorizeLAS...\n")
view(las_by_height)

# Test 3: POINT geometry with vector palette
cat("3. Testing POINT with vector palette...\n")
plot_sf_3d(ttops, dem, color_by = "Z", palette = c("blue", "green", "red"))

# Test 4: POINT geometry with function palette (should work now)
cat("4. Testing POINT with function palette...\n")
plot_sf_3d(ttops, dem, color_by = "Z", palette = heat.colors)

# Test 5: Polygon lines
cat("5. Testing polygon lines...\n")
plot_sf_3d(metrics, dem,
          color_by = "treeID",
          palette = brewer.pal(8, "Dark2"))

# Polygons colored by height (Z) with custom gradient
plot_sf_3d(metrics, dem,
          color_by = "Z",
          palette = random.colors(50),
          density = 2)

# ========================================================
# 3. POINT GEOMETRIES (Tree tops)
# ========================================================

# Display tree tops as colored points
plot_sf_3d(ttops, dem,
          color_by = "Z",
          palette = heat.colors)

# ========================================================
# 4. DEM WIREFRAME VISUALIZATION
# ========================================================

# Create DEM wireframe using triangulation (like lidR's tin())
dem_wireframe <- create_dem_wireframe(dem, sample_factor = 3)

# Display DEM wireframe only
plot_sf_3d(dem_wireframe, dem,
          default_color = "#808080")

# ========================================================
# 5. COMPREHENSIVE COMBINED VISUALIZATION
# ========================================================

# Enhanced view() function with multiple components:

# Method 1: Point cloud + tree crowns + tree tops
view(filter_poi(las, Z > 0.5),
     polygons = metrics,
     dem = dem,
     color_by = "treeID",
     palette = rainbow(50))

# Method 2: Colored point cloud + colored polygons
view(las_by_height,
     polygons = metrics,
     dem = dem,
     color_by = "Z",
     palette = c("blue", "green", "yellow", "red"))

# ========================================================
# 6. CREATING CUSTOM GEOMETRY TYPES FOR TESTING
# ========================================================

# Create sample LINESTRING geometries (e.g., transects, roads)
create_sample_lines <- function(dem) {
  ext_dem <- terra::ext(dem)

  # Create some sample transect lines
  line1 <- matrix(c(
    ext_dem[1] + 10, ext_dem[3] + 10,
    ext_dem[2] - 10, ext_dem[4] - 10
  ), ncol = 2, byrow = TRUE)

  line2 <- matrix(c(
    ext_dem[1] + 20, ext_dem[4] - 20,
    ext_dem[2] - 20, ext_dem[3] + 20
  ), ncol = 2, byrow = TRUE)

  sample_lines <- st_sf(
    id = 1:2,
    type = c("transect", "road"),
    width = c(2, 5),
    geometry = st_sfc(st_linestring(line1), st_linestring(line2))
  )

  st_crs(sample_lines) <- st_crs(dem)
  return(sample_lines)
}

# Create sample POINT geometries (e.g., sensors, landmarks)
create_sample_points <- function(dem) {
  ext_dem <- terra::ext(dem)

  # Create some sample point locations
  sample_points <- st_sf(
    id = 1:4,
    type = c("sensor", "tower", "benchmark", "sensor"),
    elevation = c(5, 15, 0, 8),
    geometry = st_sfc(
      st_point(c(ext_dem[1] + 30, ext_dem[3] + 30)),
      st_point(c(ext_dem[2] - 30, ext_dem[4] - 30)),
      st_point(c(ext_dem[1] + 50, ext_dem[4] - 50)),
      st_point(c(ext_dem[2] - 50, ext_dem[3] + 50))
    )
  )

  st_crs(sample_points) <- st_crs(dem)
  return(sample_points)
}

# Create sample geometries
sample_lines <- create_sample_lines(dem)
sample_points <- create_sample_points(dem)

# Display different geometry types with custom styling:

# LINESTRING visualization
plot_sf_3d(sample_lines, dem,
          color_by = "type",
          palette = c("blue", "red"))

# POINT visualization
plot_sf_3d(sample_points, dem,
          color_by = "type",
          palette = c("red", "blue", "green", "orange"))

# ========================================================
# 7. ADVANCED COMBINED SCENARIOS
# ========================================================

# Scenario A: Forest analysis with all components
view(las_by_trees,                    # Colored point cloud
     polygons = metrics,              # Tree crowns
     dem = dem,
     color_by = "treeID",
     palette = brewer.pal(8, "Set3"))

# Scenario B: Terrain analysis with wireframe
view(filter_poi(las, Z < 2),          # Ground points only
     polygons = dem_wireframe,        # DEM wireframe
     dem = dem,
     default_color = "#606060")

# Scenario C: Infrastructure mapping
view(las,                            # Full point cloud
     polygons = sample_lines,         # Roads/transects
     dem = dem,
     color_by = "type",
     palette = c("red", "blue"))

# ========================================================
# 8. UTILITY FUNCTIONS AND TIPS
# ========================================================

# Function to create custom color palettes
create_forest_palette <- function(n) {
  colorRampPalette(c("darkgreen", "forestgreen", "lightgreen", "yellow"))(n)
}

# Function to create elevation palette
create_elevation_palette <- function(n) {
  colorRampPalette(c("blue", "cyan", "green", "yellow", "orange", "red"))(n)
}

# Usage examples:
plot_sf_3d(metrics, dem, color_by = "Z", palette = create_elevation_palette)
plot_sf_3d(ttops, dem, color_by = "Z", palette = create_forest_palette)

# ========================================================
# 9. PERFORMANCE AND MEMORY TIPS
# ========================================================

# For large datasets:
# 1. Filter point clouds before visualization:
   las_filtered <- filter_poi(las, Z > 0.5 & Z < 30)
#
# 2. Reduce DEM resolution for wireframe:
   dem_wireframe <- create_dem_wireframe(dem, sample_factor = 5)
#
# 3. Use detached mode for non-blocking operation:
   plot_sf_3d(metrics, dem, detach = TRUE)  # <- Don't do this!
#
# 4. Limit polygon density for smoother rendering:
   view(las, polygons = metrics, dem = dem, density = 0.5)

# ========================================================
# NOTES AND FEATURES SUMMARY
# ========================================================

# NEW FEATURES:
# colorizeLAS() - Color point clouds by any attribute with custom palettes
# plot_sf_3d() - Universal sf geometry viewer (POINT, LINESTRING, POLYGON)
# create_dem_wireframe() - DEM triangulation wireframe using tin()-like algorithm
# Enhanced view() - Supports all new color-by and styling parameters
# Color-by-attribute - Any field with custom color palettes
# Custom styling - custom colors and density control
# Multiple geometry support - Mixed visualization of different geometry types

# VIEWER CONTROLS:
# - Mouse: Left = rotate, Right = pan, Wheel = zoom
# - Keyboard: r/g/b = RGB colors, z = height, i = intensity, c = classification
# - Keyboard: +/- = point size, l = lighting toggle
# - IMPORTANT: Only one viewer window at a time - close before opening another

# PARAMETERS:
# - color_by: Column name for color mapping
# - palette: Color vector or function (e.g., rainbow, brewer.pal, heat.colors)
# - density: Line sampling density for smoother curves (default: 1)
# - detach: Non-blocking viewer mode (default: FALSE)
