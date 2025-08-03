#' Display big 3D point clouds with optional sf geometry overlay
#'
#' Display arbitrary large in memory 3D point clouds from the lidR package with optional 
#' sf geometry overlay. Now supports separate styling for different geometry types.
#' Keyboard can be used to control the rendering:
#' - Rotate with left mouse button
#' - Zoom with mouse wheel
#' - Pan with right mouse button
#' - Keyboard <kbd>r</kbd> or <kbd>g</kbd> or <kbd>b</kbd> to color with RGB
#' - Keyboard <kbd>z</kbd> to color with Z
#' - Keyboard <kbd>i</kbd> to color with Intensity
#' - Keyboard <kbd>c</kbd> to color with Classification
#' - Keyboard <kbd>+</kbd> or <kbd>-</kbd> to change the point size
#' - Keyboard <kbd>l</kbd> to enable/disable eyes-dome lightning
#'
#' @param x a point cloud with minimally 3 columns named X,Y,Z, or an sf object with supported geometries
#' @param polygons optional sf object with geometries to overlay on the point cloud (legacy parameter)
#' @param dem SpatRaster object to extract Z values from (required when displaying sf geometries)
#' @param density numeric. Point density for polygon outline sampling when displaying polygons (points per unit distance). Default is 1.
#' @param color_by character. Column name to color by (legacy parameter). 
#' @param palette character vector or function. Color palette (legacy parameter).
#' @param default_color character. Default color in hex format (legacy parameter). Default is "#FFFFFF".
#' @param sf_objects named list. List of sf objects to display with separate styling. E.g., list(polygons = metrics, points = ttops, lines = roads)
#' @param sf_styling named list. Styling parameters for each sf object type. E.g., list(polygons = list(color_by = "treeID", palette = rainbow), points = list(color_by = "Z", palette = heat.colors))
#' @param ... Support detach = TRUE
#' @export
#' @importClassesFrom lidR LAS
#' @useDynLib lidRviewer, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom sf st_coordinates st_cast st_geometry st_as_sf st_crs
#' @importFrom terra extract vect
#' @examples
#' \dontrun{
#' # View LAS point cloud only
#' view(las)
#' 
#' # Legacy: View LAS point cloud with polygon overlay
#' view(las, metrics, dem)
#' 
#' # New: Multiple sf objects with separate styling
#' view(las,
#'      sf_objects = list(
#'        polygons = metrics,
#'        points = ttops,
#'        lines = sample_lines
#'      ),
#'      sf_styling = list(
#'        polygons = list(color_by = "treeID", palette = rainbow),
#'        points = list(color_by = "Z", palette = heat.colors),
#'        lines = list(color_by = "type", palette = c("red", "blue"))
#'      ),
#'      dem = dem)
#' 
#' # View sf geometries only
#' view(metrics, dem)
#' }
#' @md

# Helper function to create line segments from polygon vertices
create_line_segments <- function(coords, density = 1, min_points_per_unit = 10) {
  if (nrow(coords) < 2) return(coords)
  
  # Ensure the polygon is closed (first point = last point)
  if (!all(coords[1, ] == coords[nrow(coords), ])) {
    coords <- rbind(coords, coords[1, ])
  }
  
  all_line_points <- list()
  
  for (i in 1:(nrow(coords) - 1)) {
    x1 <- coords[i, "X"]
    y1 <- coords[i, "Y"] 
    x2 <- coords[i + 1, "X"]
    y2 <- coords[i + 1, "Y"]
    
    # Calculate distance between points
    distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
    
    # Determine number of interpolation points
    # Use density parameter but ensure minimum density for line visibility
    n_points <- max(min_points_per_unit, round(distance * density * min_points_per_unit))
    
    if (n_points > 1 && distance > 0) {
      # Create interpolated points along the line segment
      t_seq <- seq(0, 1, length.out = n_points)
      x_interp <- x1 + t_seq * (x2 - x1)
      y_interp <- y1 + t_seq * (y2 - y1)
      
      # Store all points except the last one (to avoid duplication)
      segment_points <- cbind(X = x_interp[-length(x_interp)], 
                             Y = y_interp[-length(y_interp)])
      all_line_points[[length(all_line_points) + 1]] <- segment_points
    }
  }
  
  if (length(all_line_points) > 0) {
    return(do.call(rbind, all_line_points))
  } else {
    return(coords[1:(nrow(coords)-1), , drop = FALSE])  # Return original points minus last duplicate
  }
}

# Helper function to create line segments with line IDs for proper line rendering
create_line_segments_with_ids <- function(coords, density = 1, min_points_per_unit = 2, line_id_start = 1) {
  if (nrow(coords) < 2) return(list(coords = coords, line_ids = rep(line_id_start, nrow(coords))))
  
  # Ensure the polygon is closed (first point = last point)
  if (!all(coords[1, ] == coords[nrow(coords), ])) {
    coords <- rbind(coords, coords[1, ])
  }
  
  all_line_points <- list()
  all_line_ids <- list()
  current_line_id <- line_id_start
  
  for (i in 1:(nrow(coords) - 1)) {
    x1 <- coords[i, "X"]
    y1 <- coords[i, "Y"] 
    x2 <- coords[i + 1, "X"]
    y2 <- coords[i + 1, "Y"]
    
    # Calculate distance between points
    distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
    
    # Determine number of interpolation points
    n_points <- max(min_points_per_unit, round(distance * density * min_points_per_unit))
    
    if (n_points > 1 && distance > 0) {
      # Create interpolated points along the line segment
      t_seq <- seq(0, 1, length.out = n_points)
      x_interp <- x1 + t_seq * (x2 - x1)
      y_interp <- y1 + t_seq * (y2 - y1)
      
      # Store all points for this segment
      segment_points <- cbind(X = x_interp, Y = y_interp)
      segment_ids <- rep(current_line_id, length(x_interp))
      
      all_line_points[[length(all_line_points) + 1]] <- segment_points
      all_line_ids[[length(all_line_ids) + 1]] <- segment_ids
      current_line_id <- current_line_id + 1
    }
  }
  
  if (length(all_line_points) > 0) {
    combined_coords <- do.call(rbind, all_line_points)
    combined_ids <- do.call(c, all_line_ids)
    return(list(coords = combined_coords, line_ids = combined_ids))
  } else {
    original_coords <- coords[1:(nrow(coords)-1), , drop = FALSE]
    return(list(coords = original_coords, line_ids = rep(line_id_start, nrow(original_coords))))
  }
}

#' Display big 3D point clouds with optional polygon overlay
#'
#' Enhanced view function supporting both LAS point clouds and sf objects with
#' optional color-by-attribute functionality and custom styling parameters.
#'
#' @param x a point cloud with minimally 3 columns named X,Y,Z, or an sf object
#' @param polygons optional sf object with geometries to overlay (when x is a LAS object)
#' @param dem SpatRaster object to extract Z values from (required when x is sf or for polygons)
#' @param density numeric. Point density for line sampling. Default is 1.
#' @param color_by character. Column name to color by (for sf objects)
#' @param palette color vector or function for color mapping
#' @param point_size numeric. Size for POINT geometries. Default is 3.
#' @param line_width numeric. Width for lines/polygon outlines. Default is 2.
#' @param default_color character. Default color when no color_by specified. Default "#FFFFFF".
#' @param ... Additional arguments including detach = TRUE for non-blocking mode
#' @export
#' @examples
#' \dontrun{
#' # View LAS point cloud only
#' view(las)
#' 
#' # View sf polygons with color-by-attribute
#' view(metrics, dem = dem, color_by = "treeID", palette = rainbow(10))
#' 
#' # View LAS + sf with custom colors
#' view(las, polygons = metrics, dem = dem, color_by = "Z", palette = heat.colors)
#' }
#' @md
view = function(x, polygons = NULL, dem = NULL, density = 1, 
                # Legacy parameters for backward compatibility
                color_by = NULL, palette = NULL, default_color = "#FFFFFF",
                # New separate coloring parameters
                sf_objects = NULL, sf_styling = NULL,
                ...)
{
  p = list(...)
  detach = isTRUE(p$detach)
  
  # Check if x is an sf object
  if (inherits(x, "sf")) {
    # When x is sf, dem should be provided either as second positional arg or named dem parameter
    dem_param <- if (!is.null(polygons)) {
      # Check if polygons parameter is actually a raster (when called as view(sf, raster))
      if (inherits(polygons, "SpatRaster") || inherits(polygons, "RasterLayer") || inherits(polygons, "RasterStack")) {
        polygons  # dem passed as second positional argument
      } else {
        dem       # polygons is actually polygons, dem should be the third parameter
      }
    } else {
      dem       # dem passed as named parameter (view(sf, dem = raster))
    }
    
    if (is.null(dem_param)) {
      stop("A SpatRaster DEM must be provided when viewing sf objects")
    }
    
    # Use the new plot_sf_3d function
    plot_sf_3d(x, dem_param, density = density, color_by = color_by, palette = palette,
               default_color = default_color, detach = detach)
    return()
  }
  
  # x is a LAS object - check for different SF object approaches
  las_data <- x@data
  
  # New approach: Multiple SF objects with separate styling
  if (!is.null(sf_objects)) {
    if (is.null(dem)) {
      stop("A SpatRaster DEM must be provided when displaying sf objects with LAS data")
    }
    
    all_sf_data <- list()
    current_line_id <- 1
    
    # Process each SF object type
    for (obj_name in names(sf_objects)) {
      sf_obj <- sf_objects[[obj_name]]
      if (!inherits(sf_obj, "sf")) {
        warning(paste("Object", obj_name, "is not an sf object, skipping"))
        next
      }
      
      # Get styling for this object
      obj_styling <- if (!is.null(sf_styling) && obj_name %in% names(sf_styling)) {
        sf_styling[[obj_name]]
      } else {
        list()  # Use defaults
      }
      
      # Extract styling parameters with defaults
      obj_color_by <- obj_styling$color_by
      obj_palette <- obj_styling$palette
      obj_default_color <- if (!is.null(obj_styling$default_color)) obj_styling$default_color else "#FFFFFF"
      obj_density <- if (!is.null(obj_styling$density)) obj_styling$density else density
      
      # Process this SF object
      sf_result <- plot_sf_3d(sf_obj, dem, density = obj_density, 
                             color_by = obj_color_by, palette = obj_palette,
                             default_color = obj_default_color, 
                             detach = FALSE, return_data = TRUE)
      
      if (!is.null(sf_result)) {
        # Adjust line IDs to avoid conflicts
        if (!is.null(sf_result$line_ids)) {
          sf_result$line_ids <- sf_result$line_ids + current_line_id - 1
          current_line_id <- max(sf_result$line_ids, na.rm = TRUE) + 1
        }
        all_sf_data[[obj_name]] <- sf_result
      }
    }
    
    # Combine all SF data
    if (length(all_sf_data) > 0) {
      combined_sf_coords <- do.call(rbind, lapply(all_sf_data, function(x) x$coords))
      combined_sf_colors <- do.call(rbind, lapply(all_sf_data, function(x) x$colors))
      combined_sf_line_ids <- do.call(c, lapply(all_sf_data, function(x) x$line_ids))
      
      # Handle missing RGB columns in LAS data
      las_r <- if ("R" %in% names(las_data)) las_data$R else rep(255, nrow(las_data))
      las_g <- if ("G" %in% names(las_data)) las_data$G else rep(255, nrow(las_data))
      las_b <- if ("B" %in% names(las_data)) las_data$B else rep(255, nrow(las_data))
      
      combined_coords <- rbind(
        cbind(las_data[, c("X", "Y", "Z")], R = las_r, G = las_g, B = las_b, LINE_ID = NA),
        cbind(combined_sf_coords, combined_sf_colors, LINE_ID = combined_sf_line_ids)
      )
      
      viewer(combined_coords, detach, "")
      return()
    }
  }
  
  # Legacy approach: Single SF object (for backward compatibility)
  if (!is.null(polygons) && inherits(polygons, "sf")) {
    # Mixed LAS + sf display
    if (is.null(dem)) {
      stop("A SpatRaster DEM must be provided when displaying polygons with LAS data")
    }
    
    # Process sf geometry to get coordinates
    sf_result <- plot_sf_3d(polygons, dem, density = density, color_by = color_by, 
                           palette = palette, default_color = default_color, 
                           detach = FALSE, return_data = TRUE)
    
    if (!is.null(sf_result)) {
      # Combine LAS data with sf line data
      # Handle missing RGB columns in LAS data
      las_r <- if ("R" %in% names(las_data)) las_data$R else rep(255, nrow(las_data))
      las_g <- if ("G" %in% names(las_data)) las_data$G else rep(255, nrow(las_data))
      las_b <- if ("B" %in% names(las_data)) las_data$B else rep(255, nrow(las_data))
      
      combined_coords <- rbind(
        cbind(las_data[, c("X", "Y", "Z")], R = las_r, G = las_g, B = las_b, LINE_ID = NA),
        cbind(sf_result$coords, sf_result$colors, LINE_ID = sf_result$line_ids)
      )
      
      viewer(combined_coords, detach, "")
    } else {
      # Fall back to LAS only
      viewer(las_data, detach, "")
    }
  } else {
    # Original LAS-only behavior
    viewer(las_data, detach, "")
  }
}

render = function(f)
{
  x = normalizePath(x)
  las = lidR::readLAS(x)
  hnof = paste0(substr(x, 1, nchar(x) - 3), "hno")
  f = if (file.exists(hnof)) hnof else x
  viewer(las@data, FALSE, f)
}


#' Deprecated backward compatible function
#'
#' @param x,y,z numeric vector
#' @param r,g,b integer vector
#' @param id index
#' @param size not used
#' @export
plot_xyzrgb <- function(x, y, z, r, g, b, id = NULL, size = 4)
{
  xtxt = deparse(substitute(x))
  ytxt = deparse(substitute(y))
  ztxt = deparse(substitute(z))
  rtxt = deparse(substitute(r))
  gtxt = deparse(substitute(g))
  btxt = deparse(substitute(b))
  itxt = deparse(substitute(id))

  if (!is.vector(x)) {stop(paste(xtxt, "is not a vector"))}
  if (!is.vector(y)) {stop(paste(ytxt, "is not a vector"))}
  if (!is.vector(z)) {stop(paste(ztxt, "is not a vector"))}
  if (!is.vector(r)) {stop(paste(rtxt, "is not a vector"))}
  if (!is.vector(g)) {stop(paste(gtxt, "is not a vector"))}
  if (!is.vector(b)) {stop(paste(btxt, "is not a vector"))}

  if (!is.integer(r)) {stop(paste(rtxt, "must contain integers"))}
  if (!is.integer(g)) {stop(paste(gtxt, "must contain integers"))}
  if (!is.integer(b)) {stop(paste(btxt, "must contain integers"))}

  if (length(x) != length(y)) {stop(paste(xtxt, "is not same length as", ytxt))}
  if (length(x) != length(z)) {stop(paste(xtxt, "is not same length as", ztxt))}
  if (length(r) != length(g)) {stop(paste(rtxt, "is not same length as", gtxt))}
  if (length(r) != length(b)) {stop(paste(rtxt, "is not same length as", btxt))}

  if (length(x) != length(r) & is.null(id)) {stop(paste(xtxt, "is not same length as", rtxt))}
  if (length(x) != length(g) & is.null(id)) {stop(paste(xtxt, "is not same length as", gtxt))}
  if (length(x) != length(b) & is.null(id)) {stop(paste(xtxt, "is not same length as", btxt))}

  if (!is.null(id))
  {
    if (!is.integer(id))
      stop("'id' must contain integers")

    if (max(id) > length(x))
      stop("Index out of bound. 'id' contains wrong ids.")

    if (min(id) <= 0)
      stop("Index out of bound. 'id' contains negative or 0 values.")

    if (length(x) != length(id))
        stop(paste(xtxt, "is not same length as", itxt))
  }
  else
  {
    id = 0
  }

  message("Point cloud viewer must be closed before to run other R code")

  df = data.frame(X = x, Y = y, Z = z, R = r, G = g, B = b)
  viewer(df)
}

#' Display 3D sf geometries with Z values from SpatRaster
#'
#' Display sf geometries (POINT, LINESTRING, POLYGON) in 3D using Z values 
#' extracted from a SpatRaster DEM. Supports styling and coloring by attributes.
#'
#' @param sf_object sf object containing any supported geometry type
#' @param dem SpatRaster object to extract Z values from
#' @param density numeric. Point density for line/polygon outline sampling (points per unit distance). Default is 1.
#' @param color_by character. Name of column to color by. If NULL, uses default colors.
#' @param palette character vector. Color palette to use when color_by is specified. If NULL, uses default palette.
#' @param default_color character. Default color in hex format. Default is "#FFFFFF" (white).
#' @param detach logical. If TRUE, viewer runs in background. Default is FALSE.
#' @param return_data logical. If TRUE, returns data instead of displaying viewer. Used internally for mixed displays.
#' @export
#' @importFrom sf st_coordinates st_cast st_geometry st_geometry_type
#' @importFrom terra extract vect
#' @importFrom grDevices colorRampPalette
#' @md
plot_sf_3d <- function(sf_object, dem, density = 1, color_by = NULL, palette = NULL, 
                       default_color = "#FFFFFF", detach = FALSE, return_data = FALSE)
{
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for this function")
  }
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for this function")
  }
  
  # Determine geometry types
  geom_types <- unique(sf::st_geometry_type(sf_object))
  
  # Setup color palette if color_by is specified
  colors_rgb <- NULL
  if (!is.null(color_by)) {
    if (!color_by %in% names(sf_object)) {
      stop(paste("Column", color_by, "not found in sf object"))
    }
    
    unique_values <- unique(sf_object[[color_by]])
    n_colors <- length(unique_values)
    
    if (is.null(palette)) {
      # Default palette
      if (n_colors <= 8) {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          palette <- RColorBrewer::brewer.pal(max(3, n_colors), "Dark2")[1:n_colors]
        } else {
          palette <- rainbow(n_colors)
        }
      } else {
        palette <- rainbow(n_colors)
      }
    } else if (is.function(palette)) {
      # Handle palette functions (e.g., heat.colors, rainbow)
      palette <- palette(n_colors)
    } else if (length(palette) < n_colors) {
      # Extend palette if needed
      color_func <- grDevices::colorRampPalette(palette)
      palette <- color_func(n_colors)
    }
    
    # Create color mapping
    color_map <- setNames(palette, unique_values)
    colors_rgb <- lapply(sf_object[[color_by]], function(val) {
      col <- color_map[as.character(val)]
      rgb_vals <- col2rgb(col)
      list(r = rgb_vals[1], g = rgb_vals[2], b = rgb_vals[3])
    })
  } else {
    # Use default color
    default_rgb <- col2rgb(default_color)
    colors_rgb <- list(list(r = default_rgb[1], g = default_rgb[2], b = default_rgb[3]))
  }
  
  all_coords <- list()
  all_line_ids <- list()
  all_colors <- list()
  current_line_id <- 1
  
  for (i in 1:nrow(sf_object)) {
    geom <- sf_object[i, ]
    geom_type <- sf::st_geometry_type(geom)
    
    # Get color for this feature
    if (!is.null(color_by)) {
      feature_color <- colors_rgb[[i]]
    } else {
      feature_color <- colors_rgb[[1]]
    }
    
    if (geom_type == "POINT") {
      # Handle POINT geometries
      coords <- sf::st_coordinates(geom)
      # For single points, coords might be a vector, so ensure it's a matrix
      if (is.vector(coords)) {
        coords <- matrix(coords, nrow = 1)
        colnames(coords) <- c("X", "Y", "Z")[1:length(coords)]
      }
      # Ensure we have column names
      if (is.null(colnames(coords))) {
        colnames(coords) <- c("X", "Y", "Z")[1:ncol(coords)]
      }
      
      # Extract only X, Y coordinates
      xy_coords <- coords[, c("X", "Y"), drop = FALSE]
      colnames(xy_coords) <- c("X", "Y")
      
      all_coords[[length(all_coords) + 1]] <- xy_coords
      all_line_ids[[length(all_line_ids) + 1]] <- rep(NA, nrow(xy_coords))  # Points don't have line IDs
      all_colors[[length(all_colors) + 1]] <- feature_color
      
    } else if (geom_type == "LINESTRING") {
      # Handle LINESTRING geometries
      coords <- sf::st_coordinates(geom)
      line_coords <- coords[, c("X", "Y"), drop = FALSE]
      
      # Create line segments with IDs
      line_result <- create_line_segments_with_ids(line_coords, density, line_id_start = current_line_id)
      all_coords[[length(all_coords) + 1]] <- line_result$coords
      all_line_ids[[length(all_line_ids) + 1]] <- line_result$line_ids
      all_colors[[length(all_colors) + 1]] <- feature_color
      current_line_id <- max(line_result$line_ids) + 1
      
    } else if (geom_type == "POLYGON") {
      # Handle POLYGON geometries (convert to linestrings for outline)
      linestrings <- sf::st_cast(geom, "LINESTRING")
      coords <- sf::st_coordinates(linestrings)
      
      if ("L1" %in% colnames(coords)) {
        for (feature_id in unique(coords[, "L1"])) {
          feature_coords <- coords[coords[, "L1"] == feature_id, , drop = FALSE]
          
          if ("L2" %in% colnames(feature_coords)) {
            # Handle multiple rings (exterior + holes)
            for (ring_id in unique(feature_coords[, "L2"])) {
              ring_coords <- feature_coords[feature_coords[, "L2"] == ring_id, c("X", "Y"), drop = FALSE]
              ring_result <- create_line_segments_with_ids(ring_coords, density, line_id_start = current_line_id)
              all_coords[[length(all_coords) + 1]] <- ring_result$coords
              all_line_ids[[length(all_line_ids) + 1]] <- ring_result$line_ids
              all_colors[[length(all_colors) + 1]] <- feature_color
              current_line_id <- max(ring_result$line_ids) + 1
            }
          } else {
            # Single ring
            ring_coords <- feature_coords[, c("X", "Y"), drop = FALSE]
            ring_result <- create_line_segments_with_ids(ring_coords, density, line_id_start = current_line_id)
            all_coords[[length(all_coords) + 1]] <- ring_result$coords
            all_line_ids[[length(all_line_ids) + 1]] <- ring_result$line_ids
            all_colors[[length(all_colors) + 1]] <- feature_color
            current_line_id <- max(ring_result$line_ids) + 1
          }
        }
      } else {
        # Simple case - single linestring
        ring_coords <- coords[, c("X", "Y"), drop = FALSE]
        ring_result <- create_line_segments_with_ids(ring_coords, density, line_id_start = current_line_id)
        all_coords[[length(all_coords) + 1]] <- ring_result$coords
        all_line_ids[[length(all_line_ids) + 1]] <- ring_result$line_ids
        all_colors[[length(all_colors) + 1]] <- feature_color
        current_line_id <- max(ring_result$line_ids) + 1
      }
    }
  }
  
  # Combine all coordinates
  combined_coords <- do.call(rbind, all_coords)
  # Ensure column names are set
  if (is.null(colnames(combined_coords))) {
    colnames(combined_coords) <- c("X", "Y")
  }
  combined_line_ids <- do.call(c, all_line_ids)
  
  # Replicate colors for all points
  combined_colors <- list()
  for (i in 1:length(all_coords)) {
    n_points <- nrow(all_coords[[i]])
    color <- all_colors[[i]]
    combined_colors <- c(combined_colors, rep(list(color), n_points))
  }
  
  # Extract Z values from DEM
  points_sf <- sf::st_as_sf(data.frame(X = combined_coords[, "X"], Y = combined_coords[, "Y"]), 
                           coords = c("X", "Y"), crs = sf::st_crs(sf_object))
  z_values <- terra::extract(dem, terra::vect(points_sf))[, 2]
  
  # Handle NA values
  valid_idx <- !is.na(z_values)
  if (sum(valid_idx) == 0) {
    stop("No valid Z values could be extracted from the DEM")
  }
  
  x <- combined_coords[valid_idx, "X"]
  y <- combined_coords[valid_idx, "Y"]
  z <- z_values[valid_idx]
  line_ids_valid <- combined_line_ids[valid_idx]
  colors_valid <- combined_colors[valid_idx]
  
  # Extract RGB values
  r <- sapply(colors_valid, function(col) col$r)
  g <- sapply(colors_valid, function(col) col$g)
  b <- sapply(colors_valid, function(col) col$b)
  
  df <- data.frame(X = x, Y = y, Z = z, R = r, G = g, B = b, LINE_ID = line_ids_valid)
  
  if (return_data) {
    # Return data for mixed display instead of showing viewer
    return(list(
      coords = df[, c("X", "Y", "Z")],
      colors = df[, c("R", "G", "B")],
      line_ids = df$LINE_ID
    ))
  } else {
    message("SF geometry viewer must be closed before to run other R code")
    viewer(df, detach, "")
  }
}

#' Colorize LAS point cloud by attribute using custom palette
#'
#' This function assigns RGB colors to a LAS point cloud based on a specified attribute
#' and color palette function. It normalizes the attribute values and maps them to colors.
#'
#' @param las LAS object to colorize
#' @param attribute character. Name of the attribute column to color by
#' @param palette_func function. Color palette function (e.g., rainbow, heat.colors, terrain.colors)
#' @param n integer. Number of colors in the palette. Default is 100.
#' @return LAS object with R, G, B columns added/updated
#' @export
#' @importFrom grDevices col2rgb
#' @examples
#' \dontrun{
#' # Colorize by Z values using rainbow palette
#' las_colored <- colorizeLAS(las, "Z", rainbow)
#' 
#' # Colorize by intensity using heat colors
#' las_colored <- colorizeLAS(las, "Intensity", heat.colors, n = 50)
#' 
#' # Use with RColorBrewer
#' library(RColorBrewer)
#' las_colored <- colorizeLAS(las, "Classification", 
#'                           function(n) brewer.pal(min(n, 9), "Set1"))
#' }
#' @md
colorizeLAS <- function(las, attribute, palette_func, n = 100) {
  # Ensure the attribute exists in the LAS data
  if (!attribute %in% names(las@data)) {
    stop(paste("Attribute", attribute, "not found in LAS data."))
  }
  
  attr_values <- las@data[[attribute]]
  
  # Handle NA values by replacing them with the minimum value
  if (any(is.na(attr_values))) {
    min_val <- min(attr_values, na.rm = TRUE)
    attr_values[is.na(attr_values)] <- min_val
  }
  
  # Normalize attribute values to 0–1
  rng <- range(attr_values, na.rm = TRUE)
  if (rng[1] == rng[2]) {
    norm_attr <- rep(0.5, length(attr_values))
  } else {
    norm_attr <- (attr_values - rng[1]) / (rng[2] - rng[1])
  }
  
  # Create the color palette
  colors <- palette_func(n)
  
  # Map normalized values to color indices
  color_indices <- as.integer(norm_attr * (n - 1)) + 1
  mapped_colors <- colors[color_indices]
  
  # Convert colors to RGB
  rgb_values <- grDevices::col2rgb(mapped_colors)
  
  # Assign RGB to the LAS object
  las@data$R <- rgb_values[1, ]
  las@data$G <- rgb_values[2, ]
  las@data$B <- rgb_values[3, ]
  
  return(las)
}

#' Create wireframe visualization of DEM using triangulation
#'
#' This function converts a SpatRaster DEM to a wireframe representation using 
#' Delaunay triangulation, similar to lidR's tin() algorithm. The wireframe can
#' be displayed in the 3D viewer as line segments connecting triangulated points.
#'
#' @param dem SpatRaster object (DEM)
#' @param sample_factor numeric. Factor to reduce DEM resolution for triangulation (default: 4)
#' @param line_color character. Color for wireframe lines (default: "#808080" gray)
#' @return sf object with LINESTRING geometries representing the wireframe
#' @export
#' @importFrom terra as.points crds values ext
#' @importFrom sf st_sf st_sfc st_linestring st_crs
#' @examples
#' \dontrun{
#' # Create wireframe from DEM
#' dem_wireframe <- create_dem_wireframe(dem, sample_factor = 2)
#' 
#' # Display wireframe only
#' plot_sf_3d(dem_wireframe, dem, line_width = 1, default_color = "#808080")
#' 
#' # Combine with point cloud and polygons
#' view(las, polygons = metrics, dem_wireframe = dem_wireframe, dem = dem)
#' }
#' @md
create_dem_wireframe <- function(dem, sample_factor = 4, line_color = "#808080") {
  # Sample the DEM to reduce computational load
  if (sample_factor > 1) {
    # Get every nth pixel based on sample_factor
    dem_matrix <- as.matrix(dem, wide = TRUE)
    nrows <- nrow(dem_matrix)
    ncols <- ncol(dem_matrix)
    
    row_indices <- seq(1, nrows, by = sample_factor)
    col_indices <- seq(1, ncols, by = sample_factor)
    
    # Create coordinate matrices
    ext_dem <- terra::ext(dem)
    x_coords <- seq(ext_dem[1], ext_dem[2], length.out = ncols)[col_indices]
    y_coords <- seq(ext_dem[4], ext_dem[3], length.out = nrows)[row_indices]  # Note: y is flipped
    
    # Sample the matrix
    sampled_matrix <- dem_matrix[row_indices, col_indices]
    
    # Create data frame of points
    points_df <- expand.grid(X = x_coords, Y = rev(y_coords))  # Rev to match raster orientation
    points_df$Z <- as.vector(sampled_matrix)
    
    # Remove NA values
    points_df <- points_df[!is.na(points_df$Z), ]
  } else {
    # Use all points
    points_terra <- terra::as.points(dem, values = TRUE, na.rm = TRUE)
    coords <- terra::crds(points_terra)
    values <- terra::values(points_terra)[, 1]
    points_df <- data.frame(X = coords[, 1], Y = coords[, 2], Z = values)
  }
  
  if (nrow(points_df) < 3) {
    stop("Not enough valid DEM points for triangulation")
  }
  
  # Perform Delaunay triangulation using geometry package if available
  tryCatch({
    # Use geometry package for triangulation
    if (!requireNamespace("geometry", quietly = TRUE)) {
      stop("Package 'geometry' is required for DEM wireframe generation. Install with: install.packages('geometry')")
    }
    
    # Create triangulation
    tri <- geometry::delaunayn(points_df[, 1:2])
    
    # Create line segments from triangles
    line_geometries <- list()
    for (i in 1:nrow(tri)) {
      # Get triangle vertices
      p1 <- tri[i, 1]
      p2 <- tri[i, 2] 
      p3 <- tri[i, 3]
      
      # Create three line segments for each triangle
      # Line 1: p1 -> p2
      line1 <- matrix(c(points_df$X[p1], points_df$Y[p1],
                       points_df$X[p2], points_df$Y[p2]), ncol = 2, byrow = TRUE)
      
      # Line 2: p2 -> p3
      line2 <- matrix(c(points_df$X[p2], points_df$Y[p2],
                       points_df$X[p3], points_df$Y[p3]), ncol = 2, byrow = TRUE)
      
      # Line 3: p3 -> p1
      line3 <- matrix(c(points_df$X[p3], points_df$Y[p3],
                       points_df$X[p1], points_df$Y[p1]), ncol = 2, byrow = TRUE)
      
      line_geometries[[length(line_geometries) + 1]] <- sf::st_linestring(line1)
      line_geometries[[length(line_geometries) + 1]] <- sf::st_linestring(line2)
      line_geometries[[length(line_geometries) + 1]] <- sf::st_linestring(line3)
    }
    
    # Create sf object
    wireframe_sf <- sf::st_sf(
      id = 1:length(line_geometries),
      type = rep("wireframe", length(line_geometries)),
      color = rep(line_color, length(line_geometries)),
      geometry = sf::st_sfc(line_geometries)
    )
    
    # Set CRS to match DEM
    sf::st_crs(wireframe_sf) <- sf::st_crs(dem)
    
    return(wireframe_sf)
    
  }, error = function(e) {
    stop(paste("Triangulation failed:", e$message, 
               "\nTry reducing sample_factor or check DEM validity."))
  })
}
