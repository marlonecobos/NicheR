#' This function applies one or more bias surfaces (sampling probability layers) to an
#' ecological suitability layer to determine the final realized suitability.
#'
#' It handles multiple bias inputs, aligns them to a common template, standardizes them
#' individually, applies user-defined directionality (e.g., $1 - \text{layer}$), and
#' multiplies them together before applying the pooled bias to the ecological suitability.
#'
#' @param bias_surface A single bias layer or a collection of bias layers. This can be:
#'   \itemize{
#'     \item A \code{terra::SpatRaster} (single layer or stack).
#'     \item A \code{terra::SpatVector} (which will be rasterized).
#'     \item A path to a single raster/vector file.
#'     \item A \code{list} of any of the above (paths, SpatRasters, SpatVectors).
#'   }
#'   Values should ideally be scaled from 0 to 1.
#' @param ecological_suitability Optional. A \code{terra::SpatRaster} (values 0/1) representing the
#'   ecological suitability layer. If \code{NULL} (default), the function returns only the
#'   \strong{Pooled Bias Surface}.
#' @param bias_direction Numeric vector. Defines the directionality for each bias layer.
#'   \itemize{
#'     \item \code{1} (or \code{-1}) applies to all layers if a single value is provided.
#'     \item A vector (e.g., \code{c(1, -1, 1)}) must match the number of \code{bias_surface} layers.
#'   }
#'   A value of \code{1} indicates higher raw layer values mean higher bias. A value of \code{-1}
#'   indicates higher raw layer values mean \strong{lower} bias (applying the $1 - \text{layer}$ reversal).
#' @param raster_resolution Numeric. The resolution to use when rasterizing vector inputs
#'   or creating a template if no \code{ecological_suitability} raster is provided. Defaults to 0.5.
#' @param output Character string, controlling the returned layers.
#'   \itemize{
#'     \item \code{"default"} (or \code{NULL}): Returns the \code{final_suitability} (or Pooled Bias) and the \code{pooled_bias} layer.
#'     \item \code{"all"}: Returns all layers: \code{final_suitability} (or Pooled Bias), \code{pooled_bias}, and the stack of \code{directional_bias_stack}.
#'   }
#'
#' @return A list of class \code{biased_suitable_env} containing:
#'   \item{pooled_bias_sp}{The pooled (multiplied) and standardized bias layer (\code{SpatRaster}).}
#'   \item{pooled_bias_df}{Data frame of coordinates/values for non-zero cells of \code{pooled_bias_sp}.}
#'   \item{final_suitability_sp}{*Only when \code{ecological_suitability} is provided.* The final \code{terra::SpatRaster} result (Niche * Pooled Bias).}
#'   \item{final_suitability_df}{*Only when \code{ecological_suitability} is provided.* Data frame of coordinates/values for non-zero cells of \code{final_suitability_sp}.}
#'   \item{directional_bias_stack}{*Only with \code{output="all"}.* \code{SpatRaster} stack of all individual, standardized, and direction-corrected bias layers.}
#'   \item{directional_bias_df}{*Only with \code{output="all"}.* Data frame of coordinates/values for non-zero cells of \code{directional_bias_stack}.}
#'   \item{combination_formula}{*Internal*: A string detailing how the layers were combined.}
#'   \item{is_niche_applied}{*Internal*: Boolean indicating if \code{ecological_suitability} was used.}
#'
#' @family niche functions
#' @importFrom terra rast resample values as.data.frame vect rasterize app nlyr as.list
#' @export
apply_bias_surface <- function(bias_surface, ecological_suitability = NULL, bias_direction = 1, raster_resolution = 0.5, output = "default") {

  # --- 1. Load, Validate, and Determine Template ---

  if (inherits(bias_surface, "list")) {
    bias_list <- bias_surface
  } else if (inherits(bias_surface, "SpatRaster")) {
    # If a SpatRaster stack (like from c(r1, r2)) is provided, split it into individual layers
    if (terra::nlyr(bias_surface) > 1) {
      bias_list <- terra::as.list(bias_surface)
    } else {
      bias_list <- list(bias_surface)
    }
  } else if (inherits(bias_surface, c("SpatVector", "character"))) {
    bias_list <- list(bias_surface)
  } else {
    stop("'bias_surface' must be a SpatRaster, SpatVector, path, or a list of these objects.")
  }

  if (length(bias_list) == 0) {
    stop("No bias layers provided.")
  }

  # Template raster selection
  is_niche_applied <- !is.null(ecological_suitability)

  if (is_niche_applied) {
    if (!inherits(ecological_suitability, "SpatRaster")) {
      stop("'ecological_suitability' must be a SpatRaster if provided.")
    }
    template_raster <- ecological_suitability
  } else {
    # If no suitability, use the first rasterized bias layer as the template
    # Find the first SpatRaster or SpatVector to set template extent
    first_bias_obj <- bias_list[[1]]
    if (inherits(first_bias_obj, "character")) {
      first_bias_obj <- tryCatch(
        terra::rast(first_bias_obj),
        error = function(e) tryCatch(terra::vect(first_bias_obj), error = function(e2) stop(paste("Could not load first bias layer:", e2$message)))
      )
    }

    if (inherits(first_bias_obj, "SpatRaster")) {
      template_raster <- first_bias_obj
    } else if (inherits(first_bias_obj, "SpatVector")) {
      template_raster <- terra::rast(first_bias_obj, resolution = raster_resolution)
      message(paste0("No template raster provided; using first vector extent and custom resolution (", raster_resolution, ") for template creation."))
    } else {
      stop("Could not determine template from first bias layer.")
    }
  }

  # Validate bias_direction length
  if (length(bias_direction) == 1) {
    bias_direction <- rep(bias_direction, length(bias_list))
  } else if (length(bias_direction) != length(bias_list)) {
    stop("The length of 'bias_direction' (", length(bias_direction), ") must match the number of bias layers (", length(bias_list), ") or be 1.")
  }
  if (!all(bias_direction %in% c(1, -1))) {
    stop("'bias_direction' must only contain values of 1 or -1.")
  }

  # --- 2. Process, Align, and Standardize Layers ---
  directional_bias_list <- list()
  dir_message_parts <- c()

  for (i in seq_along(bias_list)) {
    bias_obj <- bias_list[[i]]
    current_direction <- bias_direction[i]
    layer_name <- paste0(names(bias_list[[i]]))

    # Load/Convert to SpatRaster
    if (inherits(bias_obj, "character")) {
      bias_obj <- tryCatch(
        terra::rast(bias_obj),
        error = function(e) tryCatch(terra::vect(bias_obj), error = function(e2) stop(paste("Error loading layer", i, ":", e2$message)))
      )
    }

    if (inherits(bias_obj, "SpatVector")) {
      # Rasterize vector input
      field_names <- names(bias_obj)
      field_to_use <- if (length(field_names) == 0) 1 else field_names[1]
      if (length(field_names) == 0) {
        warning(paste0("Layer ", i, " (vector) has no attribute fields. Rasterizing to binary presence (value=1)."), call. = FALSE)
      }
      bias_raster_raw <- terra::rasterize(bias_obj, template_raster, fun = "max", background = 0, field = field_to_use)
    } else if (inherits(bias_obj, "SpatRaster")) {
      bias_raster_raw <- bias_obj
    } else {
      stop(paste("Layer", i, "must be a SpatRaster, SpatVector, or a valid path."))
    }

    # Align to template
    if (!identical(terra::res(bias_raster_raw), terra::res(template_raster)) || !identical(terra::ext(bias_raster_raw), terra::ext(template_raster))) {
      bias_raster_aligned <- terra::resample(bias_raster_raw, template_raster, method = "near")
    } else {
      bias_raster_aligned <- bias_raster_raw
    }

    # Normalize to [0, 1]
    min_val <- min(terra::values(bias_raster_aligned), na.rm = TRUE)
    max_val <- max(terra::values(bias_raster_aligned), na.rm = TRUE)
    range_val <- max_val - min_val

    if (range_val == 0) {
      warning(paste0("Layer ", i, " contains only a single value. Setting all non-NA values to 1."), call. = FALSE)
      scaled_layer <- bias_raster_aligned
      scaled_layer[!is.na(scaled_layer)] <- 1
    } else {
      # Apply min-max scaling: (X - min) / range
      scaled_layer <- (bias_raster_aligned - min_val) / range_val
    }

    # Apply directionality
    if (current_direction == -1) {
      directional_layer <- 1 - scaled_layer
      dir_message_parts[i] <- paste0("(1-", layer_name, ")")
    } else {
      directional_layer <- scaled_layer
      dir_message_parts[i] <- layer_name
    }
    names(directional_layer) <- layer_name
    directional_bias_list[[i]] <- directional_layer
  }

  # Stack the directional layers
  directional_bias_stack <- if (length(directional_bias_list) == 1) directional_bias_list[[1]] else do.call(c, directional_bias_list)

  # --- 3. Pool Bias Layers and Generate Message ---

  # Multiply all layers in the stack (pool them)
  # Use terra::app with multiplication for robustness
  if (terra::nlyr(directional_bias_stack) > 1) {
    pooled_bias_sp <- terra::app(directional_bias_stack, fun = "prod")
  } else {
    pooled_bias_sp <- directional_bias_stack
  }
  names(pooled_bias_sp) <- "pooled_bias"

  # Calculate pooled_bias_df (always returned)
  pooled_bias_df <- terra::as.data.frame(pooled_bias_sp, xy = TRUE, na.rm = TRUE)
  pooled_bias_df <- pooled_bias_df[pooled_bias_df[["pooled_bias"]] > 0, ]


  # Generate explanatory message (for storage and printing)
    combination_formula <- paste0("Bias layers were combined: ", paste(dir_message_parts, collapse = " * "))

  # --- 4. Calculate Final Suitability (Conditional) ---
  final_suitability_sp <- NULL
  final_suitability_df <- NULL

  if (is_niche_applied) {
    # If niche is applied, calculate final suitability and generate its data frame
    final_suitability_sp <- ecological_suitability * pooled_bias_sp
    names(final_suitability_sp) <- "final_suitability"

    final_suitability_df <- terra::as.data.frame(final_suitability_sp, xy = TRUE, na.rm = TRUE)
    final_suitability_df <- final_suitability_df[final_suitability_df[["final_suitability"]] > 0, ]

    if (nrow(final_suitability_df) == 0) {
      warning("No area remains after processing and filtering.", call. = FALSE)
    }
  }


  # --- 5. Prepare Output List (Conditional return keys) ---
  res <- list(
    pooled_bias_sp = pooled_bias_sp,
    pooled_bias_df = pooled_bias_df,
    combination_formula = combination_formula, # Store formula
    is_niche_applied = is_niche_applied       # Store status
  )

  # Only add final suitability layers if niche was applied, based on user request
  if (is_niche_applied) {
    res$final_suitability_sp <- final_suitability_sp
    res$final_suitability_df <- final_suitability_df
  }

  # Add individual layers if output="all"
  if (tolower(output) == "all") {
    directional_bias_df <- terra::as.data.frame(directional_bias_stack, xy = TRUE, na.rm = TRUE)
    # Filter out rows where ALL bias values are zero for efficiency
    bias_value_cols <- names(directional_bias_df)[!names(directional_bias_df) %in% c("x", "y")]
    directional_bias_df <- directional_bias_df[rowSums(directional_bias_df[bias_value_cols], na.rm = TRUE) > 0, ]

    res$directional_bias_stack <- directional_bias_stack
    res$directional_bias_df <- directional_bias_df
  }

  class(res) <- "biased_suitable_env"
  return(res)
}

#' @export
print.biased_suitable_env <- function(x, ...) {
  # 1. Print the combination formula
  if (!is.null(x$combination_formula)) {
    cat(x$combination_formula, "\n")
  }

  # 2. Print suitability calculation status
  if (isTRUE(x$is_niche_applied)) {
    cat("Ecological suitability provided: (Niche * Pooled Bias).\n")
    # 3. Print layer summary for final suitability
    if (!is.null(x$final_suitability_df)) {
    }
  } else {
    # Custom message when ecological suitability is NULL, as requested by the user
    cat("Ecological suitability is NULL.\nReturning the Pooled Bias layer only.\n")
  }

  invisible(x)
}
