#' Combine one or more bias surfaces into a pooled bias raster
#'
#' This helper function combines one or more sampling bias layers into a single
#' pooled bias surface. Each bias layer is:
#' \itemize{
#'   \item aligned to a common template grid,
#'   \item scaled to \eqn{[0, 1]},
#'   \item optionally reversed via \eqn{1 - x},
#'   \item multiplied with other bias layers to form a pooled bias surface.
#' }
#'
#' The resulting pooled bias raster can then be used elsewhere in NicheR to
#' modulate sampling or suitability. Optionally, the function can also return
#' the standardized bias layers as a raster stack.
#'
#' @param bias_surface One or more bias layers which must be:
#'   \itemize{
#'     \item a \code{terra::SpatRaster} (single or multi-layer), or
#'     \item a \code{list} of \code{terra::SpatRaster} objects.
#'   }
#'   Any multi-layer \code{SpatRaster} is internally split into single-layer
#'   rasters before standardization.
#' @param bias_dir Numeric values of \code{1} or \code{-1} controlling
#'   directionality. A length 1 value is recycled across all layers. A value
#'   of \code{-1} applies \code{1 - layer} after scaling to \eqn{[0, 1]}.
#' @param suitable_env Optional \code{terra::SpatRaster} or related object used
#'   as a mask and template. Typically this is a binary or continuous
#'   suitability raster or a \code{suitable_env} object. Masking is applied
#'   after standardization so that bias is only defined where the species is
#'   suitable. Required when \code{out.bias} is \code{"biased"} or \code{"both"}.
#' @param out.bias Character, one of \code{"biased"}, \code{"standardized"},
#'   or \code{"both"}:
#'   \itemize{
#'     \item \code{"standardized"}: returns only the standardized,
#'       direction-adjusted bias stack.
#'     \item \code{"biased"}: returns only the pooled bias surface, masked by
#'       \code{suitable_env}. Requires \code{suitable_env}.
#'     \item \code{"both"}: returns both the pooled bias surface and the
#'       standardized stack. Requires \code{suitable_env}.
#'   }
#' @param verbose Logical. If \code{TRUE}, prints basic progress messages.
#'
#' @return A list of class \code{"nicheR_bias_surface"} containing some or all
#'   of:
#'   \item{pooled_bias_sp}{Pooled bias raster (scaled 0–1; \code{SpatRaster}),
#'         masked by \code{suitable_env} when supplied. Present only when
#'         \code{out.bias} is \code{"biased"} or \code{"both"}.}
#'   \item{directional_bias_stack}{Stack of standardized, direction-corrected
#'         bias layers (\code{SpatRaster}). Present only when
#'         \code{out.bias} is \code{"standardized"} or \code{"both"}.}
#'   \item{combination_formula}{String describing how layers were combined
#'         (e.g., \code{"bias_1 * (1-bias_2)"}).}
#'
#' @importFrom terra rast resample values app nlyr as.list ext res mask
#' @export
set_bias_surface <- function(bias_surface,
                             bias_dir = 1,
                             suitable_env = NULL,
                             out.bias = c("biased", "standardized", "both"),
                             verbose = TRUE) {

  gc()
  out.bias <- match.arg(out.bias)

  if (isTRUE(verbose)) {
    message("Starting set_bias_surface()...")
  }

  # ---- 0. Basic checks ------------------------------------------------------

  if (missing(bias_surface) || is.null(bias_surface)) {
    stop("'bias_surface' is required and must be a SpatRaster or a list of SpatRasters.")
  }

  # Enforce raster inputs -----------------------------------------------------
  if (inherits(bias_surface, "SpatRaster")) {
    # split multi-layer rasters into single-layer list
    if (isTRUE(verbose)) {
      message("bias_surface provided as SpatRaster. Splitting into single layers...")
    }
    bias_list <- terra::as.list(bias_surface)
  } else if (is.list(bias_surface) &&
             length(bias_surface) > 0 &&
             all(vapply(bias_surface, inherits, logical(1), "SpatRaster"))) {
    if (isTRUE(verbose)) {
      message("bias_surface provided as list of SpatRasters. Flattening layers...")
    }
    # flatten list of SpatRasters into single-layer rasters
    bias_list <- unlist(lapply(bias_surface, terra::as.list), recursive = FALSE)
  } else {
    stop(
      "'bias_surface' must be either:\n",
      "  * a terra::SpatRaster (single or multi-layer), or\n",
      "  * a list of terra::SpatRaster objects."
    )
  }

  if (length(bias_list) == 0) {
    stop("No bias layers provided.")
  }

  # Suitable env / mask checks -----------------------------------------------

  if (out.bias %in% c("biased", "both") && is.null(suitable_env)) {
    stop(
      "out.bias = '", out.bias, "' requires 'suitable_env' (SpatRaster or suitable_env object) ",
      "to be provided so it can be used as a mask."
    )
  }

  mask_ras <- NULL
  if (!is.null(suitable_env)) {

    # Helper: extract first SpatRaster found inside a variety of objects
    get_first_raster <- function(obj) {
      # Direct SpatRaster
      if (inherits(obj, "SpatRaster")) {
        return(obj[[1]])
      }

      # "suitable_env" object from get_suitable_env(out.suit = "both")
      if (inherits(obj, "suitable_env") && "suitable_env_sp" %in% names(obj)) {
        return(get_first_raster(obj$suitable_env_sp)) # tiny recursive call
      }

      # List: look for first SpatRaster in elements
      if (is.list(obj)) {
        # If named, try "suitable" first
        if (!is.null(names(obj)) && "suitable" %in% names(obj)) {
          if (inherits(obj[["suitable"]], "SpatRaster")) {
            return(obj[["suitable"]][[1]])
          }
        }
        # Otherwise, scan through elements
        for (el in obj) {
          if (inherits(el, "SpatRaster")) {
            return(el[[1]])
          }
        }
      }

      # Nothing usable found
      NULL
    }

    mask_ras <- get_first_raster(suitable_env)

    if (is.null(mask_ras)) {
      stop(
        "'suitable_env' must contain at least one terra::SpatRaster to be used as a mask.\n",
        "Examples:\n",
        "  - a SpatRaster directly (e.g. suitable_ras),\n",
        "  - a 'suitable_env' object from get_suitable_env(out.suit = 'both'),\n",
        "  - a list with a 'suitable' SpatRaster element ",
        "    (e.g. suitable_env$suitable_env_sp$suitable)."
      )
    }

    if (isTRUE(verbose)) {
      message("Using 'suitable_env' raster as template and mask.")
    }
  }

  # Template raster: use mask if available, else the first bias layer --------

  if (!is.null(mask_ras)) {
    template_raster <- mask_ras
  } else {
    first_bias <- bias_list[[1]]

    if (!inherits(first_bias, "SpatRaster")) {
      stop("First bias layer is not a SpatRaster. All bias layers must be SpatRasters.")
    }

    template_raster <- first_bias

    if (isTRUE(verbose)) {
      message("Using first bias raster as template.")
    }
  }

  # Validate bias_dir ---------------------------------------------------------

  if (length(bias_dir) == 1) {
    bias_dir <- rep(bias_dir, length(bias_list))
  } else if (length(bias_dir) != length(bias_list)) {
    stop("The length of 'bias_dir' (", length(bias_dir),
         ") must match the number of bias layers (", length(bias_list),
         ") or be 1.")
  }
  if (!all(bias_dir %in% c(1, -1))) {
    stop("'bias_dir' must only contain values of 1 or -1.")
  }

  # ---- 1. Process, align, and standardize each bias layer -------------------

  if (isTRUE(verbose)) {
    message("Processing and standardizing ", length(bias_list),
            " bias layer(s) on the template grid...")
  }

  directional_bias_list <- vector("list", length(bias_list))
  dir_message_parts     <- character(length(bias_list))

  for (i in seq_along(bias_list)) {

    bias_raster_raw <- bias_list[[i]]
    this_dir        <- bias_dir[i]

    # Try to derive a readable name
    layer_name <- names(bias_raster_raw)
    if (is.null(layer_name) || length(layer_name) == 0 || layer_name[1] == "") {
      layer_name <- paste0("bias_", i)
    } else {
      layer_name <- layer_name[1]
    }

    # Align resolution + extent to template
    if (!identical(terra::res(bias_raster_raw), terra::res(template_raster)) ||
        !identical(terra::ext(bias_raster_raw), terra::ext(template_raster))) {

      if (isTRUE(verbose)) {
        message("  - Aligning bias layer ", i, " (", layer_name, ") to template grid (resample).")
      }
      bias_raster_aligned <- terra::resample(bias_raster_raw, template_raster, method = "near")
    } else {
      bias_raster_aligned <- bias_raster_raw
    }

    # Normalize to [0, 1] (preserving NA)
    vals      <- terra::values(bias_raster_aligned)
    min_val   <- min(vals, na.rm = TRUE)
    max_val   <- max(vals, na.rm = TRUE)
    range_val <- max_val - min_val

    if (range_val == 0) {
      warning("Bias layer ", i,
              " ('", layer_name, "') contains a single unique value. ",
              "Setting all non-NA values to 1.",
              call. = FALSE)
      scaled_layer <- bias_raster_aligned
      scaled_layer[!is.na(scaled_layer)] <- 1
    } else {
      scaled_layer <- (bias_raster_aligned - min_val) / range_val
    }

    # Apply directionality
    if (this_dir == -1) {
      directional_layer    <- 1 - scaled_layer
      dir_message_parts[i] <- paste0("(1-", layer_name, ")")
    } else {
      directional_layer    <- scaled_layer
      dir_message_parts[i] <- layer_name
    }

    names(directional_layer)   <- layer_name
    directional_bias_list[[i]] <- directional_layer
  }

  if (isTRUE(verbose)) {
    message("Finished processing bias layers.")
  }

  # Stack directional layers ---------------------------------------------------

  directional_bias_stack <- if (length(directional_bias_list) == 1) {
    directional_bias_list[[1]]
  } else {
    do.call(c, directional_bias_list)
  }

  # ---- 2. Pool bias layers into a single surface ----------------------------

  pooled_bias_sp <- NULL

  if (out.bias %in% c("biased", "both")) {

    if (isTRUE(verbose)) {
      message("Pooling ", terra::nlyr(directional_bias_stack),
              " bias layer(s) into a single surface...")
    }

    if (terra::nlyr(directional_bias_stack) > 1) {
      pooled_bias_sp <- terra::app(
        directional_bias_stack,
        fun = function(x) {
          if (all(is.na(x))) {
            NA_real_           # keep fully-missing cells as NA
          } else {
            prod(x, na.rm = TRUE)
          }
        }
      )
    } else {
      pooled_bias_sp <- directional_bias_stack
    }

    names(pooled_bias_sp) <- "pooled_bias"

    # Mask to suitable area if provided
    if (!is.null(mask_ras)) {
      pooled_bias_sp <- terra::mask(pooled_bias_sp, mask_ras)
    }

    if (isTRUE(verbose)) {
      message("Finished pooling bias layers.")
    }
  }

  # Also mask the standardized stack if a mask was provided -------------------

  if (!is.null(mask_ras)) {
    directional_bias_stack <- terra::mask(directional_bias_stack, mask_ras)
  }

  combination_formula <- paste(dir_message_parts, collapse = " * ")

  # ---- 3. Build result object -----------------------------------------------

  res <- list(
    pooled_bias_sp         = NULL,
    directional_bias_stack = NULL,
    combination_formula    = combination_formula
  )

  if (out.bias %in% c("biased", "both")) {
    res$pooled_bias_sp <- pooled_bias_sp
  }

  if (out.bias %in% c("standardized", "both")) {
    res$directional_bias_stack <- directional_bias_stack
  }

  class(res) <- "nicheR_bias_surface"

  gc()

  if (isTRUE(verbose)) {
    message("set_bias_surface() completed successfully.")
  }

  return(res)
}

#' @export
print.nicheR_bias_surface <- function(x, ...) {

  cat("NicheR bias surface object\n")

  # Combination formula (how layers were combined)
  if (!is.null(x$combination_formula)) {
    cat("  Combination:\n")
    cat("    ", x$combination_formula, "\n\n")
  }

  # pooled_bias_sp summary
  if (!is.null(x$pooled_bias_sp) && inherits(x$pooled_bias_sp, "SpatRaster")) {
    r <- x$pooled_bias_sp
    cat("  pooled_bias_sp (terra::SpatRaster)\n")
    cat("    layers :", terra::nlyr(r), "\n")
    cat("    names  :", paste(names(r), collapse = ", "), "\n\n")
  } else {
    cat("  pooled_bias_sp : NULL\n\n")
  }

  # directional_bias_stack summary (if present)
  if (!is.null(x$directional_bias_stack) &&
      inherits(x$directional_bias_stack, "SpatRaster")) {

    r <- x$directional_bias_stack
    cat("  directional_bias_stack (terra::SpatRaster)\n")
    cat("    layers :", terra::nlyr(r), "\n")
    cat("    names  :", paste(names(r), collapse = ", "), "\n")
  }

  invisible(x)
}
