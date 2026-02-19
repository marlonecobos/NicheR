#' Apply Sampling Bias to a Prediction Surface
#'
#' Modifies a prediction surface by applying a spatial bias layer.
#' The bias layer is dynamically cropped and masked to match the prediction
#' layer, adjusted for effect direction, and combined multiplicatively.
#'
#' @param prepared_bias A \code{nicheR_bias_surface} object (output from
#'   \code{prepare_bias}) or a \code{terra::SpatRaster}.
#' @param pred A \code{terra::SpatRaster} representing the prediction surface
#'   (e.g., habitat suitability).
#' @param effect_direction Numeric vector. \code{1} for a direct effect (\code{pred * bias}),
#'   or \code{-1} for an indirect/inverse effect (\code{pred * (1 - bias)}).
#'   Can be a single value applied to all layers, or a vector (e.g., \code{c(1, -1)})
#'   to apply different effects to different layers. Default is \code{1}.
#' @param layers Character vector. Names of the layers in \code{pred} to which
#'   the bias should be applied. Defaults to \code{c("suitability", "suitability_trunc")}.
#'
#' @return A \code{terra::SpatRaster} containing only the bias-adjusted prediction layers.
#' @export
apply_bias <- function(prepared_bias,
                       pred,
                       effect_direction = 1,
                       layers = c("suitability", "suitability_trunc")) {

  # 1. Input Validation
  if (!inherits(pred, "SpatRaster")) {
    stop("'pred' must be a terra::SpatRaster.")
  }

  # Use all() to allow vectors of length > 1
  if (!all(effect_direction %in% c(1, -1))) {
    stop("'effect_direction' must only contain 1 (direct) and/or -1 (indirect).")
  }

  # 2. Subset prediction to target layers dynamically
  available_layers <- names(pred)
  missing_layers <- setdiff(layers, available_layers)

  if (length(missing_layers) > 0) {
    warning(sprintf("The following layers were not found in 'pred' and will be skipped: %s",
                    paste(missing_layers, collapse = ", ")))
  }

  valid_layers <- intersect(layers, available_layers)

  if (length(valid_layers) == 0) {
    stop("None of the target 'layers' were found in 'pred'.")
  }

  pred_sub <- pred[[valid_layers]]

  # 3. Synchronize effect_direction length with valid layers
  if (length(effect_direction) != length(valid_layers)) {
    if (length(effect_direction) > 1 && length(effect_direction) < length(valid_layers)) {
      warning("Length of 'effect_direction' is > 1 but doesn't match the number of valid layers. It will be recycled.")
    }
    # Safely recycle or truncate the vector to match the number of layers
    effect_direction <- rep_len(effect_direction, length(valid_layers))
  }

  # 4. Extract the bias raster
  if (inherits(prepared_bias, "nicheR_bias_surface")) {
    if (!is.null(prepared_bias$pooled_bias)) {
      bias_ras <- prepared_bias$pooled_bias
    } else {
      stop("'prepared_bias' is missing the pooled layer. Did you use out_bias = 'standardized'?")
    }
  } else if (inherits(prepared_bias, "SpatRaster")) {
    bias_ras <- prepared_bias
  } else {
    stop("'prepared_bias' must be a 'nicheR_bias_surface' list or a 'terra::SpatRaster'.")
  }

  # 5. Crop, Mask, Apply Effect, and Combine
  result_list <- list()

  for (i in seq_along(valid_layers)) {
    layer_name <- valid_layers[i]
    current_pred <- pred_sub[[layer_name]]

    # Safely select the corresponding bias layer
    bias_idx <- min(i, terra::nlyr(bias_ras))
    current_bias <- bias_ras[[bias_idx]]

    # Crop and Mask
    b_cropped <- terra::crop(current_bias, current_pred, mask = TRUE)

    # Apply Effect Direction for THIS specific layer in the loop
    if (effect_direction[i] == -1) {
      b_cropped <- 1 - b_cropped
    }

    # Apply the bias
    applied <- current_pred * b_cropped
    names(applied) <- paste0(layer_name, "_biased")

    result_list[[i]] <- applied
  }

  # Re-combine into a single SpatRaster
  applied_bias <- do.call(c, result_list)

  return(applied_bias)
}
