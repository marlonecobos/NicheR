#' Apply sampling bias to suitability surfaces
#'
#' @description
#' Applies a prepared composite sampling bias surface to a suitability raster
#' by multiplication. The bias surface is aligned to the suitability grid when
#' needed and the result is cropped and masked to the suitability domain. The
#' output is a product of suitability and bias and is therefore no longer
#' interpretable as a probability.
#'
#' @param prepared_bias A single-layer \code{SpatRaster} composite bias
#'   surface, or the list output from \code{\link{prepare_bias}} containing
#'   a \code{composite_surface} element.
#' @param prediction A \code{SpatRaster} containing one or more suitability
#'   layers with values in \code{[0, 1]}.
#' @param prediction_layer Character. Name of the layer to extract from
#'   \code{prediction} when it contains multiple layers. If \code{NULL}
#'   (default) and \code{prediction} has a single layer, that layer is used.
#' @param effect_direction Character. How the bias surface is applied to the
#'   suitability layer. \code{"direct"} (default) multiplies suitability by
#'   the bias directly — higher bias increases sampling probability.
#'   \code{"inverse"} multiplies by \eqn{1 - \text{bias}} — higher bias
#'   decreases sampling probability.
#' @param verbose Logical. If \code{TRUE} (default), prints progress messages.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts the composite bias surface from \code{prepared_bias}.
#'   \item Verifies both bias and suitability values are within \code{[0, 1]}.
#'   \item Aligns the bias surface to the suitability grid if geometries differ,
#'   using \code{terra::resample()} with nearest-neighbor interpolation.
#'   \item Multiplies suitability by the (possibly inverted) bias surface.
#'   \item Crops and masks the output to the suitability domain.
#' }
#'
#' @return
#' A named list of class \code{"nicheR_biased_surface"} containing:
#' \itemize{
#'   \item One \code{SpatRaster} per input suitability layer, named
#'   \code{"<layer>_biased"}. The raster layer name includes the applied
#'   direction (e.g., \code{"suitability_biased_direct"}).
#'   \item \code{combination_formula}: a character string describing the
#'   operation applied (e.g., \code{"suitability * bias"} or
#'   \code{"suitability * (1-bias)"}).
#' }
#'
#' @seealso \code{\link{prepare_bias}} to build the composite bias surface,
#'   \code{\link{sample_biased_data}} to sample occurrences from the output.
#'
#' @importFrom terra compareGeom resample crop mask global nlyr
#'
#' @export
apply_bias <- function(prepared_bias,
                       prediction,
                       prediction_layer = NULL,
                       effect_direction = "direct",
                       verbose = TRUE){


  gc()
  verbose_message(verbose, "Starting: apply_bias()\n")

  # Basic Input checks --------------------------------------------------------

  if(missing(prepared_bias) || is.null(prepared_bias)){
    stop("'prepared_bias' must be provided.")
  }

  if(missing(prediction) || is.null(prediction)){
    stop("'prediction' must be provided.")
  }

  if(!inherits(prediction, "SpatRaster")){
    stop("'prediction' must be a terra::SpatRaster.")
  }

  prediction <- resolve_prediction(prediction, prediction_layer)$rast


  # 1. Extract composite bias surface ----------------------------------------

  if(inherits(prepared_bias, "SpatRaster")){
    bias_rast <- prepared_bias

  }else if(is.list(prepared_bias)){

    if(!is.null(prepared_bias$composite_surface) &&
       inherits(prepared_bias$composite_surface, "SpatRaster")){
      bias_rast <- prepared_bias$composite_surface

    }else{
      # fallback: allow a named SpatRaster layer if user passed it directly
      stop("If 'prepared_bias' is a list, it must contain a SpatRaster named 'composite_surface'.")
    }

  }else{
    stop("'prepared_bias' must be a SpatRaster or the list output from prepare_bias().")
  }

  if(terra::nlyr(bias_rast) != 1){
    stop("'prepared_bias' must contain exactly 1 layer (a composite bias surface).")
  }

  # 2. Check bias range [0, 1] ----------------------------------------------

  bias_rng <- terra::global(bias_rast,
                            fun = c("min", "max"),
                            na.rm = TRUE)

  bias_min <- as.numeric(bias_rng[1, 1])
  bias_max <- as.numeric(bias_rng[1, 2])

  if(is.finite(bias_min) && bias_min < 0){
    stop("'prepared_bias' has values < 0. Bias must be standardized to [0, 1].")
  }

  if(is.finite(bias_max) && bias_max > 1){
    stop("'prepared_bias' has values > 1. Bias must be standardized to [0, 1].")
  }

  # 3. Check prediction range [0, 1] ------------------------------

  env_rng <- terra::global(prediction,
                           fun = c("min", "max"),
                           na.rm = TRUE)

  env_min <- suppressWarnings(min(as.numeric(env_rng[1, ]), na.rm = TRUE))
  env_max <- suppressWarnings(max(as.numeric(env_rng[2, ]), na.rm = TRUE))

  if(is.finite(env_min) && env_min < 0){
    stop("'prediction' has values < 0. Expected suitability surfaces in [0, 1].")
  }

  if(is.finite(env_max) && env_max > 1){
    stop("'prediction' has values > 1. Expected suitability surfaces in [0, 1].")
  }

  # 4. Effect direction -------------------------------------------------------

  effect_direction <- match.arg(effect_direction,
                                choices = c("direct", "inverse"),
                                several.ok = TRUE)

  if(length(effect_direction) > 1) stop("'effect_direction' can only be a length of 1")

  verbose_message(verbose, "Step: applying bias with '",
                  effect_direction,
                  "' effect to to \"", names(prediction), "\" layer...\n")


  # 5. Align bias to prediction grid -------------------------------

  same_grid <- terra::compareGeom(bias_rast,
                                  prediction[[1]],
                                  stopOnError = FALSE)

  if(!isTRUE(same_grid)){
    verbose_message(verbose, "Step: resampling prepared bias to match prediction grid...\n")
    bias_rast <- terra::resample(bias_rast,
                                 prediction[[1]],
                                 method = "near")
  }

  # 6. Apply bias to each suitable layer -------------------------------------

  out_list <- vector("list", terra::nlyr(prediction))
  formula_entries <- character(terra::nlyr(prediction))


  for(i in 1:terra::nlyr(prediction)){
    s <- prediction[[i]]
    dir_i <- effect_direction[i]

      if(dir_i == "inverse"){
        bias_effect <- 1 - bias_rast
        dir_tag <- "inverse"
      }else{
        bias_effect <- bias_rast
        dir_tag <- "direct"
      }

      # Safe bias name
      bias_name <- names(bias_rast)
      if(is.null(bias_name) || length(bias_name) == 0L || !nzchar(bias_name[1]) || bias_name[1] == "lyr.1"){
        bias_name <- "bias"
      }else{
        bias_name <- bias_name[1]
      }

      # Safe suitability name
      suit_name <- names(s)
      if(is.null(suit_name) || length(suit_name) == 0L || !nzchar(suit_name[1]) || suit_name[1] == "lyr.1"){
        suit_name <- paste0("suitability_", i)
      }else{
        suit_name <- suit_name[1]
      }

      # Formula entry
      if(dir_i == "inverse"){
        formula_entries[i] <- paste0(suit_name, " * (1-", bias_name, ")")
      }else{
        formula_entries[i] <- paste0(suit_name, " * ", bias_name)
      }

      # Compute output raster
      out_r <- s * bias_effect
      out_r <- terra::crop(out_r, prediction[[1]], mask = TRUE)

      # Name list element + raster layer safely
      list_name <- paste0(suit_name, "_biased")
      if(!nzchar(list_name)){
        list_name <- paste0("suitability_", i, "_biased")
      }
      names(out_r) <- paste0(list_name, "_", dir_tag)

      out_list[[i]] <- out_r
      names(out_list)[i] <- list_name
    }


  # 8. Attach message / metadata ---------------------------------------------
  out_list$combination_formula <- formula_entries

  class(out_list) <- "nicheR_biased_surface"

  verbose_message(verbose, "Done: apply_bias(). Note: values are no longer probabilities\n")

  gc()

  out_list
}
