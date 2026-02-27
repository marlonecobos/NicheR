#' Apply Sampling Bias to Suitability Surfaces
#'
#' Applies a prepared composite sampling bias surface to one or more suitability
#' rasters by multiplication. The bias surface is aligned to the suitability grid
#' when needed, optionally inverted per suitability layer, and the result is
#' cropped and masked to the suitability domain.
#'
#' @details
#' This function expects a single-layer composite bias surface (from
#' \code{\link{prepare_bias}}) and a \code{terra::SpatRaster} of suitability values
#' (typically the raster output from a prediction step). It performs the following:
#' \enumerate{
#'   \item Extracts a single-layer composite bias surface from \code{prepared_bias}.
#'   \item Verifies both bias and suitability values are within \eqn{[0, 1]}.
#'   \item Aligns the bias surface to the suitability grid if geometries differ
#'   (\code{terra::compareGeom()} + \code{terra::resample(method = "near")}).
#'   \item Applies the bias to each suitability layer via multiplication.
#'   \item Crops and masks each output layer to the suitability domain.
#' }
#'
#' The output is a product of suitability and bias surfaces and therefore is not
#' interpretable as a probability.
#'
#' @param prepared_bias A single-layer \code{terra::SpatRaster} composite bias
#'   surface, or the list output returned by \code{\link{prepare_bias}} containing
#'   \code{composite_surface}.
#' @param prediction A \code{terra::SpatRaster} containing one or more
#'   suitability layers with values in \eqn{[0, 1]}.
#' @param effect_direction Character vector indicating how bias affects each
#'   suitability layer. Options are \code{"direct"} or \code{"inverse"}.
#'   A single value is recycled to match the number of suitability layers.
#'   \code{"inverse"} applies \eqn{1 - bias} before multiplication.
#' @param verbose Logical. If \code{TRUE}, prints progress messages.
#'
#' @return
#' A named list of class \code{"nicheR_biased_suitability"} containing:
#' \itemize{
#'   \item One element per input suitability layer, each a single-layer
#'   \code{terra::SpatRaster} of biased suitability (cropped and masked to the
#'   suitability domain). Elements are named \code{"<layer>_biased"} and the raster
#'   layer name includes the applied direction (e.g., \code{"<layer>_biased_direct"}).
#'   \item \code{combination_formula} Character vector (length = number of
#'   suitability layers) describing the operation applied to each layer
#'   (e.g., \code{"suitability * bias"} or \code{"suitability * (1-bias)"}).
#' }
#'
#' @seealso \code{\link{prepare_bias}}
#'
#' @export
apply_bias <- function(prepared_bias,
                       prediction,
                       effect_direction = "direct",
                       verbose = TRUE){

  verbose_message <- function(...) if(isTRUE(verbose)) cat(...)
  gc()


  verbose_message("Starting: apply_bias()\n")

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

  if(length(effect_direction) == 1L){
    effect_direction <- rep(effect_direction, terra::nlyr(prediction))

    verbose_message("Step: applying bias with '",
                    effect_direction[1],
                    "' effect to all suitability layer(s)...\n")

  }else if(length(effect_direction) != terra::nlyr(prediction)){
    # recycle to match layers
    effect_direction <- rep(effect_direction,
                            length.out = terra::nlyr(prediction))

    verbose_message("Step: applying layer-specific bias directions with repetition: ",
                    paste(effect_direction, collapse = ", "),
                    "\n")
  }else{
    verbose_message("Step: applying layer-specific bias directions: ",
                    paste(effect_direction, collapse = ", "),
                    "\n")
  }

  # 5. Align bias to prediction grid -------------------------------

  same_grid <- terra::compareGeom(bias_rast,
                                  prediction[[1]],
                                  stopOnError = FALSE)

  if(!isTRUE(same_grid)){
    verbose_message("Step: resampling prepared bias to match prediction grid...\n")
    bias_rast <- terra::resample(bias_rast,
                                 prediction[[1]],
                                 method = "near")
  }

  # 6. Apply bias to each suitable layer -------------------------------------

  verbose_message("Step: applying bias to", terra::nlyr(prediction), "prediction layer/s...\n")

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

  verbose_message("Done: apply_bias(). Note: values are no longer probabilities\n")

  gc()

  out_list
}
