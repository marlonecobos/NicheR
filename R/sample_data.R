
#' @export
resolve_prediction <- function(prediction, prediction_layer){

  # ---- Case 1: prediction is a data.frame -----------------------------------
  if(is.data.frame(prediction)){

    if(is.null(prediction_layer) || !nzchar(prediction_layer)){
      stop("If 'prediction' is a data.frame, 'prediction_layer' must be provided (column name).")
    }

    if(!all(c("x", "y") %in% names(prediction))){
      stop("If 'prediction' is a data.frame, it must contain columns named 'x' and 'y'.")
    }

    if(!(prediction_layer %in% names(prediction))){
      stop("Column '", prediction_layer, "' not found in data.frame prediction.")
    }

    return(list(type = "data.frame",
                df = prediction,
                pred_name = prediction_layer))
  }

  # ---- Case 2: prediction is a SpatRaster -----------------------------------
  if(inherits(prediction, "SpatRaster")){

    if(terra::nlyr(prediction) == 1L){

      if(!is.null(prediction_layer) && nzchar(prediction_layer)){
        if(is.null(names(prediction)) || !(prediction_layer %in% names(prediction))){
          stop("Layer '", prediction_layer, "' not found in SpatRaster prediction.")
        }
        prediction <- prediction[[prediction_layer]]
      }

      return(list(type = "SpatRaster", rast = prediction, pred_name = names(prediction)[1]))
    }

    # multi-layer raster: require layer name
    if(is.null(prediction_layer) || !nzchar(prediction_layer)){
      stop("If 'prediction' is a multi-layer SpatRaster, 'prediction_layer' must be provided (layer name).")
    }

    if(is.null(names(prediction)) || !(prediction_layer %in% names(prediction))){
      stop("Layer '", prediction_layer, "' not found in SpatRaster prediction. Available: ",
           paste(names(prediction), collapse = ", "))
    }

    r <- prediction[[prediction_layer]]
    return(list(type = "SpatRaster", rast = r, pred_name = prediction_layer))
  }

  # ---- Case 3: prediction is a list -----------------------------------------
  if(is.list(prediction)){

    if(length(prediction) == 0L){
      stop("'prediction' is an empty list.")
    }

    # If list has one element: repeat rules above
    if(length(prediction) == 1L){
      return(resolve_prediction(prediction[[1]], prediction_layer))
    }

    # list has multiple elements: require prediction_layer
    if(is.null(prediction_layer) || !nzchar(prediction_layer)){
      stop("If 'prediction' is a list with multiple elements, 'prediction_layer' must be provided.")
    }

    # (A) First: check if prediction_layer matches list names
    if(!is.null(names(prediction)) && prediction_layer %in% names(prediction)){
      return(resolve_prediction(prediction[[prediction_layer]], prediction_layer))
      # note: passing prediction_layer again is fine; for SpatRaster single-layer it will no-op,
      # and for multi-layer it will select inside that element if needed.
    }

    # (B) Otherwise: go element-by-element to find prediction_layer
    for(i in seq_along(prediction)){
      obj <- prediction[[i]]

      if(is.data.frame(obj)){
        if(prediction_layer %in% names(obj)){
          return(list(type = "data.frame", df = obj, pred_name = prediction_layer))
        }
      }

      if(inherits(obj, "SpatRaster")){
        if(!is.null(names(obj)) && prediction_layer %in% names(obj)){
          return(list(type = "SpatRaster", rast = obj[[prediction_layer]], pred_name = prediction_layer))
        }
      }
    }

    stop("Could not find 'prediction_layer' in list names or within any list element.")
  }

  stop("'prediction' must be a SpatRaster, a data.frame, or a list.")
}

#' @export
sample_data <- function(n_occ,
                        prediction,
                        prediction_layer = NULL,
                        sampling = "weighted",
                        biased = FALSE,
                        sampling_mask = NULL,
                        seed = 1,
                        verbose = TRUE){

  verbose_message <- function(...) if(isTRUE(verbose)) cat(...)
  gc()

  verbose_message("Starting: sample_data()\n")

  # Basic Input checks --------------------------------------------------------

  sampling <- match.arg(
    sampling,
    choices = c("weighted", "center", "edge", "random"),
    several.ok = FALSE
  )

  if(!is.numeric(n_occ) || length(n_occ) != 1L || is.na(n_occ) || n_occ <= 0){
    stop("'n_occ' must be a single positive number.")
  }
  n_occ <- as.integer(n_occ)

  if(!is.logical(biased) || length(biased) != 1L){
    stop("'biased' must be TRUE or FALSE.")
  }

  if(!is.numeric(seed) || length(seed) != 1L || is.na(seed)){
    stop("'seed' must be a single number.")
  }

  if(isTRUE(biased) && sampling != "weighted"){
    stop("If 'biased = TRUE', sampling must be 'weighted'.")
  }

  # Resolve prediction input --------------------------------------------------

  resolved <- resolve_prediction(prediction, prediction_layer)

  # Enforce sampling_mask rules ----------------------------------------------

  if(!is.null(sampling_mask) && resolved$type != "SpatRaster"){
    stop("'sampling_mask' is only supported when 'prediction' resolves to a SpatRaster.")
  }

  # Apply sampling mask (raster only) ----------------------------------------

  if(resolved$type == "SpatRaster"){
    r <- resolved$rast

    if(!is.null(sampling_mask)){
      if(inherits(sampling_mask, "SpatRaster")){
        r <- terra::mask(r, sampling_mask)
      } else if(inherits(sampling_mask, "SpatVector")){
        r <- terra::mask(r, sampling_mask)
      } else {
        stop("'sampling_mask' must be a SpatRaster or SpatVector.")
      }
    }

    df <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
    if(nrow(df) == 0L){
      stop("No non-NA cells available for sampling.")
    }

    pred_cols <- setdiff(names(df), c("x", "y"))
    if(length(pred_cols) != 1L){
      stop("Unexpected: raster extraction did not produce exactly one prediction column.")
    }

    df$pred <- df[[pred_cols[1]]]

  } else {

    df <- resolved$df
    df$pred <- df[[resolved$pred_name]]
    df <- df[is.finite(df$pred), , drop = FALSE]

    if(nrow(df) == 0L){
      stop("No finite prediction values available for sampling.")
    }
  }

  if(n_occ > nrow(df)){
    stop("'n_occ' is larger than the number of available samples.")
  }

  # Determine value type ------------------------------------------------------

  tol <- 1e-9
  eps <- 1e-12

  zmin <- suppressWarnings(min(df$pred, na.rm = TRUE))
  zmax <- suppressWarnings(max(df$pred, na.rm = TRUE))
  standardized <- is.finite(zmin) && is.finite(zmax) &&
    (zmin >= -tol) && (zmax <= 1 + tol)

  if(isTRUE(standardized)){
    df$pred <- pmax(pmin(df$pred, 1), 0)
  }

  # Weights -------------------------------------------------------------------

  if(sampling == "random"){

    w <- rep(1, nrow(df))

  } else if(sampling == "weighted"){

    w <- df$pred
    w[!is.finite(w)] <- NA_real_

    if(any(w < 0, na.rm = TRUE)){
      stop("Weighted sampling requires non-negative values. Found values < 0.")
    }

    w <- w + eps

  } else {

    if(isTRUE(standardized)){

      if(sampling == "center"){
        w <- df$pred + eps
      } else if(sampling == "edge"){
        w <- (1 - df$pred) + eps
      }

    } else {

      d <- df$pred
      d[!is.finite(d)] <- NA_real_

      if(any(d < 0, na.rm = TRUE)){
        stop("For non-standardized surfaces, expected non-negative distance-like values for center/edge.")
      }

      if(sampling == "center"){
        w <- 1/(d + eps)
      } else if(sampling == "edge"){
        w <- d + eps
      }
    }
  }

  w[!is.finite(w)] <- 0
  if(sum(w) <= 0){
    stop("Sampling weights are all zero. Check your inputs.")
  }

  # Sample --------------------------------------------------------------------

  set.seed(seed)
  idx <- sample.int(nrow(df), size = n_occ, replace = FALSE, prob = w)
  df$pred <- NULL
  out <- df[idx, , drop = FALSE]

  verbose_message("Done: sampled ", nrow(out), " occurrences.\n")
  gc()

  out
}


