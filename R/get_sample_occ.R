# Sample Occurrence Points from a Suitable Environment
library(dplyr)


get_sample_occ <- function(n_occ,
                           env_bg, # accepts Spat.Rasters and data.frames
                           niche,
                           method = c("random", "center", "edge"),
                           seed = NULL) {

  method <- tolower(match.arg(method))


  # --- 1) Input validation and coercion ---

  # 1a) Check niche structure early
  if (!inherits(niche, "ellipsoid")) {
    stop("'niche' must be an object of class 'ellipsoid' produced by build_ellipsoid().")
  }

  need <- c("center", "Sigma_inv", "dimen")
  miss <- setdiff(need, names(niche))

  if (length(miss)) {
    stop(sprintf("'niche' is missing required fields: %s", paste(miss, collapse = ", ")))
  }

  if (!(is.numeric(niche$center) && length(niche$center) == niche$dimen)) {
    stop("'niche$center' must be numeric with length equal to 'niche$dimen'.")
  }

  if (!is.matrix(niche$Sigma_inv) || any(!is.finite(niche$Sigma_inv))) {
    stop("'niche$Sigma_inv' must be a finite numeric matrix.")
  }

  # 1b) Accept tibble -> data.frame (keeps names and types)
  if (inherits(env_bg, "tbl_df")) {
    env_bg <- as.data.frame(env_bg)
  }

  # 1c) Accept raster::Raster* by converting to terra::SpatRaster
  if (inherits(env_bg, "Raster")) {
    env_bg <- terra::rast(env_bg)
  }

  # 1d) Validate env_bg type now
  if (!inherits(env_bg, c("SpatRaster", "data.frame", "matrix"))) {
    stop("'env_bg' must be a terra::SpatRaster, data.frame, or matrix.")
  }

  # 1e) Branch on desired output type
  if (out %in% c("spatial", "both")) {
    # For spatial outputs we REQUIRE a SpatRaster so geometry is preserved
    if (!inherits(env_bg, "SpatRaster")) {
      stop("For spatial output, 'env_bg' must be a terra::SpatRaster.")
    }

    # Build a data.frame with XY and layer values for lookups later
    env_bg_df <- terra::as.data.frame(env_bg, xy = TRUE, na.rm = FALSE)

    # Ensure XY names exist
    if (!all(c("x", "y") %in% names(env_bg_df))) {
      stop("Could not find 'x' and 'y' columns after converting raster to data.frame.")
    }

    # Predictor columns will be all raster layers (exclude x,y)
    niche_vars <- setdiff(names(env_bg_df), c("x", "y"))

    if (length(niche_vars) < niche$dimen) {
      stop("Raster has fewer predictor layers than 'niche$dimen'.")
    }

  } else {
    # data.frame or matrix path (non-spatial outputs)
    if (!(is.matrix(env_bg) || is.data.frame(env_bg))) {
      stop("For 'data.frame' output, 'env_bg' must be a matrix or data.frame.")
    }

    env_bg_df <- as.data.frame(env_bg)
    niche_vars <- setdiff(names(env_bg_df), c("x", "y"))

    # Ensure XY names exist
    if (!all(c("x", "y") %in% names(env_bg_df))) {
      stop("'evn_bg' must include spatial columns, x and y are necessary. Update data.frame or use Spat.Raster")
    }

  }

  # Calculate mahalanobis distance
  suitable_pool <- get_suitable_env(niche,
                                    env_bg,
                                    out = "data.frame",
                                    distances = TRUE)

  # If there is no suitable area send message
  if (nrow(suitable_pool) < 1)
    stop("No suitable environments were found in the provided niche space.")

  # If suitable are exists, continue

  # --- 3. Sampling weights by method ---
  # d in [0, 1] (distance from center to boundary in Mahalanobis metric)
  # Making sure that there are numbers not larger than 1 or less than 0?
  d <- sqrt(pmax(0, suitable_pool$dist_sq))
  d <- pmin(d, 1)

  w <- switch(
    method,
    "random" = rep(1, nrow(suitable_pool)),
    "center" = 1 - d,   # higher near center
    "edge"   = d        # higher near boundary
  )

  w[!is.finite(w)] <- 0

  if (sum(w) == 0) {
    warning("All weights are zero; falling back to uniform sampling.",
            call. = FALSE)
    w <- rep(1, length(w))
  }

  # --- 4. Sample indices (ensure we return exactly n rows) ---
  if (!is.null(seed)) set.seed(seed)
  replace_flag <- n_occ > nrow(suitable_pool)

  idx <- sample.int(n = nrow(suitable_pool),
                    size = n_occ,
                    replace = replace_flag,
                    prob = w)

  # --- 5. Return sampled rows (drop temp column) ---
  occ <- suitable_pool[idx, , drop = FALSE]

  # Remove distances
  occ$dist_sq <- NULL
  rownames(occ) <- NULL

  message(sprintf("Done sampling %d occurrences", n_occ))
  return(occ)
}
