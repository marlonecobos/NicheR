#' Extract Suitable Environmental Area from a Niche Ellipsoid
#'
#' This function identifies and extracts environmental grid cells or data
#' points that fall within a defined ellipsoid niche based on Mahalanobis
#' distance.
#'
#' @param niche An object of class \code{ellipsoid} created by
#'   \code{build_ellipsoid()}.
#' @param env_bg A \code{terra::SpatRaster}, \code{data.frame}, or \code{matrix}
#'   of environmental predictor variables. It must contain the variables
#'   referenced by the \code{niche} object. If a \code{data.frame} is used and
#'   spatial output is requested, it should also contain \code{x} and \code{y}
#'   columns for spatial referencing.
#' @param out.suit A character string specifying the desired output format.
#'   One of \code{"data.frame"}, \code{"spatial"}, or \code{"both"}:
#'   \itemize{
#'     \item \code{"data.frame"}: Returns a data frame of all suitable points.
#'     \item \code{"spatial"}: Returns one or more \code{terra::SpatRaster}
#'       objects (see Details).
#'     \item \code{"both"}: Returns a list containing both the spatial
#'       raster output and the data frame of suitable points.
#'   }
#' @param distances Logical; if \code{TRUE}, an additional column named
#'   \code{dist_sq} is added to the output data frame containing the squared
#'   Mahalanobis distance for each suitable point. For spatial output,
#'   a separate raster of \code{dist_sq} values is also returned.
#' @param verbose Logical; if \code{TRUE}, prints basic progress messages.
#'
#' @details
#' For raster inputs, a size check is performed. If the estimated size of a
#' full in-memory data.frame representation of \code{env_bg} exceeds an
#' internal threshold (in MB), the function stops and asks the user to supply
#' a data.frame instead, recommending \code{as.data.frame.nicheR()} as a
#' memory-safe helper.
#'
#' When \code{out.suit} includes spatial output and \code{distances = TRUE},
#' the function returns two separate rasters:
#' \itemize{
#'   \item \code{suitable}: binary 0/1 raster of the suitable area,
#'   \item \code{dist_sq}: raster with squared Mahalanobis distance inside the
#'         ellipsoid (NA outside).
#' }
#'
#' When \code{env_bg} is provided as a data.frame, these rasters are
#' reconstructed from the \code{x} and \code{y} coordinates using
#' \code{terra::rast(..., type = "xyz")}, assuming a regular grid.
#'
#' @return
#' Depending on \code{out.suit}:
#' \itemize{
#'   \item \code{"data.frame"}: a data.frame of suitable points (with optional
#'         \code{dist_sq} column).
#'   \item \code{"spatial"}: a named list of \code{SpatRaster} objects
#'         (e.g. \code{list(suitable = <SpatRaster>, dist_sq = <SpatRaster>)}).
#'   \item \code{"both"}: a list with elements
#'         \code{suitable_env_sp} (spatial output as above) and
#'         \code{suitable_env_df} (data.frame).
#' }
#'
#' @seealso \code{\link{build_ellipsoid}}, \code{\link{get_sample_occ}},
#'   \code{\link{as.data.frame.nicheR}}
#'
#' @export
get_suitable_env <- function(niche,
                             env_bg,
                             out.suit = c("data.frame", "spatial", "both"),
                             distances = FALSE,
                             verbose = TRUE) {

  gc()
  out.suit <- tolower(match.arg(out.suit))

  # --- 1) Check niche structure ----------------------------------------------

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

  if (!is.logical(distances) || length(distances) != 1) {
    stop("'distances' must be a single logical value.")
  }

  if (missing(env_bg) || is.null(env_bg)) {
    stop("'env_bg' is required and cannot be NULL.")
  }

  # --- 2) Coerce env_bg and build env_bg_df ----------------------------------

  # tibble -> data.frame
  if (inherits(env_bg, "tbl_df")) {
    env_bg <- as.data.frame(env_bg)
  }

  # raster::Raster* -> SpatRaster
  if (inherits(env_bg, "Raster")) {
    env_bg <- terra::rast(env_bg)
  }

  env_bg_rast <- NULL
  env_is_raster <- inherits(env_bg, "SpatRaster")

  if (!inherits(env_bg, c("SpatRaster", "data.frame", "matrix"))) {
    stop("'env_bg' must be a terra::SpatRaster, data.frame, or matrix.")
  }

  if (env_is_raster) {
    # size check: approximate memory if fully converted to df
    ncell <- terra::ncell(env_bg)
    nlyr  <- terra::nlyr(env_bg)
    est_mb <- (ncell * nlyr * 8) / 1024^2  # 8 bytes per numeric

    # you can tweak this threshold
    size_threshold_mb <- 5000

    if (est_mb > size_threshold_mb) {
      stop(
        "The provided 'env_bg' raster stack is large (estimated ~",
        round(est_mb, 1),
        " MB if converted to a full data.frame).\n",
        "For memory safety, please convert it to a data.frame yourself, e.g. using\n",
        "  as.data.frame.nicheR(env_bg, use_cache = TRUE)\n",
        "and then pass that data.frame as 'env_bg'."
      )
    }

    env_bg_rast <- env_bg
    # Use your helper for memory-aware conversion
    env_bg_df <- as.data.frame.nicheR(env_bg_rast, verbose = verbose, use_cache = TRUE)

  } else {
    env_bg_df <- as.data.frame(env_bg)
  }

  # --- 3) Determine predictor columns ----------------------------------------

  if (all(c("x", "y") %in% names(env_bg_df))) {
    niche_vars <- setdiff(names(env_bg_df), c("x", "y"))
  } else {
    niche_vars <- names(env_bg_df)
  }

  if (length(niche_vars) < niche$dimen) {
    stop("The number of predictor columns in 'env_bg' is less than 'niche$dimen'.")
  }

  # --- 4) Compute Mahalanobis distance and filter inside ---------------------

  cc <- stats::complete.cases(env_bg_df[, niche_vars, drop = FALSE])
  if (!any(cc)) {
    stop("All candidate rows contain NA in predictor columns. Provide complete predictors or impute values.")
  }

  env_bg_df_cc <- env_bg_df[cc, , drop = FALSE]

  pts   <- as.matrix(env_bg_df_cc[, niche_vars, drop = FALSE])
  diffs <- sweep(pts, 2, niche$center, "-")
  m_sq  <- rowSums((diffs %*% niche$Sigma_inv) * diffs)

  is_inside_cc <- is.finite(m_sq) & (m_sq <= 1)

  if (!any(is_inside_cc)) {
    stop("No points fall inside the ellipsoid after removing rows with NA predictors.")
  }

  inside_rows <- which(cc)[is_inside_cc]

  return_df <- env_bg_df[inside_rows, , drop = FALSE]
  if (isTRUE(distances)) {
    return_df$dist_sq <- m_sq[is_inside_cc]
  }

  # --- 5) Spatial output (binary + optional distance rasters) ----------------

  suitable_sp_list <- NULL

  if (out.suit %in% c("spatial", "both")) {

    if (!all(c("x", "y") %in% names(env_bg_df))) {
      stop("For spatial output, 'env_bg' (or its data.frame representation) must contain 'x' and 'y' columns.")
    }

    # we will always build *two separate rasters* when distances = TRUE,
    # and a single raster when distances = FALSE.

    if (env_is_raster && !is.null(env_bg_rast)) {
      # Case A: small raster – use it as template, avoid re-building from xyz

      # 5A.1 suitable raster (0/1)
      suitable_ras <- env_bg_rast[[1]]
      vals <- terra::values(suitable_ras)
      vals[!is.na(vals)] <- NA_real_
      terra::values(suitable_ras) <- vals

      # coordinates for all "inside" rows
      xy_all <- env_bg_df[inside_rows, c("x", "y"), drop = FALSE]
      ok_xy  <- stats::complete.cases(xy_all)

      if (!any(ok_xy)) {
        stop("No valid XY coordinates found for suitable environments.")
      }

      xy_ok <- as.matrix(xy_all[ok_xy, , drop = FALSE])
      inside_cells <- terra::cellFromXY(env_bg_rast, xy_ok)
      inside_cells <- inside_cells[!is.na(inside_cells)]

      if (length(inside_cells)) {
        suitable_ras[as.integer(inside_cells)] <- 1
      }
      names(suitable_ras) <- "suitable"

      suitable_sp_list <- list(suitable = suitable_ras)

      # 5A.2 distance raster (only if requested)
      if (isTRUE(distances)) {
        dist_ras <- env_bg_rast[[1]]
        vals_d <- terra::values(dist_ras)
        vals_d[] <- NA_real_
        terra::values(dist_ras) <- vals_d

        # m_sq[is_inside_cc] is aligned with inside_rows and xy_all
        m_sq_inside <- m_sq[is_inside_cc]
        dist_vals   <- m_sq_inside[ok_xy]

        if (length(dist_vals) != length(inside_cells)) {
          warning(
            "Length mismatch when assigning distance values to raster cells; ",
            "using matched subset."
          )
          len <- min(length(dist_vals), length(inside_cells))
          dist_vals    <- dist_vals[seq_len(len)]
          inside_cells <- inside_cells[seq_len(len)]
        }

        dist_ras[as.integer(inside_cells)] <- dist_vals
        names(dist_ras) <- "dist_sq"

        suitable_sp_list$dist_sq <- dist_ras
      }
    } else {
      # Case B: env_bg was provided as data.frame / matrix
      # Build rasters from full env_bg_df grid using xyz representation

      if (!all(c("x", "y") %in% names(env_bg_df))) {
        stop("For spatial output from data.frame 'env_bg', columns 'x' and 'y' are required.")
      }

      # track inside_rows indices into env_bg_df
      # Build binary suitable column on full grid
      env_xy <- env_bg_df[, c("x", "y"), drop = FALSE]
      env_xy$suitable <- NA_real_
      env_xy$suitable[inside_rows] <- 1

      # suitable raster
      suitable_ras <- terra::rast(env_xy[, c("x", "y", "suitable")],
                                  type = "xyz")
      names(suitable_ras) <- "suitable"

      suitable_sp_list <- list(suitable = suitable_ras)

      if (isTRUE(distances)) {
        env_xy$dist_sq <- NA_real_
        env_xy$dist_sq[inside_rows] <- m_sq[is_inside_cc]

        dist_ras <- terra::rast(env_xy[, c("x", "y", "dist_sq")],
                                type = "xyz")
        names(dist_ras) <- "dist_sq"

        suitable_sp_list$dist_sq <- dist_ras
      }
    }
  }

  # --- 6) Pack and return ----------------------------------------------------

  res <- switch(
    out.suit,
    "data.frame" = return_df,
    "spatial"    = suitable_sp_list,
    "both"       = {
      tmp <- list(
        suitable_env_sp = suitable_sp_list,
        suitable_env_df = return_df
      )
      class(tmp) <- c("suitable_env", class(tmp))
      tmp
    }
  )

  gc()
  return(res)
}

#' @export
print.suitable_env <- function(x, ...) {

  if (is.list(x) && all(c("suitable_env_sp", "suitable_env_df") %in% names(x))) {

    cat("Suitable environment object:\n")

    # Spatial part
    sp <- x$suitable_env_sp
    if (inherits(sp, "SpatRaster")) {
      cat("  Spatial (SpatRaster):\n")
      print(sp)
    } else if (is.list(sp) &&
               length(sp) > 0 &&
               all(vapply(sp, inherits, logical(1), "SpatRaster"))) {

      cat("  Spatial rasters (list):\n")
      for (nm in names(sp)) {
        cat("   -", nm, ":\n")
        print(sp[[nm]])
      }
    } else if (!is.null(sp)) {
      cat("  Spatial component present but not a SpatRaster or list of SpatRasters.\n")
      str(sp)
    }

    cat("\n  Data frame (showing first 6 rows):\n")
    print(utils::head(x$suitable_env_df))

  } else if (inherits(x, "SpatRaster")) {

    cat("Suitable environment raster:\n")
    print(x)

  } else if (is.list(x) &&
             length(x) > 0 &&
             all(vapply(x, inherits, logical(1), "SpatRaster"))) {

    cat("Suitable environment rasters (list):\n")
    for (nm in names(x)) {
      cat(" -", nm, ":\n")
      print(x[[nm]])
    }

  } else if (is.data.frame(x)) {

    cat("Suitable environment data frame (showing first 6 rows):\n")
    print(utils::head(x))

  } else {
    NextMethod()
  }

  invisible(x)
}
