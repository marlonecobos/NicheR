#' Sample occurrence points from a suitable environment pool
#'
#' This function samples occurrence points from a precomputed pool of
#' suitable environments. Sampling can be biased towards the center,
#' the edge, or be purely random. Optionally, a 0–1 sampling bias
#' surface can be supplied to modulate the probability of selecting
#' each suitable location.
#'
#' @param n_occ Integer; the number of occurrence points to sample.
#' @param suitable_env A precomputed suitability object describing the
#'   suitable environment pool. May be:
#'   \itemize{
#'     \item a \code{terra::SpatRaster} or \code{raster::Raster*} with
#'           at least one layer containing squared distances
#'           (\code{dist_sq}), or
#'     \item a \code{data.frame} / \code{matrix} with columns \code{x},
#'           \code{y}, and (for non-random methods) \code{dist_sq},
#'           typically produced by
#'           \code{get_suitable_env(..., out.suit = "data.frame", distances = TRUE)}, or
#'     \item a \code{suitable_env} object returned by
#'           \code{get_suitable_env(..., out.suit = "both")}.
#'   }
#' @param method Character; sampling strategy, one of:
#'   \itemize{
#'     \item \code{"random"}: uniform probability over the suitable pool.
#'     \item \code{"center"}: probability inversely proportional to
#'           Mahalanobis distance (higher near the niche center).
#'     \item \code{"edge"}: probability proportional to Mahalanobis
#'           distance (higher near the niche boundary).
#'   }
#' @param bias_surface Optional 1-layer \code{terra::SpatRaster} or
#'   \code{raster::Raster*} describing sampling bias, scaled to [0, 1]
#'   on the same grid/CRS as the suitability surface. Bias values are
#'   extracted at suitable locations and multiplied with the method-
#'   based weights. Cells with \code{NA} bias are treated as 0 (i.e.,
#'   not sampled).
#' @param seed Optional integer used to set the random number generator
#'   seed for reproducible results.
#' @param verbose Logical. If \code{TRUE}, print progress messages.
#'
#' @return A \code{data.frame} containing \code{n_occ} sampled rows from the
#'   suitable environment pool, including \code{x}, \code{y}, and any
#'   environmental predictor columns. The \code{dist_sq} column is dropped
#'   from the returned data.
#'
#' @details
#' This function assumes that the suitability / distance calculations
#' have already been performed (e.g. via \code{get_suitable_env()}).
#' It only handles the sampling step:
#' \enumerate{
#'   \item Coerce \code{suitable_env} to a data.frame with \code{x}, \code{y},
#'         and (for non-random methods) \code{dist_sq}.
#'   \item Compute method-specific weights.
#'   \item Optionally multiply weights by a user-supplied 0–1 bias surface.
#'   \item Draw \code{n_occ} samples with probability proportional to the
#'         final weights.
#' }
#'
#' @seealso \code{\link{get_suitable_env}}, \code{\link{set_bias_surface}}
#'
#' @export
get_sample_occ <- function(n_occ,
                           suitable_env,
                           method = c("random", "center", "edge"),
                           bias_surface = NULL,
                           seed = NULL,
                           verbose = TRUE) {

  gc()
  method <- tolower(match.arg(method))

  if (isTRUE(verbose)) {
    message("Starting get_sample_occ()...")
  }

  # --- 0) Basic checks -------------------------------------------------------

  if (missing(n_occ) || length(n_occ) != 1 || !is.numeric(n_occ)) {
    stop("'n_occ' must be a single numeric value.")
  }
  n_occ <- as.integer(n_occ)
  if (!is.finite(n_occ) || n_occ <= 0) {
    stop("'n_occ' must be a positive integer.")
  }

  if (is.null(seed) && isTRUE(verbose)) {
    warning("Sampling seed not set; results will differ each time this function is run.",
            call. = FALSE, immediate. = TRUE)
  }

  if (missing(suitable_env) || is.null(suitable_env)) {
    stop("'suitable_env' must be provided.")
  }

  suitable_pool <- NULL
  suitable_rast <- NULL

  if (isTRUE(verbose)) {
    message("Coercing 'suitable_env' to a data.frame of candidate environments...")
  }

  # --- 1) Coerce suitable_env to data.frame with required columns -----------

  ## 1A) Handle 'suitable_env' objects returned by get_suitable_env(...)
  if (inherits(suitable_env, "suitable_env") ||
      (is.list(suitable_env) &&
       any(c("suitable_env_df", "suitable_env_sp") %in% names(suitable_env)))) {

    # Prefer the data.frame if present
    if ("suitable_env_df" %in% names(suitable_env) &&
        is.data.frame(suitable_env$suitable_env_df)) {

      suitable_pool <- suitable_env$suitable_env_df

    } else if ("suitable_env_sp" %in% names(suitable_env)) {

      sp <- suitable_env$suitable_env_sp

      # sp can be a SpatRaster or a list of SpatRaster objects
      if (inherits(sp, "Raster")) {
        sp <- terra::rast(sp)
      }

      if (inherits(sp, "SpatRaster")) {
        suitable_rast <- sp
      } else if (is.list(sp) &&
                 length(sp) > 0 &&
                 all(vapply(sp, inherits, logical(1), "SpatRaster"))) {

        # for out.suit = "both", this is usually a named list: list(suitable = ..., dist_sq = ...)
        if (method == "random") {
          # any raster is fine if we just need x,y; prefer "suitable" if named
          if ("suitable" %in% names(sp)) {
            suitable_rast <- sp[["suitable"]]
          } else {
            suitable_rast <- sp[[1]]
          }
        } else {
          # center / edge → need dist_sq
          if ("dist_sq" %in% names(sp)) {
            suitable_rast <- sp[["dist_sq"]]
          } else {
            # try to find a layer named "dist_sq" inside any raster
            idx_found <- NA_integer_
            for (i in seq_along(sp)) {
              ly_names <- names(sp[[i]])
              if ("dist_sq" %in% ly_names) {
                idx_found <- i
                break
              }
            }
            if (is.na(idx_found)) {
              stop("For method = '", method,
                   "', 'suitable_env' does not contain a 'dist_sq' raster.\n",
                   "Run get_suitable_env(..., distances = TRUE) or pass a data.frame with 'dist_sq'.")
            }
            suitable_rast <- sp[[idx_found]]
          }
        }

      } else {
        stop(
          "'suitable_env$suitable_env_sp' must be a SpatRaster or a list of SpatRaster objects."
        )
      }

      # convert chosen raster to data.frame (x, y, layer values)
      if (!is.null(suitable_rast)) {
        suitable_pool <- as.data.frame.nicheR(suitable_rast)
      }

    } else {
      stop(
        "'suitable_env' object must contain either 'suitable_env_df' or 'suitable_env_sp'."
      )
    }

    ## 1B) Otherwise: raw SpatRaster / Raster* / data.frame / matrix ----------
  } else {

    if (inherits(suitable_env, "Raster")) {
      suitable_env <- terra::rast(suitable_env)
    }

    if (inherits(suitable_env, "SpatRaster")) {

      suitable_rast <- suitable_env
      suitable_pool <- as.data.frame.nicheR(suitable_rast)

    } else if (is.data.frame(suitable_env) || is.matrix(suitable_env)) {

      suitable_pool <- as.data.frame(suitable_env)

    } else {
      stop(
        "'suitable_env' must be a 'suitable_env' object, a terra::SpatRaster, ",
        "raster::Raster*, data.frame, or matrix."
      )
    }
  }

  # Need x, y always
  required_xy <- c("x", "y")
  missing_xy  <- setdiff(required_xy, names(suitable_pool))

  if (length(missing_xy) > 0) {
    stop(
      "'suitable_env' (after coercion) must contain columns: ",
      paste(required_xy, collapse = ", "),
      ". Missing: ", paste(missing_xy, collapse = ", "), "."
    )
  }

  # For center/edge methods, need dist_sq
  if (method %in% c("center", "edge")) {
    if (!"dist_sq" %in% names(suitable_pool)) {
      stop(
        "For method = '", method,
        "', 'suitable_env' must contain a 'dist_sq' column.\n",
        "If needed, run get_suitable_env(..., distances = TRUE)."
      )
    }
  }

  if (nrow(suitable_pool) < 1) {
    stop("No suitable environments were found in 'suitable_env'.")
  }

  if (isTRUE(verbose)) {
    message("Suitable pool contains ", nrow(suitable_pool), " candidate environments.")
  }

  # --- 2) Optional bias surface ---------------------------------------------

  bias_values <- NULL

  if (!is.null(bias_surface)) {

    if (inherits(bias_surface, "Raster")) {
      bias_surface <- terra::rast(bias_surface)
    }

    if (!inherits(bias_surface, "SpatRaster")) {
      stop("'bias_surface' must be a 1-layer terra::SpatRaster or raster::Raster*.")
    }

    if (terra::nlyr(bias_surface) != 1) {
      stop("'bias_surface' must have exactly one layer (0–1 bias values).")
    }

    if (isTRUE(verbose)) {
      message("Extracting bias values at suitable locations...")
    }

    # Check values in [0, 1]
    vals <- terra::values(bias_surface)
    rng  <- range(vals, na.rm = TRUE)

    if (rng[1] < 0 - 1e-6 || rng[2] > 1 + 1e-6) {
      stop(
        "'bias_surface' must be scaled to [0, 1]. ",
        "Observed range: [", signif(rng[1], 3), ", ", signif(rng[2], 3), "]."
      )
    }

    # Extract bias at suitable locations
    coords <- cbind(suitable_pool$x, suitable_pool$y)
    ext_df <- terra::extract(bias_surface, coords)

    bias_values <- ext_df[, 1]
    bias_values[is.na(bias_values)] <- 0

  } else if (isTRUE(verbose)) {

    message("No 'bias_surface' supplied; sampling based only on 'method'.")
  }

  # --- 3) Sampling weights ---------------------------------------------------

  if (isTRUE(verbose)) {
    message("Computing sampling weights using method = '", method, "'...")
  }

  # Base weights from method
  if (method == "random") {

    w <- rep(1, nrow(suitable_pool))

  } else {

    d <- sqrt(pmax(0, suitable_pool$dist_sq))
    # clamp distances to [0, 1] if they were normalized that way
    d <- pmin(d, 1)

    w <- switch(
      method,
      "center" = 1 - d,   # higher near center
      "edge"   = d        # higher near boundary
    )

    w[!is.finite(w)] <- 0
  }

  # Apply bias if available
  if (!is.null(bias_values)) {
    if (isTRUE(verbose)) {
      message("Applying bias surface to sampling weights...")
    }
    w <- w * bias_values
  }

  if (sum(w, na.rm = TRUE) == 0) {
    warning("All sampling weights are zero; falling back to uniform sampling.",
            call. = FALSE)
    w <- rep(1, length(w))
  }

  # --- 4) Draw samples -------------------------------------------------------

  if (!is.null(seed)) {
    if (isTRUE(verbose)) {
      message("Setting seed to ", seed, " and drawing ", n_occ, " samples...")
    }
    set.seed(seed)
  } else if (isTRUE(verbose)) {
    message("Drawing ", n_occ, " samples without a fixed seed...")
  }

  replace_flag <- n_occ > nrow(suitable_pool)

  idx <- sample.int(
    n       = nrow(suitable_pool),
    size    = n_occ,
    replace = replace_flag,
    prob    = w
  )

  occ <- suitable_pool[idx, , drop = FALSE]
  occ$dist_sq <- NULL
  rownames(occ) <- NULL

  if (isTRUE(verbose)) {
    message("Finished get_sample_occ(): sampled ", n_occ, " occurrences.")
  }

  gc()
  return(occ)
}
