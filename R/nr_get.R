#' Universal Getter for NicheR Objects
#'
#' `nr_get()` provides a unified accessor interface for all core NicheR object
#' types, including:
#' \itemize{
#'   \item `suitable_env` objects returned by \code{\link{get_suitable_env}},
#'   \item `nicheR_bias_surface` objects from \code{\link{set_bias_surface_new}},
#'   \item full `NicheR_species` objects produced by
#'         \code{\link{create_virtual_species}},
#'   \item lists of `SpatRaster` objects,
#'   \item direct `terra::SpatRaster` inputs.
#' }
#'
#' The function searches intelligently through nested objects and lists and
#' returns the requested component without the user needing to know the
#' object structure. This makes it suitable for internal use in the NicheR
#' workflow and for user-facing utilities.
#'
#' @section What You Can Extract:
#' The `what` argument controls which component is returned:
#'
#' \describe{
#'   \item{\code{"suitable"}}{A binary (0–1) raster of suitable areas, if present.}
#'
#'   \item{\code{"dist_sq"}}{Raster of squared Mahalanobis distance from the
#'   niche centroid (inside the ellipsoid), if present.}
#'
#'   \item{\code{"all_suitable"}}{The full suitability object stored in the input:
#'     for a \code{suitable_env} this returns the entire object; for a
#'     \code{NicheR_species} this returns the \code{$suitability} slot; for
#'     lists or \code{SpatRaster} inputs it returns the object itself.}
#'
#'   \item{\code{"df"}}{Data frame of suitable points, if stored.}
#'
#'   \item{\code{"env"}}{The environmental background (`env_bg`) originally
#'   supplied to the virtual species workflow.}
#'
#'   \item{\code{"niche"}}{The ellipsoid object created by
#'   \code{\link{build_ellps}}.}
#'
#'   \item{\code{"occ"}}{Occurrence data sampled by
#'   \code{\link{get_sample_occ}}.}
#'
#'   \item{\code{"bias"}}{The pooled bias surface (a single-layer SpatRaster).}
#'
#'   \item{\code{"bias_pooled"}}{Alias for `"bias"`.}
#'
#'   \item{\code{"bias_stack"}}{The stack of standardized, direction-adjusted
#'   bias layers.}
#' }
#'
#' The function returns `NULL` if the requested component cannot be located.
#'
#'
#' @param obj Any NicheR object: a `NicheR_species`, a `suitable_env`, a
#'   `nicheR_bias_surface`, a list, or a direct `SpatRaster`.
#'
#' @param what A character string specifying the item to retrieve.
#'   One of:
#'   `"suitable"`, `"dist_sq"`, `"all_suitable"`, `"df"`,
#'   `"env"`, `"niche"`, `"occ"`,
#'   `"bias"`, `"bias_pooled"`, `"bias_stack"`.
#'
#' @return
#' A wide variety of possible outputs depending on `what`:
#' \itemize{
#'   \item A \code{SpatRaster},
#'   \item A list of SpatRasters,
#'   \item A data frame,
#'   \item An ellipsoid object,
#'   \item A vector or data frame of occurrences,
#'   \item \code{NULL} if not found.
#' }
#'
#' @details
#' This function uses recursive searching to extract nested components from
#' deeply embedded lists. It allows downstream functions (e.g.,
#' \code{\link{set_bias_surface_new}}, \code{\link{plot_g_space}}) to accept
#' any NicheR object and internally resolve the required raster or data component.
#'
#' @seealso
#'   \code{\link{create_virtual_species}},
#'   \code{\link{get_suitable_env}},
#'   \code{\link{set_bias_surface_new}},
#'   \code{\link{get_sample_occ}}
#'
#' @export
nr_get <- function(obj, what = c(
  "suitable", "dist_sq", "all_suitable", "df",
  "env", "niche", "occ",
  "bias", "bias_pooled", "bias_stack")) {

  what <- match.arg(what)

  # 0. Direct SpatRaster
  if (inherits(obj, "SpatRaster")) {
    if (what %in% c("suitable", "dist_sq", "bias", "bias_pooled", "all_suitable")) {
      return(obj)
    }
  }

  # 1. Suitable environment object (class suitable_env)
  if (inherits(obj, "suitable_env")) {

    sp <- obj$suitable_env_sp
    df <- obj$suitable_env_df

    if (what == "df") return(df)

    if (what == "suitable" && "suitable" %in% names(sp)) {
      return(sp$suitable)
    }

    if (what == "dist_sq" && "dist_sq" %in% names(sp)) {
      return(sp$dist_sq)
    }

    if (what == "all_suitable") {
      # NEW: return the full suitable_env object, not just the SpatRaster stack
      return(obj)
    }
  }

  # 2. Direct spatial-only output from get_suitable_env (list of rasters)
  if (is.list(obj) && all(vapply(obj, inherits, logical(1), "SpatRaster"))) {

    if (what == "suitable" && "suitable" %in% names(obj))
      return(obj$suitable)

    if (what == "dist_sq" && "dist_sq" %in% names(obj))
      return(obj$dist_sq)

    if (what == "all_suitable")
      return(obj)
  }

  # 3. bias_surface object
  if (inherits(obj, "nicheR_bias_surface")) {

    if (what == "bias" || what == "bias_pooled")
      return(obj$pooled_bias_sp)

    if (what == "bias_stack")
      return(obj$directional_bias_stack)
  }

  # 4. NicheR_species object
  if (inherits(obj, "NicheR_species")) {

    # niche
    if (what == "niche") return(obj$niche)

    # occurrences
    if (what == "occ") return(obj$occurrences)

    # env
    if (what == "env") {
      if (!is.null(obj$call_args$env_bg)) return(obj$call_args$env_bg)
    }

    # suitability object
    suit <- obj$suitability

    if (!is.null(suit)) {

      if (what == "all_suitable")
        # NEW: return entire suitability slot (e.g. suitable_env object)
        return(suit)

      # ask recursively inside suitability for specific components
      if (what == "suitable")
        return(nr_get(suit, "suitable"))

      if (what == "dist_sq")
        return(nr_get(suit, "dist_sq"))

      if (what == "df")
        return(nr_get(suit, "df"))
    }
  }

  # 5. Try nested list search (deep recursive)
  if (is.list(obj)) {

    # direct match
    if (what %in% names(obj))
      return(obj[[what]])

    # deep search
    for (el in obj) {
      hit <- nr_get(el, what)
      if (!is.null(hit))
        return(hit)
    }
  }

  return(NULL)
}


#' @describeIn nr_get Extract the suitable-area raster.
#' @export
nr_get_suitable <- function(obj) nr_get(obj, "suitable")

#' @describeIn nr_get Extract the squared-distance raster.
#' @export
nr_get_dist_sq <- function(obj) nr_get(obj, "dist_sq")

#' @describeIn nr_get Extract all suitability rasters.
#' @export
nr_get_suitable_all <- function(obj) nr_get(obj, "all_suitable")

#' @describeIn nr_get Extract suitability as a data frame.
#' @export
nr_get_suitable_df <- function(obj) nr_get(obj, "df")

#' @describeIn nr_get Extract occurrence data.
#' @export
nr_get_occ <- function(obj) nr_get(obj, "occ")

#' @describeIn nr_get Extract the ellipsoid niche object.
#' @export
nr_get_niche <- function(obj) nr_get(obj, "niche")

#' @describeIn nr_get Extract the environmental background.
#' @export
nr_get_env <- function(obj) nr_get(obj, "env")

#' @describeIn nr_get Extract pooled bias surface.
#' @export
nr_get_bias <- function(obj) nr_get(obj, "bias")

#' @describeIn nr_get Pooled bias alias.
#' @export
nr_get_bias_pooled <- function(obj) nr_get(obj, "bias_pooled")

#' @describeIn nr_get Extract standardized bias stack.
#' @export
nr_get_bias_stack <- function(obj) nr_get(obj, "bias_stack")


# Testing -----------------------------------------------------------------

# suit <- suitable_area_t4
#
# nr_get_suitable(suit)      # returns SpatRaster
# nr_get_dist_sq(suit)       # returns SpatRaster
# nr_get_suitable_df(suit)   # returns data frame
# nr_get_suitable_all(suit)  # returns list of rasters
#
# weird <- list(a = list(b = list(suit)))
# nr_get_suitable(weird)
# nr_get_dist_sq(weird)       # returns SpatRaster
#
# nr_get_bias_pooled(bias_obj)
# nr_get_bias_stack(bias_obj)
#
# nr_get_suitable(vs)
# nr_get_dist_sq(vs)
# nr_get_occ(vs)
# nr_get_env(vs)


