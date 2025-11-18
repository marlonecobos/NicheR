#' Subset environments to a bounding box around an ellipsoid in E space
#'
#' This function uses the ellipsoid surface stored in a \code{niche} object
#' to build an axis aligned rectangular bounding box in environmental space.
#' The box is defined by the minimum and maximum values of each axis on the
#' ellipsoid surface, optionally expanded by a user specified factor. All
#' rows in \code{env_df} that fall inside this multidimensional box are kept.
#'
#' @param env_df A \code{data.frame} of background environmental values.
#'   It must contain the columns listed in \code{vars}.
#' @param niche An ellipsoid object returned by \code{build_ellps()} that
#'   contains a \code{$surface} component with coordinates in E space.
#' @param vars A character vector giving the names of the columns in
#'   \code{env_df} that correspond to the environmental axes of the ellipsoid.
#'   The order of \code{vars} must match the column order of \code{niche$surface}.
#' @param expand A numeric expansion factor applied to the width of the box
#'   along each axis. For example, \code{expand = 0.1} expands the minimum
#'   and maximum of each axis by ten percent of that axis range.
#'
#' @return A \code{data.frame} containing only the rows of \code{env_df}
#'   whose values on all \code{vars} fall inside the constructed bounding
#'   box in E space.
subset_env_to_niche_box <- function(env_df,
                                    niche,
                                    vars,
                                    expand = 0.1) {
  if (is.null(niche$surface)) {
    stop("`niche$surface` is required to build the bounding box.")
  }
  if (missing(vars) || is.null(vars)) {
    stop("`vars` must be provided and must match the dimensions of `niche$surface`.")
  }

  surf <- as.matrix(niche$surface)

  if (ncol(surf) != length(vars)) {
    stop("Number of columns in `niche$surface` (",
         ncol(surf), ") does not match length of `vars` (", length(vars), ").")
  }
  if (!all(vars %in% names(env_df))) {
    missing_cols <- vars[!vars %in% names(env_df)]
    stop("These `vars` are not columns in `env_df`: ",
         paste(missing_cols, collapse = ", "))
  }

  # 1D ranges per axis
  box <- apply(surf, 2, range, na.rm = TRUE)
  width <- box[2, ] - box[1, ]
  box[1, ] <- box[1, ] - expand * width
  box[2, ] <- box[2, ] + expand * width

  keep <- rep(TRUE, nrow(env_df))
  for (i in seq_along(vars)) {
    v <- vars[i]
    keep <- keep &
      env_df[[v]] >= box[1, i] &
      env_df[[v]] <= box[2, i]
  }

  env_df[keep, , drop = FALSE]
}
