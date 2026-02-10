#' Compute Ellipsoidal Niche Volume
#'
#' Calculates the geometric volume (or area in 2D) of a dimensions-dimensional ellipsoid
#' defined by its semi-axis lengths.
#'
#' For an ellipsoid with semi-axes \eqn{a_1, a_2, \ldots, a_p}, the volume is:
#'
#' \deqn{
#' V_p = V_{\text{n},dimensions} \prod_{i=1}^{dimensions} a_i
#' }
#'
#' where the volume of the dimensions-dimensional n is:
#'
#' \deqn{
#' V_{\text{n},dimensions} =
#' \frac{\pi^{dimensions/2}}{\Gamma\left(\frac{dimensions}{2} + 1\right)}
#' }
#'
#' In the context of probabilistic niche modeling, the semi-axis lengths are
#' typically defined as:
#'
#' \deqn{
#' a_i = \sqrt{\lambda_i \, c^2}
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{\lambda_i} are the eigenvalues of the covariance matrix
#'   \eqn{\cov_matrix}
#'   \item \eqn{c^2 = \chi^2_p(\text{level})} is the chi-square cutoff defining the
#'   probability contour
#' }
#'
#' This formulation allows ellipsoid volume to be computed directly from either
#' covariance matrices or axis lengths derived from them.
#'
#' @param n_dimensions Integer. Number of dimensions (\eqn{dimensions}).
#' @param semi_axes_lengths Numeric vector of length \code{n_dimensions}
#'   containing the ellipsoid semi-axis lengths.
#'
#' @return Numeric. Geometric volume of the ellipsoid.
#'
#' @details
#' \strong{Interpretation:}
#'
#' \itemize{
#'   \item In 2D, the returned value is area.
#'   \item In 3D, the returned value is volume.
#'   \item In higher dimensions, the result is a hypervolume measured in the
#'     product of environmental units.
#' }
#'
#' The volume depends on the selected probability contour (\code{level}) and
#' therefore increases monotonically with increasing \code{level}.
#'
#' @seealso
#' \code{\link{build_ellipsoid}},
#' \code{\link{ellipsoid_surface_points}}
#'
#' @export
ellipsoid_volume <- function(n_dimensions, semi_axes_lengths) {

  # ---- validation ----
  if (missing(n_dimensions))
    stop("Argument 'n_dimensions' must be defined.")

  if (missing(semi_axes_lengths))
    stop("Argument 'semi_axes_lengths' must be defined.")

  if (!is.numeric(n_dimensions) || length(n_dimensions) != 1)
    stop("'n_dimensions' must be a single numeric value.")

  if (!is.numeric(semi_axes_lengths))
    stop("'semi_axes_lengths' must be numeric.")

  if (length(semi_axes_lengths) != n_dimensions)
    stop("'semi_axes_lengths' must have length equal to 'n_dimensions'.")

  if (any(semi_axes_lengths <= 0))
    stop("All semi-axis lengths must be positive.")

  # ---- dimensions-dimensional n volume ----
  unit_ball_volume <- pi^(n_dimensions / 2) /
    gamma(n_dimensions / 2 + 1)

  # ---- ellipsoid volume ----
  volume <- unit_ball_volume * prod(semi_axes_lengths)

  return(volume)
}


#' Update covariances in a nicheR ellipsoid and recompute metrics
#'
#' @param x A `nicheR_ellipsoid` object from `build_ellipsoid()`.
#' @param covariance Either a single numeric (applied to all off-diagonals) or a
#'   named numeric vector with names like `"var1-var2"`.
#' @param tol Small positive number used by covariance limit helpers.
#' @param verbose Logical; print brief progress messages.
#'
#' @return A new `nicheR_ellipsoid` object with updated covariance matrix and
#'   recomputed ellipsoid metrics. Also includes remaining safe limits for any
#'   still-zero covariances.
#' @export
update_ellipsoid_covariance <- function(object,
                                        covariance,
                                        tol = 1e-6,
                                        verbose = TRUE){

  stopifnot(inherits(object, "nicheR_ellipsoid"))

  verbose_message <- function(...) if (isTRUE(verbose)) cat(...)

  verbose_message("Starting: updating covariance values...\n")

  up <- update_covariance(object$cov_matrix, covariance = covariance, tol = tol)

  # Assuming new matrix is symmetrical
  # or can check with: Sigma_new <- (up$updated_matrix + t(up$updated_matrix)) / 2
  Sigma_new <- up$updated_matrix

  verbose_message("Step: recomputing ellipsoid metrics...\n")

  # SPD + inverse
  chol_Sigma <- tryCatch(chol(Sigma_new), error = function(e) NULL)
  if (is.null(chol_Sigma)) stop("Updated covariance matrix is not SPD.")
  Sigma_inv <- chol2inv(chol_Sigma)

  # Cutoff
  chi2_cutoff <- stats::qchisq(object$cl, df = object$dimensions)

  # Eigen + axes
  eig <- eigen(Sigma_new, symmetric = TRUE)
  vals <- pmax(eig$values, 0)
  semi_axes_lengths <- sqrt(vals * chi2_cutoff)

  axis_points <- lapply(seq_len(object$dimensions), function(i) {
    list(
      neg = object$mu_vec - semi_axes_lengths[i] * eig$vectors[, i],
      pos = object$mu_vec + semi_axes_lengths[i] * eig$vectors[, i]
    )
  })

  volume <- ellipsoid_volume(
    n_dimensions = object$dimensions,
    semi_axes_lengths = semi_axes_lengths
  )

  verbose_message("Done: updated ellipsoidal niche boundary. For remianing limits see out$cov_limits_remaining\n")

  # update object
  object$cov_matrix <- Sigma_new
  object$chol_Sigma <- chol_Sigma
  object$Sigma_inv <- Sigma_inv
  object$eigen <- list(vectors = eig$vectors, values = eig$values)
  object$chi2_cutoff <- chi2_cutoff
  object$semi_axes_lengths <- as.numeric(semi_axes_lengths)
  object$axis_points <- axis_points
  object$volume <- volume

  # new fields
  object$cov_limits_remaining <- up$remaining_limits

  object
}




