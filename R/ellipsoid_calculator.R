#' @export
ellipsoid_calculator <- function(cov_matrix,
                                 centroid,
                                 cl,
                                 verbose = TRUE){

  verbose_message <- function(...) if(isTRUE(verbose)) message(...)

  verbose_message("Step: computing ellipsoid metrics...\n")

  # SPD + inverse
  chol_Sigma <- tryCatch(chol(cov_matrix), error = function(e) NULL)
  if (is.null(chol_Sigma)) stop("Updated covariance matrix is not SPD.")
  Sigma_inv <- chol2inv(chol_Sigma)

  # Cutoff
  chi2_cutoff <- stats::qchisq(cl, df = ncol(cov_matrix))

  # Eigen + axes
  eig <- eigen(cov_matrix, symmetric = TRUE)
  vals <- pmax(eig$values, 0)
  semi_axes_lengths <- sqrt(vals * chi2_cutoff)

  axis_points <- lapply(seq_len(ncol(cov_matrix)), function(i) {
    list(
      start = centroid - semi_axes_lengths[i] * eig$vectors[, i],
      end = centroid + semi_axes_lengths[i] * eig$vectors[, i]
    )
  })

  volume <- ellipsoid_volume(
    n_dimensions = ncol(cov_matrix),
    semi_axes_lengths = semi_axes_lengths
  )

  cov_limits <- covariance_limits(cov_matrix)

  verbose_message("Done: updated ellipsoidal niche boundary. For remianing limits see out$cov_limits_remaining\n")


  out <- list(
    dimensions = ncol(cov_matrix),
    var_names = colnames(cov_matrix),
    centroid = centroid,
    cov_matrix = cov_matrix,
    Sigma_inv = Sigma_inv,
    chol_Sigma = chol_Sigma,
    eigen = list(vectors = eig$vectors, values = eig$values),
    cl = cl,
    chi2_cutoff = chi2_cutoff,
    truncation_level = cl,
    semi_axes_lengths = as.numeric(semi_axes_lengths),
    axis_points = axis_points,
    volume = volume,
    cov_limits = cov_limits
  )

  class(out) <- "nicheR_ellipsoid"

  out
}
