plot_ellipsoid <- function(object,
                           background = NULL,
                           sample = 1000,
                           lty = 1,
                           lwd = 1,
                           col_ell = "#000000",
                           col_bg = "#8A8A8A",
                           pch = 1,
                           alpha_bg = 0.5,
                           alpha_ell = 0.7,
                           cex_ell = 1,
                           cex_bg = 1, ...){

#   Check for data frame


  ell_points <- ellipsoid_surface_points(mu_vec = object$centroid,
                                         cov_matrix = object$cov_matrix,
                                         chi2_cutoff = object$chi2_cutoff,
                                         n_point = 100) # to do: make sure name of vars does not dissapear

  if(!is.null(background)){

    if(nrow(background) <= sample){
      sample <- nrow(background)
    }

    sample_bg <- sample(1:nrow(background), sample)

    plot(background[sample_bg, ],
         col = adjustcolor(col_bg, alpha.f = alpha_bg),
         pch = pch,
         cex = cex_bg, ...)

    lines(ell_points,
          lty = lty,
          lwd = lwd,
          col = adjustcolor(col_ell, alpha.f = alpha_ell),
          cex = cex_ell)

  }else{
    # Basic line for elliposid
    plot(ell_points, type = "l",
         lty = lty, lwd = lwd,
         col = adjustcolor(col_ell, alpha.f = alpha_ell),
         cex = cex_ell, ...)
  }


}






pred <- predict(ell,
                newdata = bios,
                include_mahalanobis = TRUE,
                include_suitability = TRUE)

# plot_nicheR(list(ell))
plot(pred, col=rev(viridis::viridis(100)))
plot(log(pred), col=rev(viridis::viridis(100)))

occ <- sample_data(n_occ = 50,
                   suitable_env = pred,
                   sampling = "center",
                   method = "probability",
                   seed = 42)
