bios <- rast("inst/extdata/bio_layers.tif")
plot(bios)

# Define environmental ranges
range_df <- data.frame(bio1  = c(22, 28),
                       bio12 = c(1000, 1500),
                       bio15 = c(35, 100))

# Build ellipsoid niche
ell <- build_ellipsoid(range = range_df)
plot_nicheR(list(ell))

# Predict suitability
pred <- predict(ell,
                newdata = bios,
                include_mahalanobis = TRUE,
                include_suitability = TRUE,
                suitability_truncated = TRUE,
                mahalanobis_truncated = TRUE)

plot((pred), col = rev(viridis::viridis(100)))

# import bias rasters
bias_layer <- terra::rast("inst/extdata/bias_layers.tif")
plot(bias_layer)
plot(log(pred), col = rev(viridis::viridis(100)))


bias_result <- prepare_bias(bias_surface = bias_layer,
                          bias_dir = c(-1, 1),
                          pred = pred,   # Use suitability raster as a mask
                          out_bias = "both",
                          truncated = TRUE)

plot(bias_result$pooled_bias)
plot(bias_result$directional_bias)
bias_result$pooled_bias

apply_bias_result <- apply_bias(prepared_bias = bias_result,
                                pred = pred,
                                effect_direction = c(1, -1),
                                layers = c("suitability", "suitability_trunc"))
plot(apply_bias_result,
     main = c("Direct effect", "Indirect effect"))
