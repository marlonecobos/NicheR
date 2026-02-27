#' Prepare Sampling Bias Surfaces
#'
#' Standardizes, optionally inverts, and combines one or more sampling bias
#' rasters into a composite bias surface aligned to a common template.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Splits \code{bias_surface} into single-layer rasters.
#'   \item Determines a template raster:
#'   \itemize{
#'     \item If \code{template_surface} is provided, each bias layer is cropped
#'     and masked to its extent and later resampled to its grid.
#'     \item If \code{template_surface} is \code{NULL}, the first bias layer is
#'     used as the template grid.
#'   }
#'   \item Resamples bias layers to match the template grid when needed
#'   (\code{terra::compareGeom()} + \code{terra::resample(method = "near")}).
#'   \item Standardizes each layer to \eqn{[0, 1]} using min-max scaling.
#'   \item Applies the direction of effect:
#'   \itemize{
#'     \item \code{"direct"} uses standardized values as-is.
#'     \item \code{"inverse"} inverts standardized values via \eqn{1 - x}.
#'   }
#'   \item Optionally builds a composite bias surface:
#'   \itemize{
#'     \item Where multiple layers overlap, values are combined multiplicatively.
#'     \item Where only one layer is defined, that layer's value is retained.
#'     \item Where all layers are \code{NA}, the composite is \code{NA}.
#'   }
#' }
#'
#' @param bias_surface A \code{terra::SpatRaster} or a list of
#'   \code{terra::SpatRaster} objects. Multi-layer rasters are split internally
#'   into single-layer bias surfaces.
#' @param effect_direction Character vector specifying how each bias layer
#'   contributes to sampling probability. Options are \code{"direct"} or
#'   \code{"inverse"}. A single value is recycled to match the number of bias
#'   layers.
#' @param template_surface Optional. A \code{terra::SpatRaster} used as a spatial
#'   template. If provided, bias layers are cropped and masked to its extent and
#'   resampled to its grid. If \code{NULL}, the first bias layer is used as the
#'   template grid.
#' @param include_composite Logical. If \code{TRUE} (default), returns a
#'   composite bias surface.
#' @param include_processed_layers Logical. If \code{TRUE}, also returns the
#'   processed (standardized and directional) bias layers as a multi-layer
#'   \code{SpatRaster}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages.
#'
#' @return
#' A list of class \code{"nicheR_bias_surface"} with:
#' \itemize{
#'   \item \code{composite_surface} A single-layer composite bias surface
#'   (only if \code{include_composite = TRUE}).
#'   \item \code{processed_layers} A multi-layer raster of standardized and
#'   directional bias layers (returned if \code{include_processed_layers = TRUE},
#'   or if \code{include_composite = FALSE}).
#'   \item \code{combination_formula} Character string describing how layers
#'   were combined (e.g., \code{"roads * (1-pop_density) * protected"}).
#' }
#'
#' @seealso \code{\link{apply_bias}}
#'
#' @export
prepare_bias <- function(bias_surface,
                         effect_direction = c("direct", "inverse"),
                         template_layer = NULL,
                         include_composite = TRUE,
                         include_processed_layers = FALSE,
                         mask_na = TRUE,
                         verbose = TRUE){

  verbose_message <- function(...) if(isTRUE(verbose)) cat(...)

  gc()


  verbose_message("Starting: prepare_bias()\n")

  # Basic Input checks --------------------------------------------------------

  if(missing(bias_surface) || is.null(bias_surface)){
    stop("'bias_surface' must be a SpatRaster or list of SpatRasters.")
  }

  # Normalize bias_surface → list of single-layer rasters
  if(inherits(bias_surface, "SpatRaster")){
    verbose_message("Step: splitting SpatRaster into layers...\n")
    bias_list <- terra::as.list(bias_surface)

  }else if(is.list(bias_surface) &&
           all(vapply(bias_surface, inherits, logical(1), "SpatRaster"))){
    verbose_message("Step: flattening list of SpatRasters...\n")
    bias_list <- unlist(lapply(bias_surface, terra::as.list), recursive = FALSE)
  }else{
    stop("'bias_surface' must be SpatRaster or list of SpatRasters.")
  }

  if(length(bias_list) == 0){
    stop("No usable bias layers provided.")
  }

  # 1. Determine template raster ------------------------------------------

  if(!is.null(template_surface) && inherits(template_surface, "SpatRaster")){
    verbose_message("Step: using user provided template surface layer to crop and mask to extent, will also be use to resample if necessary...\n")

  }else{
    verbose_message("Step: Using first layer in bias_surface as template surface layer to resample if necessary...\n")
    template_surface <- bias_list[[1]]
  }

  # 2. Prepare effect direction ------------------------------------------

  effect_direction <- match.arg(effect_direction,
                                choices = c("direct", "inverse"),
                                several.ok = TRUE)

  if(length(effect_direction) == 1L){
    effect_direction <- rep(effect_direction, length(bias_list))
  }

  if(length(effect_direction) != length(bias_list)){
    stop("'effect_direction' must be length 1 or the same length as the number of bias layers.")
  }

  if(!all(effect_direction %in% c("direct", "inverse"))){
    stop("'effect_direction' must contain only 'direct' or 'inverse'.")
  }

  # Output logic checks ------------------------------------------------------

  if(!is.logical(include_composite) || length(include_composite) != 1L){
    stop("'include_composite' must be TRUE or FALSE.")
  }

  if(!is.logical(include_processed_layers) || length(include_processed_layers) != 1L){
    stop("'include_processed_layers' must be TRUE or FALSE.")
  }

  if(!include_composite && !include_processed_layers){
    include_composite <- TRUE

    if(isTRUE(verbose)){
      verbose_message("Step: both 'include_composite' and 'include_processed_layers' were FALSE. ", "Defaulting to include_composite = TRUE...\n")
    }
  }

  # 3. Process layers ------------------------------------------------------

  verbose_message("Step: standarizing (min/max) and applying direction of effect to ", length(bias_list), " bias layer/s...\n")

  directional_bias_list <- vector("list", length(bias_list))
  formula_entries <- character(length(bias_list))


  for(i in seq_along(bias_list)){

    raw <- bias_list[[i]]
    this_dir <- effect_direction[i]

    nm <- names(raw)
    nm <- if(is.null(nm) || nm == "" || nm == "lyr.1") paste0("bias_", i) else nm[1]

    # Resample if needed
    same_grid <- terra::compareGeom(raw,
                                    template_surface,
                                    stopOnError = FALSE)

    if(!isTRUE(same_grid)){
      verbose_message("Step: resampling bias layer ", i, " to match template surface...\n")
      aligned <- terra::resample(raw,
                                 template_surface,
                                 method = "near")
    }else{
      aligned <- raw
    }

    aligned <- terra::crop(aligned,
                           template_surface, mask = TRUE)

    # Standardize
    mm <- terra::minmax(aligned)
    min_val <- mm[1]
    max_val <- mm[2]

    if(!is.finite(min_val) || !is.finite(max_val)){
      scaled <- aligned
      scaled[!is.na(scaled)] <- 1
    }else if((max_val - min_val) == 0){
      scaled <- aligned
      scaled[!is.na(scaled)] <- 1
    }else{
      scaled <- (aligned - min_val)/(max_val - min_val)
    }

    # Apply direction
    resampled_flag <- !isTRUE(same_grid)

    resample_tag <- if(resampled_flag) "_resampled" else ""

    if(this_dir == "inverse"){
      formula_entries[i] <- paste0("(1-", nm, ")")
      directional_bias_list[[i]] <- 1 - scaled
      names(directional_bias_list[[i]]) <- paste0("standarized_", nm, "_inverse", resample_tag)
    }else{
      formula_entries[i] <- nm
      directional_bias_list[[i]] <- scaled
      names(directional_bias_list[[i]]) <- paste0("standarized_", nm, "_direct", resample_tag)
    }
  }

  # Stack layers ----------------------------------------------------------

  directional_bias_stack <- if(length(directional_bias_list) == 1L){
    directional_bias_list[[1]]
  }else{
    do.call(c, directional_bias_list)
  }

  # 4. Combine layers ----------------------------------------------------

  out_rast <- list()

  if(isTRUE(include_composite)){

    verbose_message("Step: building standarized (min/max) directional composite bias surface...\n")

    if(terra::nlyr(directional_bias_stack) > 1){
      # Bias layers are combined multiplicatively where overlapping; in areas
      # where only one layer is defined, that layer’s value is retained.
      composite_raster <- terra::app(directional_bias_stack,
                                     fun = function(x){
                                       if(all(is.na(x))){
                                         NA_real_
                                       }else if(sum(!is.na(x)) == 1){
                                         x[!is.na(x)]
                                       }else{
                                         prod(x)
                                       }
                                     })

    }else{
      composite_raster <- directional_bias_stack
    }

    names(composite_raster) <- "standarized_composite_bias_surface"
    out_rast$composite_surface <- composite_raster

    if(isTRUE(include_processed_layers)){
      out_rast$processed_layers <- directional_bias_stack
    }

    out_rast$combination_formula <- paste(formula_entries, collapse = " * ")

  }else{
    out_rast$processed_layers <- directional_bias_stack
  }


  # 5. Build result --------------------------------------------------------

  verbose_message("Done: prepare_bias()\n")

  gc()

  out_rast
}
