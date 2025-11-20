#' Plot Geographic Space with Optional Suitability/Distance and Occurrences
#'
#' Draws a basemap and overlays (a) suitable environments as tiles,
#' (b) an optional continuous distance-to-centroid surface, and/or (c)
#' occurrence points. The map extent is primarily inferred from `env_bg`
#' (if available), otherwise from the coordinates in `suitable_env` and
#' finally from `occ_pts`, always with a small buffer where needed.
#'
#' When a `NicheR_species` object (`vs`) is supplied, any of
#' `env_bg`, `suitable_env`, `occ_pts`, or `niche` that are `NULL` are filled
#' using `nr_get_env()`, `nr_get_suitable_all()`, `nr_get_occ()`, and
#' `nr_get_niche()` respectively.
#'
#' If `surface` is `NULL`, only the basemap (optionally informed by `env_bg`)
#' and occurrence points are drawn; no suitability or distance surfaces
#' are plotted.
#'
#' @export
plot_g_space <- function(env_bg = NULL,
                         suitable_env = NULL,
                         occ_pts = NULL,
                         show.occ.density = FALSE,
                         colors = NULL,
                         palette = "default",
                         surface = NULL,
                         vs = NULL,
                         niche = NULL) {

  ## ------------------------------------------------------------------------
  ## 0. Surface mode
  ## ------------------------------------------------------------------------
  surface_mode <- if (is.null(surface)) {
    "none"     # no suitability/distance, just background + occurrences
  } else {
    match.arg(surface, c("both", "suit", "dist"))
  }

  ## ------------------------------------------------------------------------
  ## 1. Pull components from NicheR_species (using nr_get_*)
  ## ------------------------------------------------------------------------
  if (!is.null(vs)) {
    if (!inherits(vs, "NicheR_species")) {
      stop("'vs' must be an object of class 'NicheR_species' created by create_virtual_species().")
    }

    # env_bg: prefer user-supplied; otherwise from vs via nr_get_env()
    if (is.null(env_bg)) {
      env_from_vs <- nr_get_env(vs)
      if (!is.null(env_from_vs)) {
        env_bg <- env_from_vs
        message("Using env_bg retrieved from 'vs' via nr_get_env().")
      }
    }

    # niche (not strictly needed anymore here, but kept for possible future use)
    if (is.null(niche)) {
      niche_from_vs <- nr_get_niche(vs)
      if (!is.null(niche_from_vs)) {
        niche <- niche_from_vs
      }
    }

    # occurrences
    if (is.null(occ_pts)) {
      occ_from_vs <- nr_get_occ(vs)
      if (!is.null(occ_from_vs)) {
        occ_pts <- occ_from_vs
      }
    }

    # suitable_env: let nr_get_* work on vs directly
    if (is.null(suitable_env) && surface_mode != "none") {
      suitable_env <- nr_get_suitable_all(vs)
      if (!is.null(suitable_env)) {
        message("Using 'vs' as suitable_env source via nr_get_*().")
      }
    }
  }

  ## ------------------------------------------------------------------------
  ## 2. Palette / color handling
  ## ------------------------------------------------------------------------
  palettes <- list(
    default = list(
      bg           = "#F0F0F0FF",
      suitable_env = "#FED789FF",
      occ_fill     = "#B4BF3AFF",
      occ_stroke   = "#B4BF3AFF",
      dist         = "YlOrRd"
    ),
    palette2 = list(
      bg           = "#E0ECF4FF",
      suitable_env = "#9ECAE1FF",
      occ_fill     = "#08519CFF",
      occ_stroke   = "#08306BFF",
      dist         = "YlGnBu"
    ),
    palette3 = list(
      bg           = "#F0F0F0FF",
      suitable_env = "#BDBDBDFF",
      occ_fill     = "#252525FF",
      occ_stroke   = "#000000FF",
      dist         = "OrRd"
    )
  )

  if (!palette %in% names(palettes)) {
    stop("Unknown palette '", palette, "'. Use one of: ",
         paste(names(palettes), collapse = ", "), ".")
  }

  base_colors <- palettes[[palette]]

  if (is.null(colors)) {
    default_colors <- base_colors
  } else {
    if (is.null(names(colors))) {
      names(colors) <- names(base_colors)[seq_along(colors)]
      if (!is.list(colors)) colors <- as.list(colors)
      message(
        "No names detected in 'colors'. Using the provided order and filling ",
        "missing entries with defaults.\n",
        "Named options are: ",
        paste(names(base_colors), collapse = ", "), "."
      )
    }
    if (!is.list(colors)) colors <- as.list(colors)
    default_colors <- utils::modifyList(base_colors, colors)
  }

  ## ------------------------------------------------------------------------
  ## 3. Occurrence points validation
  ## ------------------------------------------------------------------------
  if (!is.null(occ_pts) && !all(c("x", "y") %in% names(occ_pts))) {
    stop("`occ_pts` must have columns named 'x' and 'y'.")
  }

  ## ------------------------------------------------------------------------
  ## 4. Standardize suitable_env into a data.frame using nr_get_*
  ## ------------------------------------------------------------------------
  suitable_df   <- NULL
  has_suitable  <- FALSE
  has_distance  <- FALSE

  if (surface_mode != "none") {

    if (is.null(suitable_env)) {
      stop(
        "A suitability/distance surface was requested (surface = '", surface_mode, "'),\n",
        "but no `suitable_env` was provided and none could be resolved from `vs`.\n",
        "Either:\n",
        "  * supply suitable_env (e.g., output from get_suitable_env()), or\n",
        "  * call plot_g_space(env_bg = ..., occ_pts = ..., surface = NULL)\n",
        "    to plot only background and occurrence points."
      )
    }

    # 4A. Try to get a df directly via nr_get()
    df_try <- suppressWarnings(nr_get_suitable_df(suitable_env))

    if (!is.null(df_try)) {

      suitable_df <- as.data.frame(df_try)

    } else {

      # 4B. Otherwise, use rasters via nr_get_suitable / nr_get_dist_sq
      r_suit <- NULL
      r_dist <- NULL

      if (surface_mode %in% c("both", "suit")) {
        r_suit <- suppressWarnings(nr_get_suitable(suitable_env))
      }
      if (surface_mode %in% c("both", "dist")) {
        r_dist <- suppressWarnings(nr_get_dist_sq(suitable_env))
      }

      # Coerce any Raster* to SpatRaster
      if (inherits(r_suit, "Raster")) r_suit <- terra::rast(r_suit)
      if (inherits(r_dist, "Raster")) r_dist <- terra::rast(r_dist)

      if (surface_mode == "suit") {

        if (is.null(r_suit)) {
          stop(
            "surface = 'suit' requested, but no suitable raster could be found via nr_get_suitable().\n",
            "Check that your `suitable_env` (or `vs`) contains a 'suitable' layer or df."
          )
        }
        names(r_suit) <- "suitable"
        suitable_df <- as.data.frame.nicheR(r_suit)

      } else if (surface_mode == "dist") {

        if (is.null(r_dist)) {
          stop(
            "surface = 'dist' requested, but no distance raster could be found via nr_get_dist_sq().\n",
            "Run get_suitable_env(..., distances = TRUE) or provide an object\n",
            "from which nr_get_dist_sq() can extract a 'dist_sq' layer."
          )
        }
        names(r_dist) <- "dist_sq"
        suitable_df <- as.data.frame.nicheR(r_dist)

      } else if (surface_mode == "both") {

        if (is.null(r_suit) || is.null(r_dist)) {
          stop(
            "surface = 'both' requested, but could not find both suitability and distance rasters.\n",
            "Ensure your object contains both 'suitable' and 'dist_sq' layers or a df with those columns."
          )
        }

        names(r_suit) <- "suitable"
        names(r_dist) <- "dist_sq"
        ras_stack <- c(r_suit, r_dist)
        suitable_df <- as.data.frame.nicheR(ras_stack)
      }
    }

    if (!is.null(suitable_df)) {
      if (!all(c("x", "y") %in% names(suitable_df))) {
        stop("`suitable_env` (after coercion) must contain 'x' and 'y' columns.")
      }
      has_suitable <- "suitable" %in% names(suitable_df) || surface_mode != "dist"
      has_distance <- "dist_sq" %in% names(suitable_df) && surface_mode != "suit"
    }
  }

  ## ------------------------------------------------------------------------
  ## 5. Legend configuration (with distance note toggle)
  ## ------------------------------------------------------------------------
  opts <- list(
    background_point = TRUE,
    suitable_point   = has_suitable && surface_mode %in% c("both", "suit"),
    occurrence_point = !is.null(occ_pts),
    dist_note        = has_distance && surface_mode %in% c("both", "dist")
  )

  legend_items <- data.frame(
    id    = c("background_point", "suitable_point", "occurrence_point", "dist_note"),
    type  = c("point", "point", "point", "text"),
    label = c("Background environments",
              "Suitable environments",
              "Occurrence",
              "Distance surface: light = closer to centroid"),
    color = c(default_colors[["bg"]],
              default_colors[["suitable_env"]],
              default_colors[["occ_fill"]],
              NA_character_),
    stringsAsFactors = FALSE
  )

  active <- logical(nrow(legend_items))
  for (i in seq_len(nrow(legend_items))) {
    active[i] <- isTRUE(opts[[ legend_items$id[i] ]])
  }
  legend_items <- legend_items[active, , drop = FALSE]

  if (nrow(legend_items) == 0) {
    legend_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  } else {
    top_y   <- 2
    spacing <- 0.25
    x_point <- 0.00

    legend_items <- legend_items %>%
      dplyr::mutate(
        row     = dplyr::row_number(),
        y       = top_y - (row - 1) * spacing,
        x_point = x_point,
        x_text  = 0.01
      )

    legend_base <- ggplot2::ggplot() +
      ggplot2::coord_cartesian(xlim = c(0, 0.2), ylim = c(0, 2.5), clip = "off") +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "none")

    legend_plot <- legend_base +
      ggplot2::geom_point(
        data  = dplyr::filter(legend_items, type == "point"),
        ggplot2::aes(x = x_point, y = y, colour = color),
        size  = 2
      ) +
      ggplot2::geom_text(
        data  = legend_items,
        ggplot2::aes(x = x_text, y = y, label = label),
        hjust = 0
      ) +
      ggplot2::scale_colour_identity()
  }

  ## ------------------------------------------------------------------------
  ## 6. Compute map extent: env_bg -> suitable_env -> occ_pts
  ## ------------------------------------------------------------------------
  x_lim <- y_lim <- NULL
  extent_source <- NULL

  # 6A. Try env_bg first
  if (!is.null(env_bg)) {
    env_bg_r <- NULL

    if (inherits(env_bg, "Raster")) {
      env_bg_r <- terra::rast(env_bg)
    } else if (inherits(env_bg, "SpatRaster")) {
      env_bg_r <- env_bg
    }

    if (!is.null(env_bg_r)) {
      ex <- terra::ext(env_bg_r)
      if (all(is.finite(as.vector(ex)))) {
        x_lim <- c(ex[1], ex[2])
        y_lim <- c(ex[3], ex[4])
        extent_source <- "env_bg"
      }
    } else if (is.data.frame(env_bg) || is.matrix(env_bg)) {
      env_bg_df <- as.data.frame(env_bg)
      if (all(c("x", "y") %in% names(env_bg_df))) {
        coords <- env_bg_df[, c("x", "y")]
        coords <- coords[stats::complete.cases(coords), , drop = FALSE]
        if (nrow(coords) > 0) {
          x_rng <- range(coords$x, na.rm = TRUE)
          y_rng <- range(coords$y, na.rm = TRUE)
          dx    <- x_rng[2] - x_rng[1]
          dy    <- y_rng[2] - y_rng[1]
          if (!is.finite(dx) || dx <= 0) dx <- 1
          if (!is.finite(dy) || dy <= 0) dy <- 1
          pad_x <- 0.05 * dx
          pad_y <- 0.05 * dy
          x_lim <- c(x_rng[1] - pad_x, x_rng[2] + pad_x)
          y_lim <- c(y_rng[1] - pad_y, y_rng[2] + pad_y)
          extent_source <- "env_bg_xy"
        }
      }
    }
  }

  # 6B. If env_bg did not give an extent, use suitable_env / occ_pts coords
  if (is.null(x_lim) || is.null(y_lim)) {
    coord_sources <- list()
    if (!is.null(suitable_df) && nrow(suitable_df) > 0) {
      coord_sources[[length(coord_sources) + 1]] <- suitable_df[, c("x", "y")]
    }
    if (!is.null(occ_pts) && all(c("x", "y") %in% names(occ_pts))) {
      coord_sources[[length(coord_sources) + 1]] <- occ_pts[, c("x", "y")]
    }

    if (length(coord_sources) > 0) {
      coords <- do.call(rbind, coord_sources)
      coords <- coords[stats::complete.cases(coords), , drop = FALSE]

      if (nrow(coords) > 0) {
        x_rng <- range(coords$x, na.rm = TRUE)
        y_rng <- range(coords$y, na.rm = TRUE)
        dx    <- x_rng[2] - x_rng[1]
        dy    <- y_rng[2] - y_rng[1]
        if (!is.finite(dx) || dx <= 0) dx <- 1
        if (!is.finite(dy) || dy <= 0) dy <- 1
        pad_x <- 0.1 * dx
        pad_y <- 0.1 * dy
        x_lim <- c(x_rng[1] - pad_x, x_rng[2] + pad_x)
        y_lim <- c(y_rng[1] - pad_y, y_rng[2] + pad_y)
        extent_source <- "suitable/occ"
      }
    }
  }

  if (!is.null(extent_source)) {
    msg_src <- switch(
      extent_source,
      "env_bg"        = "Extent derived from env_bg raster.",
      "env_bg_xy"     = "Extent derived from env_bg x/y coordinates.",
      "suitable/occ"  = "Extent derived from suitable_env / occurrence coordinates."
    )
    message(msg_src)
  } else {
    message("Could not infer extent from env_bg, suitable_env, or occ_pts; using global extent.")
  }

  ## ------------------------------------------------------------------------
  ## 7. Basemap
  ## ------------------------------------------------------------------------
  world <- ggplot2::map_data("world")

  base_map <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = world,
      ggplot2::aes(x = long, y = lat, group = group),
      fill = default_colors[["bg"]]
    ) +
    ggplot2::xlab("Longitude") +
    ggplot2::ylab("Latitude") +
    ggplot2::theme_bw()

  if (!is.null(x_lim) && !is.null(y_lim)) {
    base_map <- base_map +
      ggplot2::coord_quickmap(xlim = x_lim, ylim = y_lim, expand = FALSE)
  } else {
    base_map <- base_map +
      ggplot2::coord_quickmap()
  }

  ## ------------------------------------------------------------------------
  ## 8. Occurrence points (sf)
  ## ------------------------------------------------------------------------
  occ_pts_sp <- NULL
  if (!is.null(occ_pts)) {
    occ_pts_sp <- sf::st_as_sf(
      occ_pts[, c("x", "y")],
      coords = c("x", "y"),
      crs    = 4326
    )
  }

  add_occ <- function(p) {
    if (!is.null(occ_pts_sp)) {
      p <- p +
        ggplot2::geom_sf(
          data  = occ_pts_sp,
          ggplot2::aes(geometry = geometry),
          color = default_colors[["occ_stroke"]],
          fill  = default_colors[["occ_fill"]],
          pch   = 21,
          size  = 0.75
        )
    }
    p
  }

  ## ------------------------------------------------------------------------
  ## 9. Build one or two map panels
  ## ------------------------------------------------------------------------
  map_list <- list()

  # Panel 1: binary suitable (presence-only tiles)
  if (has_suitable && surface_mode %in% c("both", "suit")) {
    p_suit <- base_map +
      ggplot2::geom_tile(
        data = suitable_df,
        ggplot2::aes(x = x, y = y),
        fill = default_colors[["suitable_env"]]
      ) +
      ggplot2::ggtitle("Suitable Environments")

    p_suit <- add_occ(p_suit)
    map_list[["suitable"]] <- p_suit
  }

  # Panel 2: distance surface
  if (has_distance && surface_mode %in% c("both", "dist")) {
    p_dist <- base_map +
      ggplot2::geom_tile(
        data = suitable_df,
        ggplot2::aes(x = x, y = y, fill = dist_sq)
      ) +
      ggplot2::scale_fill_distiller(
        name    = "Distance to centroid",
        palette = default_colors[["dist"]],
        guide   = "none"  # explained via legend note
      ) +
      ggplot2::ggtitle("Distance surface")

    p_dist <- add_occ(p_dist)
    map_list[["dist_sq"]] <- p_dist
  }

  # If surface_mode == "none" or no suitability at all: just basemap + occs
  if (surface_mode == "none" || length(map_list) == 0L) {
    p_base <- add_occ(base_map)
    return(ggpubr::ggarrange(p_base, legend_plot, widths = c(0.7, 0.3)))
  }

  ## ------------------------------------------------------------------------
  ## 10. Arrange panels + legend
  ## ------------------------------------------------------------------------
  if (length(map_list) == 1L) {
    main_panel <- map_list[[1]]

    return(ggpubr::ggarrange(
      main_panel,
      legend_plot,
      ncol    = 1,
      heights = c(0.7, 0.3)
    ))

  } else {

    main_panel <- ggpubr::ggarrange(
      plotlist = map_list,
      nrow     = length(map_list),
      labels   = NULL)

    return(ggpubr::ggarrange(
      main_panel,
      legend_plot,
      ncol    = 1,
      heights = c(0.8, 0.2)
    ))
  }
}
