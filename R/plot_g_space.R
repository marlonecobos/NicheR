#' Plot Geographic Space with Optional Suitability/Distance and Occurrences
#'
#' Draws a world basemap and overlays (a) suitable environments as tiles,
#' (b) an optional continuous distance-to-centroid surface, and/or (c)
#' occurrence points. A compact custom legend is built to match the
#' discrete layers shown.
#'
#' @param env_bg Optional background environment object used when
#'   `suitable_env` is not supplied and suitability must be computed.
#'   Typically a `terra::SpatRaster` (or `raster::Raster*`).
#' @param n_bg Integer (currently unused; reserved for future downsampling of
#'   the background grid).
#' @param niche Optional object of class `"ellipsoid"` (from [build_ellps()]).
#'   Required when suitability is computed internally (i.e., when
#'   `suitable_env` is `NULL`).
#' @param suitable_env Optional precomputed suitable environment object. Can be:
#'   \itemize{
#'     \item a `data.frame` with columns `x`, `y`, and optionally `dist_sq`,
#'     \item a `"suitable_env"` object returned by [get_suitable_env()]
#'           with a `suitable_env_df` component.
#'   }
#'   If `dist_sq` is present, a distance panel is also drawn. If you only
#'   want one surface (e.g., just suitability or just distance), supply only
#'   that information in `suitable_env`.
#' @param occ_pts Optional `data.frame` of occurrences with columns `x`, `y`
#'   (assumed longitude/latitude, WGS84). Plotted as points if supplied. If
#'   `vs` is supplied and `occ_pts` is `NULL`, the function will try to use
#'   `vs$occurrences`.
#' @param show.occ.density Logical (currently unused; reserved for future
#'   density panels).
#' @param colors Optional named list to override aesthetics. Recognized names:
#'   `bg` (basemap fill), `suitable_env` (suitable tiles), `occ_fill`,
#'   `occ_stroke`, and `dist` (RColorBrewer palette name for distance tiles).
#' @param palette Character palette key. One of `"default"`, `"palette2"`,
#'   or `"palette3"`.
#' @param vs Optional object of class `"NicheR_species"` created by
#'   [create_virtual_species()]. When provided, any of `env_bg`, `niche`,
#'   `occ_pts`, or `suitable_env` that are `NULL` will be filled from this
#'   object if possible.
#'
#' @return A `ggpubr` object produced by [ggpubr::ggarrange()], containing
#'   one or two map panels plus a matching legend panel.
#'
#' @family plotting functions
#' @seealso [build_ellps()], [get_suitable_env()], [create_virtual_species()]
#'
#' @export
plot_g_space <- function(env_bg = NULL,
                         n_bg = 10000,
                         suitable_env = NULL,
                         occ_pts = NULL,
                         show.occ.density = FALSE,
                         colors = NULL,
                         palette = "default",
                         vs = NULL,
                         niche = NULL) {

  # ---- 0. Pull components from NicheR_species if provided -----------------
  if (!is.null(vs)) {
    if (!inherits(vs, "NicheR_species")) {
      stop("'vs' must be an object of class 'NicheR_species' created by create_virtual_species().")
    }

    if (is.null(env_bg) && !is.null(vs$call_args) && "env_bg" %in% names(vs$call_args)) {
      env_bg <- vs$call_args$env_bg
    }
    if (is.null(niche) && !is.null(vs$niche)) {
      niche <- vs$niche
    }
    if (is.null(occ_pts) && !is.null(vs$occurrences)) {
      occ_pts <- vs$occurrences
    }
    if (is.null(suitable_env) && !is.null(vs$suitability)) {
      suitable_env <- vs$suitability
    }
  }

  # ---- 1. Palette / color handling ----------------------------------------
  palettes <- list(
    default = list(
      bg           = "#FED789FF",
      suitable_env = "#B4BF3AFF",
      occ_fill     = "black",
      occ_stroke   = "black",
      dist         = "YlOrRd"   # brewer palette name
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

  # ---- 2. Occurrence points validation ------------------------------------

  if (!is.null(occ_pts) && !all(c("x", "y") %in% names(occ_pts))) {
    stop("`occ_pts` must have columns named 'x' and 'y'.")
  }
  # ---- 3. Standardize suitable_env into a data.frame ----------------------

  suitable_df <- NULL

  if (!is.null(suitable_env)) {

    if (is.data.frame(suitable_env)) {

      # Case 1: plain data.frame
      suitable_df <- suitable_env

    } else if (is.list(suitable_env) &&
               "suitable_env_df" %in% names(suitable_env)) {

      # Case 2: full suitable_env object with a df component
      suitable_df <- suitable_env$suitable_env_df

    } else {

      # Cases where suitability is stored as rasters
      rast_stack <- NULL

      # 2a) suitable_env is a full suitable_env object with only rasters
      if (is.list(suitable_env) &&
          "suitable" %in% names(suitable_env)) {

        sp <- suitable_env

        if (inherits(sp, "SpatRaster")) {
          rast_stack <- sp
        } else if (is.list(sp) &&
                   length(sp) > 0 &&
                   all(vapply(sp, inherits, logical(1), "SpatRaster"))) {

          # list of rasters (e.g. suitable, dist_sq) → stack
          rast_stack <- rast(sp)

        } else {
          stop(
            "'suitable_env$suitable_env_sp' must be a SpatRaster or a list of SpatRasters ",
            "when no 'suitable_env_df' is present."
          )
        }

      } else if (inherits(suitable_env, "SpatRaster")) {

        # 2b) single SpatRaster stack
        rast_stack <- suitable_env

      } else if (is.list(suitable_env) &&
                 length(suitable_env) > 0 &&
                 all(vapply(suitable_env, inherits, logical(1), "SpatRaster"))) {

        # 2c) list of SpatRasters directly
        rast_stack <- rast(sp)

      } else {
        stop(
          "`suitable_env` must be either:\n",
          "  * a data.frame with columns 'x', 'y' (and optionally 'dist_sq'), or\n",
          "  * a 'suitable_env' object with 'suitable_env_df' or 'suitable_env_sp', or\n",
          "  * a SpatRaster / list of SpatRasters.\n"
        )
      }

      # convert stack to data.frame
      suitable_df <- as.data.frame.nicheR(rast_stack)
    }

    # sanity check for coordinates
    if (!all(c("x", "y") %in% names(suitable_df))) {
      stop("`suitable_env` (after coercion) must contain 'x' and 'y' columns.")
    }

  } else {
    # no suitable_env provided; fall back to internal computation if possible
    if (!is.null(niche) && !is.null(env_bg)) {
      suitable_df <- get_suitable_env(
        niche     = niche,
        env_bg    = env_bg,
        out.suit  = "data.frame",
        distances = TRUE,
        verbose   = FALSE
      )
      message(
        "Suitability was not supplied; computed internally via get_suitable_env().\n",
        "For repeated plots, consider calling get_suitable_env() once and passing ",
        "the result to `suitable_env`."
      )
    } else {
      suitable_df <- NULL
    }
  }

  has_suitable <- !is.null(suitable_df)
  has_distance <- has_suitable && ("dist_sq" %in% names(suitable_df))

  # ---- 4. Legend configuration --------------------------------------------

  opts <- list(
    background_point = TRUE,
    suitable_point   = has_suitable,
    occurrence_point = !is.null(occ_pts)
  )

  legend_items <- data.frame(
    id    = c("background_point","suitable_point","occurrence_point"),
    type  = c("point","point","point"),
    label = c("Background environments","Suitable environments","Occurrence"),
    color = c(
      default_colors[["bg"]],
      default_colors[["suitable_env"]],
      default_colors[["occ_fill"]]
    ),
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

  # ---- 5. Basemap ----------------------------------------------------------

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

  # ---- 6. Occurrence points (convert to sf) --------------------------------

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

  # ---- 7. Build one or two map panels -------------------------------------

  map_list <- list()

  if (has_suitable) {
    # Panel 1: binary suitable (presence-only tiles)
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

  if (has_distance) {
    # Panel 2: distance surface
    p_dist <- base_map +
      ggplot2::geom_tile(
        data = suitable_df,
        ggplot2::aes(x = x, y = y, fill = dist_sq)
      ) +
      ggplot2::scale_fill_distiller(
        name    = "Distance to centroid",
        palette = default_colors[["dist"]]
      ) +
      ggplot2::ggtitle("Distance surface") +
      ggplot2::theme(legend.position = "bottom")

    p_dist <- add_occ(p_dist)
    map_list[["dist_sq"]] <- p_dist
  }

  # If no suitability at all, just show basemap + occurrences
  if (length(map_list) == 0L) {
    p_base <- add_occ(base_map)
    return(ggpubr::ggarrange(p_base, legend_plot, widths = c(0.7, 0.3)))
  }

  # ---- 8. Arrange panels + legend -----------------------------------------

  if (length(map_list) == 1L) {
    main_panel <- map_list[[1]]
  } else {
    main_panel <- ggpubr::ggarrange(
      plotlist = map_list,
      nrow     = length(map_list),
      labels   = NULL,
      heights = c(0.45, 0.55)
    )
  }

  ggpubr::ggarrange(
    main_panel,
    legend_plot,
    ncol    = 1,
    heights  = c(0.7, 0.3)
  )

}
