#' Plot Environmental Space with Optional Ellipsoid Overlays
#'
#' Produces pairwise views of a 3D environmental space with optional overlays of a
#' virtual niche ellipsoid and occurrence points. In 2D mode, it returns a grid
#' of pairwise scatterplots with projected ellipse boundaries. In 3D mode, it
#' returns an interactive \code{plotly} scatterplot.
#'
#' @param env_bg A \code{data.frame} of background environments with at least three
#'   numeric predictor columns. These columns must contain the variables
#'   referenced by \code{x}, \code{y}, and \code{z}. If \code{NULL}, a
#'   \code{NicheR_species} object supplied via \code{vs} can be used to infer
#'   \code{env_bg} from its \code{suitability} slot (either a data.frame or a
#'   SpatRaster that will be converted).
#' @param x,y,z Column specifications for the three predictors to display. Each may
#'   be a single column name (character string) or a single 1-based integer index
#'   into \code{env_bg}.
#' @param labels Character vector of length 3 giving axis labels for the x, y,
#'   and z variables in display order. Defaults to \code{c("ENV 1", "ENV 2", "ENV 3")}.
#' @param n_bg Positive integer giving the maximum number of background rows to plot.
#'   If \code{nrow(env_bg)} is greater than \code{n_bg}, a random subset of size
#'   \code{n_bg} is drawn. Using a large \code{n_bg} may slow plotting.
#' @param niche Optional object of class \code{ellipsoid} describing the niche. If
#'   provided, its boundary and center will be plotted. For 2D plots, the
#'   object must contain \code{niche$angles}. If \code{NULL} and a
#'   \code{NicheR_species} object is supplied via \code{vs}, \code{niche} is
#'   filled from \code{vs$niche}.
#' @param suitable_env Optional suitable environment object (data.frame,
#'   \code{suitable_env} list, or SpatRaster), typically produced by
#'   \code{\link{get_suitable_env}}. When provided, points inside the ellipsoid
#'   are plotted in E space.
#' @param occ_pts Optional \code{data.frame} of occurrence points that includes the
#'   same predictor columns used for \code{x}, \code{y}, and \code{z}. These are
#'   overplotted if supplied. If \code{NULL} and a \code{NicheR_species} object
#'   is supplied via \code{vs}, this is filled from \code{vs$occurrences} when available.
#' @param rand_seed Integer used to set the random number generator seed for
#'   reproducible background downsampling.
#' @param show.occ.density Logical. If \code{TRUE} and \code{occ_pts} is provided,
#'   adds marginal density panels for each variable (2D plots only;
#'   \code{plot.3d = FALSE}).
#' @param plot.3d Logical. If \code{TRUE}, returns an interactive \code{plotly} 3D scatter
#'   plot. If \code{FALSE} (the default), returns a static \code{ggpubr} grid of 2D panels.
#' @param colors Optional named list of colors to override the default palette.
#'   Valid names are: \code{bg}, \code{ellipsoid}, \code{centroid}, \code{tolerance},
#'   \code{suitable_env}, \code{occ}.
#' @param palette Character name of the internal palette to use.
#'   One of \code{"default"}, \code{"palette2"}, ..., \code{"palette6"}.
#' @param vs Optional \code{NicheR_species} object returned by
#'   \code{\link{create_virtual_species}}. If provided, \code{niche},
#'   \code{occ_pts}, \code{env_bg}, and \code{suitable_env} are auto-filled from
#'   this object when they are not supplied explicitly.
#'
#' @return
#' If \code{plot.3d = TRUE}, a \code{plotly} object.
#' If \code{plot.3d = FALSE}, a \code{ggpubr} object containing arranged
#' \code{ggplot2} panels.
#'
#' @family plotting functions
#' @seealso \code{\link{validate_plot_e_space_args}}, \code{\link{build_ellps}},
#'   \code{\link{get_suitable_env}}, \code{\link{create_virtual_species}}
#' @import RColorBrewer ggplot2 dplyr
#' @importFrom rlang .data
#' @export
plot_e_space <- function(env_bg,
                         x, y, z,
                         labels = c("ENV 1", "ENV 2", "ENV 3"),
                         n_bg = 10000,
                         niche = NULL,
                         suitable_env = NULL,
                         occ_pts = NULL,
                         rand_seed = 1234,
                         show.occ.density = FALSE, # only for 2D plots
                         plot.3d = FALSE,
                         colors = NULL,
                         palette = "default",
                         vs = NULL) {

  # --- 0. Optionally fill from NicheR_species object -------------------------

  if (!is.null(vs)) {
    if (!inherits(vs, "NicheR_species")) {
      stop("'vs' must be a NicheR_species object created by create_virtual_species().")
    }

    # Fill niche if missing
    if (is.null(niche) && !is.null(vs$niche)) {
      niche <- vs$niche
    }

    # Fill occurrences if missing
    if (is.null(occ_pts) && !is.null(vs$occurrences)) {
      occ_pts <- vs$occurrences
    }

    # Fill env_bg if missing, using suitability slot
    if (missing(env_bg) || is.null(env_bg)) {
      if (inherits(vs$suitability, "SpatRaster")) {
        env_bg <- as.data.frame.nicheR(vs$suitability)
      } else if (is.data.frame(vs$suitability)) {
        env_bg <- vs$suitability
      }
    }

    # Fill suitable_env if missing, using suitability slot directly
    if (is.null(suitable_env) && !is.null(vs$suitability)) {
      suitable_env <- vs$suitability
    }
  }

  # Helper to resolve column names (numeric or character)
  resolve_col <- function(df, col_spec) {
    if (is.numeric(col_spec)) {
      nm <- names(df)[col_spec]
    } else {
      nm <- col_spec
    }
    if (!nm %in% names(df)) {
      stop("Column '", nm, "' not found in env_bg.")
    }
    nm
  }

  # --- 1. Validate core arguments --------------------------------------------

  validate_plot_e_space_args(env_bg, x, y, z,
                             labels, n_bg, niche,
                             occ_pts, show.occ.density)

  # Resolve column names up front
  x_col <- resolve_col(env_bg, x)
  y_col <- resolve_col(env_bg, y)
  z_col <- resolve_col(env_bg, z)

  # Coerce occ_pts and ensure it has needed columns ---------------------------

  if (!is.null(occ_pts)) {
    occ_pts <- as.data.frame(occ_pts)
    missing_occ <- setdiff(c(x_col, y_col, z_col), names(occ_pts))
    if (length(missing_occ) > 0) {
      warning(
        "occ_pts is missing the following environmental columns: ",
        paste(missing_occ, collapse = ", "),
        ".\nOccurrences will not be plotted in E-space."
      )
      occ_pts <- NULL
    }
  }

  # Coerce suitable_env to a data.frame if provided ---------------------------

  pts_in <- NULL
  if (!is.null(suitable_env)) {

    if (inherits(suitable_env, "suitable_env") ||
        (is.list(suitable_env) &&
         any(c("suitable_env_df", "suitable_env_sp") %in% names(suitable_env)))) {

      if ("suitable_env_df" %in% names(suitable_env) &&
          is.data.frame(suitable_env$suitable_env_df)) {

        pts_in <- suitable_env$suitable_env_df

      } else if ("suitable_env_sp" %in% names(suitable_env)) {

        sp <- suitable_env$suitable_env_sp
        if (inherits(sp, "Raster")) sp <- terra::rast(sp)

        if (inherits(sp, "SpatRaster")) {
          pts_in <- as.data.frame.nicheR(sp)
        } else if (is.list(sp) &&
                   length(sp) > 0 &&
                   all(vapply(sp, inherits, logical(1), "SpatRaster"))) {

          # prefer 'suitable' if present, else first raster
          if ("suitable" %in% names(sp)) {
            pts_in <- as.data.frame.nicheR(sp[["suitable"]])
          } else {
            pts_in <- as.data.frame.nicheR(sp[[1]])
          }

        }
      }

    } else if (inherits(suitable_env, "Raster") || inherits(suitable_env, "SpatRaster")) {

      if (inherits(suitable_env, "Raster")) {
        suitable_env <- terra::rast(suitable_env)
      }
      pts_in <- as.data.frame.nicheR(suitable_env)

    } else if (is.data.frame(suitable_env) || is.matrix(suitable_env)) {

      pts_in <- as.data.frame(suitable_env)
    }

    if (!is.null(pts_in)) {
      # make sure it has the required columns; if not, drop it silently
      if (!all(c(x_col, y_col, z_col) %in% names(pts_in))) {
        warning(
          "suitable_env does not contain all of x, y, z columns; ",
          "suitable points will not be plotted."
        )
        pts_in <- NULL
      }
    }
  }

  # --- 2. Colors / palette ---------------------------------------------------

  palettes <- list(
    default = list(
      bg           = "#9093A2FF",
      ellipsoid    = "#2A363BFF",
      centroid     = "#D72000FF",
      tolerance    = "#EE6100FF",
      suitable_env = "#FED789FF",
      occ          = "#B4BF3AFF"
    ),
    palette2 = list(
      bg           = "#9CA9BAFF",
      ellipsoid    = "#3D619DFF",
      centroid     = "#345084FF",
      tolerance    = "#693829FF",
      suitable_env = "#CFB267FF",
      occ          = "#A56A3EFF"
    ),
    palette3 = list(
      bg           = "#C8CCC6FF",
      ellipsoid    = "#023743FF",
      centroid     = "#72874EFF",
      tolerance    = "#476F84FF",
      suitable_env = "#FED789FF",
      occ          = "#A4BED5FF"
    ),
    palette4 = list(
      bg           = "#C0D1CEFF",
      ellipsoid    = "#859B6CFF",
      centroid     = "#B74954FF",
      tolerance    = "#A99364FF",
      suitable_env = "#C2DDB2FF",
      occ          = "#EBA49EFF"
    ),
    palette5 = list(
      bg           = "#A89F8EFF",
      ellipsoid    = "#7887A4FF",
      centroid     = "#A8CDECFF",
      tolerance    = "#682C37FF",
      suitable_env = "#F6955EFF",
      occ          = "#9B6981FF"
    ),
    palette6 = list(
      bg           = "#D3D4D8FF",
      ellipsoid    = "#731A12FF",
      centroid     = "#F2D43DFF",
      tolerance    = "#3F858CFF",
      suitable_env = "#D9814EFF",
      occ          = "#707322FF"
    )
  )

  if (!palette %in% names(palettes)) {
    stop("Unknown palette '", palette, "'.")
  }

  base_colors <- palettes[[palette]]

  if (is.null(colors)) {
    colors <- base_colors
  } else {
    if (is.null(names(colors))) {
      names(colors) <- names(base_colors)[seq_along(colors)]
      if (!is.list(colors)) colors <- as.list(colors)
      message(
        "No names detected in 'colors'. ",
        "Using the order provided and filling missing entries with defaults.\n",
        "If you want to change a specific object, use a named list.\n",
        "Available options are: ",
        paste(names(base_colors), collapse = ", ")
      )
    }
    if (!is.list(colors)) colors <- as.list(colors)
    colors <- utils::modifyList(base_colors, colors)
  }

  # Legend options ------------------------------------------------------------

  opts <- list(
    background_point     = TRUE,
    trace_line           = !is.null(niche),
    centroid_point       = !is.null(niche),
    tolerance_range_line = !is.null(niche),
    suitable_point       = !is.null(pts_in),
    occurrence_point     = !is.null(occ_pts)
  )

  legend_items <- data.frame(
    id       = c("background_point","trace_line","centroid_point",
                 "tolerance_range_line","suitable_point","occurrence_point"),
    type     = c("point","line","point","line","point","point"),
    label    = c("Background environments","Niche boundary","Niche centroid",
                 "Tolerance ranges","Suitable environments","Occurrences"),
    color    = c(colors[["bg"]], colors[["ellipsoid"]], colors[["centroid"]],
                 colors[["tolerance"]], colors[["suitable_env"]], colors[["occ"]]),
    linetype = c(NA, 1, NA, 2, NA, NA),
    stringsAsFactors = FALSE
  )

  active <- logical(nrow(legend_items))
  for (i in seq_len(nrow(legend_items))) {
    active[i] <- isTRUE(opts[[ legend_items$id[i] ]])
  }
  legend_items <- legend_items[active, , drop = FALSE]

  # Build legend plot ---------------------------------------------------------

  if (nrow(legend_items) == 0) {
    legend_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  } else {
    top_y   <- 2
    spacing <- 0.25
    x_point <- 0.00
    x_text  <- 0.1
    x0_line <- -0.2
    x1_line <-  0.00

    legend_items <- legend_items %>%
      dplyr::mutate(
        row     = dplyr::row_number(),
        y       = top_y - (row - 1) * spacing,
        x_point = x_point,
        x_text  = x_text,
        x0_line = x0_line,
        x1_line = x1_line
      )

    legend_base <- ggplot2::ggplot() +
      ggplot2::coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(0, 2.5), clip = "off") +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "none")

    legend_plot <- legend_base +
      ggplot2::geom_point(
        data  = dplyr::filter(legend_items, type == "point"),
        ggplot2::aes(x = x_point, y = y, colour = color),
        size = 2
      ) +
      ggplot2::geom_segment(
        data  = dplyr::filter(legend_items, type == "line"),
        ggplot2::aes(x = x0_line, xend = x1_line, y = y, yend = y,
                     colour = color, linetype = linetype),
        linewidth = 0.5
      ) +
      ggplot2::geom_text(
        data  = legend_items,
        ggplot2::aes(x = x_text, y = y, label = label),
        hjust = 0
      ) +
      ggplot2::scale_colour_identity() +
      ggplot2::scale_linetype_identity()
  }

  # Downsample background -----------------------------------------------------

  if (nrow(env_bg) > n_bg) {
    message(sprintf("Sampling %d of %d rows from 'env_bg' for plotting.", n_bg, nrow(env_bg)))
    set.seed(rand_seed)
    env_bg <- env_bg[sample.int(nrow(env_bg), size = n_bg, replace = FALSE), ]
  }

  # ---------------------------------------------------------------------------
  # 3D plotting branch
  # ---------------------------------------------------------------------------

  if (isTRUE(plot.3d)) {

    p3 <- plotly::plot_ly(
      data = env_bg,
      x    = env_bg[[x_col]],
      y    = env_bg[[y_col]],
      z    = env_bg[[z_col]],
      type = "scatter3d",
      mode = "markers",
      marker = list(color = colors[["bg"]], size = 2),
      name   = "Background Environments"
    ) %>%
      plotly::layout(
        title = list(text = "Background Environments (E-space)"),
        scene = list(
          xaxis = list(title = list(text = labels[1])),
          yaxis = list(title = list(text = labels[2])),
          zaxis = list(title = list(text = labels[3]))
        ),
        legend = list(x = 0.05, y = 0.95)
      )

    # Ellipsoid surface + centroid
    if (!is.null(niche) && !is.null(niche$surface)) {

      surf <- as.data.frame(niche$surface)
      if (ncol(surf) >= 3) {
        p3 <- p3 %>%
          plotly::add_trace(
            data = surf,
            x    = surf[[1]],
            y    = surf[[2]],
            z    = surf[[3]],
            type = "scatter3d",
            mode = "lines",
            line = list(color = colors[["ellipsoid"]]),
            name = "Niche Boundary",
            inherit = FALSE
          )
      }

      if (!is.null(niche$center) && length(niche$center) >= 3) {
        p3 <- p3 %>%
          plotly::add_markers(
            x = niche$center[1],
            y = niche$center[2],
            z = niche$center[3],
            marker = list(color = colors[["centroid"]], size = 5),
            name   = "Niche Centroid"
          )
      }

      if (!is.null(pts_in)) {
        p3 <- p3 %>%
          plotly::add_markers(
            data = pts_in,
            x    = pts_in[[x_col]],
            y    = pts_in[[y_col]],
            z    = pts_in[[z_col]],
            marker = list(color = colors[["suitable_env"]], size = 3),
            name   = "Suitable Environments",
            inherit = FALSE
          )
      }
    }

    # Occurrence points in 3D
    if (!is.null(occ_pts)) {
      p3 <- p3 %>%
        plotly::add_markers(
          data = occ_pts,
          x    = occ_pts[[x_col]],
          y    = occ_pts[[y_col]],
          z    = occ_pts[[z_col]],
          marker = list(color = colors[["occ"]], size = 3),
          name   = "Sampled Occurrences",
          inherit = FALSE
        ) %>%
        plotly::layout(
          title  = "Virtual Niche and Sampled Occurrences in E-space",
          legend = list(x = 0.05, y = 0.95)
        )
    }

    return(p3)
  }

  # ---------------------------------------------------------------------------
  # 2D plotting branch
  # ---------------------------------------------------------------------------

  p_main_y_x <- ggplot2::ggplot(
    env_bg,
    ggplot2::aes(x = .data[[y_col]], y = .data[[x_col]])
  ) +
    ggplot2::geom_point(alpha = 0.5, color = colors[["bg"]], pch = ".") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title = ggplot2::element_blank())

  p_main_z_x <- ggplot2::ggplot(
    env_bg,
    ggplot2::aes(x = .data[[z_col]], y = .data[[x_col]])
  ) +
    ggplot2::geom_point(alpha = 0.5, color = colors[["bg"]], pch = ".") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title = ggplot2::element_blank())

  p_main_z_y <- ggplot2::ggplot(
    env_bg,
    ggplot2::aes(x = .data[[z_col]], y = .data[[y_col]])
  ) +
    ggplot2::geom_point(alpha = 0.5, color = colors[["bg"]], pch = ".") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title = ggplot2::element_blank())

  x_name <- ggplot2::ggplot() + ggplot2::theme_void() +
    ggplot2::geom_text(ggplot2::aes(0, 0, label = labels[1]),
                       fontface = "bold")
  y_name <- ggplot2::ggplot() + ggplot2::theme_void() +
    ggplot2::geom_text(ggplot2::aes(0, 0, label = labels[2]),
                       fontface = "bold")
  z_name <- ggplot2::ggplot() + ggplot2::theme_void() +
    ggplot2::geom_text(ggplot2::aes(0, 0, label = labels[3]),
                       fontface = "bold")

  return_plot <- ggpubr::ggarrange(
    x_name,      p_main_y_x, p_main_z_x,
    legend_plot, y_name,     p_main_z_y,
    NULL,        NULL,       z_name,
    ncol = 3, nrow = 3,
    widths  = c(0.15, 0.425, 0.425),
    heights = c(0.45, 0.45, 0.15)
  )

  # --- Ellipsoid overlays and extras (2D) ------------------------------------

  if (!is.null(niche)) {

    # Build 2D ellipses from 3D center/axes/angles
    center_y_x <- c(niche$center[2], niche$center[1])
    axes_y_x   <- c(niche$axes[2],   niche$axes[1])
    angle_y_x  <- c(niche$angles[2], niche$angles[1])

    center_z_x <- c(niche$center[3], niche$center[1])
    axes_z_x   <- c(niche$axes[3],   niche$axes[1])
    angle_z_x  <- c(niche$angles[3], niche$angles[1])

    center_z_y <- c(niche$center[3], niche$center[2])
    axes_z_y   <- c(niche$axes[3],   niche$axes[2])
    angle_z_y  <- c(niche$angles[3], niche$angles[2])

    ell2d_y_x <- build_ellps(center = center_y_x, axes = axes_y_x, angles = angle_y_x)
    ell2d_z_x <- build_ellps(center = center_z_x, axes = axes_z_x, angles = angle_z_x)
    ell2d_z_y <- build_ellps(center = center_z_y, axes = axes_z_y, angles = angle_z_y)

    ell_y_x <- p_main_y_x
    ell_z_x <- p_main_z_x
    ell_z_y <- p_main_z_y

    if (any(niche$angles != 0)) {
      angle_warn <- ggplot2::ggplot() + ggplot2::theme_void() +
        ggplot2::geom_text(
          ggplot2::aes(0, 0,
                       label = "Note: The ellipsoid is angled; its shape may appear distorted and some points may fall outside due to dimensionality."),
          size = 2
        )
    } else {
      angle_warn <- NULL
    }

    # Suitable points (inside)
    if (!is.null(pts_in)) {

      ell_y_x <- ell_y_x +
        ggplot2::geom_point(
          data = pts_in,
          ggplot2::aes(
            x = .data[[y_col]],
            y = .data[[x_col]]
          ),
          color = colors[["suitable_env"]],
          size  = 0.5
        )

      ell_z_x <- ell_z_x +
        ggplot2::geom_point(
          data = pts_in,
          ggplot2::aes(
            x = .data[[z_col]],
            y = .data[[x_col]]
          ),
          color = colors[["suitable_env"]],
          size  = 0.5
        )

      ell_z_y <- ell_z_y +
        ggplot2::geom_point(
          data = pts_in,
          ggplot2::aes(
            x = .data[[z_col]],
            y = .data[[y_col]]
          ),
          color = colors[["suitable_env"]],
          size  = 0.5
        )
    }

    # Occurrence points in 2D
    if (!is.null(occ_pts)) {
      ell_y_x <- ell_y_x +
        ggplot2::geom_point(
          data = occ_pts,
          ggplot2::aes(
            x = .data[[y_col]],
            y = .data[[x_col]]
          ),
          color = colors[["occ"]],
          size  = 0.5
        )

      ell_z_x <- ell_z_x +
        ggplot2::geom_point(
          data = occ_pts,
          ggplot2::aes(
            x = .data[[z_col]],
            y = .data[[x_col]]
          ),
          color = colors[["occ"]],
          size  = 0.5
        )

      ell_z_y <- ell_z_y +
        ggplot2::geom_point(
          data = occ_pts,
          ggplot2::aes(
            x = .data[[z_col]],
            y = .data[[y_col]]
          ),
          color = colors[["occ"]],
          size  = 0.5
        )
    }

    # Ellipse boundaries + axes + centroid -----------------------------------

    ell_y_x <- ell_y_x +
      ggplot2::geom_path(
        data   = ell2d_y_x$surface,
        ggplot2::aes(x, y),
        color = colors[["ellipsoid"]],
        linewidth = 0.5
      ) +
      ggplot2::annotate(
        "segment",
        x    = ell2d_y_x$center[1] - ell2d_y_x$axes[1],
        xend = ell2d_y_x$center[1] + ell2d_y_x$axes[1],
        y    = ell2d_y_x$center[2],
        yend = ell2d_y_x$center[2],
        color = colors[["tolerance"]],
        linetype = "dashed"
      ) +
      ggplot2::annotate(
        "segment",
        y    = ell2d_y_x$center[2] - ell2d_y_x$axes[2],
        yend = ell2d_y_x$center[2] + ell2d_y_x$axes[2],
        x    = ell2d_y_x$center[1],
        xend = ell2d_y_x$center[1],
        color = colors[["tolerance"]],
        linetype = "dashed"
      ) +
      ggplot2::annotate(
        "point",
        x = ell2d_y_x$center[1],
        y = ell2d_y_x$center[2],
        color = colors[["centroid"]],
        size  = 2
      )

    ell_z_x <- ell_z_x +
      ggplot2::geom_path(
        data   = ell2d_z_x$surface,
        ggplot2::aes(x, y),
        color = colors[["ellipsoid"]],
        linewidth = 0.5
      ) +
      ggplot2::annotate(
        "segment",
        x    = ell2d_z_x$center[1] - ell2d_z_x$axes[1],
        xend = ell2d_z_x$center[1] + ell2d_z_x$axes[1],
        y    = ell2d_z_x$center[2],
        yend = ell2d_z_x$center[2],
        color = colors[["tolerance"]],
        linetype = "dashed"
      ) +
      ggplot2::annotate(
        "segment",
        y    = ell2d_z_x$center[2] - ell2d_z_x$axes[2],
        yend = ell2d_z_x$center[2] + ell2d_z_x$axes[2],
        x    = ell2d_z_x$center[1],
        xend = ell2d_z_x$center[1],
        color = colors[["tolerance"]],
        linetype = "dashed"
      ) +
      ggplot2::annotate(
        "point",
        x = ell2d_z_x$center[1],
        y = ell2d_z_x$center[2],
        color = colors[["centroid"]],
        size  = 2
      )

    ell_z_y <- ell_z_y +
      ggplot2::geom_path(
        data   = ell2d_z_y$surface,
        ggplot2::aes(x, y),
        color = colors[["ellipsoid"]],
        linewidth = 0.5
      ) +
      ggplot2::annotate(
        "segment",
        x    = ell2d_z_y$center[1] - ell2d_z_y$axes[1],
        xend = ell2d_z_y$center[1] + ell2d_z_y$axes[1],
        y    = ell2d_z_y$center[2],
        yend = ell2d_z_y$center[2],
        color = colors[["tolerance"]],
        linetype = "dashed"
      ) +
      ggplot2::annotate(
        "segment",
        y    = ell2d_z_y$center[2] - ell2d_z_y$axes[2],
        yend = ell2d_z_y$center[2] + ell2d_z_y$axes[2],
        x    = ell2d_z_y$center[1],
        xend = ell2d_z_y$center[1],
        color = colors[["tolerance"]],
        linetype = "dashed"
      ) +
      ggplot2::annotate(
        "point",
        x = ell2d_z_y$center[1],
        y = ell2d_z_y$center[2],
        color = colors[["centroid"]],
        size  = 2
      )

    # Re-assemble main grid with ellipsoid overlays --------------------------

    return_plot <- ggpubr::ggarrange(
      x_name,      ell_y_x,   ell_z_x,
      legend_plot, y_name,    ell_z_y,
      NULL,        angle_warn, z_name,
      ncol = 3, nrow = 3,
      widths  = c(0.15, 0.425, 0.425),
      heights = c(0.45, 0.45, 0.15)
    )

    # Optional occurrence density panels -------------------------------------

    if (isTRUE(show.occ.density) && !is.null(occ_pts)) {

      rng_z <- range(env_bg[[z_col]], na.rm = TRUE)
      rng_y <- range(env_bg[[y_col]], na.rm = TRUE)
      rng_x <- range(env_bg[[x_col]], na.rm = TRUE)

      env_z_top <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[z_col]])) +
        ggplot2::geom_density(fill = colors[["occ"]], alpha = 0.6) +
        ggplot2::scale_x_continuous(limits = rng_z) +
        ggplot2::scale_y_continuous(n.breaks = 3) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.text.y  = ggplot2::element_text(size = 5)
        )

      env_y_top <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[y_col]])) +
        ggplot2::geom_density(fill = colors[["occ"]], alpha = 0.6) +
        ggplot2::scale_x_continuous(limits = rng_y) +
        ggplot2::scale_y_continuous(n.breaks = 3) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.text.y  = ggplot2::element_text(size = 5)
        )

      env_x_right <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[x_col]])) +
        ggplot2::geom_density(fill = colors[["occ"]], alpha = 0.6) +
        ggplot2::scale_x_continuous(limits = rng_x) +
        ggplot2::coord_flip() +
        ggplot2::scale_y_continuous(n.breaks = 3) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.y  = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.text.x  = ggplot2::element_text(size = 5)
        )

      env_y_right <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[y_col]])) +
        ggplot2::geom_density(fill = colors[["occ"]], alpha = 0.6) +
        ggplot2::scale_x_continuous(limits = rng_y) +
        ggplot2::coord_flip() +
        ggplot2::scale_y_continuous(n.breaks = 3) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.y  = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.text.x  = ggplot2::element_text(size = 5)
        )

      return_plot <- ggpubr::ggarrange(
        NULL,        env_y_top, env_z_top, NULL,
        x_name,      ell_y_x,   ell_z_x,   env_x_right,
        legend_plot, y_name,    ell_z_y,   env_y_right,
        NULL,        angle_warn, z_name,   NULL,
        ncol = 4, nrow = 4,
        widths  = c(0.1, 0.4, 0.4, 0.1),
        heights = c(0.1, 0.4, 0.4, 0.1)
      )
    }
  }

  return(return_plot)
}
