#' Internal imports for NicheR
#'
#' This file collects roxygen2 @import directives so the NAMESPACE
#' remains clean. Nothing here is exported.
#'
#' @keywords internal
#' @name NicheR-internal-imports
NULL


# Core tidyverse-style utilities
# pipe
#' @importFrom magrittr %>%
NULL

# data pronoun for non-standard evaluation (ggplot2/dplyr)
#' @importFrom rlang .data
NULL

# dplyr verbs used unqualified in helpers / plotting
#' @importFrom dplyr mutate filter row_number
NULL


# Plotting: ggplot2 + ggpubr + sf + plotly
# ggplot2 (only the bits we actually use in NicheR)
#' @importFrom ggplot2 ggplot aes
#' @importFrom ggplot2 geom_point geom_tile geom_path geom_segment geom_text geom_density
#' @importFrom ggplot2 theme_bw theme_void theme_minimal theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 coord_quickmap coord_cartesian coord_flip
#' @importFrom ggplot2 xlab ylab ggtitle
#' @importFrom ggplot2 scale_colour_identity scale_linetype_identity
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 map_data
NULL

# panel layout for plot_e_space() and plot_g_space()
#' @importFrom ggpubr ggarrange
NULL

# simple features for occurrence plotting in geographic space
#' @importFrom sf st_as_sf
NULL

# interactive 3D plotting in plot_e_space(plot.3d = TRUE)
# (You can drop this if you always use plotly:: explicitly.)
#' @importFrom plotly plot_ly add_trace add_markers layout
NULL


# Spatial / raster handling
# terra used throughout for raster handling and conversion
#' @importFrom terra ncell nlyr ext rast as.data.frame xmin xmax ymin ymax
NULL

# raster is only used for class detection / legacy objects
# (if you no longer call raster::ncell / raster::extent, you can remove these)
#' @importFrom raster ncell extent
NULL


# DB / large raster utilities (used in as.data.frame.nicheR() etc.)
#' @importFrom DBI dbConnect dbDisconnect dbWriteTable dbReadTable
#' @importFrom RSQLite SQLite
NULL


# Base / stats / utils helpers


# stats helpers used without namespace (e.g. runif(), complete.cases())
#' @importFrom stats runif complete.cases
NULL

# utils helpers (e.g. modifyList() in palette merging)
#' @importFrom utils modifyList
NULL
