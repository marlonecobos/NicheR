#' Internal imports for NicheR
#'
#' These imports are centralized here to keep roxygen/NAMESPACE tidy.
#' Nothing in this file is exported.
#'
#' @keywords internal
#'
#' @import ggplot2
#' @importFrom dplyr mutate filter row_number
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggpubr ggarrange
#' @importFrom plotly plot_ly add_trace add_markers layout
#' @importFrom sf st_as_sf
#' @importFrom terra ncell as.data.frame ext nlyr xmin xmax ymin ymax crop
#' @importFrom raster ncell extent
#'
#' # Database and SQL utilities used in convert_large_raster()
#' @importFrom DBI dbConnect dbDisconnect dbWriteTable dbReadTable
#' @importFrom RSQLite SQLite
NULL
