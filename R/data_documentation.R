#' Background environmental data for examples
#'
#' A dataset containing geographic coordinates and two bioclimatic variables
#' used as a background point cloud for generating and testing niche ellipsoids.
#'
#' @format A data frame with 12,396 rows and 4 variables:
#' \describe{
#'   \item{x}{Longitude in decimal degrees (WGS84).}
#'   \item{y}{Latitude in decimal degrees (WGS84).}
#'   \item{bio1}{Annual Mean Temperature (units: degrees Celsius).}
#'   \item{bio12}{Annual Precipitation (units: mm).}
#' }
#'
#' @details
#' This dataset represents an irregular point cloud typical of environmental
#' background data used in ecological niche modeling (ENM). It is primarily
#' used in the `nicheR` package to provide the environmental space for
#' functions like \code{\link{conserved_ellipses}}.
#'
#' @source \url{https://www.worldclim.org}
#' @examples
#' data(back_data)
#' head(back_data)
"back_data"