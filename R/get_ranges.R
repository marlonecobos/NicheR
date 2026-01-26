#' Computes ranges from a data.frame of values with optional expansion
#'
#' @param data A data.frame of at least one column. Each column should
#' contain numeric values.
#' @param expand_min A named vector or list of percentages (e.g., 0.1 for 10%)
#' defining how much to expand the minimum value of specific columns.
#' @param expand_max A named vector or list of percentages defining how
#' much to expand the maximum value of specific columns.
#'
#' @return
#' A data.frame with the (potentially expanded) minimum and maximum
#' values of each variable.
#'
#' @export
#'
#' @examples
#' data <- data.frame(var1 = c(0, 10), var2 = c(50, 100))
#' # Expand var1 min by 10% and var2 max by 20%
#' ranges_from_data(data, expand_min = list(var1 = 0.1),
#'                  expand_max = list(var2 = 0.2))

ranges_from_data <- function(data, expand_min = NULL, expand_max = NULL) {
  if (!is.data.frame(data) || ncol(data) < 1) {
    stop("'data' must be a 'data.frame' with at least one column.")
  }

  # Initial range computation: returns a matrix where row 1 is min, row 2 is max
  ranges <- sapply(data, range, na.rm = TRUE)
  
  # Iterate through columns to apply expansion
  for (col_name in colnames(ranges)) {
    current_min <- ranges[1, col_name]
    current_max <- ranges[2, col_name]
    abs_range <- current_max - current_min
    
    # Expand Minimum
    if (!is.null(expand_min) && col_name %in% names(expand_min)) {
      ranges[1, col_name] <- current_min - (abs_range * expand_min[[col_name]])
    }
    
    # Expand Maximum
    if (!is.null(expand_max) && col_name %in% names(expand_max)) {
      ranges[2, col_name] <- current_max + (abs_range * expand_max[[col_name]])
    }
  }

  ranges_df <- as.data.frame(ranges)
  rownames(ranges_df) <- c("min", "max")

  return(ranges_df)
}



#' Computes ranges from mean and standard deviation with optional expansion
#'
#' @param mean A named numeric vector of mean values for each variable.
#' @param sd A named numeric vector of standard deviation values for each variable.
#' @param cl A numeric value indicating the confidence level (default 95).
#' @param expand_min A named vector or list of percentages (e.g., 0.1 for 10%) 
#' defining how much to expand the minimum value of specific columns.
#' @param expand_max A named vector or list of percentages defining how 
#' much to expand the maximum value of specific columns.
#'
#' @return
#' A data.frame with the (potentially expanded) minimum and maximum 
#' values based on the confidence limits of a Normal distribution.
#'
#' @export
#'
#' @examples
#' m <- c(var1 = 10, var2 = 100)
#' s <- c(var1 = 2, var2 = 15)
#' # Compute 95% CL and expand var1 min by 10%
#' ranges_from_stats(mean = m, sd = s, cl = 95, expand_min = list(var1 = 0.1))

ranges_from_stats <- function(mean, sd, cl = 95, expand_min = NULL, expand_max = NULL) {
  if (!all(names(mean) %in% names(sd)) || length(mean) != length(sd)) {
    stop("mean and sd must be named vectors of the same length and share the same names.")
  }

  # Calculate tail probabilities
  alpha <- 1 - (cl / 100)
  p_lower <- alpha / 2
  p_upper <- 1 - (alpha / 2)

  # Compute confidence limits using Normal distribution
  # These serve as our base min and max
  mins <- qnorm(p_lower, mean = mean, sd = sd)
  maxs <- qnorm(p_upper, mean = mean, sd = sd)
  
  # Initialize results matrix
  var_names <- names(mean)
  ranges <- matrix(c(mins, maxs), nrow = 2, byrow = TRUE)
  colnames(ranges) <- var_names
  
  # Iterate through variables to apply expansion
  for (var in var_names) {
    current_min <- ranges[1, var]
    current_max <- ranges[2, var]
    abs_range <- current_max - current_min
    
    # Expand Minimum
    if (!is.null(expand_min) && var %in% names(expand_min)) {
      ranges[1, var] <- current_min - (abs_range * expand_min[[var]])
    }
    
    # Expand Maximum
    if (!is.null(expand_max) && var %in% names(expand_max)) {
      ranges[2, var] <- current_max + (abs_range * expand_max[[var]])
    }
  }

  ranges_df <- as.data.frame(ranges)
  rownames(ranges_df) <- c("min", "max")

  return(ranges_df)
}