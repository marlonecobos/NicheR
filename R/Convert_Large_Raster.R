#' Convert a Large Raster Stack to a Disk-Backed Table (Chunked) + RData
#'
#' Converts a (potentially very large) raster stack into a SQLite table by
#' iterating over row-chunks, appending each chunk as a data frame of cell
#' values with `x`/`y` coordinates. After processing, the full table is read
#' back into R and saved as an `.RData` file. Optionally keeps the SQLite file
#' for later DB-backed workflows.
#'
#' @param raster_stack A multi-layer raster object. Supports
#'   \code{terra::SpatRaster} or \code{raster::Raster*}. Must have valid
#'   georeferencing so that \code{xmin}, \code{xmax}, \code{ymin}, and
#'   row-to-\eqn{y} conversions are defined.
#' @param out_filename Character base name (without extension) for outputs.
#'   Files \code{<out_filename>.sqlite} and \code{<out_filename>.RData} are
#'   written into the current working directory.
#' @param chunk_height Positive integer. Number of raster rows per chunk to
#'   process and append to the SQL table. Smaller values reduce memory usage
#'   but increase the number of DB writes.
#' @param keep_sql Logical. If \code{TRUE}, the temporary SQLite file is kept
#'   on disk after the `.RData` is saved. If \code{FALSE} (default), it is
#'   deleted at the end.
#' @param verbose Logical. If \code{TRUE}, prints progress and file paths.
#'
#' @return (Invisibly) returns a \code{data.frame} with columns \code{x}, \code{y},
#'   and one column per raster layer (matching \code{names(raster_stack)}).
#'   Also writes \code{<out_filename>.RData} containing an object named
#'   \code{env_df}. If \code{keep_sql = TRUE}, also leaves
#'   \code{<out_filename>.sqlite} on disk with table \code{raster_table}.
#'
#' @details
#' The function:
#' \enumerate{
#' \item Creates (or overwrites) a SQLite DB and an empty table \code{raster_table}
#'       with schema \code{x, y, layer1, ..., layerK}.
#' \item Iterates over raster rows in blocks of \code{chunk_height}, crops each
#'       block by geographic extent, converts it to a data frame with
#'       \code{xy = TRUE}, and appends it to \code{raster_table}.
#' \item Reads the full table back into R as \code{env_df}, disconnects from
#'       the DB, saves \code{env_df} to \code{<out_filename>.RData}, and
#'       optionally removes the SQLite file.
#' }
#'
#' This design avoids loading the entire raster into memory and is suitable for
#' stacks with many layers or very large dimensions.
#'
#' @section Notes:
#' \itemize{
#' \item Output files are written in \code{getwd()} unless absolute paths are provided
#'       in \code{out_filename}.
#' \item For reproducibility and resilience, consider wrapping the DB connection in
#'       \code{on.exit(DBI::dbDisconnect(conn), add = TRUE)}.
#' }
#'
#' @seealso \code{\link[DBI]{dbConnect}}, \code{\link[RSQLite]{SQLite}},
#'   \code{\link[terra]{as.data.frame}}, \code{\link[terra]{crop}},
#'   \code{\link[terra]{ext}}
#'
#' @importFrom DBI dbConnect dbDisconnect dbWriteTable dbReadTable
#' @importFrom RSQLite SQLite
#' @importFrom terra nlyr xmin xmax ymin ext crop yFromRow as.data.frame
#' @export
convert_large_raster <- function(raster_stack,
                                 out_filename = "raster_db",
                                 chunk_height = 500,
                                 keep_sql = FALSE,
                                 verbose = TRUE) {

  # Create database connection
  sqlite_file <- paste0(out_filename, ".sqlite")
  rdata_file <- paste0(out_filename, ".RData")

  conn <- dbConnect(SQLite(), dbname = sqlite_file)

  # Create empty table
  empty_df <- as.data.frame(matrix(ncol = (nlyr(raster_stack) + 2), nrow = 0))
  colnames(empty_df) <- c("x", "y", names(raster_stack))
  dbWriteTable(conn, "raster_table", empty_df, overwrite = TRUE)

  # Get total rows for chunking
  total_rows <- nrow(raster_stack)

  # Loop over chunks
  for (start_row in seq(1, total_rows, by = chunk_height)) {
    end_row <- min(start_row + chunk_height - 1, total_rows)
    if(verbose) message("Processing rows ", start_row, " to ", end_row, " of ", total_rows)

    # Convert row numbers to y-coordinates
    y_top <- yFromRow(raster_stack, start_row)

    # Handle last chunk - use ymin instead of going beyond raster extent
    if (end_row == total_rows) {
      y_bottom <- ymin(raster_stack)
    } else {
      y_bottom <- yFromRow(raster_stack, end_row + 1)
    }

    # Build extent in spatial coordinates
    e <- ext(xmin(raster_stack), xmax(raster_stack), y_bottom, y_top)

    # Crop raster for this block
    r_block <- crop(raster_stack, e)

    # Convert to dataframe
    df_block <- as.data.frame(r_block, xy = TRUE, na.rm = TRUE)

    # Append to database
    if (nrow(df_block) > 0) {
      dbWriteTable(conn, "raster_table", df_block, append = TRUE)
    }
  }

  # Read final table from database
  env_df <- dbReadTable(conn, "raster_table")

  # Disconnect from database
  dbDisconnect(conn)

  # Save as RData
  if(verbose) message("Saving RData file: ", rdata_file)
  save(env_df, file = rdata_file)

  # Delete SQL file if requested
  if (!keep_sql) {
    if(verbose) message("Deleting SQLite file: ", sqlite_file)
    file.remove(sqlite_file)
  } else {
    if(verbose) message("SQLite file kept: ", sqlite_file)
  }

  if(verbose) message("Processing complete. Data frame has ", nrow(env_df), " rows and ", ncol(env_df), " columns")

  return(invisible(env_df))
}
