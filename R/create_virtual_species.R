#' Create a virtual species end to end
#'
#' Orchestrates the full NicheR workflow by routing `...` to
#' [build_ellps()], [get_suitable_env()], optionally [set_bias_surface()],
#' and [get_sample_occ()]. Returns an S3 object of class **NicheR_species**
#' containing the niche object, suitability, (optional) bias surface, and
#' sampled occurrences.
#'
#' @section Workflow:
#' 1. Calls [build_ellps()] to define the ellipsoid in E-space.
#' 2. Calls [get_suitable_env()] to compute suitability in G-space.
#' 3. Optionally calls [set_bias_surface()] to build a pooled bias raster
#'    when bias-construction arguments are supplied (e.g., `bias_dir`,
#'    `template`, `res`, `ext`, `out.bias`).
#' 4. Calls [get_sample_occ()] to sample occurrences from the suitable
#'    (and optionally biased) environment.
#'
#' @param ... Named arguments passed to the component functions.
#'   Arguments are routed by name to the matching formal parameters of
#'   [build_ellps()], [get_suitable_env()], [set_bias_surface()], and
#'   [get_sample_occ()]. Unknown names are ignored with a warning.
#'   \itemize{
#'     \item \strong{build_ellps()}: e.g. `center`, `axes`, `angles`, etc.
#'     \item \strong{get_suitable_env()}: e.g. `env_bg`, `out`, `distances`.
#'     \item \strong{set_bias_surface()}: e.g. `bias_surface`, `bias_dir`,
#'           `template`, `res`, `ext`, `out.bias`.
#'     \item \strong{get_sample_occ()}: e.g. `n_occ`, `method`,
#'           `bias_surface`, `suitable_env`, `seed`.
#'   }
#' @param out.file Logical. If `TRUE`, save the returned object to an `.rds`
#'   file in the working directory.
#' @param out.file.name Optional base name (without extension) for the saved
#'   file. If `NULL` or empty, a timestamped name is used.
#' @param verbose Logical. If `TRUE`, print progress messages.
#'
#' @details
#' If `env_bg` is supplied only at the top level, it is forwarded to both
#' [get_suitable_env()] and [get_sample_occ()] (if those functions accept an
#' `env_bg` argument). If a required argument (such as `env_bg`) is missing
#' for a component, the function stops with an informative error.
#'
#' \strong{Bias handling:}
#' \itemize{
#'   \item If the user passes only `bias_surface` (and no `bias_dir`,
#'         `template`, `res`, `ext`, or `out.bias`), it is assumed to be a
#'         \emph{precomputed single-layer raster} and is passed directly to
#'         [get_sample_occ()] (which expects a 0–1 layer).
#'   \item If the user passes `bias_surface` \emph{plus} any of the bias
#'         construction arguments (`bias_dir`, `template`, `res`, `ext`,
#'         `out.bias`), then [set_bias_surface()] is called once, and its
#'         pooled bias raster is passed to [get_sample_occ()] as
#'         `bias_surface`. The full bias object is stored in the output
#'         as `bias_surface`.
#' }
#'
#' @return
#' A list of class **NicheR_species** with elements:
#' \itemize{
#'   \item \code{niche}: the ellipsoid object returned by [build_ellps()].
#'   \item \code{suitability}: suitability object from [get_suitable_env()].
#'   \item \code{bias_surface}: either `NULL`, a \code{nicheR_bias_surface}
#'         object (if [set_bias_surface()] was used), or a user-supplied
#'         single-layer bias raster.
#'   \item \code{occurrences}: sampled occurrences from [get_sample_occ()].
#'   \item \code{call_args}: the original `...` captured as a named list.
#'   \item \code{routed_args}: a list showing which args went to each function.
#'   \item \code{save_path}: file path if \code{out.file = TRUE}, otherwise `NULL`.
#' }
#'
#' @seealso [build_ellps()], [get_suitable_env()], [set_bias_surface()],
#'   [get_sample_occ()]
#'
#' @export
create_virtual_species <- function(...,
                                   out.file = FALSE,
                                   out.file.name = NULL,
                                   verbose = TRUE) {

  # capture all args once
  args <- list(...)

  # ensure required functions exist
  needed_funs <- c("build_ellps", "get_suitable_env",
                   "set_bias_surface", "get_sample_occ")
  missing_funs <- needed_funs[!vapply(needed_funs, exists, logical(1), mode = "function")]

  if (length(missing_funs)) {
    stop("These functions are not available in the current environment: ",
         paste(missing_funs, collapse = ", "))
  }

  # pull formals for routing
  f_build <- names(formals(build_ellps))
  f_suit  <- names(formals(get_suitable_env))
  f_bias  <- names(formals(set_bias_surface))   # includes out.bias
  f_occ   <- names(formals(get_sample_occ))

  # split args by target function
  args_build <- args[names(args) %in% f_build]
  args_suit  <- args[names(args) %in% f_suit]
  args_bias  <- args[names(args) %in% f_bias]
  args_occ   <- args[names(args) %in% f_occ]

  # pass env_bg through if user supplied it at top level but it was not routed
  if ("env_bg" %in% names(args)) {
    if (!("env_bg" %in% names(args_suit)) && ("env_bg" %in% f_suit)) {
      args_suit$env_bg <- args$env_bg
    }
    if (!("env_bg" %in% names(args_occ)) && ("env_bg" %in% f_occ)) {
      args_occ$env_bg <- args$env_bg
    }
  }

  # check that env_bg will be available when needed for suitability
  if (("env_bg" %in% f_suit) && !("env_bg" %in% names(args_suit))) {
    stop("env_bg is required by get_suitable_env but was not supplied.")
  }

  # ---------------------------------------------------------------------------
  # 1) build ellipsoid niche
  # ---------------------------------------------------------------------------
  niche_obj <- tryCatch(
    do.call(build_ellps, args_build),
    error = function(e) stop("build_ellps failed: ", e$message)
  )

  if (verbose) message("Built niche object.")

  # ---------------------------------------------------------------------------
  # 2) suitability in G space
  # ---------------------------------------------------------------------------
  suit_env <- tryCatch(
    do.call(get_suitable_env, c(list(niche = niche_obj), args_suit)),
    error = function(e) stop("get_suitable_env failed: ", e$message)
  )

  if (verbose) message("Computed suitable environments.")

  # ---------------------------------------------------------------------------
  # 3) Optional bias construction
  # ---------------------------------------------------------------------------
  bias_obj <- NULL

  # Does the user supply ANY bias argument at all?
  has_any_bias_arg <- "bias_surface" %in% names(args) ||
    any(c("bias_dir", "template", "res", "ext", "out.bias") %in% names(args))

  if (has_any_bias_arg) {

    # Decide if we should CALL set_bias_surface() or just forward bias_surface.
    # We call set_bias_surface() only if at least one "construction" arg
    # is present in args_bias (beyond bias_surface/verbose/out.bias).
    bias_construction_args <- setdiff(names(args_bias),
                                      c("bias_surface", "verbose", "out.bias"))

    wants_bias_build <- length(bias_construction_args) > 0

    if (wants_bias_build) {
      # We will build a bias surface with set_bias_surface()

      # assemble bias args
      bias_call_args <- args_bias

      # If no template given, try to infer from env_bg (used in get_suitable_env)
      if (!("template" %in% names(bias_call_args)) ||
          is.null(bias_call_args$template)) {

        if ("env_bg" %in% names(args_suit) &&
            inherits(args_suit$env_bg, c("SpatRaster", "Raster"))) {

          bias_call_args$template <- if (inherits(args_suit$env_bg, "Raster")) {
            terra::rast(args_suit$env_bg)
          } else {
            args_suit$env_bg
          }
        }
      }

      if (!"bias_surface" %in% names(bias_call_args) ||
          is.null(bias_call_args$bias_surface)) {
        warning("Bias-construction arguments were detected but 'bias_surface' is missing; ",
                "skipping bias surface construction.", call. = FALSE)
      } else {

        bias_obj <- tryCatch(
          do.call(set_bias_surface, bias_call_args),
          error = function(e) stop("set_bias_surface failed: ", e$message)
        )

        if (verbose) message("Constructed pooled bias surface via set_bias_surface().")

        # get_sample_occ() expects a single-layer bias raster;
        # pass the pooled layer.
        if (!is.null(bias_obj$pooled_bias_sp)) {
          args_occ$bias_surface <- bias_obj$pooled_bias_sp
        } else {
          warning("set_bias_surface() returned no 'pooled_bias_sp'; ",
                  "bias will not be applied in get_sample_occ().",
                  call. = FALSE)
        }
      }

    } else {
      # No construction args: assume the user passed a precomputed 1-layer raster
      # as 'bias_surface' that should go straight into get_sample_occ().
      if ("bias_surface" %in% names(args)) {
        args_occ$bias_surface <- args$bias_surface
        bias_obj <- args$bias_surface  # store in output for transparency
        if (verbose) {
          message("Using user-supplied 'bias_surface' directly in get_sample_occ().")
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 4) sample occurrences
  # ---------------------------------------------------------------------------
  occ <- tryCatch(
    do.call(get_sample_occ, c(list(niche = niche_obj), args_occ)),
    error = function(e) stop("get_sample_occ failed: ", e$message)
  )

  if (verbose) message("Sampled occurrences.")

  # ---------------------------------------------------------------------------
  # 5) warn on unused args
  # ---------------------------------------------------------------------------
  used_names <- union(
    names(args_build),
    union(names(args_suit), union(names(args_occ), names(args_bias)))
  )
  unused <- setdiff(names(args), used_names)
  if (length(unused) && verbose) {
    warning("These arguments did not match any target function and were ignored: ",
            paste(unused, collapse = ", "))
  }

  # ---------------------------------------------------------------------------
  # 6) assemble S3 object
  # ---------------------------------------------------------------------------
  out <- list(
    niche        = niche_obj,
    suitability  = suit_env,
    bias_surface = bias_obj,
    occurrences  = occ,
    call_args    = args,
    routed_args  = list(
      build_ellps      = args_build,
      get_suitable_env = args_suit,
      set_bias_surface = args_bias,
      get_sample_occ   = args_occ
    ),
    save_path    = NULL
  )

  class(out) <- c("NicheR_species", class(out))

  # ---- Saving logic ----
  if (isTRUE(out.file)) {

    dir_path <- getwd()

    if (is.null(out.file.name) || out.file.name == "") {
      timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
      out.file.name <- paste0("NicheR_vs_", timestamp)
    }

    save_path <- file.path(dir_path, paste0(out.file.name, ".rds"))
    saveRDS(out, save_path)

    if (verbose) message("Virtual species saved to: ", normalizePath(save_path))

    out$save_path <- save_path

  } else if (verbose) {
    message("`out.file = FALSE`: object not saved to disk.")
  }

  return(out)
}

#' Print a NicheR virtual species
#'
#' @param x A \code{NicheR_species} object.
#' @param ... Not used.
#' @return \code{x}, invisibly.
#' @method print NicheR_species
#' @export
print.NicheR_species <- function(x, ...) {
  cat("NicheR virtual species components:\n")
  cat("  niche:        ", paste(class(x$niche), collapse = "/"), "\n", sep = "")
  cat("  suitability:  ", paste(class(x$suitability), collapse = "/"), "\n", sep = "")
  cat("  bias_surface: ",
      if (is.null(x$bias_surface)) "NULL" else paste(class(x$bias_surface), collapse = "/"),
      "\n", sep = "")
  cat("  occurrences:  ", paste(class(x$occurrences), collapse = "/"), "\n", sep = "")
  invisible(x)
}
