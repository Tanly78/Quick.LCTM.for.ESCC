#' Z-score Normalize Within Cohort by Timepoint
#'
#' Standardizes 7 markers to z-scores within the cohort, separately
#' for each timepoint. This ensures comparability across timepoints
#' and centers the data for the LCTM.
#'
#' @param data A data.frame with id, timepoint, and 7 markers.
#' @param by Character. Column to group by. Default "timepoint".
#' @param markers Character vector. Marker names to standardize.
#'   Default \code{c("ALB", "TC", "LYM", "CRP", "SMI", "SFI", "VFI")}.
#'
#' @return A data.frame with original columns + 7 z-scored columns
#'   (suffix \code{_z}).
#'
#' @export
zscore_within_cohort <- function(data,
                                 by = "timepoint",
                                 markers = c("ALB", "TC", "LYM", "CRP",
                                             "SMI", "SFI", "VFI")) {
  for (m in markers) {
    if (!m %in% names(data)) stop("Missing marker: ", m)
    z_col <- paste0(m, "_z")
    data[[z_col]] <- ave(data[[m]], data[[by]],
                         FUN = function(x) {
                           s <- sd(x, na.rm = TRUE)
                           if (is.na(s) || s == 0) return(rep(0, length(x)))
                           (x - mean(x, na.rm = TRUE)) / s
                         })
  }
  data
}

#' Sign-flip CRP for LCTM Convention
#'
#' In the LCTM convention, CRP is sign-flipped so that higher
#' values indicate better nutritional status (since elevated CRP
#' indicates inflammation, which is worse).
#'
#' @param data A data.frame with a \code{CRP_z} column.
#'
#' @return A data.frame with added \code{CRP_z_flip} column.
#'
#' @export
flip_crp <- function(data) {
  if (!"CRP_z" %in% names(data)) stop("Missing CRP_z column")
  data$CRP_z_flip <- -data$CRP_z
  data
}