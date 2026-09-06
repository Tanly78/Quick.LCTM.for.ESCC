#' Predict Class Membership for New Patients
#'
#' Applies a pre-trained LCTM model to new patient data. Either uses
#' the embedded v8 posterior parameters (default) or a freshly fitted
#' model object.
#'
#' @param new_data A data.frame with columns \code{id}, \code{timepoint},
#'   and the 7 raw markers \code{ALB, TC, LYM, CRP, SMI, SFI, VFI}.
#'   \code{timepoint} accepts either short form (\code{"T1","T2","T4","T5","T6"})
#'   or long form (\code{"T1_pre","T2_post","T4_3M","T5_6M","T6_1Y"}).
#'   At least 2 distinct timepoints per patient are required.
#' @param model A fitted \code{lctm7m_fit} object (from
#'   \code{\link{fit_lctm_full}}). If NULL, the embedded v8 parameters
#'   are used.
#' @param use_pretrained Logical. If TRUE and model is NULL, use the
#'   embedded v8 posterior means. Default TRUE.
#'
#' @return A data.frame with columns \code{id}, \code{p_VC}, and
#'   \code{cluster} (factor, levels RC < VC).
#'
#' @details
#' When \code{use_pretrained = TRUE}, prediction uses the v8 posterior
#' means (intercept, slope, sigma, mixing proportions) plus the v8
#' z-score reference (per-marker mean/SD) bundled in
#' \code{\link{lctm_v8_params}}. Aggregation rule follows
#' \code{_bayesian_lctm.py} lines 36-53: per-marker z-score using the
#' v8 reference, CRP sign-flipped, then mean across 7 markers at each
#' timepoint. Posterior probabilities are computed by Bayes rule with
#' the class-conditional normal likelihood and the v8 class-mixing
#' proportions.
#'
#' @examples
#' \dontrun{
#' data(example_lctm)
#' pred <- predict_lctm7m(example_lctm, use_pretrained = TRUE)
#' head(pred)
#' }
#'
#' @export
predict_lctm7m <- function(new_data,
                          model = NULL,
                          use_pretrained = TRUE) {
  if (!is.null(model) && inherits(model, "lctm7m_fit")) {
    pp  <- loadNamespace("rstanarm")$posterior_linpred(model@stan_fit, newdata = new_data)
    p_VC <- 1 - colMeans(pp)
  } else if (isTRUE(use_pretrained)) {
    params <- quick.lctm.escc::lctm_v8_params
    pre    <- predict_with_pretrained(new_data, params)
    p_VC   <- pre$p_VC
  } else {
    stop("Provide a fitted model OR set use_pretrained = TRUE.")
  }
  if (isTRUE(use_pretrained) && is.null(model)) {
    out_id <- pre$id_keys
  } else {
    out_id <- unique(new_data$id)
  }
  data.frame(
    id      = out_id,
    p_VC    = round(p_VC, 4),
    cluster = factor(ifelse(p_VC > 0.5, "VC", "RC"),
                     levels = c("RC", "VC")),
    stringsAsFactors = FALSE
  )
}

# Internal: predict using bundled v8 posterior parameters
# Returns a numeric vector of p_VC aligned with unique(new_data$id).
predict_with_pretrained <- function(new_data, params) {
  composite <- compute_composite_z(new_data, params)
  pid_keys  <- rownames(composite)         # patient ids as character
  times     <- params$times
  mu_RC     <- params$class_means["RC", "intercept"] +
               params$class_means["RC", "slope"]     * times
  mu_VC     <- params$class_means["VC", "intercept"] +
               params$class_means["VC", "slope"]     * times
  sigma     <- params$sigma
  prop_RC   <- params$class_props["RC"]
  prop_VC   <- params$class_props["VC"]
  log_prior_RC <- unname(log(prop_RC))
  log_prior_VC <- unname(log(prop_VC))

  p_VC <- vapply(pid_keys, function(pid) {
    y    <- composite[pid, ]
    obs  <- !is.na(y)
    if (sum(obs) < 2L) return(NA_real_)
    ll <- c(
      RC = sum(dnorm(y[obs], mean = mu_RC[obs], sd = sigma, log = TRUE)) +
           log_prior_RC,
      VC = sum(dnorm(y[obs], mean = mu_VC[obs], sd = sigma, log = TRUE)) +
           log_prior_VC
    )
    ll_max <- max(ll)
    exp(ll["VC"] - ll_max) / sum(exp(ll - ll_max))
  }, numeric(1))
  list(p_VC = p_VC, id_keys = pid_keys)
}

# Internal: aggregate 7 raw markers to a per-timepoint composite z
# using the v8 reference mean/SD. Output is a 219-or-N x 5 matrix
# of composite z, rows keyed by patient id (as character), columns
# keyed by timepoint position 1..5.
compute_composite_z <- function(new_data, params) {
  if (!"id" %in% names(new_data)) stop("new_data must have an 'id' column")
  if (!"timepoint" %in% names(new_data)) stop("new_data must have a 'timepoint' column")

  MARKERS <- names(params$ref_mean)
  TP_SHORT <- c("T1","T2","T4","T5","T6")
  TP_LONG  <- params$timepoint_lab

  ## Normalize timepoint to long form
  tp <- as.character(new_data$timepoint)
  hit_short <- tp %in% TP_SHORT
  hit_long  <- tp %in% TP_LONG
  if (any(!hit_short & !hit_long)) {
    stop("timepoint must be one of: ",
         paste(TP_SHORT, collapse = ","), " or ",
         paste(TP_LONG,  collapse = ","))
  }
  tp_long <- ifelse(hit_short, TP_LONG[match(tp, TP_SHORT)], tp)

  new_data$timepoint_long <- tp_long
  ids <- sort(unique(new_data$id))
  composite <- matrix(NA_real_, nrow = length(ids), ncol = 5L,
                      dimnames = list(as.character(ids), NULL))

  for (j in seq_along(TP_LONG)) {
    tp_name <- TP_LONG[j]
    sub_tp  <- new_data[new_data$timepoint_long == tp_name, , drop = FALSE]
    if (nrow(sub_tp) == 0L) next
    z_block <- matrix(NA_real_, nrow = nrow(sub_tp), ncol = length(MARKERS))
    for (k in seq_along(MARKERS)) {
      m   <- MARKERS[k]
      if (!m %in% names(sub_tp)) {
        stop("new_data is missing marker column: ", m)
      }
      raw <- sub_tp[[m]]
      z   <- (raw - params$ref_mean[[m]]) / params$ref_sd[[m]]
      if (m %in% params$negative) z <- -z
      z_block[, k] <- z
    }
    rowMeans_z <- rowMeans(z_block, na.rm = TRUE)
    ## Per-id mean of rowMeans_z (a patient with multiple rows at one
    ## timepoint is unusual but tolerated by averaging).
    agg <- tapply(rowMeans_z, as.character(sub_tp$id), mean, na.rm = TRUE)
    composite[as.character(names(agg)), j] <- as.numeric(agg)
  }
  composite
}

#' Compute 7-marker composite z at each timepoint
#'
#' Internal helper. Computes the per-timepoint composite z-score
#' using a reference mean/SD (e.g., from the v8 training cohort).
#'
#' @param data A data.frame with id, timepoint, and 7 raw markers.
#' @param ref Named numeric vector of length 7: per-marker mean.
#' @param ref_sd Named numeric vector of length 7: per-marker SD.
#' @param negative Character vector of marker names to sign-flip.
#' @param timepoint_lab Character vector of 5 timepoint labels in
#'   the canonical long form (e.g., \code{c("T1_pre",...)}).
#'
#' @return A matrix (length(unique(id)) x 5) of composite z values.
#'
#' @keywords internal
#' @export
compute_7m_score <- function(data,
                             ref,
                             ref_sd,
                             negative = c("CRP"),
                             timepoint_lab = c("T1_pre","T2_post",
                                               "T4_3M","T5_6M","T6_1Y")) {
  params <- list(ref_mean = ref,
                 ref_sd   = ref_sd,
                 negative = negative,
                 timepoint_lab = timepoint_lab)
  compute_composite_z(data, params)
}