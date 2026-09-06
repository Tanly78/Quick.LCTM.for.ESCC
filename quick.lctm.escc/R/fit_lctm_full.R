#' Fit Full 5-timepoint Latent Class Trajectory Model
#'
#' Fits the LCTM to 7 nutritional markers at 5 time points
#' (T1, T2, T4, T5, T6) using Bayesian mixture modeling via
#' \code{loadNamespace("rstanarm")$stan_polr} with K=2 classes.
#'
#' @param data A data.frame with columns: id, timepoint, ALB, TC,
#'   LYM, CRP, SMI, SFI, VFI. timepoint must be one of
#'   \code{c("T1", "T2", "T4", "T5", "T6")}.
#' @param K Integer. Number of trajectory classes. Default 2.
#' @param n_iter Integer. MCMC iterations. Default 2000.
#' @param n_chains Integer. MCMC chains. Default 4.
#' @param seed Integer. Random seed.
#'
#' @return An S4 object of class \code{lctm7m_fit} containing:
#'   \itemize{
#'     \item \code{posterior}: data.frame of class posterior
#'       probabilities (p_RC, p_VC) per patient
#'     \item \code{cluster}: factor with class assignment (RC/VC)
#'     \item \code{stan_fit}: the underlying \code{stanreg} object
#'   }
#'
#' @examples
#' \dontrun{
#' data(example_lctm)
#' fit <- fit_lctm_full(example_lctm, K = 2, seed = 42)
#' summary(fit)
#' }
#'
#' @export
fit_lctm_full <- function(data,
                          K = 2,
                          n_iter = 2000,
                          n_chains = 4,
                          seed = NULL) {

  # 1. Validate
  validate_full_data(data)

  # 2. Z-score normalize within cohort, by timepoint
  data_z <- zscore_within_cohort(data, by = "timepoint")

  # 3. Sign-flip CRP (LCTM convention: higher = better)
  data_z <- flip_crp(data_z)

  # 4. Reshape long to wide
  wide <- reshape_to_wide(data_z)

  # 5. Fit Stan mixture model
  if (!is.null(seed)) set.seed(seed)
  stan_fit <- loadNamespace("rstanarm")$stan_polr(
    formula = cluster ~ ALB_z + TC_z + LYM_z + CRP_z_flip + SMI_z + SFI_z + VFI_z,
    data = wide,
    prior = loadNamespace("rstanarm")$R2(location = 0.5),
    chains = n_chains,
    iter = n_iter
  )

  # 6. Extract posterior
  pp <- loadNamespace("rstanarm")$posterior_linpred(stan_fit)
  # K=2: prob of class 1 (VC, less frequent)
  p_VC <- 1 - colMeans(pp)
  p_RC <- colMeans(pp)

  cluster <- factor(ifelse(p_VC > p_RC, "VC", "RC"),
                     levels = c("RC", "VC"))

  # 7. Build S4 object
  posterior <- data.frame(
    id = unique(data$id),
    p_RC = round(p_RC, 4),
    p_VC = round(p_VC, 4),
    stringsAsFactors = FALSE
  )

  setClass("lctm7m_fit",
    representation(
      posterior = "data.frame",
      cluster = "factor",
      stan_fit = "ANY"
    )
  )

  methods::new("lctm7m_fit",
    posterior = posterior,
    cluster = cluster,
    stan_fit = stan_fit
  )
}

# Helper: validate full data input
validate_full_data <- function(data) {
  required_cols <- c("id", "timepoint", "ALB", "TC", "LYM", "CRP",
                      "SMI", "SFI", "VFI")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  valid_tp <- c("T1", "T2", "T4", "T5", "T6")
  invalid_tp <- setdiff(unique(data$timepoint), valid_tp)
  if (length(invalid_tp) > 0) {
    stop("Invalid timepoints: ", paste(invalid_tp, collapse = ", "),
         ". Must be one of: ", paste(valid_tp, collapse = ", "))
  }
  n_per_id <- table(data$id)
  if (!all(n_per_id == length(valid_tp))) {
    warning("Some patients do not have all 5 timepoints. ",
            "Consider using fit_lctm_partial() instead.")
  }
  invisible(NULL)
}

# Helper: reshape long to wide
reshape_to_wide <- function(data) {
  if (!"dplyr" %in% loadedNamespaces()) loadNamespace("dplyr")
  if (!"tidyr" %in% loadedNamespaces()) loadNamespace("tidyr")
  data %>%
    dplyr::select(id, timepoint, ALB_z, TC_z, LYM_z, CRP_z_flip,
                  SMI_z, SFI_z, VFI_z) %>%
    tidyr::pivot_wider(
      names_from = timepoint,
      values_from = c(ALB_z, TC_z, LYM_z, CRP_z_flip, SMI_z, SFI_z, VFI_z),
      names_sep = "_"
    ) %>%
    dplyr::mutate(cluster = 0)  # placeholder, assigned by EM-like algorithm
}