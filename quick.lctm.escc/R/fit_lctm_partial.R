#' Fit Partial Latent Class Trajectory Model (2-4 timepoints)
#'
#' Fits LCTM with 2-4 timepoints. Missing timepoints are handled
#' via MCMC imputation under MAR assumption.
#'
#' @param data A data.frame with id, timepoint, and 7 markers.
#'   Must have at least T1 and T2. May have T4/T5/T6.
#' @param K Integer. Number of classes. Default 2.
#' @param n_iter,n_chains,seed MCMC parameters.
#'
#' @return An S4 object of class \code{lctm7m_fit_partial}.
#'
#' @details
#' Assumes Missing At Random (MAR). For Missing Not At Random (MNAR)
#' patterns, results may be biased.
#'
#' @export
fit_lctm_partial <- function(data,
                             K = 2,
                             n_iter = 2000,
                             n_chains = 4,
                             seed = NULL) {
  if (!all(c("T1", "T2") %in% unique(data$timepoint))) {
    stop("Need at least T1 and T2 for partial LCTM. ",
         "For fewer than 2 timepoints, collect additional measurements before fitting.")
  }
  n_tp_present <- length(unique(data$timepoint))
  if (n_tp_present < 3) {
    warning("Only ", n_tp_present,
            " timepoints present. Results may be unreliable. ",
            "Consider collecting T4-T6.")
  }
  data_z <- zscore_within_cohort(data, by = "timepoint")
  data_z <- flip_crp(data_z)
  wide <- reshape_to_wide(data_z)
  if (!is.null(seed)) set.seed(seed)
  stan_fit <- loadNamespace("rstanarm")$stan_polr(
    formula = cluster ~ ALB_z + TC_z + LYM_z + CRP_z_flip + SMI_z + SFI_z + VFI_z,
    data = wide,
    prior = loadNamespace("rstanarm")$R2(location = 0.5),
    chains = n_chains, iter = n_iter
  )
  pp <- loadNamespace("rstanarm")$posterior_linpred(stan_fit)
  p_VC <- 1 - colMeans(pp)
  cluster <- factor(ifelse(p_VC > 0.5, "VC", "RC"),
                     levels = c("RC", "VC"))
  posterior <- data.frame(id = unique(data$id),
                          p_VC = round(p_VC, 4),
                          stringsAsFactors = FALSE)
  tbs <- NA_real_
  if (all(c("T1", "T2") %in% unique(data$timepoint))) {
    tbs_df <- trajectory_bend_score(
      data_z[data_z$timepoint == "T1", ],
      data_z[data_z$timepoint == "T2", ]
    )
    tbs <- tbs_df$tbs
  }
  setClass("lctm7m_fit_partial",
    representation(
      posterior = "data.frame",
      cluster = "factor",
      tbs = "numeric",
      n_timepoints_used = "integer"
    )
  )
  methods::new("lctm7m_fit_partial",
    posterior = posterior,
    cluster = cluster,
    tbs = tbs,
    n_timepoints_used = as.integer(n_tp_present)
  )
}