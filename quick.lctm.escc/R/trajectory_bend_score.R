#' Compute Trajectory Bend Score (TBS)
#'
#' TBS = standardized change in 7-marker from T1 to T2.
#'
#' @param data_t1,data_t2 Data frames with id and 7 markers at T1, T2.
#' @param reference_sd Optional reference SD. Default uses pooled SD.
#'
#' @return A data.frame with id, tbs, and category.
#'
#' @export
trajectory_bend_score <- function(data_t1, data_t2, reference_sd = NULL) {
  z_cols <- c("ALB_z", "TC_z", "LYM_z", "CRP_z_flip",
              "SMI_z", "SFI_z", "VFI_z")
  if (!all(z_cols %in% names(data_t1))) {
    data_t1 <- zscore_within_cohort(data_t1)
    data_t1 <- flip_crp(data_t1)
  }
  if (!all(z_cols %in% names(data_t2))) {
    data_t2 <- zscore_within_cohort(data_t2)
    data_t2 <- flip_crp(data_t2)
  }
  s1 <- rowMeans(data_t1[, z_cols], na.rm = TRUE)
  s2 <- rowMeans(data_t2[, z_cols], na.rm = TRUE)
  delta <- s2 - s1
  if (is.null(reference_sd)) {
    reference_sd <- sd(c(s1, s2), na.rm = TRUE)
    if (is.na(reference_sd) || reference_sd == 0) reference_sd <- 1
  }
  tbs <- delta / reference_sd
  category <- cut(tbs,
                  breaks = c(-Inf, -0.3, 0, 0.3, Inf),
                  labels = c("worsening", "stable",
                             "weak_improvement", "strong_improvement"))
  data.frame(id = data_t1$id, tbs = round(tbs, 3), category = category)
}