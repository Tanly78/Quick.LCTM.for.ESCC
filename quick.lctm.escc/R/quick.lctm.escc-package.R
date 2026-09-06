#' quick.lctm.escc: Latent Class Trajectory Modeling for 7 Nutritional Markers
#'
#' Fits Bayesian latent class trajectory models (LCTM) on 7 markers
#' (ALB, TC, LYM, CRP, SMI, SFI, VFI) at multiple time points to
#' identify nutritional trajectory phenotypes in cancer patients.
#'
#' Two-level input design (7 markers only; no baseline covariates):
#' \itemize{
#'   \item L1 - \code{\link{fit_lctm_full}}: Full 5-timepoint LCTM
#'   \item L2 - \code{\link{fit_lctm_partial}}: 2-4 timepoints with MCMC imputation
#' }
#'
#' Includes pre-trained v8 posterior parameters (\code{\link{lctm_v8_params}})
#' for offline prediction of new patients without re-running MCMC.
#'
#' @import methods
#' @importFrom methods setClass representation new is
#' @importFrom stats model.frame model.matrix predict coef sd density
#'   setNames reformulate qlogis plogis complete.cases ave dnorm as.formula
#' @importFrom utils packageVersion globalVariables
#' @importFrom dplyr filter mutate group_by summarise left_join bind_rows `%>%`
#' @importFrom tidyr pivot_wider pivot_longer
#' @author Anonymous \email{package@@example.com} (ORCID: 0009-0000-6575-5277)
#' @keywords package
#' @docType package
#' @name quick.lctm.escc-package
NULL