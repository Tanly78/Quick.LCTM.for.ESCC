# Suppress R CMD check NOTEs about non-standard evaluation variables
# in NSE-heavy functions (reshape_to_wide uses dplyr/tidyr NSE).
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "id", "timepoint",
    "ALB_z", "TC_z", "LYM_z", "CRP_z_flip",
    "SMI_z", "SFI_z", "VFI_z"
  ))
}