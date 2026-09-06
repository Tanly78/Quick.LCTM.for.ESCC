## ============================================================
##  quick.lctm.escc End-to-End User Workflow
## ============================================================
##
##  This script demonstrates the full user journey:
##    1.  Install / load the package
##    2.  Prepare patient data in long format
##    3.  Apply the pre-trained v8 LCTM
##    4.  Inspect the per-patient posterior probabilities
##    5.  Stratify patients into RC / VC groups
##    6.  (Optional) Apply Trajectory Bend Score (TBS) for T1->T2 change
##
##  Required input columns (per row = one measurement at one time point):
##    id          integer / character  patient identifier
##    timepoint   T1 | T2 | T4 | T5 | T6     (or T1_pre / T2_post / T4_3M / T5_6M / T6_1Y)
##    ALB         numeric   g/L           (albumin)
##    TC          numeric   mmol/L        (total cholesterol)
##    LYM         numeric   10^9/L        (lymphocyte count)
##    CRP         numeric   mg/L          (C-reactive protein)
##    SMI         numeric   cm^2/m^2      (skeletal muscle index, CT)
##    SFI         numeric   cm^2/m^2      (subcutaneous fat index, CT)
##    VFI         numeric   cm^2/m^2      (visceral fat index, CT)
##
##  At least 2 distinct timepoints per patient are required.
## ============================================================

library(quick.lctm.escc)

## ---- 1. Load patient data --------------------------------------
##  Replace this with your own file path / database query.
##  The CSV must have the 9 columns listed above.
patient_csv <- system.file("extdata", "example_lctm.csv",
                           package = "quick.lctm.escc")
patients    <- read.csv(patient_csv, stringsAsFactors = FALSE)
cat("Loaded", nrow(patients), "rows covering",
    length(unique(patients$id)), "patients and",
    length(unique(patients$timepoint)), "timepoints.\n")

## ---- 2. Classify every patient with the pre-trained v8 model ----
##  This call needs no MCMC: it applies the v8 posterior means
##  (intercept, slope, sigma, class-mixing proportions) bundled in
##  the package as `lctm_v8_params`. Wall-clock time < 1 s per
##  1000 patients.
pred <- predict_lctm7m(patients, use_pretrained = TRUE)
print(head(pred, 10))
#   id     p_VC  cluster
#   550193 0.80   VC
#   495730 0.02   RC
#   ...

## ---- 3. Stratify into Resilient vs Vulnerable groups ------------
n_rc <- sum(pred$cluster == "RC")
n_vc <- sum(pred$cluster == "VC")
cat(sprintf("Stratified: %d RC, %d VC\n", n_rc, n_vc))

## ---- 4. (Optional) Compute Trajectory Bend Score (TBS) ---------
##  TBS = standardized change in the 7-marker composite from T1 to T2.
##  Higher TBS = larger early improvement (favourable).
t1 <- patients[patients$timepoint %in% c("T1", "T1_pre"), ]
t2 <- patients[patients$timepoint %in% c("T2", "T2_post"), ]
tbs <- trajectory_bend_score(t1, t2, reference_sd = NULL)
print(head(tbs))
#   id     tbs  category
#   550193 0.42 strong_improvement
#   ...

## ---- 5. (Optional) Inspect the bundled v8 parameters ----------
data(lctm_v8_params)
cat("v8 training cohort n =", lctm_v8_params$training_n, "\n")
cat("Class-conditional means (RC / VC):\n")
print(round(lctm_v8_params$class_means, 4))
cat("Class-mixing proportions:\n")
print(round(lctm_v8_params$class_props, 4))
cat("Residual SD: ", round(lctm_v8_params$sigma, 4), "\n", sep = "")

## ---- 6. (Optional) Export results to CSV -----------------------
write.csv(pred, "patient_classification.csv", row.names = FALSE)
cat("Wrote patient_classification.csv\n")

## ============================================================
##  End-to-end workflow complete.
## ============================================================