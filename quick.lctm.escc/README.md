# quick.lctm.escc

Latent Class Trajectory Modeling for 7 Nutritional Markers (ALB, TC, LYM, CRP, SMI, SFI, VFI).

## What it does

Fits a Bayesian Latent Class Trajectory Model (LCTM) on the 7 nutritional /
inflammatory markers at 5 time points and classifies patients into
Resilient Class (RC) vs Vulnerable Class (VC).

The package uses only the 7 markers.

## Installation

```r
# install.packages("devtools")
devtools::install_github("Tanly78/quick.lctm.escc")
```

## Quick start (7-marker only)

```r
library(quick.lctm.escc)

# Input: long-format data.frame with columns
#   id, timepoint (T1/T2/T4/T5/T6 or T1_pre/T2_post/...),
#   ALB, TC, LYM, CRP, SMI, SFI, VFI
str(new_patient_data)

# Apply pre-trained v8 parameters to classify new patients
# (no MCMC required; classification in <1 second)
pred <- predict_lctm7m(new_patient_data, use_pretrained = TRUE)
head(pred)
#   id     p_VC cluster
#  550193 0.80  VC
#  ...

# Fit a fresh LCTM on a training cohort with all 5 time points
fit <- fit_lctm_full(my_5tp_data, K = 2, n_iter = 2000, seed = 42)
table(fit@cluster)   # RC vs VC

# Or fit a partial LCTM when some time points are missing
fit_p <- fit_lctm_partial(my_partial_data, K = 2, seed = 42)
```

## Included pre-trained parameters

`lctm_v8_params` bundles the v8 posterior means from the n=219 ESCC nCIT
training cohort (R-hat < 1.01, ESS > 1400). Use `use_pretrained = TRUE` in
`predict_lctm7m()` to apply these without re-running MCMC.

```r
data(lctm_v8_params)
names(lctm_v8_params)
# "ref_mean" "ref_sd" "negative" "times" "timepoint_lab"
# "class_means" "class_props" "sigma" "composite_train"
# "p_VC_train" "version" "training_n" "source_rds"
# "source_csv" "source_commit"
```

## See also

- `vignette("L1_full_analysis")` -- full 5-timepoint LCTM workflow
- `?predict_lctm7m` -- classification of new patients
- `?lctm_v8_params` -- bundled posterior parameters