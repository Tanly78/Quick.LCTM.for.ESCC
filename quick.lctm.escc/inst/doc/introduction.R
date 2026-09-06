## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(quick.lctm.escc)

## ----eval = FALSE-------------------------------------------------------------
# # install.packages("devtools")
# devtools::install_github("Tanly78/quick.lctm.escc")

## -----------------------------------------------------------------------------
library(quick.lctm.escc)
patients <- read.csv(
  system.file("extdata", "example_lctm.csv", package = "quick.lctm.escc"),
  stringsAsFactors = FALSE
)
str(patients)

## -----------------------------------------------------------------------------
pred <- predict_lctm7m(patients, use_pretrained = TRUE)
head(pred)

## -----------------------------------------------------------------------------
table(pred$cluster)

## -----------------------------------------------------------------------------
t1 <- patients[patients$timepoint %in% c("T1", "T1_pre"), ]
t2 <- patients[patients$timepoint %in% c("T2", "T2_post"), ]
tbs <- trajectory_bend_score(t1, t2)
head(tbs)

## ----eval = FALSE-------------------------------------------------------------
# write.csv(pred, "patient_classification.csv", row.names = FALSE)

