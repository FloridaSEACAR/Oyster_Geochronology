library(tidyverse)
library(magrittr)
library(rintcal)

# DEP AGREEMENT NO. CZ325 Deliverable 2b (Radiocarbon sample size for age and time-averaging determination)
# this script calculates a proper calibration curve for sending to OxCal. OxCal needs all fields filled in 
# so here we estimate 14C from fractions. Some of the data have 14C and fractions, so we use those fot fit a linear
# model, then use that model to impute all missing 14C values
# we also get the 14C standard deviations from the fitted values.

Curve = read.csv("HOBS_14C_CalCurve_PredictedVals_20230524.csv")
Combined14C_forFLCalCurve = read.csv("Combined14C_forFLCalCurve.csv")

Curve %<>% mutate(Year = year, D14C = NA, D14Csd = NA, F14C = estimate, F14Csd = esterror) %>% select(-contains(c("q", "est"))) %>% select(-year)

Combined14C_forFLCalCurve %>% select()
fit = lm(delta14C~F14C_calc, Combined14C_forFLCalCurve)
Curve$D14C = coef(fit)[1] + coef(fit)[2] * Curve$F14C #F14C.D14C(Curve$F14C, t = Curve$Year)

Curve$D14Csd = coef(fit)[2] * Curve$F14Csd
names(Curve)[1] = "#year" # force a comment in the first row
write.table(Curve, "/Users/johnhandley/Documents/PaleontologyResearch/Oysters/HOBS/Contract work/Comparisons/Kowalewski_NP/Curves/HOBS_14C_CalCurve_PredictedVals_20230524.14c", sep = "\t", row.names = FALSE, quote = FALSE)
curvename = "HOBS_14C_CalCurve_PredictedVals_20230524"
curvepath = "/Users/johnhandley/Documents/PaleontologyResearch/Oysters/HOBS/Contract work/Comparisons/Kowalewski_NP/Curves/HOBS_14C_CalCurve_PredictedVals_20230524.14c"
