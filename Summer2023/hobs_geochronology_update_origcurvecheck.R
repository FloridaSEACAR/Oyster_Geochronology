#Update all of the HOBS geochronology manuscript analyses and figures: 
#   - using new calibration curve, 
#   - including additional radiocarbon dates from CZ325, 
#   - using scripts put together by John Handley for project CZ325,
#   - re-incorporating Lemon Bay data now that we have updated HP radiocarbon analyses of live-caught oysters from LB

library(sf)
library(data.table)
library(ggplot2)
library(tidyverse)
library(maps)
library(ggspatial)
library(rnaturalearth)
library(ggrepel)
library(lmodel2)
library(kableExtra)
library(rmarkdown)
library(gtable)
library(grid)
library(cowplot)
library(lemon)
library(egg)
library(ggtext)
library(colorspace)
library(brms)
library(patchwork)
library(ragg)
library(openxlsx)
library(Hmisc)
library(tidyverse)
library(magrittr)
library(reshape2)
library(HDInterval)
library(readxl)
library(oxcAAR)
library(rintcal)
library(lme4)
library(cmdstanr)
library(bayesplot)
library(ggpattern)
library(extrafont)
library(ggnewscale)
setOxcalExecutablePath("C:/Program Files/OxCal/bin/OxCalWin.exe")
options(scipen = 999) #disable scientific notation


#Load data
# predictedCurve <- fread(here::here("CZ325_Sampling_task_2/HOBS_14C_CalCurve_PredictedVals_20230524_14c.csv"))
# predictedCurve <- fread(here::here("Summer2023/HOBS_14C_CalCurve_PredictedVals_2023-07-20.csv"))
# predictedCurve <- fread(here::here("Summer2023/HOBS_14C_CalCurve_PredictedVals_2023-07-24.csv"))
predictedCurve <- fread(here::here("Summer2023/HOBS_14C_CalCurve_PredictedVals_2023-07-28.csv"))
predictedCurve <- janitor::clean_names(predictedCurve)

masterSourceFile <- "HOBS14CMaster.xlsx"

hobsf14c <- read.xlsx(here::here(paste0("Summer2023/", masterSourceFile)))
hobsf14c <- janitor::clean_names(hobsf14c)
setDT(hobsf14c)
hobsf14c[, `:=` (lat = as.numeric(lat), 
                 lon = as.numeric(lon), 
                 collection_date = as.Date(as.integer(hobsf14c$collection_date), origin = "1899-12-30"), 
                 station = as.integer(station),
                 depth_int = fcase(depth == "Surface", 0,
                                   depth %in% c("15-25cm", "15-30cm"), 1,
                                   default = 2))]
# hobsf14c[is.na(collection_date), collection_date := as_date(paste0(max(predictedCurve$number_year), "-12-31"))] #LB and BH replacement specimens were collected after the latest year on the calibration curve, so I am replacing their actual collection dates with the latest possible calibration curve date

hobslocs <- c("Big Hickory", "Goose Island/East Cove", "Guana River", "Hendry Creek/Mullock Creek", "Jack Island", "Lone Cabbage", "Little St. George Island", "Matanzas River", "New Pass", "Pellicer Creek", "Lemon Bay")

#remove suspect values
hobsf14c <- hobsf14c[remove == 0, ]


#New Calibration curve----------------------------------------------------------
#Code originally copied from "DR3_QA.R" script

comb14c <- fread(here::here("CalCurve/14Cdata_FL_Gulf_Caribbean_20230329.csv"))
comb14c[str_detect(study, "Durham"), `:=` (F14C = hobsf14c[analysis == "HP" & cat_no == samp, f14c],
                                           F14C_sd = hobsf14c[analysis == "HP" & cat_no == samp, f14c_sd],
                                           dataSource = masterSourceFile), by = samp]

comb14c[lon > 1, lon := lon * -1]
comb14c_sf <- st_as_sf(comb14c[!is.na(lon) & !is.na(lat), ], coords = c("lon", "lat"), crs = 4326)
# mapview::mapview(comb14c_sf)
# 
# colnames(comb14c)
# 
# ggplot(comb14c[!is.na(delta14C) & !is.na(year), ]) +
#   geom_point(aes(x = year, y = delta14C, color = study, shape = source)) +
#   theme_bw()

comb14c[, ind := seq(1, nrow(comb14c))]

comb14c[!is.na(`14Cage`), `:=` (F14C_calc = age.F14C(`14Cage`, sdev = `14Cage_sd`)[,1],
                                F14C_calc_sd = age.F14C(`14Cage`, sdev = `14Cage_sd`)[,2])]

comb14c[!is.na(pMC), `:=` (C14age_calc = pMC.age(pMC, sdev = pMC_sd)[[1]][1], 
                           C14age_sd_calc = pMC.age(pMC, sdev = pMC_sd)[[1]][2]), by = ind] #comb14c[!is.na(pMC), ind]

comb14c[, year_bp := 2022 - year]

comb14c[!is.na(delta14C) & !is.na(year_bp) & is.na(F14C_calc), F14C_calc := D14C.F14C(delta14C, year_bp)]

comb14c[!is.na(C14age_calc) & is.na(F14C_calc), `:=` (F14C_calc = age.F14C(C14age_calc, sdev = C14age_sd_calc)[,1],
                                                      F14C_calc_sd = age.F14C(C14age_calc, sdev = C14age_sd_calc)[,2])]
comb14c[!is.na(year) & !is.na(F14C) & is.na(F14C_calc), `:=` (F14C_calc = F14C, F14C_calc_sd = F14C_sd)]

comb14c[!is.na(year) & !is.na(F14C) & !is.na(F14C_calc), delta14C_calc := F14C.D14C(F14C_calc, year_bp)]

# ggplot(comb14c[!is.na(year) & !is.na(F14C_calc), ]) +
#   geom_point(aes(x = year, y = F14C_calc, color = study, shape = source)) +
#   geom_point(data = comb14c[!is.na(year) & !is.na(F14C_calc) & str_detect(species, "Mercenaria sp"), ], aes(x = year, y = F14C_calc), shape = 1, color = "black", stroke = 1, size = 1.5) +
#   theme_bw()


comb14c[!is.na(year), decade := fcase(year < 1900, "pre-1900",
                                      year >= 1900 & year < 1910, "1900s",
                                      year >= 1910 & year < 1920, "1910s",
                                      year >= 1920 & year < 1930, "1920s",
                                      year >= 1930 & year < 1940, "1930s",
                                      year >= 1940 & year < 1950, "1940s",
                                      year >= 1950 & year < 1960, "1950s",
                                      year >= 1960 & year < 1970, "1960s",
                                      year >= 1970 & year < 1980, "1970s",
                                      year >= 1980 & year < 1990, "1980s",
                                      year >= 1990 & year < 2000, "1990s",
                                      year >= 2000 & year < 2010, "2000s",
                                      year >= 2010 & year < 2020, "2010s",
                                      year >= 2020 & year < 2030, "2020s")]

comb14c_sf2 <- st_as_sf(comb14c[!is.na(year) & !is.na(lon) & !is.na(lat), ], coords = c("lon", "lat"), crs = 4326)

currents <- st_read(here::here("CalCurve/MajorOceanCurrents/Major_Ocean_Currents.shp"))
currents <- st_make_valid(currents)
currents <- st_transform(currents, crs = 4326)
# mapview::mapview(currents) +
#   mapview::mapview(subset(comb14c_sf2, !is.na(comb14c_sf2$year)), zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), label = "sta", layer.name = "F14C_calc", alpha = 0, alpha.regions = 0, legend = TRUE, popup = FALSE, highlight = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year < 1900), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "pre-1900", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1900 & year < 1910), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1900s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1910 & year < 1920), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1910s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1920 & year < 1930), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1920s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1930 & year < 1940), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1930s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1940 & year < 1950), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1940s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1950 & year < 1960), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1950s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1960 & year < 1970), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1960s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1970 & year < 1980), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1970s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1980 & year < 1990), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1980s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1990 & year < 2000), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1990s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 2000 & year < 2010), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "2000s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 2010 & year < 2020), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "2010s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 2020 & year < 2030), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "2020s", legend = FALSE)


#Remove points north of the southern Bahamas and east of the Caribbean Sea from the data, based on the ocean currents
comb14c2 <- comb14c[!is.na(year) & lon < -77.87 | lat < 26.43 & lon < -57, ]
comb14c2_sf <- st_as_sf(comb14c2[!is.na(year) & !is.na(lon) & !is.na(lat), ], coords = c("lon", "lat"), crs = 4326)
minF14C <- plyr::round_any(range(comb14c2_sf$F14C_calc)[1], 0.05, f = floor)
maxF14C <- plyr::round_any(range(comb14c2_sf$F14C_calc)[2], 0.05, f = ceiling)

# mapview::mapview(currents) +
#   mapview::mapview(subset(comb14c2_sf, !is.na(comb14c2_sf$year)), zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), label = "sta", layer.name = "F14C_calc", alpha = 0, alpha.regions = 0, legend = TRUE, popup = FALSE, highlight = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year < 1900), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "pre-1900", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1900 & year < 1910), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1900s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1910 & year < 1920), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1910s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1920 & year < 1930), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1920s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1930 & year < 1940), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1930s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1940 & year < 1950), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1940s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1950 & year < 1960), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1950s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1960 & year < 1970), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1960s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1970 & year < 1980), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1970s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1980 & year < 1990), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1980s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1990 & year < 2000), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1990s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 2000 & year < 2010), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "2000s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 2010 & year < 2020), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "2010s", legend = FALSE) +
#   mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 2020 & year < 2030), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "2020s", legend = FALSE)

comb14c2_sf$current <- st_nearest_feature(comb14c2_sf, currents)
comb14c2_sf$current2 <- lapply(comb14c2_sf$current, function(x) currents$NAME[x])
comb14c2_sf$current2 <- unlist(comb14c2_sf$current2)
comb14c2_sf$current2[which(is.na(comb14c2_sf$current2))] <- "Gulf Stream"


#Calibration curve model--------------------------------------------------------
gamtest7 <- brm(bf(F14C_calc ~ s(year, k = 20) + current + source), data = subset(comb14c2_sf, !is.na(comb14c2_sf$year)), family = gaussian(), 
                cores = 4, seed = 456, chains = 4, iter = 5000, warmup = 1500, thin = 3, 
                control = list(adapt_delta = 0.99, max_treedepth = 15), backend = "cmdstanr", threads = threading(2), file = here::here("gamtest7b.rds")) #Not sure what happened to "gamtest7.rds" but it started plotting weird so I re-ran the model and saved it as "gamtest7b.rds" and the plotting issues were fixed

gamtest7 <- update(gamtest7, newdata = subset(comb14c2_sf, !is.na(comb14c2_sf$year)), family = gaussian(), cores = 4, seed = 456, 
                   chains = 4, iter = 5000, warmup = 1500, thin = 3, control = list(adapt_delta = 0.99, max_treedepth = 15), 
                   backend = "cmdstanr", threads = threading(2), file = here::here(paste0("gamtest7b_update", Sys.Date(), ".rds")))

#Note that the HP analyses specimens from our study are included in the calibration data set, but not the LP-only analyses specimens
comb14c2_sf_mod <- comb14c2_sf
setDT(comb14c2_sf_mod)
# 
# #find the values that fell within the IQR for their year
# for(y in sort(unique(round(comb14c2_sf_mod$year)))){
#   comb14c2_sf_mod[year == y, `:=` (nobs = .N,
#                                    an05ci = quantile(F14C_calc, 0.25),
#                                    an95ci = quantile(F14C_calc, 0.75))]
# }
# 
# #find the values that fell within the IQR for their year and the year before and after, combined
# for(y in sort(unique(round(comb14c2_sf_mod$year)))){
#   comb14c2_sf_mod[year %in% c(y - 1, y, y + 1), `:=` (nobs_3yr = .N,
#                                                       an05ci_3yr = quantile(F14C_calc, 0.25),
#                                                       an95ci_3yr = quantile(F14C_calc, 0.75))]
# }
# 
# comb14c2_sf_mod[nobs == 1 | (F14C_calc > an05ci & F14C_calc < an95ci), inIQR_1yr := TRUE][is.na(inIQR_1yr), inIQR_1yr := FALSE]
# comb14c2_sf_mod[nobs_3yr == 1 | (F14C_calc > an05ci_3yr & F14C_calc < an95ci_3yr), inIQR_3yr := TRUE][is.na(inIQR_3yr), inIQR_3yr := FALSE]
# 
# #try the model again, this time only using observations that fell within the IQR for its year
# gamtest7c <- brm(bf(F14C_calc ~ s(year, k = 20) + current + source), data = comb14c2_sf_mod[inIQR_1yr == TRUE, ], family = gaussian(), 
#                  cores = 4, seed = 456, chains = 4, iter = 5000, warmup = 1500, thin = 3, 
#                  control = list(adapt_delta = 0.99, max_treedepth = 15), backend = "cmdstanr", threads = threading(2), file = here::here("Summer2023/gamtest7c.rds"))
# 
# #try dropping current and source variables since their effects look small
# gamtest7d <- brm(bf(F14C_calc ~ s(year, k = 20)), data = comb14c2_sf_mod[inIQR_1yr == TRUE, ], family = gaussian(), 
#                  cores = 4, seed = 456, chains = 4, iter = 5000, warmup = 1500, thin = 3, 
#                  control = list(adapt_delta = 0.99, max_treedepth = 15), backend = "cmdstanr", threads = threading(2), file = here::here("Summer2023/gamtest7d.rds"))
# 
# loo(gamtest7c, gamtest7d) #model 7c is a little better than 7d
# 
# #try the model again, this time only using observations that fell within the IQR for its year, the year before, and the year after
# gamtest7e <- brm(bf(F14C_calc ~ s(year, k = 20) + current + source), data = comb14c2_sf_mod[inIQR_3yr == TRUE, ], family = gaussian(), 
#                  cores = 4, seed = 456, chains = 4, iter = 5000, warmup = 1500, thin = 3, 
#                  control = list(adapt_delta = 0.99, max_treedepth = 15), backend = "cmdstanr", threads = threading(2), file = here::here("Summer2023/gamtest7e.rds"))
# 
# #try dropping current and source variables since their effects look small
# gamtest7f <- brm(bf(F14C_calc ~ s(year, k = 20)), data = comb14c2_sf_mod[inIQR_3yr == TRUE, ], family = gaussian(), 
#                  cores = 4, seed = 456, chains = 4, iter = 5000, warmup = 1500, thin = 3, 
#                  control = list(adapt_delta = 0.99, max_treedepth = 15), backend = "cmdstanr", threads = threading(2), file = here::here("Summer2023/gamtest7f.rds"))
# 
# loo(gamtest7e, gamtest7f) #model 7e is a little better than 7f
# 
# 
# newdat <- data.table(year = seq(min(gamtest7c$data$year), 2022, 0.5),
#                      current = 65,
#                      source = "otolith")
# predvals <- posterior_predict(gamtest7c, newdata = newdat, ndraws = 1000, allow_new_levels = TRUE)
# newcurvepreds <- newdat[, .(year)]
# for(i in seq(1, nrow(newcurvepreds))){
#   newcurvepreds[i, `:=` (f14c_mn = mean(predvals[, i]),
#                          f14c_sd = sd(predvals[, i]))]
# }
# 
# newdat2 <- data.table(year = seq(min(gamtest7e$data$year), 2022, 0.5),
#                       current = 65,
#                       source = "otolith")
# predvals <- posterior_predict(gamtest7e, newdata = newdat2, ndraws = 1000, allow_new_levels = TRUE)
# newcurvepreds2 <- newdat2[, .(year)]
# for(i in seq(1, nrow(newcurvepreds2))){
#   newcurvepreds2[i, `:=` (f14c_mn = mean(predvals[, i]),
#                           f14c_sd = sd(predvals[, i]))]
# }
# 
# #Decided to go with the 3yr running IQR criteria
# fwrite(newcurvepreds2, here::here(paste0("Summer2023/HOBS_14C_CalCurve_PredictedVals_", Sys.Date(), ".csv")))

#Talked with John Handley about these new IQR-based curves on 7/24/2023 and he persuaded me that the unknown cost of reducing 
#variability in the data that the curve is based on is probably greater than the cost of including some data points that don't 
#represent a marine condition, so I went back to our original "new" curve.
set.seed(888)
newdat3 <- data.table(year = seq(min(gamtest7$data$year), 2022, 0.5),
                      current = 65,
                      source = "shell")
predvals <- posterior_predict(gamtest7, newdata = newdat3, ndraws = 1000, allow_new_levels = TRUE)
newcurvepreds3 <- newdat3[, .(year)]
for(i in seq(1, nrow(newcurvepreds3))){
  newcurvepreds3[i, `:=` (f14c_mn = mean(predvals[, i]),
                          f14c_sd = sd(predvals[, i]))]
}

newdat3_oto <- data.table(year = seq(min(gamtest7$data$year), 2022, 0.5),
                          current = 65,
                          source = "otolith")
predvals_oto <- posterior_predict(gamtest7, newdata = newdat3_oto, ndraws = 1000, allow_new_levels = TRUE)
newcurvepreds3_oto <- newdat3_oto[, .(year)]
for(i in seq(1, nrow(newcurvepreds3_oto))){
  newcurvepreds3_oto[i, `:=` (f14c_mn = mean(predvals_oto[, i]),
                              f14c_sd = sd(predvals_oto[, i]))]
}

newdat3_cor <- data.table(year = seq(min(gamtest7$data$year), 2022, 0.5),
                          current = 65,
                          source = "otolith")
predvals_cor <- posterior_predict(gamtest7, newdata = newdat3_cor, ndraws = 1000, allow_new_levels = TRUE)
newcurvepreds3_cor <- newdat3_cor[, .(year)]
for(i in seq(1, nrow(newcurvepreds3_cor))){
  newcurvepreds3_cor[i, `:=` (f14c_mn = mean(predvals_cor[, i]),
                              f14c_sd = sd(predvals_cor[, i]))]
}

fwrite(newcurvepreds3, here::here(paste0("Summer2023/HOBS_14C_CalCurve_PredictedVals_", Sys.Date(), ".csv")))

f14c_quan <- fread(here::here("CalCurve/QuanBombCurve_20230112.csv"), sep = "|")

studycolslist <- qualitative_hcl(28)
studycols <- setNames(studycolslist, sort(unique(c(comb14c2_sf$study, f14c_quan$study))))
f14c_quan[, study := factor(study, levels = names(studycols))]
comb14c2_sf_mod[, study := factor(study, levels = names(studycols))]
f14c_quan[, source := factor(source, levels = c("coral", "otolith", "shell", "surface water"))]
comb14c2_sf_mod[, source := factor(source, levels = c("coral", "otolith", "shell", "surface water"))]

#New curve with Marine20 values for years <= 1950
ncm20_dat <- fread(here::here("Summer2023/NewCurve_wMarine20_vals.csv"))
ncm20_dat[, `:=` (f14c_calc = age.F14C(`14Cage`, sdev = `14Cage_sd`)[,1],
                  f14c_calc_sd = age.F14C(`14Cage`, sdev = `14Cage_sd`)[,2])]

#Save tab-delimited file that can be appended to a .14c file to update the curve (e.g., NNewcurveto10000calBP_Marine20.14c)
write.table(formerge[, .(`#Year` = format(round(`#year`, digits = 1), nsmall = 1), 
                         D14C = format(round(D14C, digits = 1), nsmall = 1),
                         `1sigma` = format(round(D14Csd, digits = 1), nsmall = 1),
                         F14C = format(round(F14C, digits = 4), nsmall = 4),
                         `1sigma` = format(round(F14Csd, digits = 4), nsmall = 4))],
            here::here("Summer2023/HOBS_14C_CalCurve_PredictedVals_20230728_14c_formerge.txt"), 
            sep = '\t', row.names = FALSE, quote = FALSE)

#Update predictedCurve object with Marine20 values
predictedCurve <- predictedCurve[year > 1950, ]
m20 <- ncm20_dat[!is.na(f14c_calc), .(year = Year_ce, f14c_mn = f14c_calc, f14c_sd = f14c_calc_sd)]
predictedCurve <- rbind(m20, predictedCurve)

# qcurve <- ggplot(data = f14c_quan, aes(x = Year, y = F14C, color = study, shape = source)) +
#   geom_line(aes(x = Year, y = F14C), lty = "solid", color = "firebrick", lwd = 0.75, inherit.aes = FALSE) +
#   geom_point() +
#   theme_bw() +
#   labs(x = NULL, y = "F14C", color = "Study", shape = "Source") +
#   coord_cartesian(xlim = c(1950, 2023), ylim = c(0.75, 1.2))


#Calculate reference values for dead C calculations (refval calculations described in "Weighted mean calculation (002).doc" from Quan Hua)
# refval <- hobsf14c[locality == "Alligator Harbor", .(f14c_ref = sum(f14c/f14c_sd^2)/sum(1/f14c_sd^2), 
#                                                      f14c_ref_sd = max(c(sqrt(1/sum(1/f14c_sd^2)),
#                                                                          sqrt((.N * sum(((f14c - sum(f14c/f14c_sd^2)/sum(1/f14c_sd^2))^2/f14c_sd^2)))/((.N - 1) * sum(1/f14c_sd^2))))))]

yrs <- unique(sort(unlist(lapply(unique(year(hobsf14c$collection_date)), function(x) seq(x, x - 2, -0.5)))))
yrs <- unique(sort(unlist(lapply(yrs, function(x) ifelse(x <= 1950, plyr::round_any(x, 10), x))))) #round pre-1950 yrs to nearest decade
refvals <- predictedCurve[year %in% yrs, .(year, f14c_mn, f14c_sd)]
refvals[, `:=` (f14c_high = f14c_mn + f14c_sd, f14c_low = f14c_mn - f14c_sd)]
# refval_orig <- data.table(f14c_mn = 1.01059284984934,
#                           f14c_sd = 0.00170241365804406,
#                           f14c_high = 1.01059284984934 + 0.00170241365804406,
#                           f14c_low = 1.01059284984934 - 0.00170241365804406)

#Calculate dead carbon values for each locality & year
# deadC <- hobsf14c[locality != "Alligator Harbor" & analysis == "HP" & live_or_dead == "Live-caught", .(n = .N,
#                                                                                                        f14c_mn = sum(f14c/f14c_sd^2)/sum(1/f14c_sd^2),
#                                                                                                        f14c_sd = max(c(sqrt(1/sum(1/f14c_sd^2)),
#                                                                                                                        sqrt((.N * sum((f14c - sum(f14c/f14c_sd^2)/sum(1/f14c_sd^2))^2/f14c_sd^2))/((.N - 1) * sum(1/f14c_sd^2)))))), by = list(year(collection_date), locality)]

deadC <- hobsf14c[analysis == "HP" & live_or_dead == "Live-caught", .(n = .N,
                                                                      f14c_mn = sum(f14c/f14c_sd^2)/sum(1/f14c_sd^2),
                                                                      f14c_sd = max(c(sqrt(1/sum(1/f14c_sd^2)),
                                                                                      sqrt((.N * sum((f14c - sum(f14c/f14c_sd^2)/sum(1/f14c_sd^2))^2/f14c_sd^2))/((.N - 1) * sum(1/f14c_sd^2)))))), by = list(year(collection_date), locality)]
setnames(deadC, c("locality", "year"), c("loc", "yr"))

#reconstruct original Quan corrections
deadC <- hobsf14c[analysis == "HP" & live_or_dead == "Live-caught", .(n = .N,
                                                                      f14c_mn = sum(f14c/f14c_sd^2)/sum(1/f14c_sd^2),
                                                                      f14c_sd = max(c(sqrt(1/sum(1/f14c_sd^2)),
                                                                                      sqrt((.N * sum((f14c - sum(f14c/f14c_sd^2)/sum(1/f14c_sd^2))^2/f14c_sd^2))/((.N - 1) * sum(1/f14c_sd^2)))))), by = locality]

setnames(deadC, c("locality"), c("loc"))

# deadC[, `:=` (loc = locality, yr = year)]
deadC[n == 1, f14c_sd := hobsf14c[year(collection_date) == yr & locality == loc, f14c_sd], by = list(yr, loc)]
deadC[, yr2 := ifelse(yr <= 1950, plyr::round_any(yr, 10), yr)] #round years <= 1950 to nearest decade to match Marine20 resolution
deadC[, `:=` (ref_mn = refvals[year %in% seq(yr2, yr2 - 2, -0.5), sum(f14c_mn/f14c_sd^2)/sum(1/f14c_sd^2)],
              ref_sd = refvals[year %in% seq(yr2, yr2 - 2, -0.5), max(c(sqrt(1/sum(1/f14c_sd^2)),
                                                                      sqrt((.N * sum((f14c_mn - sum(f14c_mn/f14c_sd^2)/sum(1/f14c_sd^2))^2/f14c_sd^2))/((.N - 1) * sum(1/f14c_sd^2))))[which(!is.nan(c(sqrt(1/sum(1/f14c_sd^2)),
                                                                                                                                                                                                       sqrt((.N * sum((f14c_mn - sum(f14c_mn/f14c_sd^2)/sum(1/f14c_sd^2))^2/f14c_sd^2))/((.N - 1) * sum(1/f14c_sd^2))))))])],
              ref_mn_l = refvals[year %in% seq(yr2, yr2 - 2, -0.5), sum(f14c_low/f14c_sd^2)/sum(1/f14c_sd^2)],
              ref_sd_l = refvals[year %in% seq(yr2, yr2 - 2, -0.5), max(c(sqrt(1/sum(1/f14c_sd^2)),
                                                                      sqrt((.N * sum((f14c_low - sum(f14c_low/f14c_sd^2)/sum(1/f14c_sd^2))^2/f14c_sd^2))/((.N - 1) * sum(1/f14c_sd^2))))[which(!is.nan(c(sqrt(1/sum(1/f14c_sd^2)),
                                                                                                                                                                                                         sqrt((.N * sum((f14c_low - sum(f14c_low/f14c_sd^2)/sum(1/f14c_sd^2))^2/f14c_sd^2))/((.N - 1) * sum(1/f14c_sd^2))))))])],
              ref_mn_h = refvals[year %in% seq(yr2, yr2 - 2, -0.5), sum(f14c_high/f14c_sd^2)/sum(1/f14c_sd^2)],
              ref_sd_h = refvals[year %in% seq(yr2, yr2 - 2, -0.5), max(c(sqrt(1/sum(1/f14c_sd^2)),
                                                                        sqrt((.N * sum((f14c_high - sum(f14c_high/f14c_sd^2)/sum(1/f14c_sd^2))^2/f14c_sd^2))/((.N - 1) * sum(1/f14c_sd^2))))[which(!is.nan(c(sqrt(1/sum(1/f14c_sd^2)),
                                                                                                                                                                                                             sqrt((.N * sum((f14c_high - sum(f14c_high/f14c_sd^2)/sum(1/f14c_sd^2))^2/f14c_sd^2))/((.N - 1) * sum(1/f14c_sd^2))))))])]), by = yr2]
#reconstruct original Quan corrections
deadC[n == 1, f14c_sd := hobsf14c[locality == loc, f14c_sd], by = loc]
deadC[, `:=` (ref_mn_orig = refval_orig$f14c_mn,
              ref_sd_orig = refval_orig$f14c_sd,
              ref_mn_l_orig = refval_orig$f14c_low,
              ref_sd_l_orig = refval_orig$f14c_sd,
              ref_mn_h_orig = refval_orig$f14c_high,
              ref_sd_h_orig = refval_orig$f14c_sd)]


dc18 <- deadC[yr %in% c(2019, 2020, 2022), .(n = .N,
                                             yr = 2018,
                                             yr2 = 2018,
                                             f14c_mn = sum(f14c_mn/f14c_sd^2)/sum(1/f14c_sd^2),
                                             f14c_sd = max(c(sqrt(1/sum(1/f14c_sd^2)),
                                                             sqrt((.N * sum((f14c_mn - sum(f14c_mn/f14c_sd^2)/sum(1/f14c_sd^2))^2/f14c_sd^2))/((.N - 1) * sum(1/f14c_sd^2))))),
                                             ref_mn = sum(ref_mn/ref_sd^2)/sum(1/ref_sd^2),
                                             ref_sd = max(c(sqrt(1/sum(1/ref_sd^2)),
                                                            sqrt((.N * sum((ref_mn - sum(ref_mn/ref_sd^2)/sum(1/ref_sd^2))^2/ref_sd^2))/((.N - 1) * sum(1/ref_sd^2))))),
                                             ref_mn_l = sum(ref_mn_l/ref_sd_l^2)/sum(1/ref_sd_l^2),
                                             ref_sd_l = max(c(sqrt(1/sum(1/ref_sd_l^2)),
                                                            sqrt((.N * sum((ref_mn_l - sum(ref_mn_l/ref_sd_l^2)/sum(1/ref_sd_l^2))^2/ref_sd_l^2))/((.N - 1) * sum(1/ref_sd_l^2))))),
                                             ref_mn_h = sum(ref_mn_h/ref_sd_h^2)/sum(1/ref_sd_h^2),
                                             ref_sd_h = max(c(sqrt(1/sum(1/ref_sd_h^2)),
                                                            sqrt((.N * sum((ref_mn_h - sum(ref_mn_h/ref_sd_h^2)/sum(1/ref_sd_h^2))^2/ref_sd_h^2))/((.N - 1) * sum(1/ref_sd_h^2)))))), by = loc]
setnames(dc18, c("yr", "loc"), c("year", "locality"))
dc18[n == 1, `:=` (f14c_sd = deadC[loc == locality, unique(f14c_sd)],
                   ref_sd = deadC[loc == locality, unique(ref_sd)],
                   ref_sd_l = deadC[loc == locality, unique(ref_sd_l)],
                   ref_sd_h = deadC[loc == locality, unique(ref_sd_h)]), by = locality][, n := 2]
deadC <- rbind(deadC, dc18[, .(yr = year, yr2, loc = locality, n, f14c_mn, f14c_sd, ref_mn, ref_sd, ref_mn_l, ref_sd_l, ref_mn_h, ref_sd_h)])
deadC[, `:=` (deadc = 1 - f14c_mn/ref_mn,
              deadc_sd = sqrt((f14c_sd/ref_mn)^2 + (ref_sd * f14c_mn/ref_mn^2)^2),
              deadc_l = 1 - f14c_mn/ref_mn_l,
              deadc_sd_l = sqrt((f14c_sd/ref_mn_l)^2 + (ref_sd_l * f14c_mn/ref_mn_l^2)^2),
              deadc_h = 1 - f14c_mn/ref_mn_h,
              deadc_sd_h = sqrt((f14c_sd/ref_mn_h)^2 + (ref_sd_h * f14c_mn/ref_mn_h^2)^2)), by = row.names(deadC)]
deadC[yr %in% c(1932, 1938, 1940, 1952, 1979), loc_rectest := fcase(yr == 1932, "Big Hickory",
                                                                    yr == 1938, "Little St. George Island",
                                                                    yr == 1940, "Lone Cabbage",
                                                                    yr == 1952, "Lone Cabbage",
                                                                    yr == 1979, "Lone Cabbage")]
deadC[is.na(loc_rectest), loc_rectest := loc]

#reconstruct original Quan corrections
deadC[, `:=` (deadc = 1 - f14c_mn/ref_mn_orig,
              deadc_sd = sqrt((f14c_sd/ref_mn_orig)^2 + (ref_sd_orig * f14c_mn/ref_mn_orig^2)^2),
              deadc_l = 1 - f14c_mn/ref_mn_l_orig,
              deadc_sd_l = sqrt((f14c_sd/ref_mn_l_orig)^2 + (ref_sd_l_orig * f14c_mn/ref_mn_l_orig^2)^2),
              deadc_h = 1 - f14c_mn/ref_mn_h_orig,
              deadc_sd_h = sqrt((f14c_sd/ref_mn_h_orig)^2 + (ref_sd_h_orig * f14c_mn/ref_mn_h_orig^2)^2)), by = row.names(deadC)]


# setnames(deadC, "locality", "loc")


# ##Aside - build table comparing weighted means and uncertainties with original values from Quan Hua
# refval2 <- hobsf14c[locality == "Alligator Harbor", .(n = .N,
#                                                       f14c_mn = sum(f14c/f14c_sd^2)/sum(1/f14c_sd^2), 
#                                                       f14c_sd = c(sqrt(1/sum(1/f14c_sd^2)),
#                                                                   sqrt((.N * sum(((f14c - sum(f14c/f14c_sd^2)/sum(1/f14c_sd^2))^2/f14c_sd^2)))/((.N - 1) * sum(1/f14c_sd^2)))))]
# refval2[, locality := "Alligator Harbor (reference)"]
# 
# deadC2 <- hobsf14c[locality != "Alligator Harbor" & analysis == "HP" & live_or_dead == "Live-caught", .(n = .N,
#                                                                                                         f14c_mn = sum(f14c/f14c_sd^2)/sum(1/f14c_sd^2),
#                                                                                                         f14c_sd = c(sqrt(1/sum(1/f14c_sd^2)),
#                                                                                                                     sqrt((.N * sum((f14c - refval$f14c_ref)^2/f14c_sd^2))/((.N - 1) * sum(1/f14c_sd^2))))), by = locality]
# deadC2 <- rbind(deadC2, refval2)
# deadC2[, `:=` (deadc = 1 - f14c_mn/refval$f14c_ref,
#                deadc_sd = sqrt((f14c_sd/refval$f14c_ref)^2 + (refval$f14c_ref_sd * f14c_mn/refval$f14c_ref^2)^2)), by = row.names(deadC2)]
# deadC2[, sd_method := rep(c("Error of the mean", "Standard deviation"), 15)]
# deadC2 <- deadC2[n == 2, ]
# 
# qhvals <- data.table(locality = unique(deadC2$locality))
# qhvals[, orig_f14c := fcase(str_detect(locality, "Alligator"), 1.01059284984934,
#                             str_detect(locality, "New"), 0.997242957605182,
#                             str_detect(locality, "Hendry"), 0.976840620995996,
#                             str_detect(locality, "Goose"), 1.01203166770136,
#                             str_detect(locality, "Lone"), 0.922489182101228,
#                             str_detect(locality, "Little"), 1.00259065847587,
#                             str_detect(locality, "Guana"), 0.997442866257221,
#                             str_detect(locality, "Pellicer"), 1.01088972848189,
#                             str_detect(locality, "Matanzas"), 1.02288020544133,
#                             str_detect(locality, "Jack"), 1.00820441264155)]
# 
# qhvals[, orig_f14c_sd := fcase(str_detect(locality, "Alligator"), 0.00170241365804406,
#                                str_detect(locality, "New"), 0.00419556564165041,
#                                str_detect(locality, "Hendry"), 0.00302793707820219,
#                                str_detect(locality, "Goose"), 0.00136572458906164,
#                                str_detect(locality, "Lone"), 0.00388633006535513,
#                                str_detect(locality, "Little"), 0.00692501483225965,
#                                str_detect(locality, "Guana"), 0.00377890564371161,
#                                str_detect(locality, "Pellicer"), 0.00389536046654143,
#                                str_detect(locality, "Matanzas"), 0.00140742080211584,
#                                str_detect(locality, "Jack"), 0.00210247016092719)]
# 
# qhvals[, `:=` (orig_deadc = 1 - orig_f14c/qhvals[str_detect(locality, "Alligator"), orig_f14c],
#                orig_deadc_sd = sqrt((orig_f14c_sd/qhvals[str_detect(locality, "Alligator"), orig_f14c])^2 + (qhvals[str_detect(locality, "Alligator"), orig_f14c_sd] * orig_f14c/qhvals[str_detect(locality, "Alligator"), orig_f14c]^2)^2))]
# 
# deadC3 <- merge(deadC2, qhvals, by = "locality", all = TRUE)
# setcolorder(deadC3, c("locality", "n", "f14c_mn", "orig_f14c", "sd_method", "f14c_sd", "orig_f14c_sd", "deadc", "orig_deadc", "deadc_sd", "orig_deadc_sd"))
# fwrite(deadC3, here::here("DeadC_CalcCompare.csv"))


#Correct the f14c values using the dead carbon value for the relevant locality
hobsf14c[!is.na(locality), `:=` (f14c_corr = f14c/(1 - deadC[loc == locality & yr == year(collection_date), deadc]),
                                 f14c_corr_sd = sqrt((f14c_sd/(1 - deadC[loc == locality & yr == year(collection_date), deadc]))^2 + (deadC[loc == locality & yr == year(collection_date), deadc_sd] * f14c/(1 - deadC[loc == locality & yr == year(collection_date), deadc])^2)^2),
                                 f14c_corr_l = f14c/(1 - deadC[loc == locality & yr == year(collection_date), deadc_l]),
                                 f14c_corr_sd_l = sqrt((f14c_sd/(1 - deadC[loc == locality & yr == year(collection_date), deadc_l]))^2 + (deadC[loc == locality & yr == year(collection_date), deadc_sd_l] * f14c/(1 - deadC[loc == locality & yr == year(collection_date), deadc_l])^2)^2),
                                 f14c_corr_h = f14c/(1 - deadC[loc == locality & yr == year(collection_date), deadc_h]),
                                 f14c_corr_sd_h = sqrt((f14c_sd/(1 - deadC[loc == locality & yr == year(collection_date), deadc_h]))^2 + (deadC[loc == locality & yr == year(collection_date), deadc_sd_h] * f14c/(1 - deadC[loc == locality & yr == year(collection_date), deadc_h])^2)^2)), by = row.names(hobsf14c)]

#reconstruct original Quan corrections
hobsf14c[!is.na(locality), `:=` (f14c_corr = f14c/(1 - deadC[loc == locality, deadc]),
                                 f14c_corr_sd = sqrt((f14c_sd/(1 - deadC[loc == locality, deadc]))^2 + (deadC[loc == locality, deadc_sd] * f14c/(1 - deadC[loc == locality, deadc])^2)^2),
                                 f14c_corr_l = f14c/(1 - deadC[loc == locality, deadc_l]),
                                 f14c_corr_sd_l = sqrt((f14c_sd/(1 - deadC[loc == locality, deadc_l]))^2 + (deadC[loc == locality, deadc_sd_l] * f14c/(1 - deadC[loc == locality, deadc_l])^2)^2),
                                 f14c_corr_h = f14c/(1 - deadC[loc == locality, deadc_h]),
                                 f14c_corr_sd_h = sqrt((f14c_sd/(1 - deadC[loc == locality, deadc_h]))^2 + (deadC[loc == locality, deadc_sd_h] * f14c/(1 - deadC[loc == locality, deadc_h])^2)^2)), by = row.names(hobsf14c)]


hobsf14c[, locality_rectest := fcase(year(collection_date) == 1932, "Big Hickory",
                                     year(collection_date) == 1938, "Little St. George Island",
                                     year(collection_date) == 1940, "Lone Cabbage",
                                     year(collection_date) == 1952, "Lone Cabbage",
                                     year(collection_date) == 1979, "Lone Cabbage",
                                     default = NA)]
hobsf14c[is.na(locality_rectest), locality_rectest := locality]
hobsf14c[, `:=` (f14c_corr_rectest = f14c/(1 - deadC[loc == locality_rectest & yr == 2018, deadc]),
                 f14c_corr_sd_rectest = sqrt((f14c_sd/(1 - deadC[loc == locality_rectest & yr == 2018, deadc]))^2 + (deadC[loc == locality & yr == 2018, deadc_sd] * f14c/(1 - deadC[loc == locality & yr == 2018, deadc])^2)^2),
                 f14c_corr_l_rectest = f14c/(1 - deadC[loc == locality_rectest & yr == 2018, deadc_l]),
                 f14c_corr_sd_l_rectest = sqrt((f14c_sd/(1 - deadC[loc == locality_rectest & yr == 2018, deadc_l]))^2 + (deadC[loc == locality & yr == 2018, deadc_sd_l] * f14c/(1 - deadC[loc == locality & yr == 2018, deadc_l])^2)^2),
                 f14c_corr_h_rectest = f14c/(1 - deadC[loc == locality_rectest & yr == 2018, deadc_h]),
                 f14c_corr_sd_h_rectest = sqrt((f14c_sd/(1 - deadC[loc == locality_rectest & yr == 2018, deadc_h]))^2 + (deadC[loc == locality & yr == 2018, deadc_sd_h] * f14c/(1 - deadC[loc == locality & yr == 2018, deadc_h])^2)^2)), by = analysis_id]

#reconstruct original Quan corrections
hobsf14c[, `:=` (f14c_corr_rectest = f14c/(1 - deadC[loc == locality_rectest, deadc]),
                 f14c_corr_sd_rectest = sqrt((f14c_sd/(1 - deadC[loc == locality_rectest, deadc]))^2 + (deadC[loc == locality, deadc_sd] * f14c/(1 - deadC[loc == locality, deadc])^2)^2),
                 f14c_corr_l_rectest = f14c/(1 - deadC[loc == locality_rectest, deadc_l]),
                 f14c_corr_sd_l_rectest = sqrt((f14c_sd/(1 - deadC[loc == locality_rectest, deadc_l]))^2 + (deadC[loc == locality, deadc_sd_l] * f14c/(1 - deadC[loc == locality, deadc_l])^2)^2),
                 f14c_corr_h_rectest = f14c/(1 - deadC[loc == locality_rectest, deadc_h]),
                 f14c_corr_sd_h_rectest = sqrt((f14c_sd/(1 - deadC[loc == locality_rectest, deadc_h]))^2 + (deadC[loc == locality, deadc_sd_h] * f14c/(1 - deadC[loc == locality, deadc_h])^2)^2)), by = analysis_id]


#Calibrate radiocarbon results and generate sample-level estimates of total age variation (TAV)
source(here::here("CZ325_Sampling_task_2/Utilities 2.R"))
source(here::here("CZ325_Sampling_task_2/write_oxcal_sum_script 2.R"))
source(here::here("Summer2023/SetUpNewCurve_20230720.R")) #Don't forget to update curve data and path names in the "SetUpNewCurve" script if new prediction data has been generated since the last time it was sourced.

# oxcARR_Calibration = function(fractions, sample_name = "TEST", curvename = "HOBS_14C_CalCurve_PredictedVals_20230724", curvepath = here::here("Summer2023/HOBS_14C_CalCurve_PredictedVals_20230724.14c")){
# oxcARR_Calibration = function(fractions, sample_name = "TEST", curvename = "bombBahamasto10000calBP_Marine20", curvepath = here::here("CZ325_Sampling_task_2/bombBahamasto10000calBP_Marine20.14c")){
oxcARR_Calibration = function(fractions, sample_name = "TEST", curvename = "NNewCurveto10000calBP_Marine20", curvepath = here::here("Summer2023/NNewCurveto10000calBP_Marine20.14c")){
    
  shell_names = unique(fractions$name)
  code = write_oxcal_sum_script(fractions, path = "OxCal_Sum_scripts", filename = sample_name, curvename = curvename, curvepath = curvepath)
  
  my_result_file <- executeOxcalScript(code) # looks like it stores it in JSON
  my_result_text <- readOxcalOutput(my_result_file)
  
  #Or you get the whole output of Oxcal as object:
  
  #my_result_data <- parseFullOxcalOutput(my_result_text)
  my_result_data = parseOxcalOutput(my_result_text, only.R_Date = F)
  
  local_posterior_distribution = data.frame() # one per object
  for( p in my_result_data) {
    
    if( p$name %in% c(shell_names,"Sum"))
    {
      df = data.frame(year = p$raw_probabilities$dates,
                      probability = p$raw_probabilities$probabilities/sum( p$raw_probabilities$probabilities),
                      like_det = p$bp, #save likelihood determinations for inclusion in appendices
                      like_det_sd = p$std) #save likelihood determinations for inclusion in appendices
      df %<>% mutate(name = p$name)
      local_posterior_distribution = rbind(local_posterior_distribution, df)
    }
    
  }
  
  return(local_posterior_distribution)
}

all_posteriors <- data.table(igsn = character(),
                             name = character(),
                             analysis = character(),
                             deadcref = character(),
                             year = numeric(),
                             probability = numeric())
all_post_sum <- data.table(igsn = character(),
                           analysis = character(),
                           n = integer(),
                           n_oor = integer(),
                           n_tot = integer(),
                           tav_median = numeric(),
                           tav_min = numeric(),
                           tav_max = numeric(),
                           tav_range = numeric(),
                           tav_iqr = numeric())
all_post_specsum <- data.table(igsn = character(),
                               analysis = character(),
                               spec_field_id = character(),
                               oor = logical())


# for(s in hobsf14c[!is.na(igsn) & !is.na(f14c_corr), unique(igsn)]){
#   for(a in hobsf14c[igsn == s & !is.na(f14c_corr), unique(analysis)]){

set.seed(452)
for(a in hobsf14c[!is.na(f14c_corr), unique(analysis)]){
    #Generate and save posteriors
    # posteriors_s = oxcARR_Calibration(hobsf14c[igsn == s & analysis == a & !is.na(f14c_corr), .(name = spec_field_id, F14C_corr = f14c_corr, F14C_corr_sd = f14c_corr_sd, Depth = depth, Station = station)], sample_name = s)
    posteriors_s = oxcARR_Calibration(hobsf14c[analysis == a & !is.na(f14c_corr), .(name = spec_field_id, F14C_corr = f14c_corr, F14C_corr_sd = f14c_corr_sd, Depth = depth, Station = station)], sample_name = "All specimens_newcurveorigrefs")
    setDT(posteriors_s)
    posteriors_s <- posteriors_s[name != "Sum", ] # do not use "Sum" results; un-comment the for loops above if you want "Sum" values by sample and analysis
    posteriors_s_l = oxcARR_Calibration(hobsf14c[analysis == a & !is.na(f14c_corr_l), .(name = spec_field_id, F14C_corr = f14c_corr_l, F14C_corr_sd = f14c_corr_sd_l, Depth = depth, Station = station)], sample_name = "All specimens_lowdc_newcurveorigrefs")
    setDT(posteriors_s_l)
    posteriors_s_l <- posteriors_s_l[name != "Sum", ] # do not use "Sum" results; un-comment the for loops above if you want "Sum" values by sample and analysis
    posteriors_s_h = oxcARR_Calibration(hobsf14c[analysis == a & !is.na(f14c_corr_h), .(name = spec_field_id, F14C_corr = f14c_corr_h, F14C_corr_sd = f14c_corr_sd_h, Depth = depth, Station = station)], sample_name = "All specimens_highdc_newcurveorigrefs")
    setDT(posteriors_s_h)
    posteriors_s_h <- posteriors_s_h[name != "Sum", ] # do not use "Sum" results; un-comment the for loops above if you want "Sum" values by sample and analysis
    
    # posteriors_s[, `:=` (igsn = s, analysis = a)]
    posteriors_s[, `:=` (analysis = a, deadcref = "meanref")]
    posteriors_s_l[, `:=` (analysis = a, deadcref = "lowref")]
    posteriors_s_h[, `:=` (analysis = a, deadcref = "highref")]
    
    out_of_range <- posteriors_s[, .(oor = min(year) > 2022), by = list(analysis, name)]
    out_of_range_l <- posteriors_s_l[, .(oor_l = min(year) > 2022), by = list(analysis, name)]
    out_of_range_h <- posteriors_s_h[, .(oor_h = min(year) > 2022), by = list(analysis, name)]
    
    posteriors_s <- merge(posteriors_s, distinct(hobsf14c[!is.na(igsn) & analysis == a, .(name = spec_field_id, igsn, analysis)]), by = c("name", "analysis"), all.x = TRUE)
    all_posteriors <- bind_rows(all_posteriors, posteriors_s)
    posteriors_s_l <- merge(posteriors_s_l, distinct(hobsf14c[!is.na(igsn) & analysis == a, .(name = spec_field_id, igsn, analysis)]), by = c("name", "analysis"), all.x = TRUE)
    all_posteriors <- bind_rows(all_posteriors, posteriors_s_l)
    posteriors_s_h <- merge(posteriors_s_h, distinct(hobsf14c[!is.na(igsn) & analysis == a, .(name = spec_field_id, igsn, analysis)]), by = c("name", "analysis"), all.x = TRUE)
    all_posteriors <- bind_rows(all_posteriors, posteriors_s_h)
    
    #Summarize specimen-level posteriors and save results
    post_s_specsum <- posteriors_s[probability > 0, .(cal_age_med = unique(like_det),
                                                      cal_age_sd = unique(like_det_sd),
                                                      spec_median = wtd.quantile(year, probability * 10000, 0.5),
                                                      spec_min = wtd.quantile(year, probability * 10000, 0),
                                                      spec_max = wtd.quantile(year, probability * 10000, 1),
                                                      spec_range = wtd.quantile(year, probability * 10000, 1) - wtd.quantile(year, probability * 10000, 0),
                                                      spec_iqr = wtd.quantile(year, probability * 10000, 0.75) - wtd.quantile(year, probability * 10000, 0.25)), by = list(igsn, analysis, name)]
    post_s_specsum_l <- posteriors_s_l[probability > 0, .(cal_age_med_l = unique(like_det),
                                                          cal_age_sd_l = unique(like_det_sd),
                                                          spec_median_l = wtd.quantile(year, probability * 10000, 0.5),
                                                          spec_min_l = wtd.quantile(year, probability * 10000, 0),
                                                          spec_max_l = wtd.quantile(year, probability * 10000, 1),
                                                          spec_range_l = wtd.quantile(year, probability * 10000, 1) - wtd.quantile(year, probability * 10000, 0),
                                                          spec_iqr_l = wtd.quantile(year, probability * 10000, 0.75) - wtd.quantile(year, probability * 10000, 0.25)), by = list(igsn, analysis, name)]
    post_s_specsum_h <- posteriors_s_h[probability > 0, .(cal_age_med_h = unique(like_det),
                                                          cal_age_sd_h = unique(like_det_sd),
                                                          spec_median_h = wtd.quantile(year, probability * 10000, 0.5),
                                                          spec_min_h = wtd.quantile(year, probability * 10000, 0),
                                                          spec_max_h = wtd.quantile(year, probability * 10000, 1),
                                                          spec_range_h = wtd.quantile(year, probability * 10000, 1) - wtd.quantile(year, probability * 10000, 0),
                                                          spec_iqr_h = wtd.quantile(year, probability * 10000, 0.75) - wtd.quantile(year, probability * 10000, 0.25)), by = list(igsn, analysis, name)]
   
    all_post_specsum_a <- data.table(igsn = character(),
                                     analysis = character(),
                                     spec_field_id = character(),
                                     oor = logical()) 

    #note which specimen results were out of range
    post_s_specsum <- merge(post_s_specsum, out_of_range, by = c("name", "analysis"), all = TRUE)
    setnames(post_s_specsum, "name", "spec_field_id")
    all_post_specsum_a <- bind_rows(all_post_specsum_a, post_s_specsum)
    # all_post_specsum <- merge(all_post_specsum, post_s_specsum, by = c("igsn", "spec_field_id", "analysis"))
    
    post_s_specsum_l <- merge(post_s_specsum_l, out_of_range_l, by = c("name", "analysis"), all = TRUE)
    setnames(post_s_specsum_l, "name", "spec_field_id")
    # if(a == "HP"){
    all_post_specsum_a <- merge(all_post_specsum_a, post_s_specsum_l, by = c("igsn", "spec_field_id", "analysis"))
    # } else{
    #   all_post_specsum <- bind_rows(all_post_specsum, post_s_specsum_l)
    # }
    
    post_s_specsum_h <- merge(post_s_specsum_h, out_of_range_h, by = c("name", "analysis"), all = TRUE)
    setnames(post_s_specsum_h, "name", "spec_field_id")
    # if(a == "HP"){
    all_post_specsum_a <- merge(all_post_specsum_a, post_s_specsum_h, by = c("igsn", "spec_field_id", "analysis"))
    # } else{
    #   all_post_specsum <- bind_rows(all_post_specsum, post_s_specsum_h)
    # }
    
    all_post_specsum <- bind_rows(all_post_specsum, all_post_specsum_a)
    
    #If more than one specimen was calibrated successfully, calculate sample-level summaries and save results.
    if(a == "LP" & nrow(out_of_range[oor == FALSE, ]) > 1){
      tav <- TAV2(posteriors_s[probability > 0 & name %in% out_of_range[oor == FALSE, name], ])
      setDT(tav)
      tav[, igsn2 := igsn]
      tav[, `:=` (analysis = a,
                  n = posteriors_s[igsn == igsn2 & probability > 0 & name %in% out_of_range[oor == FALSE, name], length(unique(name))],
                  n_oor = posteriors_s[igsn == igsn2 & probability > 0 & name %in% out_of_range[oor == TRUE, name], length(unique(name))],
                  n_tot = posteriors_s[igsn == igsn2 & probability > 0, length(unique(name))]), by = igsn2]
      
      tav_l <- TAV2(posteriors_s_l[probability > 0 & name %in% out_of_range_l[oor_l == FALSE, name], ])
      setDT(tav_l)
      tav_l[, igsn2 := igsn]
      tav_l[, `:=` (analysis = a,
                    n_l = posteriors_s_l[igsn == igsn2 & probability > 0 & name %in% out_of_range_l[oor_l == FALSE, name], length(unique(name))],
                    n_l_oor = posteriors_s_l[igsn == igsn2 & probability > 0 & name %in% out_of_range_l[oor_l == TRUE, name], length(unique(name))],
                    n_l_tot = posteriors_s_l[igsn == igsn2 & probability > 0, length(unique(name))]), by = igsn2]
      
      tav_h <- TAV2(posteriors_s_h[probability > 0 & name %in% out_of_range_h[oor_h == FALSE, name], ])
      setDT(tav_h)
      tav_h[, igsn2 := igsn]
      tav_h[, `:=` (analysis = a,
                    n_h = posteriors_s_h[igsn == igsn2 & probability > 0 & name %in% out_of_range_h[oor_h == FALSE, name], length(unique(name))],
                    n_h_oor = posteriors_s_h[igsn == igsn2 & probability > 0 & name %in% out_of_range_h[oor_h == TRUE, name], length(unique(name))],
                    n_h_tot = posteriors_s_h[igsn == igsn2 & probability > 0, length(unique(name))]), by = igsn2]
      
      # post_sum_s <- tav[, .(igsn = unique(igsn),
      #                       analysis = unique(analysis),
      #                       n = unique(n),
      #                       n_oor = unique(n_oor),
      #                       n_tot = unique(n_tot),
      #                       tav_median = wtd.quantile(year, probability * 10000, 0.5),
      #                       tav_min = wtd.quantile(year, probability * 10000, 0),
      #                       tav_max = wtd.quantile(year, probability * 10000, 1),
      #                       tav_range = wtd.quantile(year, probability * 10000, 1) - wtd.quantile(year, probability * 10000, 0),
      #                       tav_iqr = wtd.quantile(year, probability * 10000, 0.75) - wtd.quantile(year, probability * 10000, 0.25))]
      post_sum_s <- tav[, .(tav_median = wtd.quantile(year, probability * 10000, 0.5),
                            tav_min = wtd.quantile(year, probability * 10000, 0),
                            tav_max = wtd.quantile(year, probability * 10000, 1),
                            tav_range = wtd.quantile(year, probability * 10000, 1) - wtd.quantile(year, probability * 10000, 0),
                            tav_iqr = wtd.quantile(year, probability * 10000, 0.75) - wtd.quantile(year, probability * 10000, 0.25)), by = list(igsn, analysis, n, n_oor, n_tot)]
      post_sum_s_l <- tav_l[, .(tav_median_l = wtd.quantile(year, probability * 10000, 0.5),
                                tav_min_l = wtd.quantile(year, probability * 10000, 0),
                                tav_max_l = wtd.quantile(year, probability * 10000, 1),
                                tav_range_l = wtd.quantile(year, probability * 10000, 1) - wtd.quantile(year, probability * 10000, 0),
                                tav_iqr_l = wtd.quantile(year, probability * 10000, 0.75) - wtd.quantile(year, probability * 10000, 0.25)), by = list(igsn, analysis, n_l, n_l_oor, n_l_tot)]
      post_sum_s_h <- tav_h[, .(tav_median_h = wtd.quantile(year, probability * 10000, 0.5),
                                tav_min_h = wtd.quantile(year, probability * 10000, 0),
                                tav_max_h = wtd.quantile(year, probability * 10000, 1),
                                tav_range_h = wtd.quantile(year, probability * 10000, 1) - wtd.quantile(year, probability * 10000, 0),
                                tav_iqr_h = wtd.quantile(year, probability * 10000, 0.75) - wtd.quantile(year, probability * 10000, 0.25)), by = list(igsn, analysis, n_h, n_h_oor, n_h_tot)]
      
      
      all_post_sum <- bind_rows(all_post_sum, post_sum_s)
      # all_post_sum <- merge(all_post_sum, post_sum_s, by = c("igsn", "analysis"), all = TRUE)
      all_post_sum <- merge(all_post_sum, post_sum_s_l, by = c("igsn", "analysis"), all = TRUE)
      all_post_sum <- merge(all_post_sum, post_sum_s_h, by = c("igsn", "analysis"), all = TRUE)
    }
    
    #Re-calibrate HP specimens using alternative F14C_corr using HOBS locality Dead C (i.e., the calculated "2018" values).
    if(a == "HP"){
      posteriors_nocorr = oxcARR_Calibration(hobsf14c[!is.na(f14c_corr_rectest) & analysis == a, .(name = spec_field_id, F14C_corr = f14c, F14C_corr_sd = f14c_sd, Depth = depth, Station = station)], sample_name = "NoCorr_newcurveorigrefs")
      setDT(posteriors_nocorr)
      posteriors_nocorr <- posteriors_nocorr[name != "Sum", ] # do not use "Sum" results; un-comment the for loops above if you want "Sum" values by sample and analysis
      posteriors_nocorr[, analysis := paste0(a, "_nocorr")]
      out_of_range_nocorr <- posteriors_nocorr[, .(oor_nocorr = min(year) > 2020), by = list(analysis, name)]
      posteriors_nocorr <- merge(posteriors_nocorr, distinct(hobsf14c[!is.na(f14c_corr_rectest) & analysis == a, .(name = spec_field_id, igsn)]), by = "name", all.x = TRUE)
      all_posteriors <- bind_rows(all_posteriors, posteriors_nocorr)
      # all_posteriors <- merge(all_posteriors, posteriors_nocorr, by = "name", all.x = TRUE)
      
      posteriors_rectest = oxcARR_Calibration(hobsf14c[!is.na(f14c_corr_rectest) & analysis == a, .(name = spec_field_id, F14C_corr = f14c_corr_rectest, F14C_corr_sd = f14c_corr_sd_rectest, Depth = depth, Station = station)], sample_name = "RecTest_newcurveorigrefs")
      setDT(posteriors_rectest)
      posteriors_rectest <- posteriors_rectest[name != "Sum", ] # do not use "Sum" results; un-comment the for loops above if you want "Sum" values by sample and analysis
      posteriors_rectest[, `:=` (analysis = paste0(a, "_rectest"), deadcref = "meanref")]
      out_of_range_rectest <- posteriors_rectest[, .(oor_rectest = min(year) > 2020), by = list(analysis, name)]
      posteriors_rectest <- merge(posteriors_rectest, distinct(hobsf14c[!is.na(f14c_corr_rectest) & analysis == a, .(name = spec_field_id, igsn)]), by = "name", all.x = TRUE)
      all_posteriors <- bind_rows(all_posteriors, posteriors_rectest)
      # all_posteriors <- merge(all_posteriors, posteriors_rectest, by = "name", all.x = TRUE)
      
      posteriors_rectest_l = oxcARR_Calibration(hobsf14c[!is.na(f14c_corr_l_rectest) & analysis == a, .(name = spec_field_id, F14C_corr = f14c_corr_l_rectest, F14C_corr_sd = f14c_corr_sd_l_rectest, Depth = depth, Station = station)], sample_name = "RecTest_l_newcurveorigrefs")
      setDT(posteriors_rectest_l)
      posteriors_rectest_l <- posteriors_rectest_l[name != "Sum", ] # do not use "Sum" results; un-comment the for loops above if you want "Sum" values by sample and analysis
      posteriors_rectest_l[, `:=` (analysis = paste0(a, "_rectest"), deadcref = "lowref")]
      out_of_range_rectest_l <- posteriors_rectest_l[, .(oor_rectest_l = min(year) > 2020), by = list(analysis, name)]
      posteriors_rectest_l <- merge(posteriors_rectest_l, distinct(hobsf14c[!is.na(f14c_corr_l_rectest) & analysis == a, .(name = spec_field_id, igsn)]), by = "name", all.x = TRUE)
      all_posteriors <- bind_rows(all_posteriors, posteriors_rectest_l)
      
      posteriors_rectest_h = oxcARR_Calibration(hobsf14c[!is.na(f14c_corr_h_rectest) & analysis == a, .(name = spec_field_id, F14C_corr = f14c_corr_h_rectest, F14C_corr_sd = f14c_corr_sd_h_rectest, Depth = depth, Station = station)], sample_name = "RecTest_h_newcurveorigrefs")
      setDT(posteriors_rectest_h)
      posteriors_rectest_h <- posteriors_rectest_h[name != "Sum", ] # do not use "Sum" results; un-comment the for loops above if you want "Sum" values by sample and analysis
      posteriors_rectest_h[, `:=` (analysis = paste0(a, "_rectest"), deadcref = "highref")]
      out_of_range_rectest_h <- posteriors_rectest_h[, .(oor_rectest_h = min(year) > 2020), by = list(analysis, name)]
      posteriors_rectest_h <- merge(posteriors_rectest_h, distinct(hobsf14c[!is.na(f14c_corr_h_rectest) & analysis == a, .(name = spec_field_id, igsn)]), by = "name", all.x = TRUE)
      all_posteriors <- bind_rows(all_posteriors, posteriors_rectest_h)
      
      #Summarize specimen-level posteriors and save results
      posteriors_nocorr[, analysis := "HP"]
      post_nocorr_specsum <- posteriors_nocorr[probability > 0, .(cal_age_med_nocorr = unique(like_det),
                                                                  cal_age_sd_nocorr = unique(like_det_sd),
                                                                  spec_median_nocorr = wtd.quantile(year, probability * 10000, 0.5),
                                                                  spec_min_nocorr = wtd.quantile(year, probability * 10000, 0),
                                                                  spec_max_nocorr = wtd.quantile(year, probability * 10000, 1),
                                                                  spec_range_nocorr = wtd.quantile(year, probability * 10000, 1) - wtd.quantile(year, probability * 10000, 0),
                                                                  spec_iqr_nocorr = wtd.quantile(year, probability * 10000, 0.75) - wtd.quantile(year, probability * 10000, 0.25)), by = list(igsn, analysis, name)]
    
      posteriors_rectest[, analysis := "HP"]
      post_rectest_specsum <- posteriors_rectest[probability > 0, .(cal_age_med_rectest = unique(like_det),
                                                                    cal_age_sd_rectest = unique(like_det_sd),
                                                                    spec_median_rectest = wtd.quantile(year, probability * 10000, 0.5),
                                                                    spec_min_rectest = wtd.quantile(year, probability * 10000, 0),
                                                                    spec_max_rectest = wtd.quantile(year, probability * 10000, 1),
                                                                    spec_range_rectest = wtd.quantile(year, probability * 10000, 1) - wtd.quantile(year, probability * 10000, 0),
                                                                    spec_iqr_rectest = wtd.quantile(year, probability * 10000, 0.75) - wtd.quantile(year, probability * 10000, 0.25)), by = list(igsn, analysis, name)]
      
      posteriors_rectest_l[, analysis := "HP"]
      post_rectest_specsum_l <- posteriors_rectest_l[probability > 0, .(cal_age_med_l_rectest = unique(like_det),
                                                                        cal_age_sd_l_rectest = unique(like_det_sd),
                                                                        spec_median_l_rectest = wtd.quantile(year, probability * 10000, 0.5),
                                                                        spec_min_l_rectest = wtd.quantile(year, probability * 10000, 0),
                                                                        spec_max_l_rectest = wtd.quantile(year, probability * 10000, 1),
                                                                        spec_range_l_rectest = wtd.quantile(year, probability * 10000, 1) - wtd.quantile(year, probability * 10000, 0),
                                                                        spec_iqr_l_rectest = wtd.quantile(year, probability * 10000, 0.75) - wtd.quantile(year, probability * 10000, 0.25)), by = list(igsn, analysis, name)]
      
      posteriors_rectest_h[, analysis := "HP"]
      post_rectest_specsum_h <- posteriors_rectest_h[probability > 0, .(cal_age_med_h_rectest = unique(like_det),
                                                                        cal_age_sd_h_rectest = unique(like_det_sd),
                                                                        spec_median_h_rectest = wtd.quantile(year, probability * 10000, 0.5),
                                                                        spec_min_h_rectest = wtd.quantile(year, probability * 10000, 0),
                                                                        spec_max_h_rectest = wtd.quantile(year, probability * 10000, 1),
                                                                        spec_range_h_rectest = wtd.quantile(year, probability * 10000, 1) - wtd.quantile(year, probability * 10000, 0),
                                                                        spec_iqr_h_rectest = wtd.quantile(year, probability * 10000, 0.75) - wtd.quantile(year, probability * 10000, 0.25)), by = list(igsn, analysis, name)]
      
      #note which specimen results were out of range
      out_of_range_nocorr[, analysis := "HP"]
      post_nocorr_specsum <- merge(post_nocorr_specsum, out_of_range_nocorr, by = c("name", "analysis"), all = TRUE)
      setnames(post_nocorr_specsum, "name", "spec_field_id")
      # all_post_specsum <- bind_rows(all_post_specsum, post_nocorr_specsum)
      all_post_specsum <- merge(all_post_specsum, post_nocorr_specsum, by = c("igsn", "spec_field_id", "analysis"), all = TRUE)
      
      out_of_range_rectest[, analysis := "HP"]
      post_rectest_specsum <- merge(post_rectest_specsum, out_of_range_rectest, by = c("name", "analysis"), all = TRUE)
      setnames(post_rectest_specsum, "name", "spec_field_id")
      # all_post_specsum <- bind_rows(all_post_specsum, post_rectest_specsum)
      all_post_specsum <- merge(all_post_specsum, post_rectest_specsum, by = c("igsn", "spec_field_id", "analysis"), all = TRUE)
      
      out_of_range_rectest_l[, analysis := "HP"]
      post_rectest_specsum_l <- merge(post_rectest_specsum_l, out_of_range_rectest_l, by = c("name", "analysis"), all = TRUE)
      setnames(post_rectest_specsum_l, "name", "spec_field_id")
      all_post_specsum <- merge(all_post_specsum, post_rectest_specsum_l, by = c("igsn", "spec_field_id", "analysis"), all = TRUE)
      
      out_of_range_rectest_h[, analysis := "HP"]
      post_rectest_specsum_h <- merge(post_rectest_specsum_h, out_of_range_rectest_h, by = c("name", "analysis"), all = TRUE)
      setnames(post_rectest_specsum_h, "name", "spec_field_id")
      all_post_specsum <- merge(all_post_specsum, post_rectest_specsum_h, by = c("igsn", "spec_field_id", "analysis"), all = TRUE)
    }
    
  #  print(paste0("Done with ", s, ", ", a))
  # }
  print(paste0("Done with ", a, " specimens"))
}

hobsf14c2 <- merge(hobsf14c, all_post_sum, by = c("igsn", "analysis"), all = TRUE)
hobsf14c2 <- merge(hobsf14c2, all_post_specsum, by = c("igsn", "analysis", "spec_field_id"), all = TRUE)
hobsf14c2[, `:=` (median_age_bp = 2022 - tav_median)]
hobsf14c2[!is.na(igsn), reef := str_sub(unique(spec_field_id), 1, str_locate(unique(spec_field_id), "_R")[2] + 1), by = spec_field_id]
hobsf14c2[!is.na(igsn), hole := str_sub(unique(spec_field_id), str_locate(unique(spec_field_id), unique(reef))[2] + 1, str_locate(unique(spec_field_id), unique(reef))[2] + 2), by = spec_field_id]
hobsf14c2[!is.na(reef), `:=` (loc = str_split_fixed(hobsf14c2[!is.na(reef), reef], "_R", 2)[, 1],
                              reef_num = as.integer(str_split_fixed(hobsf14c2[!is.na(reef), reef], "_R", 2)[, 2]))]
reef_ind <- distinct(hobsf14c2[!is.na(reef), .(loc, reef_num)])
setorder(reef_ind, loc, reef_num)
reef_ind[, reef_ind := seq(1, .N), by = loc]
hobsf14c2 <- merge(hobsf14c2, reef_ind, by = c("loc", "reef_num"), all = TRUE)
hobsf14c2[loc == "JI-WC", loc := "JI"]
hobsf14c2[, loc := factor(loc, levels = c("LSG", "GI-EC", "LC", "LB", "HC-MC", "NP", "BH", "JI", "PC", "MR", "GR"), ordered = TRUE)]
hobsf14c2 <- distinct(hobsf14c2[!is.na(analysis_id), ])

out_of_range <- hobsf14c2[oor == TRUE, ] #63 specimens had impossible median calibrated ages; range(Curve$F14C) is 0.9058292 to 1.1407439, max F14C value that was calibrated was 1.159851, min F14C value that failed was 1.156391.
                                         #Calculating dead C based on the calibration curve f14c for each live-caught specimen's year of collection reduced the number of out of range specimens to 50.
                                         #Combining the new bomb-pulse curve with Marine20 reduced the number oor slightly to 48 specimens.

#rectest
ggplot() +
  geom_ribbon(data = Curve, aes(x = `#year`, ymin = F14C - F14Csd, ymax = F14C + F14Csd), fill = "mistyrose") +
  geom_line(data = Curve, aes(x = `#year`, y = F14C), color = "lightpink") +
  geom_ribbon(data = newcurvepreds, aes(x = year, ymin = f14c_mn - f14c_sd, ymax = f14c_mn + f14c_sd), fill = "grey90") +
  geom_line(data = newcurvepreds, aes(x = year, y = f14c_mn)) + 
  geom_point(data = hobsf14c2[analysis_id %in% c("UAL19523", "UAL19520", "UAL19521", "UAL19522"), ], aes(x = year(collection_date), y = f14c, shape = spec_field_id, color = "no corr."), size = 2) + 
  geom_point(data = hobsf14c2[analysis_id %in% c("UAL19523", "UAL19520", "UAL19521", "UAL19522"), ], aes(x = year(collection_date), y = f14c_corr, shape = spec_field_id, color = "coll. yr. corr."), size = 2) + 
  geom_point(data = hobsf14c2[analysis_id %in% c("UAL19523", "UAL19520", "UAL19521", "UAL19522"), ], aes(x = year(collection_date), y = f14c_corr_rectest, shape = spec_field_id, color = "rec.-only corr."), size = 2) + 
  geom_segment(data = hobsf14c2[analysis_id %in% c("UAL19523", "UAL19520", "UAL19521", "UAL19522"), ], aes(x = year(collection_date), y = f14c, xend = spec_median_nocorr, yend = f14c, color = "no corr."), lwd = 0.5, arrow = arrow(length = unit(4, "points"), type = "closed", angle = 30)) + 
  geom_segment(data = hobsf14c2[analysis_id %in% c("UAL19523", "UAL19520", "UAL19521", "UAL19522"), ], aes(x = year(collection_date), y = f14c_corr, xend = spec_median, yend = f14c_corr, color = "coll. yr. corr."), lwd = 0.5, arrow = arrow(length = unit(4, "points"), type = "closed", angle = 30)) + 
  geom_segment(data = hobsf14c2[analysis_id %in% c("UAL19523", "UAL19520", "UAL19521", "UAL19522"), ], aes(x = year(collection_date), y = f14c_corr_rectest, xend = spec_median_rectest, yend = f14c_corr_rectest, color = "rec.-only corr."), lwd = 0.5, arrow = arrow(length = unit(4, "points"), type = "closed", angle = 30)) + 
  theme_classic() +
  scale_color_manual(values = c("no corr." = "green", "coll. yr. corr." = "firebrick", "rec.-only corr." = "dodgerblue")) +
  coord_cartesian(xlim = c(1750, 2023))

ggplot() +
  geom_ribbon(data = Curve, aes(x = `#year`, ymin = F14C - F14Csd, ymax = F14C + F14Csd), fill = "mistyrose", alpha = 0.5) +
  geom_line(data = Curve, aes(x = `#year`, y = F14C), color = "lightpink") +
  geom_ribbon(data = newcurvepreds, aes(x = year, ymin = f14c_mn - f14c_sd, ymax = f14c_mn + f14c_sd), fill = "grey90", alpha = 0.5) +
  geom_line(data = newcurvepreds, aes(x = year, y = f14c_mn)) + 
  geom_segment(data = hobsf14c2[!is.na(locality_rectest) & is.na(igsn.y) & is.na(igsn), ], aes(x = year(collection_date), y = f14c, xend = spec_median_nocorr, yend = f14c), lwd = 0.5, arrow = arrow(length = unit(4, "points"), type = "closed", angle = 30)) + 
  geom_segment(data = hobsf14c2[!is.na(locality_rectest) & is.na(igsn.y) & is.na(igsn), ], aes(x = year(collection_date), y = f14c_corr, xend = spec_median, yend = f14c_corr), lwd = 0.5, arrow = arrow(length = unit(4, "points"), type = "closed", angle = 30)) + 
  geom_segment(data = hobsf14c2[!is.na(locality_rectest) & is.na(igsn.y) & is.na(igsn), ], aes(x = year(collection_date), y = f14c_corr_rectest, xend = spec_median_rectest, yend = f14c_corr_rectest), lwd = 0.5, arrow = arrow(length = unit(4, "points"), type = "closed", angle = 30)) + 
  geom_point(data = hobsf14c2[!is.na(locality_rectest) & is.na(igsn.y) & is.na(igsn), ], aes(x = year(collection_date), y = f14c, shape = "no corr.", size = "no corr.", fill = spec_field_id), color = "black") + 
  geom_point(data = hobsf14c2[!is.na(locality_rectest) & is.na(igsn.y) & is.na(igsn), ], aes(x = year(collection_date), y = f14c_corr, shape = "coll. yr. corr.", size = "coll. yr. corr.", fill = spec_field_id), color = "black") + 
  geom_point(data = hobsf14c2[!is.na(locality_rectest) & is.na(igsn.y) & is.na(igsn), ], aes(x = year(collection_date), y = f14c_corr_rectest, shape = "rec.-only corr.", size = "rec.-only corr.", fill = spec_field_id), color = "black") + 
  theme_classic() +
  # scale_color_manual(values = c("no corr." = "green", "coll. yr. corr." = "firebrick", "rec.-only corr." = "dodgerblue")) +
  scale_shape_manual(values = c("no corr." = 21, "coll. yr. corr." = 22, "rec.-only corr." = 23)) +
  scale_size_manual(values = c("no corr." = 4, "coll. yr. corr." = 3, "rec.-only corr." = 2)) +
  coord_cartesian(xlim = c(1950, 2023)) +
  facet_wrap(~locality, scales = "free")


lcplotdat <- hobsf14c2[locality_rectest == "Lone Cabbage" & is.na(igsn), .(locality, collection_date, spec_field_id, f14c, f14c_sd, f14c_corr, f14c_corr_sd, f14c_corr_l, 
                                                                           f14c_corr_sd_l, f14c_corr_h, f14c_corr_sd_h, f14c_corr_rectest, f14c_corr_sd_rectest, f14c_corr_l_rectest, 
                                                                           f14c_corr_sd_l_rectest, f14c_corr_h_rectest, f14c_corr_sd_h_rectest, spec_median, spec_median_l, 
                                                                           spec_median_h, spec_median_nocorr, spec_median_rectest, spec_median_l_rectest, spec_median_h_rectest)]

ggplot() +
  geom_ribbon(data = Curve, aes(x = `#year`, ymin = F14C - F14Csd, ymax = F14C + F14Csd), fill = "mistyrose", alpha = 0.5) +
  geom_line(data = Curve, aes(x = `#year`, y = F14C), color = "lightpink") +
  geom_ribbon(data = newcurvepreds, aes(x = year, ymin = f14c_mn - f14c_sd, ymax = f14c_mn + f14c_sd), fill = "grey90", alpha = 0.5) +
  geom_line(data = newcurvepreds, aes(x = year, y = f14c_mn)) + 
  geom_segment(data = lcplotdat, aes(x = year(collection_date), y = f14c, xend = spec_median_nocorr, yend = f14c, color = "no corr."), lwd = 1, arrow = arrow(length = unit(6, "points"), type = "closed", angle = 30)) + 
  geom_segment(data = lcplotdat, aes(x = year(collection_date), y = f14c_corr, xend = spec_median, yend = f14c_corr, color = "coll. yr. corr. mn."), lwd = 0.75, arrow = arrow(length = unit(5, "points"), type = "closed", angle = 30)) + 
  geom_segment(data = lcplotdat, aes(x = year(collection_date), y = f14c_corr_l, xend = spec_median_l, yend = f14c_corr_l, color = "coll. yr. corr. l."), lwd = 0.75, arrow = arrow(length = unit(5, "points"), type = "closed", angle = 30)) + 
  geom_segment(data = lcplotdat, aes(x = year(collection_date), y = f14c_corr_h, xend = spec_median_h, yend = f14c_corr_h, color = "coll. yr. corr. h."), lwd = 0.75, arrow = arrow(length = unit(5, "points"), type = "closed", angle = 30)) + 
  geom_segment(data = lcplotdat, aes(x = year(collection_date), y = f14c_corr_l_rectest, xend = spec_median_l_rectest, yend = f14c_corr_l_rectest, color = "rec.-only corr. l."), lwd = 0.5, arrow = arrow(length = unit(4, "points"), type = "closed", angle = 30)) + 
  geom_segment(data = lcplotdat, aes(x = year(collection_date), y = f14c_corr_h_rectest, xend = spec_median_h_rectest, yend = f14c_corr_h_rectest, color = "rec.-only corr. h."), lwd = 0.5, arrow = arrow(length = unit(4, "points"), type = "closed", angle = 30)) + 
  geom_segment(data = lcplotdat, aes(x = year(collection_date), y = f14c_corr_rectest, xend = spec_median_rectest, yend = f14c_corr_rectest, color = "rec.-only corr. mn."), lwd = 0.5, arrow = arrow(length = unit(4, "points"), type = "closed", angle = 30)) + 
  # geom_errorbar(data = lcplotdat, aes(x = year(collection_date), ymin = f14c - f14c_sd, ymax = f14c + f14c_sd, color = "no corr."), lwd = 1, width = 3) +
  # geom_errorbar(data = lcplotdat, aes(x = year(collection_date), ymin = f14c_corr - f14c_corr_sd, ymax = f14c_corr + f14c_corr_sd, color = "coll. yr. corr. mn."), lwd = 0.75, width = 2) + 
  # geom_errorbar(data = lcplotdat, aes(x = year(collection_date), ymin = f14c_corr_l - f14c_corr_sd_l, ymax = f14c_corr_l + f14c_corr_sd_l, color = "coll. yr. corr. l."), lwd = 0.75, width = 2) + 
  # geom_errorbar(data = lcplotdat, aes(x = year(collection_date), ymin = f14c_corr_h - f14c_corr_sd_h, ymax = f14c_corr_h + f14c_corr_sd_h, color = "coll. yr. corr. h."), lwd = 0.75, width = 2) + 
  # geom_errorbar(data = lcplotdat, aes(x = year(collection_date), ymin = f14c_corr_rectest - f14c_corr_sd_rectest, ymax = f14c_corr_rectest + f14c_corr_sd_rectest, color = "rec.-only corr. mn."), lwd = 0.5, width = 1) + 
  # geom_errorbar(data = lcplotdat, aes(x = year(collection_date), ymin = f14c_corr_l_rectest - f14c_corr_sd_l_rectest, ymax = f14c_corr_l_rectest + f14c_corr_sd_l_rectest, color = "rec.-only corr. l."), lwd = 0.5, width = 1) + 
  # geom_errorbar(data = lcplotdat, aes(x = year(collection_date), ymin = f14c_corr_h_rectest - f14c_corr_sd_h_rectest, ymax = f14c_corr_h_rectest + f14c_corr_sd_h_rectest, color = "rec.-only corr. h."), lwd = 0.5, width = 1) + 
  geom_point(data = lcplotdat, aes(x = year(collection_date), y = f14c, color = "no corr.", size = "no corr.", shape = spec_field_id)) + 
  geom_point(data = lcplotdat, aes(x = year(collection_date), y = f14c_corr, color = "coll. yr. corr. mn.", size = "coll. yr. corr.", shape = spec_field_id)) + 
  geom_point(data = lcplotdat, aes(x = year(collection_date), y = f14c_corr_l, color = "coll. yr. corr. l.", size = "coll. yr. corr.", shape = spec_field_id)) + 
  geom_point(data = lcplotdat, aes(x = year(collection_date), y = f14c_corr_h, color = "coll. yr. corr. h.", size = "coll. yr. corr.", shape = spec_field_id)) + 
  geom_point(data = lcplotdat, aes(x = year(collection_date), y = f14c_corr_rectest, color = "rec.-only corr. mn.", size = "rec.-only corr.", shape = spec_field_id)) + 
  geom_point(data = lcplotdat, aes(x = year(collection_date), y = f14c_corr_l_rectest, color = "rec.-only corr. l.", size = "rec.-only corr.", shape = spec_field_id)) + 
  geom_point(data = lcplotdat, aes(x = year(collection_date), y = f14c_corr_h_rectest, color = "rec.-only corr. h.", size = "rec.-only corr.", shape = spec_field_id)) + 
  theme_classic() +
  scale_color_manual(values = c("no corr." = "green", "coll. yr. corr. mn." = "firebrick3", "coll. yr. corr. l." = "indianred1", "coll. yr. corr. h." = "darkred", "rec.-only corr. mn." = "dodgerblue3", "rec.-only corr. l." = "lightskyblue", "rec.-only corr. h." = "darkblue")) +
  # scale_shape_manual(values = c("no corr." = 21, "coll. yr. corr." = 22, "rec.-only corr." = 23)) +
  scale_size_manual(values = c("no corr." = 4, "coll. yr. corr." = 3, "rec.-only corr." = 2)) +
  coord_cartesian(xlim = c(1800, 2023))

lcmap <- st_as_sf(hobsf14c2[locality_rectest == "Lone Cabbage" & is.na(igsn), ], coords = c("lon", "lat"), crs = 4326)
mapview::mapview(lcmap, zcol = "spec_field_id")



hobsplotdat <- hobsf14c2[live_or_dead == "Live-caught", .(locality, collection_date, spec_field_id, f14c, f14c_sd, f14c_corr, f14c_corr_sd, f14c_corr_l, 
                                                          f14c_corr_sd_l, f14c_corr_h, f14c_corr_sd_h, f14c_corr_rectest, f14c_corr_sd_rectest, f14c_corr_l_rectest, 
                                                          f14c_corr_sd_l_rectest, f14c_corr_h_rectest, f14c_corr_sd_h_rectest, spec_median, spec_median_l, 
                                                          spec_median_h, spec_median_nocorr, spec_median_rectest, spec_median_l_rectest, spec_median_h_rectest)]

for(ll in unique(hobsplotdat$locality)){
  for(ss in hobsplotdat[locality == ll, unique(spec_field_id)]){
    hobsplotdat[locality == ll & spec_field_id == ss, spec := as.factor(which(hobsplotdat[locality == ll, sort(unique(spec_field_id))] == ss))]
  }
}

ggplot() +
  geom_ribbon(data = Curve, aes(x = `#year`, ymin = F14C - F14Csd, ymax = F14C + F14Csd), fill = "mistyrose", alpha = 0.5) +
  geom_line(data = Curve, aes(x = `#year`, y = F14C), color = "lightpink") +
  geom_ribbon(data = newcurvepreds, aes(x = year, ymin = f14c_mn - f14c_sd, ymax = f14c_mn + f14c_sd), fill = "grey90", alpha = 0.5) +
  geom_line(data = newcurvepreds, aes(x = year, y = f14c_mn)) + 
  geom_segment(data = hobsplotdat, aes(x = year(collection_date), y = f14c, xend = spec_median_nocorr, yend = f14c, color = "no corr."), lwd = 1, arrow = arrow(length = unit(6, "points"), type = "closed", angle = 30)) + 
  geom_segment(data = hobsplotdat, aes(x = year(collection_date), y = f14c_corr, xend = spec_median, yend = f14c_corr, color = "coll. yr. corr. mn."), lwd = 0.75, arrow = arrow(length = unit(5, "points"), type = "closed", angle = 30)) + 
  geom_segment(data = hobsplotdat, aes(x = year(collection_date), y = f14c_corr_l, xend = spec_median_l, yend = f14c_corr_l, color = "coll. yr. corr. l."), lwd = 0.75, arrow = arrow(length = unit(5, "points"), type = "closed", angle = 30)) + 
  geom_segment(data = hobsplotdat, aes(x = year(collection_date), y = f14c_corr_h, xend = spec_median_h, yend = f14c_corr_h, color = "coll. yr. corr. h."), lwd = 0.75, arrow = arrow(length = unit(5, "points"), type = "closed", angle = 30)) + 
  geom_segment(data = hobsplotdat, aes(x = year(collection_date), y = f14c_corr_l_rectest, xend = spec_median_l_rectest, yend = f14c_corr_l_rectest, color = "rec.-only corr. l."), lwd = 0.5, arrow = arrow(length = unit(4, "points"), type = "closed", angle = 30)) + 
  geom_segment(data = hobsplotdat, aes(x = year(collection_date), y = f14c_corr_h_rectest, xend = spec_median_h_rectest, yend = f14c_corr_h_rectest, color = "rec.-only corr. h."), lwd = 0.5, arrow = arrow(length = unit(4, "points"), type = "closed", angle = 30)) + 
  geom_segment(data = hobsplotdat, aes(x = year(collection_date), y = f14c_corr_rectest, xend = spec_median_rectest, yend = f14c_corr_rectest, color = "rec.-only corr. mn."), lwd = 0.5, arrow = arrow(length = unit(4, "points"), type = "closed", angle = 30)) + 
  # geom_errorbar(data = hobsplotdat, aes(x = year(collection_date), ymin = f14c - f14c_sd, ymax = f14c + f14c_sd, color = "no corr."), lwd = 1, width = 3) +
  # geom_errorbar(data = hobsplotdat, aes(x = year(collection_date), ymin = f14c_corr - f14c_corr_sd, ymax = f14c_corr + f14c_corr_sd, color = "coll. yr. corr. mn."), lwd = 0.75, width = 2) + 
  # geom_errorbar(data = hobsplotdat, aes(x = year(collection_date), ymin = f14c_corr_l - f14c_corr_sd_l, ymax = f14c_corr_l + f14c_corr_sd_l, color = "coll. yr. corr. l."), lwd = 0.75, width = 2) + 
  # geom_errorbar(data = hobsplotdat, aes(x = year(collection_date), ymin = f14c_corr_h - f14c_corr_sd_h, ymax = f14c_corr_h + f14c_corr_sd_h, color = "coll. yr. corr. h."), lwd = 0.75, width = 2) + 
  # geom_errorbar(data = hobsplotdat, aes(x = year(collection_date), ymin = f14c_corr_rectest - f14c_corr_sd_rectest, ymax = f14c_corr_rectest + f14c_corr_sd_rectest, color = "rec.-only corr. mn."), lwd = 0.5, width = 1) + 
  # geom_errorbar(data = hobsplotdat, aes(x = year(collection_date), ymin = f14c_corr_l_rectest - f14c_corr_sd_l_rectest, ymax = f14c_corr_l_rectest + f14c_corr_sd_l_rectest, color = "rec.-only corr. l."), lwd = 0.5, width = 1) + 
  # geom_errorbar(data = hobsplotdat, aes(x = year(collection_date), ymin = f14c_corr_h_rectest - f14c_corr_sd_h_rectest, ymax = f14c_corr_h_rectest + f14c_corr_sd_h_rectest, color = "rec.-only corr. h."), lwd = 0.5, width = 1) + 
  geom_point(data = hobsplotdat, aes(x = year(collection_date), y = f14c, color = "no corr.", size = "no corr.", shape = spec)) + 
  geom_point(data = hobsplotdat, aes(x = year(collection_date), y = f14c_corr, color = "coll. yr. corr. mn.", size = "coll. yr. corr.", shape = spec)) + 
  geom_point(data = hobsplotdat, aes(x = year(collection_date), y = f14c_corr_l, color = "coll. yr. corr. l.", size = "coll. yr. corr.", shape = spec)) + 
  geom_point(data = hobsplotdat, aes(x = year(collection_date), y = f14c_corr_h, color = "coll. yr. corr. h.", size = "coll. yr. corr.", shape = spec)) + 
  geom_point(data = hobsplotdat, aes(x = year(collection_date), y = f14c_corr_rectest, color = "rec.-only corr. mn.", size = "rec.-only corr.", shape = spec)) + 
  geom_point(data = hobsplotdat, aes(x = year(collection_date), y = f14c_corr_l_rectest, color = "rec.-only corr. l.", size = "rec.-only corr.", shape = spec)) + 
  geom_point(data = hobsplotdat, aes(x = year(collection_date), y = f14c_corr_h_rectest, color = "rec.-only corr. h.", size = "rec.-only corr.", shape = spec)) + 
  theme_classic() +
  scale_color_manual(values = c("no corr." = "green", "coll. yr. corr. mn." = "firebrick3", "coll. yr. corr. l." = "indianred1", "coll. yr. corr. h." = "darkred", "rec.-only corr. mn." = "dodgerblue3", "rec.-only corr. l." = "lightskyblue", "rec.-only corr. h." = "darkblue")) +
  # scale_shape_manual(values = c("no corr." = 21, "coll. yr. corr." = 22, "rec.-only corr." = 23)) +
  scale_size_manual(values = c("no corr." = 4, "coll. yr. corr." = 3, "rec.-only corr." = 2)) +
  coord_cartesian(xlim = c(1800, 2023)) +
  facet_wrap(~locality, scales = "free")



#New Figure S4 to replace Table S5-----------------------------------------------------------
hobsplotdat <- hobsf14c2[live_or_dead == "Live-caught" & year(collection_date) < 2019, .(locality, collection_date, genus, spec_field_id, cat_no, oor, f14c, f14c_sd, f14c_corr, f14c_corr_sd, f14c_corr_l, 
                                                                                         f14c_corr_sd_l, f14c_corr_h, f14c_corr_sd_h, f14c_corr_rectest, f14c_corr_sd_rectest, f14c_corr_l_rectest, 
                                                                                         f14c_corr_sd_l_rectest, f14c_corr_h_rectest, f14c_corr_sd_h_rectest, spec_median, spec_median_l, 
                                                                                         spec_median_h, spec_median_nocorr, spec_median_rectest, spec_median_l_rectest, spec_median_h_rectest)]

hobsplotdat[, `:=` (locality2 = ifelse(str_detect(locality, "NA;"), fcase(str_detect(locality, "George"), "Little St. George Island",
                                                                          str_detect(locality, "Hickory"), "Big Hickory",
                                                                          str_detect(locality, "Cabbage"), "Lone Cabbage"), locality),
                    nearby = factor(ifelse(str_detect(locality, "NA;"), "Nearby", "Locality"), levels = c("Nearby", "Locality")))]

hobsplotdat[, `:=` (locality = factor(locality, levels = c("Little St. George Island", "NA; nearest is Little St. George Island", "Goose Island/East Cove", "Alligator Harbor", "Lone Cabbage", "NA; nearest is Lone Cabbage", "Lemon Bay", "Hendry Creek/Mullock Creek", "New Pass", "Big Hickory", "NA; nearest is Big Hickory", "Jack Island", "Pellicer Creek", "Matanzas River", "Guana River")),
                    locality2 = factor(locality2, levels = c("Little St. George Island", "Goose Island/East Cove", "Alligator Harbor", "Lone Cabbage", "Lemon Bay", "Hendry Creek/Mullock Creek", "New Pass", "Big Hickory", "Jack Island", "Pellicer Creek", "Matanzas River", "Guana River")))]

setorder(hobsplotdat, locality2, nearby)

for(ll in unique(hobsplotdat$locality2)){
  for(ss in hobsplotdat[locality2 == ll, unique(spec_field_id)]){
    hobsplotdat[locality2 == ll & spec_field_id == ss, spec := as.character(which(hobsplotdat[locality2 == ll, unique(spec_field_id)] == ss))]
  }
}



# #By locality
# s4_plots <- list()
# for(ll in unique(hobsplotdat$locality2)){
#   futurespecs <- hobsplotdat[locality2 == ll & (spec_median > 2023 | spec_median_nocorr > 2023 | spec_median_rectest > 2023), .(y_pos = ifelse(spec_median > 2023, f14c_corr - 0.009,
#                                                                                                                                              ifelse(spec_median_nocorr > 2023, f14c - 0.009, f14c_corr_rectest - 0.009)),
#                                                                                                                                label = ifelse(spec_median > 2023, spec_median,
#                                                                                                                                               ifelse(spec_median_nocorr > 2023, spec_median_nocorr, spec_median_rectest))), by = list(locality2, spec_field_id)]
#   futurespecs[, x_pos := 2033]
#   
#   pastspecs <- hobsplotdat[locality2 == ll & (spec_median < 1751 | spec_median_nocorr < 1751 | spec_median_rectest < 1751), .(y_pos = ifelse(spec_median < 1751, f14c_corr - 0.009,
#                                                                                                                                                ifelse(spec_median_nocorr < 1751, f14c - 0.009, f14c_corr_rectest - 0.009)),
#                                                                                                                               label = ifelse(spec_median < 1751, spec_median,
#                                                                                                                                              ifelse(spec_median_nocorr < 1751, spec_median_nocorr, spec_median_rectest))), by = list(locality2, spec_field_id)]
#   pastspecs[, `:=` (x_pos = 1751,
#                     y_pos = y_pos - (0.0091 * (which(pastspecs$y_pos == y_pos) - 1))), by = y_pos]
#   
#   
#   oorspec_med <- hobsplotdat[locality2 == ll & is.na(spec_median), .(y_pos = unique(f14c_corr),
#                                                                      x_pos = year(collection_date) + 5)]
#   oorspec_nc <- hobsplotdat[locality2 == ll & is.na(spec_median_nocorr), .(y_pos = unique(f14c),
#                                                                            x_pos = year(collection_date) + 5)]
#   oorspec_rt <- hobsplotdat[locality2 == ll & is.na(spec_median_rectest), .(y_pos = unique(f14c_corr_rectest),
#                                                                             x_pos = year(collection_date) + 5)]
#   oorspec <- bind_rows(oorspec_med, oorspec_nc, oorspec_rt)
#   
#   nspec <- hobsplotdat[locality2 == ll, max(spec)]
#   
#   minyr <- min(unlist(hobsplotdat[locality2 == ll, .(min(year(collection_date)), min(spec_median), min(spec_median_nocorr), min(spec_median_rectest))]), na.rm = TRUE)
#   minyr <- ifelse(minyr < 1751, 1751,
#                   ifelse(minyr > 2000, 2000, minyr))
#   
#   plot_ll <- ggplot() +
#     geom_ribbon(data = predictedCurve[year >= minyr, ], 
#                 aes(x = year, ymin = f14c_mn - f14c_sd, ymax = f14c_mn + f14c_sd), fill = "grey50", alpha = 0.5) +
#     geom_line(data = predictedCurve[year >= minyr, ], 
#               aes(x = year, y = f14c_mn), lwd = 0.25) + 
#     geom_segment(data = hobsplotdat[locality2 == ll, ], aes(x = year(collection_date), y = f14c, xend = spec_median_nocorr, yend = f14c, color = "None"), lwd = 0.75, arrow = arrow(length = unit(5, "points"), type = "closed", angle = 30), show.legend = FALSE) + 
#     geom_segment(data = hobsplotdat[locality2 == ll, ], aes(x = year(collection_date), y = f14c_corr, xend = spec_median, yend = f14c_corr, color = "Collection\nyear ref."), lwd = 0.5, arrow = arrow(length = unit(4, "points"), type = "closed", angle = 30), show.legend = FALSE) + 
#     geom_segment(data = hobsplotdat[locality2 == ll, ], aes(x = year(collection_date), y = f14c_corr_rectest, xend = spec_median_rectest, yend = f14c_corr_rectest, color = "Recent-only ref."), lwd = 0.25, arrow = arrow(length = unit(3, "points"), type = "closed", angle = 30), show.legend = FALSE) + 
#     # geom_errorbar(data = hobsplotdat[locality2 == ll, ], aes(x = year(collection_date), ymin = f14c - f14c_sd, ymax = f14c + f14c_sd, fill = "None"), lwd = 0.75, width = 3) +
#     # geom_errorbar(data = hobsplotdat[locality2 == ll, ], aes(x = year(collection_date), ymin = f14c_corr - f14c_corr_sd, ymax = f14c_corr + f14c_corr_sd, fill = "Collection\nyear ref."), lwd = 0.5, width = 2) + 
#     # geom_errorbar(data = hobsplotdat[locality2 == ll, ], aes(x = year(collection_date), ymin = f14c_corr_rectest - f14c_corr_sd_rectest, ymax = f14c_corr_rectest + f14c_corr_sd_rectest, fill = "Recent-only ref."), lwd = 0.25, width = 1) + 
#     scale_color_manual(values = c("None" = "gold2", "Collection\nyear ref." = "firebrick3", "Recent-only ref." = "dodgerblue3")) +
#     # scale_fill_manual(values = c("None" = "gold2", "Collection\nyear ref." = "firebrick3", "Recent-only ref." = "dodgerblue3")) +
#     # scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25, "6" = 3)) +
#     # scale_size_manual(values = c("None" = 4, "Collection\nyear ref." = 3, "Recent-only ref." = 2)) +
#     new_scale_color() +
#     geom_point(data = hobsplotdat[locality2 == ll, ], aes(x = year(collection_date), y = f14c, fill = "None", color = nearby, size = "None", shape = spec, stroke = 0.8)) + 
#     geom_point(data = hobsplotdat[locality2 == ll, ], aes(x = year(collection_date), y = f14c_corr, fill = "Collection\nyear ref.", color = nearby, size = "Collection\nyear ref.", shape = spec, stroke = 0.7)) + 
#     geom_point(data = hobsplotdat[locality2 == ll, ], aes(x = year(collection_date), y = f14c_corr_rectest, fill = "Recent-only ref.", color = nearby, size = "Recent-only ref.", shape = spec, stroke = 0.6)) +
#     {if(nrow(futurespecs) > 0){
#       geom_text(data = futurespecs, aes(x = x_pos, y = y_pos, label = label), hjust = "right", color = "black", family = "Arial", size = 1.5)
#     }} +
#     {if(nrow(pastspecs) > 0){
#       geom_text(data = pastspecs, aes(x = x_pos, y = y_pos, label = label), hjust = "left", color = "black", family = "Arial", size = 1.5)
#     }} +
#     {if(nrow(oorspec) > 0){
#       geom_text(data = oorspec, aes(x = x_pos, y = y_pos, label = "oor"), hjust = "left", vjust = "center", color = "black", family = "Arial", size = 1.5)
#     }} +
#     {if(nspec == "6"){
#       list(geom_point(data = hobsplotdat[locality2 == ll & spec == "6", ], aes(x = year(collection_date), y = f14c, size = "None", shape = spec, stroke = 0.5), color = "gold2"),
#            geom_point(data = hobsplotdat[locality2 == ll & spec == "6", ], aes(x = year(collection_date), y = f14c_corr, size = "Collection\nyear ref.", shape = spec, stroke = 0.4), color = "firebrick3"),
#            geom_point(data = hobsplotdat[locality2 == ll & spec == "6", ], aes(x = year(collection_date), y = f14c_corr_rectest, size = "Recent-only ref.", shape = spec, stroke = 0.3), color = "dodgerblue3"))
#     }} +
#     theme_classic() +
#     scale_color_manual(values = c("Nearby" = "grey90", "Locality" = "black")) +
#     scale_fill_manual(values = c("None" = "gold2", "Collection\nyear ref." = "firebrick3", "Recent-only ref." = "dodgerblue3")) +
#     scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25, "6" = 3)) +
#     scale_size_manual(values = c("None" = 3, "Collection\nyear ref." = 2, "Recent-only ref." = 1)) +
#     scale_x_continuous(n.breaks = 4) +
#     scale_y_continuous(n.breaks = 4) +
#     labs(x = "Year", y = expression(paste(F^14, "C")), color = "Dead C\nCorrection", shape = "Specimen", title = ll) +
#     {if(nspec == "6"){
#       guides(size = "none", 
#              color = "none", 
#              fill = guide_legend(title = "Dead C\nCorrection", override.aes = list(shape = 15, size = 4, color = c("firebrick3", "gold2", "dodgerblue3"))),
#              shape = guide_legend(title = "Specimen", override.aes = list(color = "black")))
#     } else{
#       guides(size = "none", 
#              color = "none", 
#              fill = guide_legend(title = "Dead C\nCorrection", override.aes = list(shape = 15, size = 4, color = c("firebrick3", "gold2", "dodgerblue3"))),
#              shape = "none")
#     }} +
#     coord_cartesian(xlim = c(minyr, 2023)) +
#     theme(axis.title = element_text(size = 7, colour = "black", family = "Arial"), 
#           axis.text = element_text(size = 7, colour = "black", family = "Arial"), 
#           legend.text = element_text(size = 7, colour = "black", family = "Arial"), 
#           legend.title = element_text(size = 7, colour = "black", family = "Arial"),
#           text = element_text(size = 7, colour = "black", family = "Arial"))
#   
#   s4_plots[[which(unique(hobsplotdat$locality2) == ll)]] <- plot_ll
# }

#By specimen
s4_plots <- list()
for(ll in unique(hobsplotdat$cat_no)){
  futurespecs <- hobsplotdat[cat_no == ll & (spec_median > 2023 | spec_median_nocorr > 2023 | spec_median_rectest > 2023), .(y_pos = ifelse(spec_median > 2023, f14c_corr - 0.009,
                                                                                                                                            ifelse(spec_median_nocorr > 2023, f14c - 0.009, f14c_corr_rectest - 0.009)),
                                                                                                                             label = ifelse(spec_median > 2023, spec_median,
                                                                                                                                            ifelse(spec_median_nocorr > 2023, spec_median_nocorr, spec_median_rectest))), by = list(locality2, spec_field_id)]
  futurespecs[, x_pos := 2023]
  
  # pastspecs <- hobsplotdat[cat_no == ll & (spec_median < 1751 | spec_median_nocorr < 1751 | spec_median_rectest < 1751), .(y_pos = ifelse(spec_median < 1751, f14c_corr - 0.009,
  #                                                                                                                                            ifelse(spec_median_nocorr < 1751, f14c - 0.009, f14c_corr_rectest - 0.009)),
  #                                                                                                                             label = ifelse(spec_median < 1751, spec_median,
  #                                                                                                                                            ifelse(spec_median_nocorr < 1751, spec_median_nocorr, spec_median_rectest))), by = list(locality2, spec_field_id)]
  # pastspecs[, `:=` (x_pos = 1751,
  #                   y_pos = y_pos - (0.0091 * (which(pastspecs$y_pos == y_pos) - 1))), by = y_pos]
  
  
  oorspec_med <- hobsplotdat[cat_no == ll & is.na(spec_median), .(y_pos = unique(f14c_corr),
                                                                  x_pos = year(collection_date) + 2)]
  oorspec_nc <- hobsplotdat[cat_no == ll & is.na(spec_median_nocorr), .(y_pos = unique(f14c),
                                                                        x_pos = year(collection_date) + 2)]
  oorspec_rt <- hobsplotdat[cat_no == ll & is.na(spec_median_rectest), .(y_pos = unique(f14c_corr_rectest),
                                                                         x_pos = year(collection_date) + 2)]
  oorspec <- bind_rows(oorspec_med, oorspec_nc, oorspec_rt)
  
  nspec <- hobsplotdat[cat_no == ll, max(spec)]
  
  minyr <- min(unlist(hobsplotdat[cat_no == ll, .(min(year(collection_date)), min(spec_median), min(spec_median_nocorr), min(spec_median_rectest))]), na.rm = TRUE)
  minyr <- ifelse(minyr > 1970, 1970, minyr)
  # minyr <- ifelse(minyr < 1751, 1751,
  #                 ifelse(minyr > 2000, 2000, minyr))
  
  plot_ll <- ggplot(data = hobsplotdat[cat_no == ll, ]) +
    geom_ribbon(data = ncm20_dat[Year_ce >= 1950, ], #[year >= minyr, ]
                aes(x = Year_ce, ymin = F14C - F14C_sd, ymax = F14C + F14C_sd), fill = "grey50", alpha = 0.5) +
    geom_line(data = ncm20_dat[Year_ce >= 1950, ], #
              aes(x = Year_ce, y = F14C), lwd = 0.25) + 
    geom_ribbon(data = ncm20_dat[Year_ce >= minyr & Year_ce <= 1950, ], #Year_ce >= minyr & 
                aes(x = Year_ce, ymin = f14c_calc - f14c_calc_sd, ymax = f14c_calc + f14c_calc_sd), fill = "grey50", alpha = 0.5) +
    geom_line(data = ncm20_dat[Year_ce >= minyr & Year_ce <= 1950, ], # 
              aes(x = Year_ce, y = f14c_calc), lwd = 0.25) +
    geom_segment(aes(x = year(collection_date), y = f14c, xend = spec_median_nocorr, yend = f14c, color = "None"), lwd = 0.75, arrow = arrow(length = unit(5, "points"), type = "closed", angle = 30), show.legend = FALSE) + 
    geom_segment(aes(x = year(collection_date), y = f14c_corr, xend = spec_median, yend = f14c_corr, color = "Collection\nyear ref."), lwd = 0.5, arrow = arrow(length = unit(4, "points"), type = "closed", angle = 30), show.legend = FALSE) + 
    geom_segment(aes(x = year(collection_date), y = f14c_corr_rectest, xend = spec_median_rectest, yend = f14c_corr_rectest, color = "Recent-only ref."), lwd = 0.25, arrow = arrow(length = unit(3, "points"), type = "closed", angle = 30), show.legend = FALSE) + 
    # geom_errorbar(data = hobsplotdat[cat_no == ll, ], aes(x = year(collection_date), ymin = f14c - f14c_sd, ymax = f14c + f14c_sd, fill = "None"), lwd = 0.75, width = 3) +
    # geom_errorbar(data = hobsplotdat[cat_no == ll, ], aes(x = year(collection_date), ymin = f14c_corr - f14c_corr_sd, ymax = f14c_corr + f14c_corr_sd, fill = "Collection\nyear ref."), lwd = 0.5, width = 2) + 
    # geom_errorbar(data = hobsplotdat[cat_no == ll, ], aes(x = year(collection_date), ymin = f14c_corr_rectest - f14c_corr_sd_rectest, ymax = f14c_corr_rectest + f14c_corr_sd_rectest, fill = "Recent-only ref."), lwd = 0.25, width = 1) + 
    scale_color_manual(values = c("None" = "gold2", "Collection\nyear ref." = "firebrick3", "Recent-only ref." = "dodgerblue3")) +
    # scale_fill_manual(values = c("None" = "gold2", "Collection\nyear ref." = "firebrick3", "Recent-only ref." = "dodgerblue3")) +
    # scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25, "6" = 3)) +
    # scale_size_manual(values = c("None" = 4, "Collection\nyear ref." = 3, "Recent-only ref." = 2)) +
    # new_scale_color() +
    geom_point(aes(x = year(collection_date), y = f14c, fill = "None", size = "None", stroke = 0.8), shape = 21, color = "grey10") + 
    geom_point(aes(x = year(collection_date), y = f14c_corr, fill = "Collection\nyear ref.", size = "Collection\nyear ref.", stroke = 0.7), shape = 21, color = "grey10") + 
    geom_point(aes(x = year(collection_date), y = f14c_corr_rectest, fill = "Recent-only ref.", size = "Recent-only ref.", stroke = 0.6), shape = 21, color = "grey10") +
    {if(nrow(futurespecs) > 0){
      geom_text(data = futurespecs, aes(x = x_pos, y = y_pos, label = label), hjust = "right", color = "black", family = "Arial", size = 1.5)
    }} +
    # {if(nrow(pastspecs) > 0){
    #   geom_text(data = pastspecs, aes(x = x_pos, y = y_pos, label = label), hjust = "left", color = "black", family = "Arial", size = 1.5)
    # }} +
    {if(nrow(oorspec) > 0){
      geom_text(data = oorspec, aes(x = x_pos, y = y_pos, label = "oor"), hjust = "left", vjust = "center", color = "black", family = "Arial", size = 1.5)
    }} +
    # {if(nspec == "6"){
    #   list(geom_point(data = hobsplotdat[cat_no == ll & spec == "6", ], aes(x = year(collection_date), y = f14c, size = "None", stroke = 0.5), color = "gold2"),
    #        geom_point(data = hobsplotdat[cat_no == ll & spec == "6", ], aes(x = year(collection_date), y = f14c_corr, size = "Collection\nyear ref.", stroke = 0.4), color = "firebrick3"),
    #        geom_point(data = hobsplotdat[cat_no == ll & spec == "6", ], aes(x = year(collection_date), y = f14c_corr_rectest, size = "Recent-only ref.", stroke = 0.3), color = "dodgerblue3"))
    # }} +
    theme_classic() +
    # scale_color_manual(values = c("Nearby" = "grey90", "Locality" = "black")) +
    scale_fill_manual(values = c("None" = "gold2", "Collection\nyear ref." = "firebrick3", "Recent-only ref." = "dodgerblue3")) +
    # scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25, "6" = 3)) +
    scale_size_manual(values = c("None" = 3, "Collection\nyear ref." = 2, "Recent-only ref." = 1)) +
    scale_x_continuous(n.breaks = 4) +
    scale_y_continuous(n.breaks = 4) +
    labs(x = "Year", y = expression(paste(F^14, "C")), color = "Dead C\nCorrection", shape = "Specimen", title = ll, subtitle = paste0("Nearest locality: ", hobsplotdat[cat_no == ll, locality2], 
                                                                                                                                       "<br> Year collected: ", hobsplotdat[cat_no == ll, year(collection_date)],
                                                                                                                                       "<br> Genus: *", hobsplotdat[cat_no == ll, genus], "*")) + 
    # {if(nspec == "6"){
    #   guides(size = "none", 
    #          color = "none", 
    #          fill = guide_legend(title = "Dead C\nCorrection", override.aes = list(shape = 15, size = 4, color = c("firebrick3", "gold2", "dodgerblue3"))),
    #          shape = guide_legend(title = "Specimen", override.aes = list(color = "black")))
    # } else{
      guides(size = "none",
             color = "none",
             fill = guide_legend(title = "Dead C\nCorrection", override.aes = list(shape = 21, size = c(2, 3, 1.25), color = "grey10", fill = c("firebrick3", "gold2", "dodgerblue3"))),
             shape = "none") +
    # }} +
    coord_cartesian(xlim = c(minyr, 2023)) +
    theme(axis.title = element_text(size = 7, colour = "black", family = "Arial"), 
          axis.text = element_text(size = 7, colour = "black", family = "Arial"), 
          legend.text = element_text(size = 7, colour = "black", family = "Arial"), 
          legend.title = element_text(size = 7, colour = "black", family = "Arial"),
          text = element_text(size = 7, colour = "black", family = "Arial"),
          plot.subtitle = ggtext::element_markdown())
  
  s4_plots[[which(unique(hobsplotdat$cat_no) == ll)]] <- plot_ll
}

# FigS4 <- wrap_plots(s4_plots[1:6], ncol = 2, guides = "collect") + plot_annotation(tag_levels = "A")
# FigS5 <- wrap_plots(s4_plots[7:12], ncol = 2, guides = "collect") + plot_annotation(tag_levels = "A")
FigS4 <- wrap_plots(s4_plots, ncol = 2, guides = "collect") + plot_annotation(tag_levels = "A")
# FigS4

# ggsave(here::here("Summer2023/Figures/FigS4.png"),
ggsave(here::here("Summer2023/Figures/FigS4.png"),
       FigS4,
       width = 40*(1/6), # 1 pica = ~1/6 inch; 3.5, #column width for 2-column page layout GSA Pub spec
       height = 6,
       units = "in")

# # ggsave(here::here("Summer2023/Figures/FigS5.png"),
# ggsave(here::here("Summer2023/Figures/FigS5.png"),
#        FigS5,
#        width = 42*(1/6), # 1 pica = ~1/6 inch; 3.5, #column width for 2-column page layout GSA Pub spec
#        height = 6,
#        units = "in")

#Geology paper figure 1 --------------------------------------------------------------------------
HOBSdat <- hobsf14c2[analysis == "LP" & !is.na(igsn), ]
HOBSdat[ , locality2 := ifelse(str_detect(locality, "Hendry"), "Hendry/\nMullock Creeks", locality)]
HOBSdat[ , locality2 := str_replace_all(locality2, "Island", "Is.")]
HOBSdat[ , locality2 := str_replace_all(locality2, "Creek", "Cr.")]
HOBSdat[ , locality2 := str_replace_all(locality2, "River", "R.")]
HOBSdat[ , locality2 := str_replace_all(locality2, "Cr.s", "Cr.")]
HOBSdat_sites <- HOBSdat[, .(locality = locality2, lat, lon)]
HOBSdat_sites <- unique(HOBSdat_sites)
HOBSdat_sites <- st_as_sf(HOBSdat_sites, coords = c("lon", "lat"), remove = F, crs = 4326)
fwcoy <- st_read(here::here("Oyster_Beds_in_Florida_05-25-2021/Oyster_Beds_in_Florida.shp"))
fwcoy_m <- st_transform(fwcoy, 4326)
rcp <- st_read(here::here("orcp_all_sites/ORCP_Managed_Areas.shp"))
rcp_m <- st_transform(rcp, 4326)
FL <- map_data("county", "florida") %>%
  select(lon = long, lat, group, id = subregion)
FL2 <- st_as_sf(FL, coords = c("lon", "lat"), remove = F, crs = 4326)

HOBSdat_loc <- HOBSdat_sites[0,]
for(i in unique(HOBSdat_sites$locality)){
  loc_i <- subset(HOBSdat_sites, HOBSdat_sites$locality == i)[1,]
  HOBSdat_loc <- rbind(HOBSdat_loc, loc_i)
}

world_map_data <- ne_countries(scale = "medium", returnclass = "sf")
world_map_data_m <- st_transform(world_map_data, 4326)
state_map_data <- map('state', fill = TRUE, plot = FALSE) %>% st_as_sf()
state_map_data_m <- st_transform(state_map_data, 4326)

HOBSdat_loc %>% 
  mutate(x_nudge = case_when( locality == "Goose Is./East Cove" ~ -0.3
                              ,locality == "Little St. George Is." ~ -1.7
                              ,locality == "Lone Cabbage" ~ -1.4
                              ,locality == "Lemon Bay" ~ -1.3
                              ,locality == "Hendry/\nMullock Cr." ~ -1.4
                              ,locality == "New Pass" ~ -1.4
                              ,locality == "Big Hickory" ~ -1.1
                              ,locality == "Jack Is." ~ 1.2
                              ,locality == "Pellicer Cr." ~ 1.5
                              ,locality == "Matanzas R." ~ 1.6
                              ,locality == "Guana R." ~ 1.6
                              ,TRUE ~ 0)
         ,y_nudge = case_when(locality == "Goose Is./East Cove" ~ -0.8
                              ,locality == "Little St. George Is." ~ -0.4
                              ,locality == "Lone Cabbage" ~ -1
                              ,locality == "Lemon Bay" ~ 0.2
                              ,locality == "Hendry/\nMullock Cr." ~ 0.2
                              ,locality == "New Pass" ~ -0.5
                              ,locality == "Big Hickory" ~ -1
                              ,locality == "Jack Is." ~ 0.2
                              ,locality == "Pellicer Cr." ~ -0.2
                              ,locality == "Matanzas R." ~ 0.1
                              ,locality == "Guana R." ~ 0.2
                              ,TRUE ~ 0)
  ) ->
  HOBSdat_loc


# HOBSmap <- ggplot(data = FL2, aes(lon, lat, group = group)) +
HOBSmap <- ggplot() +
  geom_sf(data = world_map_data, lwd = 0.25, inherit.aes = FALSE) +
  geom_sf(data = state_map_data, lwd = 0.25, inherit.aes = FALSE) +
  geom_sf(data = subset(state_map_data, state_map_data$ID == "florida"), fill = "antiquewhite1", lwd = 0.25, inherit.aes = FALSE) +
  geom_sf(data = rcp_m, aes(fill = "dodgerblue"), color = "dodgerblue", alpha = 0.5, lwd = 0.25, inherit.aes = FALSE) +
  geom_sf(data = fwcoy_m, aes(fill = "firebrick"), color = "firebrick", lwd = 0.5, inherit.aes = FALSE) +
  geom_text_repel(data = HOBSdat_loc, aes(lon, lat, label = locality), size = 7*0.35, nudge_x = HOBSdat_loc$x_nudge, nudge_y = HOBSdat_loc$y_nudge, segment.size = 0.25, min.segment.length = 0.4, inherit.aes = FALSE) +
  geom_point(data = HOBSdat_loc, aes(lon, lat), color = "black", fill = "yellow1", shape = 21, size = 1.5, stroke = 0.5, inherit.aes = FALSE) +
  coord_sf(xlim = c(-87.5, -79), ylim = c(24.5, 31)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightblue1"),
        legend.title = element_blank(),
        legend.position = c(0.17, 0.45),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(1,"line"),
        text = element_text(size = 7),
        plot.margin = unit(c(0, 0, 0, 0), "mm")) +
  labs(x = "Lat.", y = "Lon.") +
  scale_fill_manual(name = NULL,
                    labels = c("FDEP ORCP\nmanaged areas", "Mapped oyster\nreefs (not to scale)"),
                    values = c("dodgerblue", "firebrick"),
                    guide = guide_legend(override.aes = list(color = "grey10", lwd = 0.25)))

HOBSmap <- HOBSmap +
  ggspatial::annotation_scale(
    location = "br",
    text_cex = .pt/5.5,
    pad_x = unit(0.05, "in"), pad_y = unit(0.05, "in"),
    line_width = 0.5,
    height = unit(0.05, "in"),
    bar_cols = c("grey60", "white")) +
  ggspatial::annotation_north_arrow(
    location = "br", which_north = "true",
    pad_x = unit(0.1, "in"), pad_y = unit(0.15, "in"),
    height = unit(0.3, "in"), width = unit(0.3, "in"),
    style = ggspatial::north_arrow_fancy_orienteering(
      fill = c("grey40", "white"),
      line_col = "grey20"))

refbox <- data.table(lat = c(rep(HOBSmap$coordinates$limits$x[1], 2), rep(HOBSmap$coordinates$limits$x[2], 2)), lon = c(HOBSmap$coordinates$limits$y, HOBSmap$coordinates$limits$y[2], HOBSmap$coordinates$limits$y[1]))
refbox2 <- st_as_sf(refbox, coords = c("lat", "lon"), crs = st_crs(4326))
refbox3 <- st_combine(refbox2)
refbox4 <- st_cast(refbox3, "POLYGON")

outbox <- data.table(lat = c(-165, -165, -20, -20), lon = c(-2.5, 82.5, 82.5, -2.5))
outbox2 <- st_as_sf(outbox, coords = c("lat", "lon"), crs = st_crs(4326))
outbox3 <- st_combine(outbox2)
outbox4 <- st_cast(outbox3, "POLYGON")

inset <- ggplot() +
  geom_sf(data = subset(world_map_data, world_map_data$continent %in% c("North America", "South America")), lwd = 0.25, fill = "grey75", color = "grey75") +
  geom_sf(data = refbox4, lwd = 0.4, color = "black", fill = "transparent") +
  geom_sf(data = outbox4, lwd = 0.3, color = "black", fill = "transparent") +
  coord_sf(xlim = c(-159, -26.5), ylim = c(1, 78.5)) +
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.background = element_rect(color = "transparent", size = 0.5, fill = "white"),
        plot.background = element_blank())

layout = c(
  area(t = 1, l = 1, b = 100, r = 120),
  area(t = 65, l = 5, b = 100, r = 47)
)

HOBSmap2 <- HOBSmap + inset + plot_layout(design = layout)

ggsave(here::here("Summer2023/Figures/Fig1.pdf"),
       HOBSmap2,
       width = 21.5*(1/6), # 1 pica = ~1/6 inch; 3.5, #column width for 2-column page layout GSA Pub spec
       height = 3,
       units = "in",
       device = grDevices::cairo_pdf)

ggsave(here::here("Summer2023/Figures/Fig1.png"),
       HOBSmap2,
       width = 21.5*(1/6), # 1 pica = ~1/6 inch; 3.5, #column width for 2-column page layout GSA Pub spec
       height = 3,
       units = "in")




#Figure S3--------------------------------------------------------
HPspec <- hobsf14c2[analysis == "HP", unique(spec_field_id)]
LPvHP_dat <- hobsf14c2[spec_field_id %in% HPspec, ]
LPvHP_dat2 <- distinct(LPvHP_dat[spec_field_id %in% intersect(HPspec, LPvHP_dat[analysis == "LP", unique(spec_field_id)]) & !is.na(locality), .(analysis, locality, station, spec_field_id, f14c, f14c_sd, f14c_corr, f14c_corr_sd)])
LPvHP_dat2 <- pivot_wider(LPvHP_dat2, names_from = analysis, values_from = c(f14c, f14c_sd, f14c_corr, f14c_corr_sd))
setDT(LPvHP_dat2)
LPvHP_dat2 <- LPvHP_dat2[!is.na(f14c_HP), ]
LPvHP_mod <- lmodel2(f14c_LP ~ f14c_HP, data = LPvHP_dat2, range.y = "relative", range.x = "relative", nperm = 10000)

coefs <- LPvHP_mod$regression.results
names(coefs) <- c("method", "intercept", "slope", "angle", "p-value")
CIs <- LPvHP_mod$confidence.intervals
names(CIs) <- c("method", "intercept_low", "intercept_high", "slope_low", "slope_high")
m <- c("SMA")


# setDT(LPvHP_dat2)
LPvHP_dat2[, loc_lab := fcase(str_detect(locality, "George"), "Little St.\nGeorge Island", 
                              str_detect(locality, "Lemon"), "Lemon Bay", 
                              str_detect(locality, "Hendry"), "Hendry Creek/\nMullock Creek", 
                              str_detect(locality, "New"), "New Pass", 
                              str_detect(locality, "Jack"), "Jack Island", 
                              str_detect(locality, "Pellicer"), "Pellicer Creek", 
                              str_detect(locality, "Guana"), "Guana River")]
LPvHP_dat2[, loc_lab := factor(loc_lab, levels = c("Little St.\nGeorge Island", "Lemon Bay", "Hendry Creek/\nMullock Creek", "New Pass", "Jack Island", "Pellicer Creek", "Guana River"))]

# #Trying to set up the SMA regression lines as segments instead of ablines, but I will have to come back to this later if there is time. 
# setDT(CIs)
# setDT(coefs)
# seginfo <- data.table(seg = c("int_low", "int_high", "int"),
#                       x = rep(min(LPvHP_dat2$f14c_HP), 3),
#                       xend = rep(max(LPvHP_dat2$f14c_HP), 3))
# seginfo[, `:=` (y = c(CIs[method == m, slope_high] * min(LPvHP_dat2$f14c_HP) + CIs[method == m, intercept_low],
#                       CIs[method == m, slope_low] * min(LPvHP_dat2$f14c_HP) + CIs[method == m, intercept_high],
#                       coefs[method == m, slope] * min(LPvHP_dat2$f14c_HP) + coefs[method == m, intercept]),
#                 yend = c(CIs[method == m, slope_high] * max(LPvHP_dat2$f14c_HP) + CIs[method == m, intercept_low],
#                          CIs[method == m, slope_low] * max(LPvHP_dat2$f14c_HP) + CIs[method == m, intercept_high],
#                          coefs[method == m, slope] * max(LPvHP_dat2$f14c_HP) + coefs[method == m, intercept]))]
# 
# seginfo <- data.table(seg = "int",
#                       x = min(LPvHP_dat2$f14c_HP),
#                       y = coefs[method == m, slope] * min(LPvHP_dat2$f14c_HP) + coefs[method == m, intercept],
#                       xend = max(LPvHP_dat2$f14c_HP),
#                       yend = coefs[method == m, slope] * max(LPvHP_dat2$f14c_HP) + coefs[method == m, intercept])
# seginfo2 <- data.table(seg = c("int_low", "int_high"),
#                        x = c((coefs[method == m, slope] * min(LPvHP_dat2$f14c_HP) + coefs[method == m, intercept])/CIs[method == m, slope_high] - CIs[method == m, intercept_low],
#                              CIs[method == m, slope_low] * min(LPvHP_dat2$f14c_HP) + CIs[method == m, intercept_high]),
#                        y = c(CIs[method == m, slope_high] * min(LPvHP_dat2$f14c_HP) + CIs[method == m, intercept_low],
#                              (coefs[method == m, slope] * min(LPvHP_dat2$f14c_HP) + coefs[method == m, intercept])/CIs[method == m, slope_low] - CIs[method == m, intercept_high]),
#                        xend = c(CIs[method == m, slope_high] * max(LPvHP_dat2$f14c_HP) + CIs[method == m, intercept_low],
#                                 (coefs[method == m, slope] * max(LPvHP_dat2$f14c_HP) + coefs[method == m, intercept])/CIs[method == m, slope_low] - CIs[method == m, intercept_high]),
#                        yend = c((coefs[method == m, slope] * max(LPvHP_dat2$f14c_HP) + coefs[method == m, intercept])/CIs[method == m, slope_high] - CIs[method == m, intercept_low],
#                                 CIs[method == m, slope_low] * max(LPvHP_dat2$f14c_HP) + CIs[method == m, intercept_high]))
# seginfo <- rbind(seginfo, seginfo2)

fsize <- 7
LPvHPplot <- ggplot() +
  scale_color_brewer(palette = "Accent", aesthetics = "fill") +
  geom_point(data = LPvHP_dat2, aes(f14c_HP, f14c_LP, fill = loc_lab), shape = 21, color = "black", size = 2.5) +
  geom_abline(data = subset(CIs, CIs$method %in% m), aes(intercept = intercept_low, slope = slope_high), color = "grey50", lwd = 0.75, lty = 2) +
  # geom_segment(data = seginfo[seg != "int", ], aes(x = x, y = y, xend = xend, yend = yend), color = "grey50", lwd = 0.75, lty = 2) +
  geom_abline(data = subset(CIs, CIs$method %in% m), aes(intercept = intercept_high, slope = slope_low), color = "grey50", lwd = 0.75, lty = 2) +
  geom_abline(data = subset(coefs, coefs$method %in% m), aes(intercept = intercept, slope = slope), color = "black", lwd = 0.75) +
  # geom_segment(data = seginfo[seg == "int", ], aes(x = x, y = y, xend = xend, yend = yend), color = "black", lwd = 0.75) +
  coord_cartesian(xlim = c(0.97, 1.23), ylim = c(0.99, 1.21)) +
  theme_bw(base_size = fsize, base_family = "Arial") + 
  theme(axis.text = element_text(size = fsize),
        axis.title = element_text(size = fsize),
        legend.text = element_text(size = fsize), 
        legend.title = element_text(size = fsize)) +
  labs(x = expression(paste("HP ", F^14, "C")), y = expression(paste("LP ", F^14, "C")), fill = "Locality") +
  annotate(geom = "text", 
           x = 0.99, 
           y = 1.205,
           size = fsize*0.35,
           label = paste0("y = ", 
                          round(subset(coefs, coefs$method %in% m)$slope, 3), 
                          "x", 
                          ifelse(subset(coefs, coefs$method %in% m)$intercept > 0, 
                                 paste0(" + ", round(subset(coefs, coefs$method %in% m)$intercept, 3)), 
                                 paste0(" - ", round(abs(subset(coefs, coefs$method %in% m)$intercept), 3)))),
           hjust = 0) +
  annotate(geom = "text",
           x = 0.99,
           y = 1.191,
           size = fsize*0.35,
           label = paste0("R^2 == ", round(LPvHP_mod$rsquare, 3)),
           parse = TRUE,
           hjust = 0)

# ggsave(here::here("Summer2023/Figures/FigS1.pdf"),
#        LPvHPplot,
#        #width = 7.125, #full page width GSA Pub spec
#        #height = 5,
#        width = 5,
#        height = 3,
#        units = "in",
#        device = grDevices::cairo_pdf)

ggsave(here::here("Summer2023/Figures/FigS3.png"),
       LPvHPplot,
       width = 5,
       height = 3,
       units = "in")



#Proportion of samples with pre/post-bomb median calibrated ages---------------------------------------------
prebomb <- hobsf14c2[oor == FALSE & analysis == "LP" & !is.na(tav_median) & tav_median <= 1950, length(unique(igsn))]/hobsf14c2[oor == FALSE & analysis == "LP" & !is.na(tav_median), length(unique(igsn))]
postbomb <- hobsf14c2[oor == FALSE & analysis == "LP" & !is.na(tav_median) & tav_median > 1950, length(unique(igsn))]/hobsf14c2[oor == FALSE & analysis == "LP" & !is.na(tav_median), length(unique(igsn))]


#Proportion of samples with sub-decadal time averaging------------------------------------------------------
TAprop_yr <- hobsf14c2[oor == FALSE & analysis == "LP" & !is.na(tav_median) & tav_iqr < 10, length(unique(igsn))]/hobsf14c2[oor == FALSE & analysis == "LP" & !is.na(tav_median), length(unique(igsn))]
TAprop_dec <- hobsf14c2[oor == FALSE & analysis == "LP" & !is.na(tav_median) & tav_iqr >= 10 & tav_iqr < 100, length(unique(igsn))]/hobsf14c2[oor == FALSE & analysis == "LP" & !is.na(tav_median), length(unique(igsn))]
TAprop_cen <- hobsf14c2[oor == FALSE & analysis == "LP" & !is.na(tav_median) & tav_iqr >= 100 & tav_iqr < 1000, length(unique(igsn))]/hobsf14c2[oor == FALSE & analysis == "LP" & !is.na(tav_median), length(unique(igsn))]
TAprop_mil <- hobsf14c2[oor == FALSE & analysis == "LP" & !is.na(tav_median) & tav_iqr >= 1000, length(unique(igsn))]/hobsf14c2[oor == FALSE & analysis == "LP" & !is.na(tav_median), length(unique(igsn))]


#Proportion of sample holes where the expected stratigraphic order was recovered----------------------------
leveldiffs <- data.table(station = numeric(), median_1 = numeric(), median_2 = numeric(), median_diff = numeric(), tav_1 = numeric(), tav_2 = numeric())
for(i in hobsf14c2[oor == FALSE & analysis == "LP" & depth_int != 0 & !is.na(station), unique(station)]){
  holedat <- distinct(hobsf14c2[oor == FALSE & analysis == "LP" & depth_int != 0 & station == i, .(station, depth_int, tav_median, tav_iqr)])
  if(nrow(holedat) < 2) next
  leveldiffs_i <- data.table(station = i,
                             median_1 = holedat[depth_int == 1, tav_median],
                             median_2 = holedat[depth_int == 2, tav_median],
                             median_diff = holedat[depth_int == 1, tav_median] - holedat[depth_int == 2, tav_median],
                             tav_1 = holedat[depth_int == 1, tav_iqr],
                             tav_2 = holedat[depth_int == 2, tav_iqr])
  leveldiffs <- rbind(leveldiffs, leveldiffs_i)
}

for(i in 1:nrow(leveldiffs)){
  leveldiffs$Order[i] <- ifelse(leveldiffs$median_diff[i] >= 0, "Expected", "Unexpected")
}

write.csv(leveldiffs, here::here("Summer2023/BurialDepthDiffs.csv"))
Stratprop_pos <- nrow(leveldiffs[!is.na(median_diff) & median_diff >= 0, ])/nrow(leveldiffs[!is.na(median_diff), ])
Stratprop_neg <- nrow(leveldiffs[!is.na(median_diff) & median_diff < 0, ])/nrow(leveldiffs[!is.na(median_diff), ])


#Medians BP & TAV-------------------------------------------------------------------------------------------
medbp <- distinct(hobsf14c2[analysis == "LP" & oor == FALSE, .(loc, reef_ind, depth_int, igsn, median_age_bp, tav_iqr)])

ggplot(medbp, aes(x = median_age_bp, fill = factor(depth_int))) +
  geom_histogram(alpha = 0.5) +
  theme_bw()

medbp[, .(median = median(median_age_bp), mean = mean(median_age_bp), sd = sd(median_age_bp)), by = depth_int]
medbp[, .(median = median(tav_iqr), mean = mean(tav_iqr), sd = sd(tav_iqr)), by = depth_int]
medbp[, .(median_alldepths = median(tav_iqr), mean_alldepths = mean(tav_iqr), sd_alldepths = sd(tav_iqr))]
range(medbp$tav_iqr)

minmax <- hobsf14c2[analysis == "LP" & oor == FALSE, .(loc, depth_int, reef_num, tav_median, tav_iqr)]
minmax[, `:=` (minmed = min(tav_median), maxmed = max(tav_median), mintav = min(tav_iqr), maxtav = max(tav_iqr)), by = list(loc, reef_num)]
minmax <- distinct(minmax)
minmax[, `:=` (meddiff = maxmed - minmed, tavdiff = maxtav - mintav)]
minmax[, .(mean_meddiff = mean(meddiff), sd_meddif = sd(meddiff), mean_tavdiff = mean(tavdiff), sd_tavdiff = sd(tavdiff))]

# #Based on data for Appendix DR1
# dr1_dat2 <- distinct(hobsf14c2[!is.na(tav_iqr), .(locality, reef, station, lat, lon, igsn, depth_int, n, tav_median = plyr::round_any(tav_median, 0.1), tav_min = plyr::round_any(tav_min, 0.1), tav_max = plyr::round_any(tav_max, 0.1), tav_iqr = plyr::round_any(tav_iqr, 0.01))])
# diffs <- dr1_dat2[, `:=` (meddiff = range(tav_median)[2] - range(tav_median)[1], tavdiff = range(tav_iqr)[2] - range(tav_iqr)[1]), by = list(locality, reef, depth_int)]
# # diffs[, .(meanmed = mean(meddiff), sdmed = sd(meddiff), meantav = mean(tavdiff), sdtav = sd(tavdiff)), by = depth_int]
# diffs[, .(meanmed = mean(meddiff), sdmed = sd(meddiff), meantav = mean(tavdiff), sdtav = sd(tavdiff))]


#Geology paper figure 3-------------------------------------------------------------------------------
vlines2 <- data.frame(depth = c('15-25cm', '25-35cm'), 
                      xint_mean = c(hobsf14c2[!is.na(igsn) & analysis == "LP" & oor == FALSE & depth_int == 1, mean(median_age_bp)], 
                                    hobsf14c2[!is.na(igsn) & analysis == "LP" & oor == FALSE & depth_int == 2, mean(median_age_bp)]),
                      xint_med = c(hobsf14c2[!is.na(igsn) & analysis == "LP" & oor == FALSE & depth_int == 1, median(median_age_bp)], 
                                   hobsf14c2[!is.na(igsn) & analysis == "LP" & oor == FALSE & depth_int == 2, median(median_age_bp)]))

Ketal <- data.frame(median_age_bp = 1058, tav_iqr = 1910)
label1 <- c(expression(atop(NA, atop(textstyle('Kowalewski et al. 2018'), 
                                     textstyle(paste('(', italic('Tucetona pectinata'), ')'))))))

#Changed my mind - I think all ages in this figure should be relative to the sample collection date
# #Correct Ketal ages for difference in "BP" - paper is relative to 2016, but our study is relative to 2019
# setDT(Ketal)
# Ketal[, median_age_bp := median_age_bp + 3]


Detal <- data.frame(median_age_bp = c(123, 283, 146, 148, 160, 2043), 
                    tav_iqr = c(130, 2766, 84, 96, 89, 2204))
label2 <- c(expression(atop(NA, atop(textstyle('Dominguez et al. 2016'), 
                                     textstyle(paste('(', italic('Fulvia tenuicostata'), ')'))))))

#Changed my mind - I think all ages in this figure should be relative to the sample collection date
# #Correct Detal ages for difference in "BP" - paper is relative to 2013, but our study is relative to 2019
# setDT(Detal)
# Detal[, median_age_bp := median_age_bp + 6]

Retal <- data.frame(Site = c(1, 2, 3),
                    median_age_bp = c(2658, 1318, 230),
                    tav_iqr = c(1914, 921, 200))
label3 <- c(expression(atop(NA, atop(textstyle('Ritter et al. 2017'), 
                                     textstyle(paste('(', italic('Mactra sp.'), ')'))))))

#Changed my mind - I think all ages in this figure should be relative to the sample collection date, but Retal still needs to reflect collection dates
# #Correct Retal ages for difference in "BP" - paper is relative to 1950, with Site 1 collected in 2013 and Sites 2 and 3 collected in 1974, but our study is relative to 2022. Correct dates relative to collection date.
setDT(Retal)
Retal[, median_age_bp := ifelse(Site == 1, median_age_bp + 63, median_age_bp + 24)]


Tetal1 <- readRDS(here::here("Tomasovych2018_medageandta.rds"))
label4 <- c(expression(atop(NA, atop(textstyle('Toma\u161ov\uFDch et al. 2018'),
                                     textstyle(paste0('(', italic('Corbula gibba'), ')'))))))
#Note that Tomasovych et al. (2018) use the 2.5% and 97.5% quantiles of the median specimen ages as their "TAV" estimate and the IQR as 
#their "TAV_IQR" estimate. Using our samples with >20 specimens dated, I show that the estimation methods are fairly close, with the 
#Tomasovych TAV method overestimating our TAV_IQR by about 1 to 50 years and their TAV_IQR method underestimating our TAV_IQR by about
#12 to 67 years. I think the comparison will probably work for our purposes.
specmedsiqr <- hobsf14c2[igsn %in% c("IEPRI005E", "IEPRI005Q", "IEPRI006M","IEPRI00BE") & analysis == "LP" & oor == FALSE, .(specmedstav = quantile(spec_median, 0.975) - quantile(spec_median, 0.025),
                                                                                                                             specmedsiqr = quantile(spec_median, 0.75) - quantile(spec_median, 0.25)), by = igsn]
specmedsiqr <- hobsf14c2[!is.na(igsn) & analysis == "LP" & oor == FALSE, .(n = length(unique(spec_field_id)),
                                                                           specmedstav = quantile(spec_median, 0.975) - quantile(spec_median, 0.025),
                                                                           specmedsiqr = quantile(spec_median, 0.75) - quantile(spec_median, 0.25)), by = igsn]
taviqr <- hobsf14c2[igsn %in% c("IEPRI005E", "IEPRI005Q", "IEPRI006M","IEPRI00BE") & analysis == "LP" & oor == FALSE, .(tav_iqr = unique(tav_iqr)), by = igsn]
taviqr <- hobsf14c2[!is.na(igsn) & analysis == "LP" & oor == FALSE, .(n = length(unique(spec_field_id)),
                                                                      tav_iqr = unique(tav_iqr)), by = igsn]
tavcompare <- merge(specmedsiqr, taviqr, by = c("igsn", "n"), all = TRUE)
setorder(tavcompare, -n)

ggplot(tavcompare) +
  geom_point(aes(x = n, y = specmedsiqr - tav_iqr, fill = igsn, shape = "specmeds_iqr"), color = "black", alpha = 0.6) +
  geom_point(aes(x = n, y = specmedstav - tav_iqr, fill = igsn, shape = "specmeds_tav"), color = "black", alpha = 0.6) +
  theme_bw() +
  labs(x = "n specimens", y = "years difference (Tomasovych method - our TAV_IQR)", shape = "Tomasovych\nmethod") +
  scale_shape_manual(values = c("specmeds_iqr" = 21, "specmeds_tav" = 24)) +
  guides(fill = "none")

ggplot(tavcompare) +
  geom_histogram(aes(x = specmedsiqr - tav_iqr, fill = "specmeds_iqr"), color = "black", alpha = 0.6) +
  geom_histogram(aes(x = specmedstav - tav_iqr, fill = "specmeds_tav"), color = "black", alpha = 0.6) +
  theme_bw() +
  labs(y = "n samples", x = "years difference (Tomasovych method - our TAV_IQR)", fill = "Tomasovych\nmethod") +
  scale_fill_manual(values = c("specmeds_iqr" = "firebrick", "specmeds_tav" = "dodgerblue"))

#Changed my mind - I think all ages in this figure should be relative to the sample collection date
# #Correct Tetal ages for difference in "BP" - paper is relative to 2013, but our study is relative to 2019
# setDT(Tetal)
# Tetal[, medage := medage + 6]

#USING Tomasovych et al. (2019) data isn't straightforward; they report the specimen ages in their supplementary materials, but not the age estimate IQRs or the core interval-level IQRs, so I did not incorporate them into this figure.
Tetal2 <- openxlsx::read.xlsx(here::here("TomasovychA-etal_2019/Supplementary Table 01.xlsx"), sheet = "AAR-Parvilucina-replicates_aver")
setDT(Tetal2)
Tetal2[, `:=` (trueupperdepth = as.numeric(`True.upper.depth.(cm)`), truelowerdepth = as.numeric(`True.lower.depth.(cm)`), agebpcoll = `Calibrated.shell.age.(years.before.AAR.dating)` - (`Year.of.amino.acid.dating.(BC/AD)` - `Year.of.shell.collection.(BC/AD)`))]
label5 <- c(expression(atop(NA, atop(textstyle('Toma\u161ov\uFDch et al. 2019'),
                                     textstyle(paste0('(', italic('Parvilucina tenuisculpta \\& Nuculana taphria'), ')'))))))


log10_rev_trans <- scales::trans_new(
  "log10_rev",
  transform = function(x) -1 * (log10(x)),
  inverse = function(x) 10 ^ abs(x),
  scales::log_breaks(5),
  domain = c(1e-100, Inf)
)

nonhobsfill <- "antiquewhite"
nonhobscol <- "black"

TAscalecompare <- ggplot(hobsf14c2[!is.na(igsn) & analysis == "LP" & oor == FALSE, ], aes(x = median_age_bp, y = tav_iqr)) +
  theme_classic(base_size = 7, base_family = "Arial") +
  scale_fill_discrete_qualitative(palette = "Harmonic", labels = c("15-25", "25-35")) +
  scale_color_discrete_qualitative(palette = "Harmonic", labels = c("15-25", "25-35")) +
  geom_vline(xintercept = c(seq(100, 900, 100), 1500, 2500), lwd = 0.25, color = "grey80") +
  # geom_vline(xintercept = c(seq(1, 9, 1), seq(10, 90, 10), seq(100, 900, 100), seq(1000, 5000, 1000)), lwd = 0.25, color = "grey70") +
  # geom_vline(xintercept = c(1, 10, 100, 1000), lwd = 0.25, color = "grey40") +
  geom_vline(xintercept = c(0, 500, 1000, 2000, 3000), lwd = 0.25, color = "grey40") +
  geom_hline(yintercept = c(1, 10, 100, 1000), lwd = 0.25) +
  geom_vline(data = vlines2, aes(xintercept = xint_med, color = depth), lwd = 0.5, show.legend = FALSE) +
  # geom_point(aes(fill = ordered(depth_int), shape = "Durham et al. (2023)"), color = "black", size = 1.5) +
  geom_point(aes(fill = ordered(depth_int), shape = "This study"), color = "black", size = 1.5) + #during copy editing the journal requested we change this legend label back to "This study"
  geom_point(data = Ketal, aes(shape = "Kowalewski et al. (2018)"), fill = nonhobsfill, color = nonhobscol, size = 1.5) +
  #annotate('text', x = Ketal$median_age_bp[1] - 250, y = Ketal$Sample_corrected_posterior_age_estimate[1] + 1500, 
  #         label = label1) +
  geom_point(data = Detal, aes(shape = "Dominguez et al. (2016)"), fill = nonhobsfill, color = nonhobscol, size = 1.5) +
  #annotate('text', x = Detal$median_age_bp[6] - 250, y = Detal$Sample_corrected_posterior_age_estimate[6] + 1700, 
  #         label = label2) +
  geom_point(data = Retal, aes(shape = "Ritter et al. (2017)"), fill = nonhobsfill, color = nonhobscol, size = 1.5) +
  geom_point(data = Tetal1[depth >= 10 & depth <= 35, ], aes(x = medage, y = IQR, shape = "Toma\u{161}ov\u{FD}ch et al. (2018)"), fill = nonhobsfill, color = nonhobscol, size = 1.5) +
  scale_y_continuous(breaks = c(1, 10, 100, 1000), trans='log10') +
  # scale_x_continuous(breaks = c(1, 10, 100, 1000), trans = c("log10", "reverse")) +
  scale_x_reverse(breaks = c(3000, 2000, 1000, 500, 0)) +
  coord_cartesian(xlim = c(3000, 1)) +
  labs(x = 'Median age (years before collection year)',
       y = 'Time-averaging estimate (yr)',
       fill = 'Burial depth (cm)',
       shape = "") + 
  # scale_shape_manual(values = c("Durham et al. (2023)" = 21, "Kowalewski et al. (2018)" = 22, "Dominguez et al. (2016)" = 23, "Ritter et al. (2017)" = 24, "Toma\u{161}ov\u{FD}ch et al. (2018)" = 25),
  #                    breaks = c("Durham et al. (2023)", "Dominguez et al. (2016)", "Ritter et al. (2017)", "Kowalewski et al. (2018)", "Toma\u{161}ov\u{FD}ch et al. (2018)")) +
  scale_shape_manual(values = c("This study" = 21, "Kowalewski et al. (2018)" = 22, "Dominguez et al. (2016)" = 23, "Ritter et al. (2017)" = 24, "Toma\u{161}ov\u{FD}ch et al. (2018)" = 25),
                     breaks = c("This study", "Dominguez et al. (2016)", "Ritter et al. (2017)", "Kowalewski et al. (2018)", "Toma\u{161}ov\u{FD}ch et al. (2018)")) +
  guides(fill = guide_legend(override.aes = list(shape = 15, size = 5, color = hcl.colors(palette = "Harmonic", n = 2)))) +
  theme(axis.title = element_text(size = 7), 
        axis.text = element_text(size = 7, colour = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7),
        legend.key.size = unit(1,"line"),
        panel.grid.major.x = element_blank(), #element_line(size = 0.25, color = "grey60"),
        panel.grid.minor.x = element_blank(), #element_line(size = 0.25, color = "grey80"),
        panel.grid.major.y = element_blank(), #element_line(size = 0.25, color = "grey70"),
        #panel.grid.minor.y = element_line(size = 0.25),
        axis.line = element_line(size = 0.25))

ggsave(here::here("Summer2023/Figures/Fig3.pdf"),
       TAscalecompare,
       width = 44*(1/6), # 7.125, #full page width GSA Pub spec
       height = 3, #3.5,
       #width = 3.5, #column width for 2-column page layout GSA Pub spec
       #height = 3,
       units = "in",
       device = grDevices::cairo_pdf)

ggsave(here::here("Summer2023/Figures/Fig3.png"),
       TAscalecompare,
       width = 44*(1/6), # 7.125, #full page width GSA Pub spec
       height = 3, #3.5,
       #width = 3.5, #column width for 2-column page layout GSA Pub spec
       #height = 3,
       units = "in")


# Geology paper figure 2----------------------------------------------------------------------------------

#Create a reverse log scale for the plot x-axis (https://stackoverflow.com/questions/63165462/how-to-reverse-axis-y-predetermined-from-scale-y-log10-ggplot-repeated)
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(paste0("reverselog-", format(base)), trans, inv, 
                    scales::log_breaks(base = base), 
                    domain = c(0, Inf))
}

hobsf14c2[!is.na(reef_ind), `:=` (max_reef_lab = max(unique(reef_ind)),
                                  pointsize = length(unique(spec_field_id)) * 0.3), by = locality]
hobsf14c2[, loc_reef_lab := paste0(loc, reef_ind)]
hobsf14c2[, hole_lab := str_replace(hole, "H", "Hole ")]

hobsf14c2[ , `:=` (hliney = vector('numeric', nrow(hobsf14c2)), nlaby = vector('numeric', nrow(hobsf14c2)), ylim_high = vector('numeric', nrow(hobsf14c2)), ylim_low = vector('numeric', nrow(hobsf14c2)))]

for(x in unique(hobsf14c2$reef)){
  hobsf14c2[!is.na(igsn) & oor == FALSE & analysis == "LP" & reef == x, hliney := min(ifelse(median_age_bp - (tav_iqr/2) < 0, 0, 
                                                                          median_age_bp - (tav_iqr/2)))]
  hobsf14c2[!is.na(igsn) & oor == FALSE & analysis == "LP" & reef == x, maxval := max(ifelse(median_age_bp - (tav_iqr/2) < 0, 
                                                                          tav_iqr, 
                                                                          median_age_bp + (tav_iqr/2)))]
  hobsf14c2[!is.na(igsn) & oor == FALSE & analysis == "LP" & reef == x, laby := hliney - ((maxval - hliney)/6.625)]
  hobsf14c2[!is.na(igsn) & oor == FALSE & analysis == "LP" & reef == x, ylim_high := maxval + (((maxval - hliney)/6.625)/15)]
  hobsf14c2[!is.na(igsn) & oor == FALSE & analysis == "LP" & reef == x, ylim_low := laby - (((maxval - hliney)/6.625)/2.5)]
}


hgridline <- distinct(hobsf14c2[!is.na(igsn) & oor == FALSE & analysis == "LP", .(loc, max_reef_lab)])
hgridline2 <- data.table(yint = hgridline[, seq(1.5, max_reef_lab - 0.5), by = loc])
setnames(hgridline2, c("yint.loc", "yint.V1"), c("loc", "yint"))

tiley <- hobsf14c2[oor == FALSE & analysis == "LP" & loc_reef_lab %in% c("BH1", "BH3", "GI-EC2", "GR1", "GR3", "HC-MC1", "JI1", "LB1", "LB3", "LC2", "LC4", "LSG1", "LSG3", "MR2", "NP2", "PC2"), .(loc, reef_ind, loc_reef_lab, max_reef_lab)]
tiley[, linewidth := 29/max_reef_lab, by = loc]
tiley <- distinct(tiley[, `:=` (reef_ind = as.numeric(reef_ind), loc_reef_lab = NULL, max_reef_lab = NULL)])
tiley[, lwind := fcase(round(linewidth, 5) == 7.25000, "3",
                       round(linewidth, 6) == 9.666667, "2",
                       round(linewidth, 5) == 14.50000, "1"), by = row.names(tiley)]

#Dummy data to try to fix the spacing of y-axis
ydum <- data.table(loc = rep(hobsf14c2[oor == FALSE & analysis == "LP" & !is.na(loc), unique(loc)], each = 10),
                   ypt = rep(seq(0.5, 5, 0.5), 11))
ydum2 <- ydum[0]
for(l in hobsf14c2[!is.na(igsn) & oor == FALSE & analysis == "LP" & !is.na(loc), unique(loc)]){
  ydum_i <- ydum[loc == l & ypt <= hobsf14c2[loc == l, unique(max_reef_lab)] + 0.5]
  ydum2 <- rbind(ydum2, ydum_i)
}

# Links I found useful for creating two separate size scales:
# 1. ggnewscale() function: https://stackoverflow.com/questions/17642190/how-to-set-multiple-legends-scales-for-the-same-aesthetic-in-ggplot2/17642257#17642257
# 2. re-purposing a different, unused scale: https://stackoverflow.com/questions/34893760/how-to-scale-the-size-of-line-and-point-separately-in-ggplot2

#get a look at the distributions of ages and time-averaging in the data
(ggplot(hobsf14c2[!is.na(igsn) & analysis == "LP" & depth_int > 0 & oor == FALSE, ], aes(x = median_age_bp, group = as.factor(depth_int), fill = as.factor(depth_int))) + geom_histogram(color = "black", alpha = 0.5) +
 ggplot(hobsf14c2[!is.na(igsn) & analysis == "LP" & depth_int > 0 & oor == FALSE, ], aes(x = median_age_bp, group = as.factor(depth_int), fill = as.factor(depth_int))) + geom_density(color = "black", alpha = 0.5)) /
(ggplot(hobsf14c2[!is.na(igsn) & analysis == "LP" & depth_int > 0 & oor == FALSE, ], aes(x = tav_iqr, group = as.factor(depth_int), fill = as.factor(depth_int))) + geom_histogram(color = "black", alpha = 0.5) +
 ggplot(hobsf14c2[!is.na(igsn) & analysis == "LP" & depth_int > 0 & oor == FALSE, ], aes(x = tav_iqr, group = as.factor(depth_int), fill = as.factor(depth_int))) + geom_density(color = "black", alpha = 0.5)) +
plot_layout(guides = 'collect') & labs(fill = "Depth")


medagebysample <- ggplot(hobsf14c2[!is.na(igsn) & analysis == "LP" & depth_int > 0 & oor == FALSE, ], aes(x = median_age_bp, y = as.character(reef_ind), shape = hole_lab, fill = as.ordered(depth_int))) +
  theme_bw(base_size = 7, base_family = "Arial") +
  geom_point(data = ydum2, aes(x = NA, y = ypt), alpha = 1, inherit.aes = FALSE) +
  geom_hline(data = tiley, aes(yintercept = reef_ind, size = lwind), color = "grey95") +
  geom_vline(xintercept = c(1, 10, 100, 1000), lwd = 0.25, color = "grey70") +
  geom_hline(data = hgridline2, aes(yintercept = yint), lwd = 0.25, color = "black") +
  scale_size_manual(values = c(14.50000, 9.666667, 7.25000), breaks = c(1:3),
                    guide = "none") +
  ggnewscale::new_scale("size") +
  geom_errorbar(aes(xmin = ifelse(median_age_bp - (tav_iqr/2) < 0, 0, median_age_bp - (tav_iqr/2)), 
                    xmax = ifelse(median_age_bp - (tav_iqr/2) < 0, tav_iqr, median_age_bp + (tav_iqr/2)), color = as.ordered(depth_int), group = igsn),
                position = position_dodge(0.5), width = 0.35, key_glyph = draw_key_rect, lwd = 0.3) +
  geom_point(aes(size = as.ordered(n), group = igsn), position = position_dodge(0.5), color = "black") +
  geom_point(aes(y = NA, x = ylim_low, group = igsn), alpha = 0) +
  geom_point(aes(y = NA, x = ylim_high, group = igsn), alpha = 0) +
  scale_x_continuous(breaks = c(1, 10, 100, 1000, 10000), trans = reverselog_trans(10)) +
  scale_y_discrete(limits = rev, breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5), labels = c("", "1", "", "2", "", "3", "", "4", "", "4", ""), expand = c(0,0)) +
  scale_shape_manual(values = c("Hole 1" = 21, "Hole 2" = 22, "Hole 3" = 24)) +
  # scale_size_manual(values = c(0.6, 0.9, 1.2, 1.5, 1.8, 2.1), breaks = c(2:7),
  scale_size_manual(values = c(0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 4, 4.5), breaks = hobsf14c2[!is.na(n), sort(unique(n))],
                    guide = guide_legend(override.aes = list(shape = 21,
                                                             color = "black"))) +
  scale_fill_discrete_qualitative(palette = "Harmonic", guide = "none") +
  scale_color_discrete_qualitative(palette = "Harmonic", labels = c("15-25", "25-35")) +
  labs(y = "Reef",
       x = "Median age (years before 2020)",
       shape = "Sample hole",
       color = "Burial depth (cm)",
       size = "Sample size") +
  theme(panel.grid.major.x = element_line(linetype = "blank"), 
        panel.grid.minor.x = element_line(linetype = "blank"), 
        panel.grid.major.y = element_line(linetype = "blank"), 
        panel.grid.minor.y = element_line(linetype = "blank"),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(1,"line"),
        panel.background = element_rect(fill = NA)) +
  facet_grid(rows = vars(loc), scales = "free_y", switch = "y")
# medagebysample

ggsave(here::here("Summer2023/Figures/Fig2.pdf"),
       medagebysample,
       #width = 7.125, #full page width GSA Pub spec
       #height = 5,
       width = 29 * (1/6), #3.5, #column width for 2-column page layout GSA Pub spec
       height = 53 * (1/6), #9.5, #full page depth GSA Pub spec
       units = "in",
       device = grDevices::cairo_pdf)

ggsave(here::here("Summer2023/Figures/Fig2.png"),
       medagebysample,
       #width = 7.125, #full page width GSA Pub spec
       #height = 5,
       width = 29 * (1/6), #3.5, #column width for 2-column page layout GSA Pub spec
       height = 53 * (1/6), #9.5, #full page depth GSA Pub spec
       units = "in")


#Space v. depth models from John-----------------------------------------
#I've put the code here to have all of the analysis steps in one file; John's original markdown script is "GeoChron.Rmd" and 
#the post-processing script I wrote to create the plots is called "ProcessModelResultsPost-Script.R".


#median estimated sample ages and corresponding time-averaging estimates from Dominguez et al. (2016)
Detal <- data.table(median_age_bp = c(123, 283, 146, 148, 160, 2043), 
                    tav_iqr = c(130, 2766, 84, 96, 89, 2204),
                    Locality = c(1, 1, 1, 1, 1, 1),
                    Reef = c(1, 2, 3, 4, 5, 6),
                    Reef2Locality = c(1, 1, 1, 1, 1, 1),
                    nL = 1,
                    nR = 6)

reeforder <- setorder(distinct(hobsf14c2[!is.na(igsn) & analysis == "LP" & oor == FALSE, .(reef), by = loc]), loc)
hobsf14c2[!is.na(igsn) & analysis == "LP" & oor == FALSE, loc_ind := as.integer(loc)]
hobsf14c2[!is.na(igsn) & analysis == "LP" & oor == FALSE, reef_ind2 := as.integer(factor(reef, levels = reeforder$reef))]
hobsf14c2[!is.na(igsn) & analysis == "LP" & oor == FALSE, reef2Locality := loc_ind, by = c("loc_ind", "reef_ind2")]
hobsf14c2[!is.na(igsn) & analysis == "LP" & oor == FALSE, holeID := paste0(reef, hole)]
hobsf14c2[!is.na(igsn) & analysis == "LP" & oor == FALSE, sampsperhole := length(unique(igsn)), by = holeID]

sampsum_r2l <- setorder(distinct(hobsf14c2[!is.na(igsn) & analysis == "LP" & oor == FALSE, .(reef = reef_num, 
                                                                                             locality = loc, 
                                                                                             reefID = reef,
                                                                                             Locality = loc_ind,
                                                                                             reef3 = reef_ind2,
                                                                                             reef2Locality)]), Locality, reef3)

sampsum_15to25cm <- setorder(distinct(hobsf14c2[!is.na(igsn) & analysis == "LP" & depth_int == 1 & oor == FALSE, .(sample = igsn,
                                                                                                                   Sample_median_age = tav_median,
                                                                                                                   median_age_bp,
                                                                                                                   tav_iqr,
                                                                                                                   reef = reef_num,
                                                                                                                   reef_lab = reef_ind,
                                                                                                                   locality = loc,
                                                                                                                   reefID = reef,
                                                                                                                   holeID,
                                                                                                                   Locality = loc_ind,
                                                                                                                   reef3 = reef_ind2,
                                                                                                                   reef2Locality)]), Locality, reef3)
reeford <- setorder(sampsum_15to25cm[, .(reefID = unique(reefID)), by = reef3], reef3)
holeord <- setorder(sampsum_15to25cm[, .(holeID = unique(holeID)), by = reef3], reef3, holeID)
sampsum_15to25cm[, `:=` (reefID = factor(reefID, levels = reeford$reefID),
                         holeID = factor(holeID, levels = holeord$holeID))]

sampsum_25to35cm <- setorder(distinct(hobsf14c2[!is.na(igsn) & analysis == "LP" & depth_int == 2 & oor == FALSE, .(sample = igsn,
                                                                                                                   Sample_median_age = tav_median,
                                                                                                                   median_age_bp,
                                                                                                                   tav_iqr,
                                                                                                                   reef = reef_num,
                                                                                                                   reef_lab = reef_ind,
                                                                                                                   locality = loc,
                                                                                                                   reefID = reef,
                                                                                                                   holeID,
                                                                                                                   Locality = loc_ind,
                                                                                                                   reef3 = reef_ind2,
                                                                                                                   reef2Locality)]), Locality, reef3)
sampsum_25to35cm[, `:=` (reefID = factor(reefID, levels = reeford$reefID),
                         holeID = factor(holeID, levels = holeord$holeID))]

depthdiffs <- setorder(distinct(hobsf14c2[!is.na(igsn) & analysis == "LP" & oor == FALSE & sampsperhole == 2, .(sample = igsn,
                                                                                                                Sample_median_age = tav_median,
                                                                                                                median_age_bp,
                                                                                                                tav_iqr,
                                                                                                                depth = depth_int,
                                                                                                                reef = reef_num,
                                                                                                                reef_lab = reef_ind,
                                                                                                                locality = loc,
                                                                                                                reefID = reef,
                                                                                                                holeID,
                                                                                                                Locality = loc_ind,
                                                                                                                reef3 = reef_ind2,
                                                                                                                reef2Locality)]), Locality, reef3)
depthdiffs[, `:=` (reefID = factor(reefID, levels = reeford$reefID),
                   holeID = factor(holeID, levels = holeord$holeID),
                   depthdiff_mabp = median_age_bp[depth == 2] -  median_age_bp[depth == 1], 
                   depthdiff_tav = tav_iqr[depth == 2] - tav_iqr[depth == 1]), by = holeID]

depthdiffs2 <- distinct(depthdiffs[, .(reef, locality, reefID, holeID, Locality, reef3, reef2Locality, depthdiff_mabp, depthdiff_tav)])

locs <- distinct(sampsum_15to25cm[, .(locality, Locality)])
reefs <- distinct(sampsum_15to25cm[, .(locality, reefID, reef, reef_lab, reef3)])


#for median age over space (15-25cm burial depth); scaled version produces cleaner MCMC diagnostics, but is difficult to compare with Dominguez et al. (2016) results
HOBS_medagexspace_15to25cm <- list(N = nrow(sampsum_15to25cm),
                                   Locality = sampsum_15to25cm$Locality,
                                   Reef = sampsum_15to25cm$reef3,
                                   Reef2Locality = sampsum_r2l$reef2Locality,
                                   x = sampsum_15to25cm$median_age_bp,
                                   nL = length(unique(sampsum_15to25cm$Locality)),
                                   nR = length(unique(sampsum_15to25cm$reef3)))

HOBS_medagexspace_15to25cm_scaled <- list(N = nrow(sampsum_15to25cm),
                                          Locality = sampsum_15to25cm$Locality,
                                          Reef = sampsum_15to25cm$reef3,
                                          Reef2Locality = sampsum_r2l$reef2Locality,
                                          x = as.numeric(scale(sampsum_15to25cm$median_age_bp)),
                                          nL = length(unique(sampsum_15to25cm$Locality)),
                                          nR = length(unique(sampsum_15to25cm$reef3)))

#for TAV over space (15-25cm burial depth); scaled version produces cleaner MCMC diagnostics, but is difficult to compare with Dominguez et al. (2016) results
HOBS_tavxspace_15to25cm <- list(N = nrow(sampsum_15to25cm),
                                Locality = sampsum_15to25cm$Locality,
                                Reef = sampsum_15to25cm$reef3,
                                Reef2Locality = sampsum_r2l$reef2Locality,
                                x = sampsum_15to25cm$tav_iqr,
                                nL = length(unique(sampsum_15to25cm$Locality)),
                                nR = length(unique(sampsum_15to25cm$reef3)))

HOBS_tavxspace_15to25cm_scaled <- list(N = nrow(sampsum_15to25cm),
                                       Locality = sampsum_15to25cm$Locality,
                                       Reef = sampsum_15to25cm$reef3,
                                       Reef2Locality = sampsum_r2l$reef2Locality,
                                       x = as.numeric(scale(sampsum_15to25cm$tav_iqr)),
                                       nL = length(unique(sampsum_15to25cm$Locality)),
                                       nR = length(unique(sampsum_15to25cm$reef3)))

#for median age over space (25-35cm burial depth); scaled version produces cleaner MCMC diagnostics, but is difficult to compare with Dominguez et al. (2016) results
HOBS_medagexspace_25to35cm <- list(N = nrow(sampsum_25to35cm),
                                   Locality = sampsum_25to35cm$Locality,
                                   Reef = sampsum_25to35cm$reef3,
                                   Reef2Locality = sampsum_r2l$reef2Locality,
                                   x = sampsum_25to35cm$median_age_bp,
                                   nL = length(unique(sampsum_25to35cm$Locality)),
                                   nR = length(unique(sampsum_25to35cm$reef3)))

HOBS_medagexspace_25to35cm_scaled <- list(N = nrow(sampsum_25to35cm),
                                          Locality = sampsum_25to35cm$Locality,
                                          Reef = sampsum_25to35cm$reef3,
                                          Reef2Locality = sampsum_r2l$reef2Locality,
                                          x = as.numeric(scale(sampsum_25to35cm$median_age_bp)),
                                          nL = length(unique(sampsum_25to35cm$Locality)),
                                          nR = length(unique(sampsum_25to35cm$reef3)))

#for TAV over space (25-35cm burial depth)
HOBS_tavxspace_25to35cm <- list(N = nrow(sampsum_25to35cm),
                                Locality = sampsum_25to35cm$Locality,
                                Reef = sampsum_25to35cm$reef3,
                                Reef2Locality = sampsum_r2l$reef2Locality,
                                x = sampsum_25to35cm$tav_iqr,
                                nL = length(unique(sampsum_25to35cm$Locality)),
                                nR = length(unique(sampsum_25to35cm$reef3)))

#for median age depth differences over space
HOBS_medageddxspace <- list(N = nrow(depthdiffs2),
                            Locality = depthdiffs2$Locality,
                            Reef = depthdiffs2$reef3,
                            Reef2Locality = sampsum_r2l$reef2Locality,
                            x = depthdiffs2$depthdiff_mabp,
                            nL = length(unique(depthdiffs2$Locality)),
                            nR = length(unique(depthdiffs2$reef3)))

#for TAV depth differences over space
HOBS_tavddxspace <- list(N = nrow(depthdiffs2),
                         Locality = depthdiffs2$Locality,
                         Reef = depthdiffs2$reef3,
                         Reef2Locality = sampsum_r2l$reef2Locality,
                         x = depthdiffs2$depthdiff_tav,
                         nL = length(unique(depthdiffs2$Locality)),
                         nR = length(unique(depthdiffs2$reef3)))



#For HOBS data
model_HOBS <- cmdstan_model("GeoChron.stan")

#For Dominguez et al. (2016) data
model_Detal <- cmdstan_model("GeoChron_Detal.stan")


datalists <- c(as.name("HOBS_medagexspace_15to25cm"), 
               as.name("HOBS_tavxspace_15to25cm"),
               as.name("HOBS_medagexspace_25to35cm"),
               as.name("HOBS_tavxspace_25to35cm"),
               as.name("HOBS_medageddxspace"),
               as.name("HOBS_tavddxspace"))


for(Data in datalists){
  if(str_detect(as.character(Data), "HOBS")){
    model <- cmdstan_model("GeoChron.stan")
  } else{
    model <- cmdstan_model("GeoChron_Detal.stan")
  }
  
  fit = model$sample(data = eval(Data),
                     seed = 1234,
                     chains = 8,
                     parallel_chains = 8,
                     refresh = 1000,
                     iter_warmup = 5000, #3000,
                     iter_sampling = 15000,
                     thin =  5,
                     adapt_delta = 0.99,
                     max_treedepth = 15)
  
  
  fit$save_output_files("StanOutput")
  fitdiag <- fit$cmdstan_diagnose()
  fwrite(fitdiag, file = here::here(paste0("Summer2023/StanOutput/", Data, "_diagnostics.txt")))
  dim(fit$draws()) # accessing draws seems to fill the structure with draws
  saveRDS(fit, file = here::here(paste0("Summer2023/StanOutput/", Data, ".rds")))
  
}


# plot1 <- bayesplot::mcmc_intervals(fit$draws("mu_locality")) 
# plot1 <- plot1 +
#   geom_point(data = Data, aes(x = x, y = Locality - 12), color = "red") +
#   geom_text(data = locs, aes(x = min(plot1$data$ll) - 3, y = loclab, label = locality), inherit.aes = FALSE) +
#   xlim(min(plot1$data$ll) - 5, max(plot1$data$hh))
# 
# saveRDS(plot1, here::here(paste0("StanOutput/", Data, "locintplot.rds")))
# ggsave(filename = here::here(paste0("StanOutput/", Data, "locintplot.pdf")),
#        plot = plot1,
#        width = 6,
#        height = 6,
#        units = "in",
#        dpi = 400)
# 
# plot2 <- bayesplot::mcmc_intervals(fit$draws("mu_reef"))
# plot2 <- plot2 +
#   geom_text(data = reefs, aes(x = min(plot2$data$ll) - 10, y = reeflab, label = reefID), inherit.aes = FALSE) +
#   xlim(min(plot2$data$ll) - 12, max(plot2$data$hh))
# 
# saveRDS(plot2, here::here(paste0("StanOutput/", Data, "reefintplot.rds")))
# ggsave(filename = here::here(paste0("StanOutput/", Data, "reefintplot.pdf")),
#        plot = plot2,
#        width = 6,
#        height = 6,
#        units = "in",
#        dpi = 400)



#Space v. depth model results post-processing------------------------------------------------------

#Note: data files are created in the "GeoChron.Rmd" script.

files <- list.files(here::here("Summer2023/StanOutput/"))
files <- files[str_detect(files, "StanOutput", negate = TRUE)]
files <- files[str_detect(files, ".rds")]
files <- files[str_detect(files, "plot", negate = TRUE)]
files <- files[str_detect(files, ".txt", negate = TRUE)]
files <- files[str_detect(files, "ModelResults", negate = TRUE)]
# files <- files[-c(6, 8)]

metric_choice <- "median" # or "mean"

#Scooch axis breaks
sep_ylim_min <- -1.25
sep_ylim_max <- 1.25
sep_xlim_min <- 5
sep_xlim_max <- 14


for(file in files){
  fit <- readRDS(here::here(paste0("Summer2023/StanOutput/", file)))
  #fitdiag <- fit$cmdstan_diagnose()
  #fwrite(fitdiag, file = here::here(paste0("StanOutput/", str_sub(file, 1L, -5L), "_diagnostics.txt")))
  
  if(str_detect(file, "15to25cm")){
    Data <- copy(sampsum_15to25cm)
  } else if(str_detect(file, "25to35cm")){
    Data <- copy(sampsum_25to35cm)
  } else{
    Data <- copy(depthdiffs2)
  }
  
  param <- ifelse(str_detect(file, "tavdd"), as.name("depthdiff_tav"), 
                  ifelse(str_detect(file, "medagedd"), as.name("depthdiff_mabp"), 
                         ifelse(str_detect(file, "medage"), as.name("median_age_bp"), as.name("tav_iqr"))))
  
  Data[, mean_param_loc := mean(eval(param)), by = Locality]
  Data[, med_param_loc := median(eval(param)), by = Locality]
  Data[, mean_param_loc_scaled := as.numeric(scale(Data$mean_param_loc))]
  Data[, mean_param_reef := mean(eval(param)), by = reef3]
  Data[, med_param_reef := median(eval(param)), by = reef3]
  Data[, mean_param_reef_scaled := as.numeric(scale(Data$mean_param_reef))]
  if(as.character(param) %in% c("median_age_bp", "tav_iqr")){
    Detal[, mean_param := mean(eval(param))]
    Detal[, med_param := median(eval(param))]
  }
  
  if(str_detect(file, "dd")){
    plot1 <- bayesplot::mcmc_intervals(fit$draws("mu_locality"), prob_outer = 0.95, point_size = 2, inner_size = 1.25)
    
    plot1 <- plot1 +
      geom_vline(xintercept = 0, color = "grey90", lwd = 0.5) +
      geom_point(data = plot1$data, aes(x = m, y = parameter, shape = "Model", color = "Model", fill = "Model", size = "Model"), stroke = 0.25) +
      {if(metric_choice == "median"){
        geom_point(data = Data, aes(x = med_param_loc, y = 12 - Locality, shape = "Data", color = "Data", fill = "Data", size = "Data"), stroke = 1)
      } else{
        geom_point(data = Data, aes(x = mean_param_loc, y = 12 - Locality, shape = "Data", color = "Data", fill = "Data", size = "Data"), stroke = 1)  
      }} +
      {if(metric_choice == "median"){
        xlim(min(c(plot1$data$ll, Data$med_param_loc)), max(plot1$data$hh, Data$med_param_loc))
      } else{
        xlim(min(c(plot1$data$ll, Data$mean_param_loc)), max(plot1$data$hh, Data$mean_param_loc))
      }} +
      guides(size = "none") +
      scale_shape_manual(values = c("Data" = 3, "Model" = 21)) +
      scale_color_manual(values = c("Data" = "red", "Model" = "#011f4b")) +
      scale_fill_manual(values = c("Data" = "red", "Model" = "#d1e1ec")) +
      scale_size_manual(values = c("Data" = 1, "Model" = 2)) +
      #scale_y_discrete(breaks = plot1$data$parameter, labels = str_pad(locs$locality, 15, "left"), limits = rev) + #padding no longer necessary after adding the sigma plots
      scale_y_discrete(breaks = plot1$data$parameter, labels = locs$locality, limits = rev) +
      labs(x = ifelse(param == "median_age_bp", expression("Median age"["25-35cm"] * " - Median age"["15-25cm"] * " (years before collection year)"), expression("TAV_IQR"["25-35cm"] * " - TAV_IQR"["15-25cm"] * " (yr)")), y = "Locality", title = str_sub(file, 1L, -5L), shape = NULL, colour = NULL, fill = NULL, size = NULL) +
      theme(text = element_text(family = "sans", size = 9),
            legend.position = "right",
            legend.text = element_text(size = 9))
    
    saveRDS(plot1, here::here(paste0("Summer2023/StanOutput/", str_sub(file, 1L, -5L), "_locintplot_", metric_choice, ".rds")))
    ggsave(filename = here::here(paste0("Summer2023/StanOutput/", str_sub(file, 1L, -5L), "_locintplot_", metric_choice, ".pdf")),
           plot = plot1,
           width = 7,
           height = 4.5,
           units = "in",
           device = cairo_pdf,
           dpi = 300)
    
    ggsave(filename = here::here(paste0("Summer2023/StanOutput/", str_sub(file, 1L, -5L), "_locintplot_", metric_choice, ".png")),
           plot = plot1,
           width = 7,
           height = 4.5,
           units = "in")
    
    plot2 <- bayesplot::mcmc_intervals(fit$draws("mu_reef"), prob_outer = 0.95, point_size = 2, inner_size = 1.25)
    plot2 <- plot2 +
      geom_vline(xintercept = 0, color = "grey90", lwd = 0.5) +
      geom_point(data = plot2$data, aes(x = m, y = parameter, shape = "Model", color = "Model", fill = "Model", size = "Model"), stroke = 0.25) +
      {if(metric_choice == "median"){
        geom_point(data = Data, aes(x = med_param_reef, y = 32 - reef3, shape = "Data", color = "Data", fill = "Data", size = "Data"), stroke = 1)
      } else{
        geom_point(data = Data, aes(x = mean_param_reef, y = 32 - reef3, shape = "Data", color = "Data", fill = "Data", size = "Data"), stroke = 1)  
      }
      } +
      {if(metric_choice == "median"){
        xlim(min(c(plot1$data$ll, Data$med_param_reef)), max(plot1$data$hh, Data$med_param_reef))
      } else{
        xlim(min(c(plot1$data$ll, Data$mean_param_reef)), max(plot1$data$hh, Data$mean_param_reef))
      }
      } +
      guides(size = "none") +
      scale_shape_manual(values = c("Data" = 3, "Model" = 21)) +
      scale_color_manual(values = c("Data" = "red", "Model" = "#011f4b")) +
      scale_fill_manual(values = c("Data" = "red", "Model" = "#d1e1ec")) +
      scale_size_manual(values = c("Data" = 1, "Model" = 2)) +
      scale_y_discrete(breaks = plot2$data$parameter, labels = paste0(reefs$locality, " R", reefs$reef_lab), limits = rev) +
      scale_x_continuous(n.breaks = 7) +
      labs(x = ifelse(param == "median_age_bp", expression("Median age"["25-35cm"] * " - Median age"["15-25cm"] * " (years before collection year)"), expression("TAV_IQR"["25-35cm"] * " - TAV_IQR"["15-25cm"] * " (yr)")), y = "Reef", title = str_sub(file, 1L, -5L), shape = NULL, colour = NULL, fill = NULL, size = NULL) +
      theme(text = element_text(family = "sans", size = 9),
            legend.position = "right",
            legend.text = element_text(size = 9))
    
    saveRDS(plot2, here::here(paste0("Summer2023/StanOutput/", str_sub(file, 1L, -5L), "_reefintplot_", metric_choice, ".rds")))
    ggsave(filename = here::here(paste0("Summer2023/StanOutput/", str_sub(file, 1L, -5L), "_reefintplot_", metric_choice, ".pdf")),
           plot = plot2,
           width = 7,
           height = 4.5,
           units = "in",
           device = cairo_pdf,
           dpi = 300)
    
    ggsave(filename = here::here(paste0("Summer2023/StanOutput/", str_sub(file, 1L, -5L), "_reefintplot_", metric_choice, ".png")),
           plot = plot2,
           width = 7,
           height = 4.5,
           units = "in")
    
    #sigma plot
    Data[, `:=` (sd_medage_loc = sd(depthdiff_mabp),
                 sd_tav_loc = sd(depthdiff_tav)), by = locality]
    Data[, `:=` (sd_medage_reef = sd(depthdiff_mabp),
                 sd_tav_reef = sd(depthdiff_tav)), by = reefID]
    Data[, `:=` (sd_medage_hole = sd(depthdiff_mabp),
                 sd_tav_hole = sd(depthdiff_tav)), by = holeID]
    Data[, `:=` (sd_medage_hole2 = sd(depthdiff_mabp),
                 sd_tav_hole2 = sd(depthdiff_tav))]
    Data[, `:=` (sd_medage_reef2 = sd(Data[, mean(depthdiff_mabp), by = reefID]$V1),
                 sd_tav_reef2 = sd(Data[, mean(depthdiff_tav), by = reefID]$V1))]
    Data[, `:=` (sd_medage_loc2 = sd(Data[, mean(depthdiff_mabp), by = locality]$V1),
                 sd_tav_loc2 = sd(Data[, mean(depthdiff_tav), by = locality]$V1))]
    Data_sd <- melt(Data[, .(locality, reefID, holeID, sample, 
                             sd_medage_loc, sd_tav_loc, sd_medage_reef, sd_tav_reef, sd_medage_hole, sd_tav_hole, 
                             sd_medage_loc2, sd_tav_loc2, sd_medage_reef2, sd_tav_reef2, sd_medage_hole2, sd_tav_hole2)], 
                    id.vars = c("locality", "reefID", "holeID"), measure.vars = c("sd_medage_loc", "sd_tav_loc", "sd_medage_reef", "sd_tav_reef", "sd_medage_hole", "sd_tav_hole",
                                                                                  "sd_medage_loc2", "sd_tav_loc2", "sd_medage_reef2", "sd_tav_reef2", "sd_medage_hole2", "sd_tav_hole2"))
    # Data_sd[, `:=` (yname = fcase(str_detect(variable, "loc"), "sigma_locality",
    #                               str_detect(variable, "reef"), "sigma_reef",
    #                               str_detect(variable, "hole"), "sigma_hole"),
    #                 metric = ifelse(str_detect(variable, "medage"), "medage", "tav_iqr")), by = c("locality", "reefID", "holeID", "variable")]
    setDT(Data_sd)
    Data_sd[, `:=` (yname = fcase(str_detect(variable, "loc"), "sigma_locality",
                                  str_detect(variable, "reef"), "sigma_reef",
                                  str_detect(variable, "hole"), "sigma_hole"),
                    metric = ifelse(str_detect(variable, "medage"), "medage", "tav_iqr")), by = list(locality, reefID, holeID, variable)]
    
    
    plot3 <- mcmc_intervals(fit$draws(c("sigma_locality", "sigma_reef", "sigma_hole")), prob_outer = 0.95, point_size = 2, inner_size = 1.25)
    plot3 <- plot3 +
      geom_point(data = plot3$data, aes(x = m, y = parameter, shape = "Model", color = "Model", fill = "Model", size = "Model"), stroke = 0.25) +
      # {if(metric_choice == "median"){
      #   geom_point(data = distinct(Data_sd[str_detect(variable, "2") & metric == "medage", .(value, yname)]), aes(x = value, y = yname, shape = "Data", color = "Data", fill = "Data", size = "Data"), stroke = 1)
      # } else{
      #   geom_point(data = distinct(Data_sd[str_detect(variable, "2") & metric == "cpe", .(value, yname)]), aes(x = value, y = yname, shape = "Data", color = "Data", fill = "Data", size = "Data"), stroke = 1)
      # }} +
      guides(size = "none") +
      scale_x_continuous(trans = "log10") +
      scale_shape_manual(values = c("Data" = 3, "Model" = 21)) +
      scale_color_manual(values = c("Data" = "red", "Model" = "#011f4b")) +
      # scale_linetype_manual(values = c("Data" = NA, "Model" = NA)) +
      scale_fill_manual(values = c("Data" = "red", "Model" = "#d1e1ec")) +
      scale_size_manual(values = c("Data" = 1, "Model" = 2)) +
      # scale_y_discrete(breaks = plot2$data$parameter[order(reefs2$parameter)], labels = paste0(reefs2$locality, " R", reefs2$reef)) +
      scale_y_discrete(breaks = plot3$data$parameter, labels = c(expression(sigma[c]), expression(sigma[r]), expression(sigma[h])), limits = rev) +
      #labs(x = param, y = "Reef", title = str_sub(file, 1L, -5L), shape = "Means", colour = "Means", fill = "Means", size = "Means", lty = "Means") +
      labs(x = ifelse(param == "median_age_bp", expression("SD Median age"["25-35cm"] * " - Median age"["15-25cm"] * " (years before collection year)"), expression("SD TAV_IQR"["25-35cm"] * " - TAV_IQR"["15-25cm"] * " (yr)")), y = "Sigma", title = str_sub(file, 1L, -5L), shape = NULL, colour = NULL, fill = NULL, size = NULL, lty = NULL) +
      theme(text = element_text(family = "sans", size = 9),
            axis.text.y = element_text(size = 11),
            axis.title.y = element_text(size = 9),
            legend.text = element_text(family = "sans", size = 9),
            legend.position = "right")
    #coord_cartesian(ylim = c(1, 4.5))
    #plot3
    
    saveRDS(plot3, here::here(paste0("Summer2023/StanOutput/", str_sub(file, 1L, -5L), "_sigmaintplot_", metric_choice, ".rds")))
    ggsave(filename = here::here(paste0("Summer2023/StanOutput/", str_sub(file, 1L, -5L), "_sigmaintplot_", metric_choice, ".pdf")),
           plot = plot3,
           width = 7,
           height = 4.5,
           units = "in",
           device = cairo_pdf,
           dpi = 300)
    
    ggsave(filename = here::here(paste0("Summer2023/StanOutput/", str_sub(file, 1L, -5L), "_sigmaintplot_", metric_choice, ".png")),
           plot = plot3,
           width = 7,
           height = 4.5,
           units = "in")
    
  } else{
    plot1 <- bayesplot::mcmc_intervals(fit$draws("mu_locality"), prob_outer = 0.95, point_size = 2, inner_size = 1.25)
    
    plot1 <- plot1 +
      geom_vline(xintercept = 0, color = "grey90", lwd = 0.5) +
      # geom_point(data = Data, aes(x = mean_param_loc, y = 12 - Locality, shape = "Data", color = "Data", fill = "Data", size = "Data"), stroke = 1.5) +
      geom_point(data = plot1$data, aes(x = m, y = parameter, shape = "Model", color = "Model", fill = "Model", size = "Model"), stroke = 0.25) +
      {if(metric_choice == "median"){
        geom_vline(data = Detal, aes(xintercept = med_param, lty = "Dominguez et al. (2016)"), lwd = 0.25) #, color = "Dominguez et al. (2016)", size = "Dominguez et al. (2016)"
      } else{
        geom_vline(data = Detal, aes(xintercept = mean_param, lty = "Dominguez et al. (2016)"), lwd = 0.25) #, color = "Dominguez et al. (2016)", size = "Dominguez et al. (2016)"
      }} +
      {if(metric_choice == "median"){
        geom_point(data = Data, aes(x = med_param_loc, y = 12 - Locality, shape = "Data", color = "Data", fill = "Data", size = "Data"), stroke = 1)
      } else{
        geom_point(data = Data, aes(x = mean_param_loc, y = 12 - Locality, shape = "Data", color = "Data", fill = "Data", size = "Data"), stroke = 1)  
      }} +
      {if(metric_choice == "median"){
        xlim(c(min(plot1$data$ll, Data$med_param_loc, Detal$med_param), max(plot1$data$hh, Data$med_param_loc, Detal$med_param)))
      } else{
        xlim(c(min(plot1$data$ll, Data$mean_param_loc, Detal$mean_param), max(plot1$data$hh, Data$mean_param_loc, Detal$mean_param)))
      }} +
      guides(size = "none") +
      #xlim(min(c(plot1$data$ll, Data$mean_param_loc)), max(c(plot1$data$hh, Data$mean_param_loc))) +
      #xlim(min(c(plot1$data$ll, Data$mean_param_loc_scaled)), max(c(plot1$data$hh, Data$mean_param_loc_scaled))) +
      scale_shape_manual(values = c("Data" = 3, "Model" = 21)) + #, "Dominguez et al. (2016)" = NA
      scale_color_manual(values = c("Data" = "red", "Model" = "#011f4b")) + #, "Dominguez et al. (2016)" = NA
      scale_linetype_manual(values = c("Dominguez et al. (2016)" = "dashed")) + #"Data" = NA, "Model" = NA, 
      scale_fill_manual(values = c("Data" = "red", "Model" = "#d1e1ec")) + #, "Dominguez et al. (2016)" = NA
      scale_size_manual(values = c("Data" = 1, "Model" = 2)) + #, "Dominguez et al. (2016)" = 1
      #scale_y_discrete(breaks = plot1$data$parameter, labels = locs$locality[order(locs$Locality, decreasing = TRUE)]) +
      #scale_y_discrete(breaks = plot1$data$parameter, labels = str_pad(locs$locality, 15, "left"), limits = rev) +
      scale_y_discrete(breaks = plot1$data$parameter, labels = locs$locality, limits = rev) +
      labs(x = ifelse(param == "median_age_bp", "Median age (years before collection year)", "Total age variation"), y = "Locality", title = str_sub(file, 1L, -5L), shape = NULL, colour = NULL, fill = NULL, size = NULL, lty = NULL) +
      theme(text = element_text(family = "sans", size = 9),
            legend.position = "right")
    
    if(ifelse(metric_choice == "median", max(c(plot1$data$hh, Data$med_param_loc, Detal$med_param)) == max(Detal$med_param) & max(Detal$med_param) - max(unique(sort(c(plot1$data$hh, Data$med_param_loc, Detal$med_param)))[length(unique(c(plot1$data$hh, Data$med_param_loc, Detal$med_param))) - 1]) > 150, 
              max(c(plot1$data$hh, Data$mean_param_loc, Detal$mean_param)) == max(Detal$mean_param) & max(Detal$mean_param) - max(unique(sort(c(plot1$data$hh, Data$mean_param_loc, Detal$mean_param)))[length(unique(c(plot1$data$hh, Data$mean_param_loc, Detal$mean_param))) - 1]) > 150)){
      
      plot1z <- plot1 + 
        {if(metric_choice == "median"){
          scale_x_continuous(breaks = seq(0, #plyr::round_any(min(c(plot1$data$ll, Data$mean_param_loc)), 50, floor), 
                                          plyr::round_any(max(c(plot1$data$ll, Data$med_param_loc)), 50, ceiling), by = 50))
        } else{
          scale_x_continuous(breaks = seq(0, #plyr::round_any(min(c(plot1$data$ll, Data$mean_param_loc)), 50, floor), 
                                          plyr::round_any(max(c(plot1$data$ll, Data$mean_param_loc)), 50, ceiling), by = 50))
        }} +
        {if(metric_choice == "median"){
          coord_cartesian(xlim = c(-10, max(c(plot1$data$hh, Data$med_param_loc)) + 50))
        } else{
          coord_cartesian(xlim = c(-10, max(c(plot1$data$hh, Data$mean_param_loc)) + 50))
        }} +
        theme_classic() +
        theme(legend.position = "none",
              text = element_text(size = 9),
              plot.title = element_blank(),
              axis.title = element_blank(),
              legend.text = element_text(size = 9))
      
      plot1o <- plot1 + 
        # scale_y_discrete(breaks = plot1$data$parameter, labels = rep("/\u005c", 31)) +
        # scale_y_discrete(breaks = plot1$data$parameter, labels = rep("\u2591 ", 31)) +
        {if(metric_choice == "median"){
          scale_x_continuous(breaks = seq(plyr::round_any(max(plot1$layers[[9]]$data$med_param), 100, floor), 
                                          plyr::round_any(max(plot1$layers[[9]]$data$med_param), 100, ceiling), by = 50))
        } else{
          scale_x_continuous(breaks = seq(plyr::round_any(max(plot1$layers[[9]]$data$mean_param), 100, floor), 
                                          plyr::round_any(max(plot1$layers[[9]]$data$mean_param), 100, ceiling), by = 50))
        }} +
        {if(metric_choice == "median"){
          coord_cartesian(xlim = c(plyr::round_any(max(plot1$layers[[9]]$data$med_param), 100, floor), 
                                   plyr::round_any(max(plot1$layers[[9]]$data$med_param), 100, ceiling)))
        } else{
          coord_cartesian(xlim = c(plyr::round_any(max(plot1$layers[[9]]$data$mean_param), 100, floor), 
                                   plyr::round_any(max(plot1$layers[[9]]$data$mean_param), 100, ceiling)))
        }} +
        theme_classic() +
        theme(plot.title = element_blank(),
              plot.background = element_blank(),
              text = element_text(size = 9),
              # axis.text.y = element_text(angle = -90, color = "grey50", size = 17, vjust = -0.1, hjust = 0.5),
              axis.text.y = element_blank(),
              axis.title = element_blank(),
              axis.ticks.y = element_blank(),
              #axis.ticks.length.y = unit(0.25, "in"),
              axis.line.y = element_blank(),
              legend.margin = margin(t = 0, unit='cm'),
              legend.text = element_text(margin = margin(l = 0), hjust = 0, size = 9),
              legend.spacing.x = unit(0, "in"),
              #legend.box.spacing = margin(0),
              legend.justification = "left")
      
      # plot1c <- plot1z + plot1o + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a", 
      #                                                                               tag_prefix = "(", 
      #                                                                               tag_suffix = ")", 
      #                                                                               title = str_sub(file, 1L, -5L))
      
      ylabel <- ggplot(data.frame(l = plot1$labels$y, x = 1, y = 1)) +
        geom_text(aes(x, y, label = l), size = 2.75, angle = 90) + 
        theme_void() +
        coord_cartesian(clip = "off")
      
      xlabel <- ggplot(data.frame(l = plot1$labels$x, x = 1, y = 1)) +
        geom_text(aes(x, y, label = l), size = 2.75, angle = 0) + 
        theme_void() +
        coord_cartesian(clip = "off")
      
      separator <- ggplot() +
        #geom_tile(data = data.frame(x_min = 0, x_max = 90, y_min = min(plot1$data$parameter), y_max = max(plot1$data$parameter)), aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max), fill = "grey50") + 
        #geom_tile_pattern(aes(x = 50, y = plot1$data$parameter), fill = NA, pattern = "stripe", pattern_density = 0.5, pattern_spacing = 1, pattern_colour = "grey50", pattern_fill = "grey50") + 
        #geom_hline(yintercept = 0.45, color = "black", lwd = 0.5) +
        geom_polygon(data = data.frame(x = c(7, 8, 10, 9), y = c(0, 0, 0.9, 0.9)), aes(x = x, y = y), fill = "white", color = "white") +
        geom_segment(aes(x = c(7, 8), y = c(0, 0), xend = c(9, 10), yend = c(0.9, 0.9)), color = "black", lwd = 0.5) +
        ylim(c(sep_ylim_min, sep_ylim_max)) +
        xlim(c(sep_xlim_min, sep_xlim_max)) +
        theme_void() +
        theme(plot.background = element_blank()) +
        coord_cartesian(clip = "off")
      
      layout <- c(
        area(t = 1, l = 1, b = 30, r = 5),
        area(t = 1, l = 6, b = 30, r = 50),
        area(t = 29, l = 47, b = 31, r = 50),
        area(t = 1, l = 50, b = 30, r = 59),
        area(t = 31, l = 6, b = 32, r = 59)
      )
      plot1 <- ylabel + plot1z + separator + plot1o + xlabel + plot_layout(design = layout) + plot_annotation(title = str_sub(file, 1L, -5L), theme = theme(plot.title = element_text(hjust = 0.12),
                                                                                                                                                            text = element_text(size = 7))) #,
      #legend.margin = margin(t = 0, unit='cm'),
      #legend.text = element_text(margin = margin(l = 0), hjust = 0),
      #legend.spacing.x = unit(0, "in"),
      #legend.box.spacing = margin(0),
      #legend.justification = "left"))
      #plot1
    }
    
    # layout_b <- c(
    #   area(t = 1, l = 1, b = 30, r = 5),
    #   area(t = 1, l = 6, b = 30, r = 50),
    #   area(t = 29, l = 47, b = 31, r = 50),
    #   area(t = 1, l = 50, b = 30, r = 59)
    # )
    # plot1b <- ylabel + plot1z + separator + plot1o + plot_layout(design = layout_b) + plot_annotation(title = str_sub(file, 1L, -5L), theme = theme(plot.title = element_text(hjust = 0.12),
    # text = element_text(size = 7)))
    
    saveRDS(plot1, here::here(paste0("Summer2023/StanOutput/", str_sub(file, 1L, -5L), "_locintplot_", metric_choice, ".rds")))
    ggsave(filename = here::here(paste0("Summer2023/StanOutput/", str_sub(file, 1L, -5L), "_locintplot_", metric_choice, ".pdf")),
           plot = plot1,
           width = 7,
           height = 4.5,
           units = "in",
           device = cairo_pdf,
           dpi = 300)
    
    ggsave(filename = here::here(paste0("Summer2023/StanOutput/", str_sub(file, 1L, -5L), "_locintplot_", metric_choice, ".png")),
           plot = plot1,
           width = 7,
           height = 4.5,
           units = "in")
    
    # saveRDS(plot1b, here::here(paste0("StanOutput/", str_sub(file, 1L, -5L), "_locintplot_b_", metric_choice, ".rds")))
    # ggsave(filename = here::here(paste0("StanOutput/", str_sub(file, 1L, -5L), "_locintplot_b_", metric_choice, ".pdf")),
    #        plot = plot1b,
    #        width = 7,
    #        height = 4.5,
    #        units = "in",
    #        device = cairo_pdf,
    #        dpi = 300)
    
    
    plot2 <- bayesplot::mcmc_intervals(fit$draws("mu_reef"), prob_outer = 0.95, point_size = 2, inner_size = 1.25)
    # plot2$data$reef3 <- as.integer(str_sub(plot2$data$parameter, 9L, -2L))
    # reefs2 <- merge(reefs, plot2$data[, c(1, 10)], by = "reef3", all.x = TRUE)
    # setorder(reefs2, "locality", "reef")
    
    plot2 <- plot2 +
      geom_vline(xintercept = 0, color = "grey90", lwd = 0.5) +
      geom_point(data = plot2$data, aes(x = m, y = parameter, shape = "Model", color = "Model", fill = "Model", size = "Model"), stroke = 0.25) +
      {if(param == "median_age_bp"){
        geom_vline(data = Detal, aes(xintercept = median_age_bp, lty = "Dominguez et al. (2016)"), lwd = 0.25) #, color = "Dominguez et al. (2016)", size = "Dominguez et al. (2016)"
      } else{
        geom_vline(data = Detal, aes(xintercept = tav_iqr, lty = "Dominguez et al. (2016)"), lwd = 0.25) #, color = "Dominguez et al. (2016)", size = "Dominguez et al. (2016)"
      }} +
      {if(metric_choice == "median"){
        geom_point(data = Data, aes(x = med_param_reef, y = 32 - reef3, shape = "Data", color = "Data", fill = "Data", size = "Data"), stroke = 1)
      } else{
        geom_point(data = Data, aes(x = mean_param_reef, y = 32 - reef3, shape = "Data", color = "Data", fill = "Data", size = "Data"), stroke = 1)
      }} +
      # {if(param == "median_age_bp"){
      #    geom_vline(data = Detal, aes(xintercept = median_age_bp, color = "Dominguez et al. (2016)", size = "Dominguez et al. (2016)", lty = "Dominguez et al. (2016)"), lwd = 0.5)
      #  } else{
      #    geom_vline(data = Detal, aes(xintercept = Sample_corrected_posterior_age_estimate, color = "Dominguez et al. (2016)", size = "Dominguez et al. (2016)", lty = "Dominguez et al. (2016)"), lwd = 0.5)
      #  }} +
      {if(param == "median_age_bp"){
        if(metric_choice == "median"){
          xlim(c(min(c(plot2$data$ll, Data$med_param_reef, Detal$median_age_bp)), max(c(plot2$data$hh, Data$med_param_reef, Detal$median_age_bp))))
        } else{
          xlim(c(min(c(plot2$data$ll, Data$mean_param_reef, Detal$median_age_bp)), max(c(plot2$data$hh, Data$mean_param_reef, Detal$median_age_bp))))
        }
      } else{
        if(metric_choice == "median"){
          xlim(c(min(c(plot2$data$ll, Data$med_param_reef, Detal$tav_iqr)), max(c(plot2$data$hh, Data$med_param_reef, Detal$tav_iqr))))
        } else{
          xlim(c(min(c(plot2$data$ll, Data$mean_param_reef, Detal$tav_iqr)), max(c(plot2$data$hh, Data$mean_param_reef, Detal$tav_iqr))))
        }
      }} +
      guides(size = "none") +
      scale_shape_manual(values = c("Data" = 3, "Model" = 21)) + #, "Dominguez et al. (2016)" = NA
      scale_color_manual(values = c("Data" = "red", "Model" = "#011f4b")) + #, "Dominguez et al. (2016)" = "black"
      scale_linetype_manual(values = c("Dominguez et al. (2016)" = "dashed")) + #"Data" = NA, "Model" = NA, 
      scale_fill_manual(values = c("Data" = "red", "Model" = "#d1e1ec")) + #, "Dominguez et al. (2016)" = NA
      scale_size_manual(values = c("Data" = 1, "Model" = 2)) + #, "Dominguez et al. (2016)" = 1
      # scale_y_discrete(breaks = plot2$data$parameter[order(reefs2$parameter)], labels = paste0(reefs2$locality, " R", reefs2$reef)) +
      scale_y_discrete(breaks = plot2$data$parameter, labels = paste0(reefs$locality, " R", reefs$reef_lab), limits = rev) +
      #labs(x = param, y = "Reef", title = str_sub(file, 1L, -5L), shape = "Means", colour = "Means", fill = "Means", size = "Means", lty = "Means") +
      labs(x = ifelse(param == "median_age_bp", "Median age (years before collection year)", "Total age variation IQR"), y = "Reef", title = str_sub(file, 1L, -5L), shape = NULL, colour = NULL, fill = NULL, size = NULL, lty = NULL) +
      theme(text = element_text(family = "Arial"),
            legend.position = "right")
    
    if(ifelse(param == "median_age_bp", ifelse(metric_choice == "median", max(c(plot2$data$hh, Data$med_param_reef, Detal$median_age_bp)) == max(Detal$median_age_bp) & max(Detal$median_age_bp) - max(unique(sort(c(plot1$data$hh, Data$med_param_reef, Detal$median_age_bp)))[length(unique(c(plot1$data$hh, Data$med_param_reef, Detal$median_age_bp))) - 1]) > 150, 
                                               max(c(plot2$data$hh, Data$mean_param_reef, Detal$median_age_bp)) == max(Detal$median_age_bp)) & max(Detal$median_age_bp) - max(unique(sort(c(plot1$data$hh, Data$mean_param_reef, Detal$median_age_bp)))[length(unique(c(plot1$data$hh, Data$mean_param_reef, Detal$median_age_bp))) - 1]) > 150,
              ifelse(metric_choice == "median", max(c(plot2$data$hh, Data$med_param_reef, Detal$tav_iqr)) == max(Detal$tav_iqr) & max(Detal$tav_iqr) - max(unique(sort(c(plot1$data$hh, Data$med_param_reef, Detal$tav_iqr)))[length(unique(c(plot1$data$hh, Data$med_param_reef, Detal$tav_iqr))) - 1]) > 150, 
                     max(c(plot2$data$hh, Data$mean_param_reef, Detal$tav_iqr)) == max(Detal$tav_iqr) & max(Detal$tav_iqr) - max(unique(sort(c(plot1$data$hh, Data$mean_param_reef, Detal$tav_iqr)))[length(unique(c(plot1$data$hh, Data$mean_param_reef, Detal$tav_iqr))) - 1]) > 150))){
      
      plot2z <- plot2 + 
        {if(param == "median_age_bp"){
          if(metric_choice == "median"){
            scale_x_continuous(breaks = seq(0, #min(plyr::round_any(c(min(c(plot2$data$ll, Data$mean_param_reef, Detal$median_age_bp)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, floor)), 
                                            max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)) + 40, by = 50)) 
          } else{
            scale_x_continuous(breaks = seq(0, #min(plyr::round_any(c(min(c(plot2$data$ll, Data$mean_param_reef, Detal$median_age_bp)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, floor)), 
                                            max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)) + 40, by = 50)) 
          }
        } else{
          if(metric_choice == "median"){
            scale_x_continuous(breaks = seq(0, #min(plyr::round_any(c(min(c(plot2$data$ll, Data$mean_param_reef, Detal$median_age_bp)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, floor)), 
                                            max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)) + 40, by = 50)) 
          } else{
            scale_x_continuous(breaks = seq(0, #min(plyr::round_any(c(min(c(plot2$data$ll, Data$mean_param_reef, Detal$median_age_bp)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, floor)), 
                                            max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)) + 40, by = 50)) 
          }
        }} +
        {if(param == "median_age_bp"){
          if(metric_choice == "median"){
            coord_cartesian(xlim = c(-10, #plyr::round_any(min(c(plot2$data$ll, Data$mean_param_reef, Detal$median_age_bp)), 50, floor) - 20, 
                                     plyr::round_any(max(c(plot2$data$hh, Data$med_param_reef, sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1])), 50, ceiling) + 40)) 
          } else{
            coord_cartesian(xlim = c(-10, #plyr::round_any(min(c(plot2$data$ll, Data$mean_param_reef, Detal$median_age_bp)), 50, floor) - 20, 
                                     plyr::round_any(max(c(plot2$data$hh, Data$mean_param_reef, sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1])), 50, ceiling) + 40))
          }
          
        } else{
          if(metric_choice == "median"){
            coord_cartesian(xlim = c(-10, #plyr::round_any(min(c(plot2$data$ll, Data$mean_param_reef, Detal$Sample_corrected_posterior_age_estimate)), 50, floor) - 20, 
                                     plyr::round_any(max(c(plot2$data$hh, Data$med_param_reef, sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2])), 50, ceiling) + 40)) 
          } else{
            coord_cartesian(xlim = c(-10, #plyr::round_any(min(c(plot2$data$ll, Data$mean_param_reef, Detal$Sample_corrected_posterior_age_estimate)), 50, floor) - 20, 
                                     plyr::round_any(max(c(plot2$data$hh, Data$mean_param_reef, sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2])), 50, ceiling) + 40))
          }
        }} +
        theme_classic() +
        theme(legend.position = "none",
              text = element_text(size = 9),
              plot.title = element_blank(),
              axis.title = element_blank())
      
      if(param == "median_age_bp" & 
         length(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, ifelse(metric_choice == "median", Data$med_param_reef, Data$mean_param_reef))), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))) > 1 |
         param == "tav_iqr" & length(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, ifelse(metric_choice == "median", Data$med_param_reef, Data$mean_param_reef))), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))) > 1){
        
        plot2o <- plot2 + 
          # scale_y_discrete(breaks = plot2$data$parameter, labels = rep("/\u005c", 31)) +
          #scale_y_discrete(breaks = plot2$data$parameter, labels = rep("\u2591 ", 31)) +
          {if(param == "median_age_bp"){
            if(metric_choice == "median"){
              scale_x_continuous(breaks = seq(plyr::round_any(min(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, floor), 
                                              plyr::round_any(min(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, ceiling), by = 50)[1:2])
            } else{
              scale_x_continuous(breaks = seq(plyr::round_any(min(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, floor), 
                                              plyr::round_any(min(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, ceiling), by = 50)[1:2])
            }
          } else{
            if(metric_choice == "median"){
              scale_x_continuous(breaks = seq(plyr::round_any(min(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, floor), 
                                              plyr::round_any(min(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, ceiling), by = 50)[1:2])
            } else{
              scale_x_continuous(breaks = seq(plyr::round_any(min(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, floor), 
                                              plyr::round_any(min(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, ceiling), by = 50)[1:2])
            }
          }} +
          {if(param == "median_age_bp"){
            if(metric_choice == "median"){
              coord_cartesian(xlim = c(plyr::round_any(min(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, floor), 
                                       plyr::round_any(min(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, ceiling)))
            } else{
              coord_cartesian(xlim = c(plyr::round_any(min(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, floor), 
                                       plyr::round_any(min(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, ceiling)))
            }
          } else{
            if(metric_choice == "median"){
              coord_cartesian(xlim = c(plyr::round_any(min(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, floor), 
                                       plyr::round_any(min(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, ceiling)))
            } else{
              coord_cartesian(xlim = c(plyr::round_any(min(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, floor), 
                                       plyr::round_any(min(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, ceiling)))
            }
          }} +
          theme_classic() +
          theme(legend.position = "none",
                plot.title = element_blank(),
                plot.background = element_blank(),
                text = element_text(size = 9),
                # axis.text.y = element_text(angle = -90, color = "grey50", size = 17, vjust = -0.1, hjust = 0.5),
                axis.text.y = element_blank(),
                axis.title = element_blank(),
                axis.ticks.y = element_blank(),
                #axis.ticks.length.y = unit(0.25, "in"),
                axis.line.y = element_blank()) #,
        # legend.margin = margin(t = 0, unit='cm'),
        # legend.text = element_text(margin = margin(l = 0), hjust = 0),
        # legend.spacing.x = unit(0, "in"),
        # #legend.box.spacing = margin(0),
        # legend.justification = "left")
        
        plot2o2 <- plot2 + 
          # scale_y_discrete(breaks = plot2$data$parameter, labels = rep("/\u005c", 31)) +
          #scale_y_discrete(breaks = plot2$data$parameter, labels = rep("\u2591 ", 31)) +
          {if(param == "median_age_bp"){
            if(metric_choice == "median"){
              scale_x_continuous(breaks = seq(plyr::round_any(max(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, floor), 
                                              plyr::round_any(max(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, ceiling), by = 50)[2:3])
            } else{
              scale_x_continuous(breaks = seq(plyr::round_any(max(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, floor), 
                                              plyr::round_any(max(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, ceiling), by = 50)[2:3])
            } 
          } else{
            if(metric_choice == "median"){
              scale_x_continuous(breaks = seq(plyr::round_any(max(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, floor), 
                                              plyr::round_any(max(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, ceiling), by = 50)[2:3])
            } else{
              scale_x_continuous(breaks = seq(plyr::round_any(max(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, floor), 
                                              plyr::round_any(max(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, ceiling), by = 50)[2:3])
            }
          }} +
          {if(param == "median_age_bp"){
            if(metric_choice == "median"){
              coord_cartesian(xlim = c(plyr::round_any(max(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, floor), 
                                       plyr::round_any(max(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, ceiling)))
            } else{
              coord_cartesian(xlim = c(plyr::round_any(max(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, floor), 
                                       plyr::round_any(max(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, ceiling)))
            }
          } else{
            if(metric_choice == "median"){
              coord_cartesian(xlim = c(plyr::round_any(max(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, floor), 
                                       plyr::round_any(max(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, ceiling)))
            } else{
              coord_cartesian(xlim = c(plyr::round_any(max(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, floor), 
                                       plyr::round_any(max(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, ceiling)))
            }
          }} +
          theme_classic() +
          theme(plot.title = element_blank(),
                plot.background = element_blank(),
                panel.background = element_blank(),
                text = element_text(size = 9),
                # axis.text.y = element_text(angle = -90, color = "grey50", size = 17, vjust = -0.1, hjust = 0.5),
                axis.text.y = element_blank(),
                axis.title = element_blank(),
                axis.ticks.y = element_blank(),
                #axis.ticks.length.y = unit(0.25, "in"),
                axis.line.y = element_blank(),
                legend.margin = margin(t = 0, unit='cm'),
                legend.text = element_text(margin = margin(l = 0), hjust = 0),
                legend.spacing.x = unit(0, "in"),
                #legend.box.spacing = margin(0),
                legend.justification = "left")
        
        ylabel <- ggplot(data.frame(l = plot2$labels$y, x = 1, y = 1)) +
          geom_text(aes(x, y, label = l), size = 2.75, angle = 90) + 
          theme_void() +
          coord_cartesian(clip = "off")
        
        xlabel <- ggplot(data.frame(l = plot2$labels$x, x = 1, y = 1)) +
          geom_text(aes(x, y, label = l), size = 2.75, angle = 0) + 
          theme_void() +
          coord_cartesian(clip = "off")
        
        separator <- ggplot() +
          #geom_tile(data = data.frame(x_min = 0, x_max = 90, y_min = min(plot2$data$parameter), y_max = max(plot2$data$parameter)), aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max), fill = "grey50") + 
          #geom_tile_pattern(aes(x = 50, y = plot2$data$parameter), fill = NA, pattern = "stripe", pattern_density = 0.5, pattern_spacing = 1, pattern_colour = "grey50", pattern_fill = "grey50") + 
          #geom_hline(yintercept = 0.45, color = "black", lwd = 0.5) +
          geom_polygon(data = data.frame(x = c(7, 8, 10, 9), y = c(0, 0, 0.9, 0.9)), aes(x = x, y = y), fill = "white", color = "white") +
          geom_segment(aes(x = c(7, 8), y = c(0, 0), xend = c(9, 10), yend = c(0.9, 0.9)), color = "black", lwd = 0.5) +
          ylim(c(sep_ylim_min, sep_ylim_max)) +
          xlim(c(sep_xlim_min, sep_xlim_max)) +
          theme_void() +
          theme(plot.background = element_blank()) +
          coord_cartesian(clip = "off")
        
        layout <- c(
          area(t = 1, l = 1, b = 30, r = 5),
          area(t = 1, l = 6, b = 30, r = 50),
          area(t = 29, l = 47, b = 31, r = 50),
          area(t = 1, l = 50, b = 30, r = 59),
          area(t = 1, l = 55, b = 30, r = 64),
          area(t = 29, l = 56, b = 31, r = 59),
          area(t = 31, l = 6, b = 32, r = 64)
        )
        
        plot2 <- ylabel + plot2z + separator + plot2o + plot2o2 + separator + xlabel + plot_layout(guides = "collect", design = layout) + plot_annotation(title = str_sub(file, 1L, -5L), theme = theme(plot.title = element_text(hjust = 0.12),
                                                                                                                                                                                                        text = element_text(size = 7))) #,
        #legend.margin = margin(t = 0, unit='cm'),
        #legend.text = element_text(margin = margin(l = 0), hjust = 0),
        #legend.spacing.x = unit(0, "in"),
        #legend.box.spacing = margin(0),
        #legend.justification = "left"))
      } else{
        
        plot2o <- plot2 + 
          # scale_y_discrete(breaks = plot2$data$parameter, labels = rep("/\u005c", 31)) +
          #scale_y_discrete(breaks = plot2$data$parameter, labels = rep("\u2591 ", 31)) +
          {if(param == "median_age_bp"){
            if(metric_choice == "median"){
              scale_x_continuous(breaks = seq(plyr::round_any(min(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, floor), 
                                              plyr::round_any(max(plot2$layers[[9]]$data$median_age_bp), 100, ceiling), by = 50))
            } else{
              scale_x_continuous(breaks = seq(plyr::round_any(min(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, floor), 
                                              plyr::round_any(max(plot2$layers[[9]]$data$median_age_bp), 100, ceiling), by = 50))
            }
          } else{
            if(metric_choice == "median"){
              scale_x_continuous(breaks = seq(plyr::round_any(min(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, floor), 
                                              plyr::round_any(max(plot2$layers[[9]]$data$tav_iqr), 100, ceiling), by = 50))
            } else{
              scale_x_continuous(breaks = seq(plyr::round_any(min(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, floor), 
                                              plyr::round_any(max(plot2$layers[[9]]$data$tav_iqr), 100, ceiling), by = 50))
            }
          }} +
          {if(param == "median_age_bp"){
            if(metric_choice == "median"){
              coord_cartesian(xlim = c(plyr::round_any(min(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, floor), 
                                       plyr::round_any(max(plot2$layers[[9]]$data$median_age_bp), 100, ceiling)))
            } else{
              coord_cartesian(xlim = c(plyr::round_any(min(subset(plot2$layers[[9]]$data$median_age_bp, plot2$layers[[9]]$data$median_age_bp > max(plyr::round_any(c(max(c(plot2$data$hh, Data$med_param_reef)), sort(Detal$median_age_bp)[length(Detal$median_age_bp) - 1]), 50, ceiling)))), 100, floor), 
                                       plyr::round_any(max(plot2$layers[[9]]$data$median_age_bp), 100, ceiling)))
            }
          } else{
            if(metric_choice == "median"){
              coord_cartesian(xlim = c(plyr::round_any(min(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, floor), 
                                       plyr::round_any(max(plot2$layers[[9]]$data$tav_iqr), 100, ceiling)))
            } else{
              coord_cartesian(xlim = c(plyr::round_any(min(subset(plot2$layers[[9]]$data$tav_iqr, plot2$layers[[9]]$data$tav_iqr > max(plyr::round_any(c(max(c(plot2$data$hh, Data$mean_param_reef)), sort(Detal$tav_iqr)[length(Detal$tav_iqr) - 2]), 50, ceiling)))), 100, floor), 
                                       plyr::round_any(max(plot2$layers[[9]]$data$tav_iqr), 100, ceiling)))
            }
          }} +
          theme_classic() +
          theme(plot.title = element_blank(),
                plot.background = element_blank(),
                text = element_text(size = 9),
                # axis.text.y = element_text(angle = -90, color = "grey50", size = 17, vjust = -0.1, hjust = 0.5),
                axis.text.y = element_blank(),
                axis.title = element_blank(),
                axis.ticks.y = element_blank(),
                #axis.ticks.length.y = unit(0.25, "in"),
                axis.line.y = element_blank(),
                legend.margin = margin(t = 0, unit='cm'),
                legend.text = element_text(margin = margin(l = 0), hjust = 0),
                legend.spacing.x = unit(0, "in"),
                #legend.box.spacing = margin(0),
                legend.justification = "left")
        
        # plot2c <- plot2z + plot2o + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a", 
        #                                                                               tag_prefix = "(", 
        #                                                                               tag_suffix = ")", 
        #                                                                               title = str_sub(file, 1L, -5L))
        
        ylabel <- ggplot(data.frame(l = plot2$labels$y, x = 1, y = 1)) +
          geom_text(aes(x, y, label = l), size = 2.75, angle = 90) + 
          theme_void() +
          coord_cartesian(clip = "off")
        
        xlabel <- ggplot(data.frame(l = plot2$labels$x, x = 1, y = 1)) +
          geom_text(aes(x, y, label = l), size = 2.75, angle = 0) + 
          theme_void() +
          coord_cartesian(clip = "off")
        
        separator <- ggplot() +
          #geom_tile(data = data.frame(x_min = 0, x_max = 90, y_min = min(plot2$data$parameter), y_max = max(plot2$data$parameter)), aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max), fill = "grey50") + 
          #geom_tile_pattern(aes(x = 50, y = plot2$data$parameter), fill = NA, pattern = "stripe", pattern_density = 0.5, pattern_spacing = 1, pattern_colour = "grey50", pattern_fill = "grey50") + 
          #geom_hline(yintercept = 0.45, color = "black", lwd = 0.5) +
          geom_polygon(data = data.frame(x = c(7, 8, 10, 9), y = c(0, 0, 0.9, 0.9)), aes(x = x, y = y), fill = "white", color = "white") +
          geom_segment(aes(x = c(7, 8), y = c(0, 0), xend = c(9, 10), yend = c(0.9, 0.9)), color = "black", lwd = 0.5) +
          ylim(c(sep_ylim_min, sep_ylim_max)) +
          xlim(c(sep_xlim_min, sep_xlim_max)) +
          theme_void() +
          theme(plot.background = element_blank()) +
          coord_cartesian(clip = "off")
        
        layout <- c(
          area(t = 1, l = 1, b = 30, r = 5),
          area(t = 1, l = 6, b = 30, r = 50),
          area(t = 29, l = 47, b = 31, r = 50),
          area(t = 1, l = 50, b = 30, r = 59),
          area(t = 31, l = 6, b = 32, r = 59)
        )
        
        plot2 <- ylabel + plot2z + separator + plot2o + xlabel + plot_layout(guides = "collect", design = layout) + plot_annotation(title = str_sub(file, 1L, -5L), theme = theme(plot.title = element_text(hjust = 0.12),
                                                                                                                                                                                  text = element_text(size = 7))) #,
        #legend.margin = margin(t = 0, unit='cm'),
        #legend.text = element_text(margin = margin(l = 0), hjust = 0),
        #legend.spacing.x = unit(0, "in"),
        #legend.box.spacing = margin(0),
        #legend.justification = "left"))
      }
    }
    #plot2
    
    saveRDS(plot2, here::here(paste0("Summer2023/StanOutput/", str_sub(file, 1L, -5L), "_reefintplot_", metric_choice, ".rds")))
    ggsave(filename = here::here(paste0("Summer2023/StanOutput/", str_sub(file, 1L, -5L), "_reefintplot_", metric_choice, ".pdf")),
           plot = plot2,
           width = 7,
           height = 4.5,
           units = "in",
           device = cairo_pdf,
           dpi = 300)

    ggsave(filename = here::here(paste0("Summer2023/StanOutput/", str_sub(file, 1L, -5L), "_reefintplot_", metric_choice, ".png")),
           plot = plot2,
           width = 7,
           height = 4.5,
           units = "in")
    
    #Sigma plots
    #Sigmas versus Detal sd
    Data[, `:=` (sd_medage_loc = sd(median_age_bp),
                 sd_tav_loc = sd(tav_iqr)), by = locality]
    Data[, `:=` (sd_medage_reef = sd(median_age_bp),
                 sd_tav_reef = sd(tav_iqr)), by = reefID]
    Data[, `:=` (sd_medage_hole = sd(median_age_bp),
                 sd_tav_hole = sd(tav_iqr)), by = holeID]
    Data[, `:=` (sd_medage_hole2 = sd(median_age_bp),
                 sd_tav_hole2 = sd(tav_iqr))]
    Data[, `:=` (sd_medage_reef2 = sd(Data[, mean(median_age_bp), by = reefID]$V1),
                 sd_tav_reef2 = sd(Data[, mean(tav_iqr), by = reefID]$V1))]
    Data[, `:=` (sd_medage_loc2 = sd(Data[, mean(median_age_bp), by = locality]$V1),
                 sd_tav_loc2 = sd(Data[, mean(tav_iqr), by = locality]$V1))]
    Data_sd <- melt(Data[, .(locality, reefID, holeID, sample, 
                             sd_medage_loc, sd_tav_loc, sd_medage_reef, sd_tav_reef, sd_medage_hole, sd_tav_hole, 
                             sd_medage_loc2, sd_tav_loc2, sd_medage_reef2, sd_tav_reef2, sd_medage_hole2, sd_tav_hole2)], 
                    id.vars = c("locality", "reefID", "holeID"), measure.vars = c("sd_medage_loc", "sd_tav_loc", "sd_medage_reef", "sd_tav_reef", "sd_medage_hole", "sd_tav_hole",
                                                                                  "sd_medage_loc2", "sd_tav_loc2", "sd_medage_reef2", "sd_tav_reef2", "sd_medage_hole2", "sd_tav_hole2"))
    setDT(Data_sd)
    Data_sd[, `:=` (yname = fcase(str_detect(variable, "loc"), "sigma_locality",
                                  str_detect(variable, "reef"), "sigma_reef",
                                  str_detect(variable, "hole"), "sigma_hole"),
                    metric = ifelse(str_detect(variable, "medage"), "medage", "tav_iqr")), by = list(locality, reefID, holeID, variable)]
    
    
    plot3 <- mcmc_intervals(fit$draws(c("sigma_locality", "sigma_reef", "sigma_hole")), prob_outer = 0.95, point_size = 2, inner_size = 1.25)
    plot3 <- plot3 +
      geom_point(data = plot3$data, aes(x = m, y = parameter, shape = "Model", color = "Model", fill = "Model", size = "Model"), stroke = 0.25) +
      {if(param == "median_age_bp"){
        geom_vline(aes(xintercept = sd(Detal$median_age_bp), lty = "Dominguez et al. (2016)"), lwd = 0.25) #, color = "Dominguez et al. (2016)", size = "Dominguez et al. (2016)"
      } else{
        geom_vline(aes(xintercept = sd(Detal$tav_iqr), lty = "Dominguez et al. (2016)"), lwd = 0.25) #, color = "Dominguez et al. (2016)", size = "Dominguez et al. (2016)"
      }} +
      # {if(param == "median_age_bp"){
      #   if(metric_choice == "median"){
      #     geom_point(data = distinct(Data_sd[str_detect(variable, "2") & metric == "medage", .(value, yname)]), aes(x = value, y = yname, shape = "Data", color = "Data", fill = "Data", size = "Data"), stroke = 1)
      #   } else{
      #     geom_point(data = distinct(Data_sd[str_detect(variable, "2") & metric == "medage", .(value, yname)]), aes(x = value, y = yname, shape = "Data", color = "Data", fill = "Data", size = "Data"), stroke = 1)
      #   }
      # } else{
      #   if(metric_choice == "median"){
      #     geom_point(data = distinct(Data_sd[str_detect(variable, "2") & metric == "cpe", .(value, yname)]), aes(x = value, y = yname, shape = "Data", color = "Data", fill = "Data", size = "Data"), stroke = 1)
      #   } else{
      #     geom_point(data = distinct(Data_sd[str_detect(variable, "2") & metric == "cpe", .(value, yname)]), aes(x = value, y = yname, shape = "Data", color = "Data", fill = "Data", size = "Data"), stroke = 1)
    #   }
    # }} +
      guides(size = "none") +
      scale_x_continuous(trans = "log10") +
      scale_shape_manual(values = c("Data" = 3, "Model" = 21)) + #, "Dominguez et al. (2016)" = NA
      scale_color_manual(values = c("Data" = "red", "Model" = "#011f4b")) + #, "Dominguez et al. (2016)" = "black"
      scale_linetype_manual(values = c("Dominguez et al. (2016)" = "dashed")) + #"Data" = NA, "Model" = NA, 
      scale_fill_manual(values = c("Data" = "red", "Model" = "#d1e1ec")) + #, "Dominguez et al. (2016)" = NA
      scale_size_manual(values = c("Data" = 1, "Model" = 2)) + #, "Dominguez et al. (2016)" = 1
      # scale_y_discrete(breaks = plot2$data$parameter[order(reefs2$parameter)], labels = paste0(reefs2$locality, " R", reefs2$reef)) +
      scale_y_discrete(breaks = plot3$data$parameter, labels = c(expression(sigma[c]), expression(sigma[r]), expression(sigma[h])), limits = rev) +
      #labs(x = param, y = "Reef", title = str_sub(file, 1L, -5L), shape = "Means", colour = "Means", fill = "Means", size = "Means", lty = "Means") +
      labs(x = ifelse(param == "median_age_bp", "Median age (years before collection year)", "Total age variation IQR"), y = "Sigma", title = str_sub(file, 1L, -5L), shape = NULL, colour = NULL, fill = NULL, size = NULL, lty = NULL) +
      theme(text = element_text(family = "sans", size = 9),
            axis.text.y = element_text(size = 11),
            axis.title.y = element_text(size = 9),
            legend.text = element_text(family = "sans", size = 9),
            legend.position = "right")
    #coord_cartesian(ylim = c(1, 4.5))
    #plot3
    
    plot3b <- plot3 + labs(y = NULL) + theme(plot.title = element_blank(),
                                             axis.title = element_blank())
    
    ylabel <- ggplot(data.frame(l = plot3$labels$y, x = 1, y = 1)) +
      geom_text(aes(x, y, label = l), size = 2.75, angle = 90) + 
      theme_void() +
      coord_cartesian(clip = "off")
    
    layout3 <- c(
      area(t = 1, l = 1, b = 30, r = 3),
      area(t = 1, l = 4, b = 30, r = 50)
    )
    
    plot3 <- ylabel + plot3b + plot_layout(design = layout3) + plot_annotation(title = str_sub(file, 1L, -5L), 
                                                                               theme = theme(plot.title = element_text(hjust = 0.12)))
    
    saveRDS(plot3, here::here(paste0("Summer2023/StanOutput/", str_sub(file, 1L, -5L), "_sigmaintplot_", metric_choice, ".rds")))
    ggsave(filename = here::here(paste0("Summer2023/StanOutput/", str_sub(file, 1L, -5L), "_sigmaintplot_", metric_choice, ".pdf")),
           plot = plot3,
           width = 7,
           height = 4.5,
           units = "in",
           device = cairo_pdf,
           dpi = 300)
    
    ggsave(filename = here::here(paste0("Summer2023/StanOutput/", str_sub(file, 1L, -5L), "_sigmaintplot_", metric_choice, ".png")),
           plot = plot3,
           width = 7,
           height = 4.5,
           units = "in")
    
  }
}



#Put plots together in a logical way for the supplement
tav1525loc <- readRDS(here::here(paste0("Summer2023/StanOutput/HOBS_tavxspace_15to25cm_locintplot_", metric_choice, ".rds")))
tav1525reef <- readRDS(here::here(paste0("Summer2023/StanOutput/HOBS_tavxspace_15to25cm_reefintplot_", metric_choice, ".rds")))
tav1525sigma <- readRDS(here::here(paste0("Summer2023/StanOutput/HOBS_tavxspace_15to25cm_sigmaintplot_", metric_choice, ".rds")))

# #The commented out code is only necessary if using means, in which case the loc plot needs an x-axis break, so is a patchwork
# #instead of a single ggplot.
# cpe1525loc$patches$layout$design <- c(
#   area(t = 1, l = 1, b = 30, r = 3),
#   area(t = 1, l = 4, b = 30, r = 50),
#   area(t = 27, l = 46, b = 31, r = 50),
#   area(t = 1, l = 50, b = 30, r = 59)
# )
#
# cpe1525loc[[2]] <- cpe1525loc[[2]] + 
#   {if(metric_choice == "median"){
#     xlim(c(-15, sampsum_15to25cm[, med_param_loc := median(Sample_corrected_posterior_age_estimate), by = Locality][, max(med_param_loc)] + 50))
#   } else{
#     xlim(c(-15, sampsum_15to25cm[, mean_param_loc := mean(Sample_corrected_posterior_age_estimate), by = Locality][, max(mean_param_loc)] + 50))
#   }} +
#   {if(metric_choice == "median"){
#     coord_cartesian(xlim = c(-15, sampsum_15to25cm[, med_param_loc := median(Sample_corrected_posterior_age_estimate), by = Locality][, max(med_param_loc)] + 50))
#   } else{
#     coord_cartesian(xlim = c(-15, sampsum_15to25cm[, mean_param_loc := mean(Sample_corrected_posterior_age_estimate), by = Locality][, max(mean_param_loc)] + 50))
#   }} +
#   {if(metric_choice == "median"){
#     scale_x_continuous(breaks = seq(0, #plyr::round_any(min(c(plot1$data$ll, Data$mean_param_loc)), 50, floor), 
#                                     plyr::round_any(max(c(cpe1525loc[[2]]$data$ll, sampsum_15to25cm[, med_param_loc := median(Sample_corrected_posterior_age_estimate), by = Locality][, max(med_param_loc)])), 50, ceiling), by = 50))
#   } else{
#     scale_x_continuous(breaks = seq(0, #plyr::round_any(min(c(plot1$data$ll, Data$mean_param_loc)), 50, floor), 
#                                     plyr::round_any(max(c(cpe1525loc[[2]]$data$ll, sampsum_15to25cm[, mean_param_loc := mean(Sample_corrected_posterior_age_estimate), by = Locality][, max(mean_param_loc)])), 50, ceiling), by = 50))
#   }}
# 
# cpe1525loc[[4]] <- cpe1525loc[[4]] + theme(legend.box.margin = margin(l = 30))

tav1525sigma <- tav1525sigma + labs(title = NULL, x = NULL) + theme(legend.position = "none")

tav1525sigma$patches$layout$design <- c(
  area(t = 1, l = 1, b = 30, r = 3),
  area(t = 1, l = 4, b = 30, r = 50)
)

tav1525loc <- tav1525loc + labs(title = NULL, x = NULL) #+ theme(legend.position = "none")

tav1525reef$patches$layout$design <- c(
  area(t = 1, l = 1, b = 30, r = 2),
  area(t = 1, l = 3, b = 30, r = 50),
  area(t = 29, l = 47, b = 31, r = 50),
  area(t = 1, l = 50, b = 30, r = 59),
  area(t = 1, l = 55, b = 30, r = 64),
  area(t = 29, l = 56, b = 31, r = 59),
  area(t = 31, l = 3, b = 32, r = 64)
)
tav1525reef[[5]] <- tav1525reef[[5]] + theme(legend.position = "none")
#cpe1525reef[[5]] <- cpe1525reef[[5]] + theme(legend.box.margin = margin(l = 30))
#cpe1525reef <- cpe1525reef & theme(plot.background = element_rect(colour = "black", fill = NA, size = 1))

layout4 <- "
AABB
AABB
CCCD
CCCD
CCCD
CCCD
CCCD
"

tavxspace15to25 <- ((tav1525sigma - tav1525loc) + tav1525reef + guide_area()) + plot_layout(guides = "collect", design = layout4) + plot_annotation(tag_levels = list(c("A", rep("", length(tav1525sigma$patches$layout$design$t) - 1),
                                                                                                                                                                        "B", #rep("", length(cpe1525loc$patches$layout$design$t)), 
                                                                                                                                                                        #"C", rep("", length(cpe1525reef$patches$layout$design$t) - 1))), "") &
                                                                                                                                                                        "C", rep("", length(tav1525sigma$patches$layout$design$t)))), "") &
  theme(legend.text = element_text(size = 9))
#cpexspace15to25

ggsave(filename = here::here(paste0("Summer2023/StanOutput/tavxspace15to25_", metric_choice, ".pdf")),
       plot = tavxspace15to25,
       width = 6.75,
       height = 6,
       units = "in",
       device = cairo_pdf,
       dpi = 300)

ggsave(filename = here::here(paste0("Summer2023/StanOutput/tavxspace15to25_", metric_choice, ".png")),
       plot = tavxspace15to25,
       width = 6.75,
       height = 6,
       units = "in")


tav2535loc <- readRDS(here::here(paste0("Summer2023/StanOutput/HOBS_tavxspace_25to35cm_locintplot_", metric_choice, ".rds")))
tav2535reef <- readRDS(here::here(paste0("Summer2023/StanOutput/HOBS_tavxspace_25to35cm_reefintplot_", metric_choice, ".rds")))
tav2535sigma <-  readRDS(here::here(paste0("Summer2023/StanOutput/HOBS_tavxspace_25to35cm_sigmaintplot_", metric_choice, ".rds")))

# #The commented out code is only necessary if using means, in which case the loc plot needs an x-axis break, so is a patchwork
# cpe2535loc$patches$layout$design <- c(
#   area(t = 1, l = 1, b = 30, r = 3),
#   area(t = 1, l = 4, b = 30, r = 50),
#   area(t = 27, l = 46, b = 31, r = 50),
#   area(t = 1, l = 50, b = 30, r = 59)
# )
# 
# cpe2535loc[[2]] <- cpe2535loc[[2]] +
#   {if(metric_choice == "median"){
#     xlim(c(-15, sampsum_25to35cm[, med_param_loc := median(Sample_corrected_posterior_age_estimate), by = Locality][, max(med_param_loc)] + 50))
#   } else{
#     xlim(c(-15, sampsum_25to35cm[, mean_param_loc := mean(Sample_corrected_posterior_age_estimate), by = Locality][, max(mean_param_loc)] + 50))
#   }} +
#   {if(metric_choice == "median"){
#     coord_cartesian(xlim = c(-15, sampsum_25to35cm[, med_param_loc := median(Sample_corrected_posterior_age_estimate), by = Locality][, max(med_param_loc)] + 50))
#   } else{
#     coord_cartesian(xlim = c(-15, sampsum_25to35cm[, mean_param_loc := mean(Sample_corrected_posterior_age_estimate), by = Locality][, max(mean_param_loc)] + 50))
#   }} +
#   {if(metric_choice == "median"){
#     scale_x_continuous(breaks = seq(0, #plyr::round_any(min(c(plot1$data$ll, Data$mean_param_loc)), 50, floor), 
#                                     plyr::round_any(max(c(cpe2535loc[[2]]$data$ll, sampsum_25to35cm[, med_param_loc := median(Sample_corrected_posterior_age_estimate), by = Locality][, max(med_param_loc)])), 50, ceiling), by = 50))
#   } else{
#     scale_x_continuous(breaks = seq(0, #plyr::round_any(min(c(plot1$data$ll, Data$mean_param_loc)), 50, floor), 
#                                     plyr::round_any(max(c(cpe2535loc[[2]]$data$ll, sampsum_25to35cm[, mean_param_loc := mean(Sample_corrected_posterior_age_estimate), by = Locality][, max(mean_param_loc)])), 50, ceiling), by = 50))
#   }}

tav2535sigma <- tav2535sigma + labs(title = NULL, x = NULL) + theme(legend.position = "none")

tav2535sigma$patches$layout$design <- c(
  area(t = 1, l = 1, b = 30, r = 3),
  area(t = 1, l = 4, b = 30, r = 50)
)

tav2535loc <- tav2535loc + labs(title = NULL, x = NULL) #+ theme(legend.position = "none")

tav2535reef$patches$layout$design <- c(
  area(t = 1, l = 1, b = 30, r = 2),
  area(t = 1, l = 3, b = 30, r = 50),
  area(t = 29, l = 47, b = 31, r = 50),
  area(t = 1, l = 50, b = 30, r = 59),
  area(t = 1, l = 55, b = 30, r = 64),
  area(t = 29, l = 56, b = 31, r = 59),
  area(t = 31, l = 3, b = 32, r = 64)
)

tav2535reef[[5]] <- tav2535reef[[5]] + theme(legend.position = "none")


tavxspace25to35 <- ((tav1525sigma - tav2535loc) + tav2535reef + guide_area()) + plot_layout(guides = "collect", design = layout4) + plot_annotation(tag_levels = list(c("A", rep("", length(tav2535sigma$patches$layout$design$t) - 1), 
                                                                                                                                                                        "B",
                                                                                                                                                                        "C", rep("", length(tav2535reef$patches$layout$design$t)))), "") &
  theme(legend.text = element_text(size = 9))
#cpexspace25to35

ggsave(filename = here::here(paste0("Summer2023/StanOutput/tavxspace25to35_", metric_choice, ".pdf")),
       plot = tavxspace25to35,
       width = 6.75,
       height = 6,
       units = "in",
       device = cairo_pdf,
       dpi = 300)

ggsave(filename = here::here(paste0("Summer2023/StanOutput/tavxspace25to35_", metric_choice, ".png")),
       plot = tavxspace25to35,
       width = 6.75,
       height = 6,
       units = "in")


tavddloc <- readRDS(here::here(paste0("Summer2023/StanOutput/HOBS_tavddxspace_locintplot_", metric_choice, ".rds")))
tavddreef <- readRDS(here::here(paste0("Summer2023/StanOutput/HOBS_tavddxspace_reefintplot_", metric_choice, ".rds")))
tavddsigma <- readRDS(here::here(paste0("Summer2023/StanOutput/HOBS_tavddxspace_sigmaintplot_", metric_choice, ".rds")))

tavddsigma <- tavddsigma + labs(title = NULL, x = NULL) + theme(legend.position = "none")

# cpeddsigma$patches$layout$design <- c(
#   area(t = 1, l = 1, b = 30, r = 3),
#   area(t = 1, l = 4, b = 30, r = 50)
# )

tavddloc <- tavddloc + 
  #xlim(c(min(cpeddreef$layers[[8]]$data$mean_param_reef), max(cpeddreef$layers[[8]]$data$mean_param_reef))) +
  # {if(metric_choice == "median"){
  #   scale_x_continuous(breaks = seq(plyr::round_any(min(cpeddreef$layers[[8]]$data$med_param_reef), 100, ceiling), 
  #                                   plyr::round_any(max(cpeddreef$layers[[8]]$data$med_param_reef), 100, floor), by = 100))
  # } else{
  #   scale_x_continuous(breaks = seq(plyr::round_any(min(cpeddreef$layers[[8]]$data$mean_param_reef), 100, ceiling), 
  #                                   plyr::round_any(max(cpeddreef$layers[[8]]$data$mean_param_reef), 100, floor), by = 100))
  # }} +
  # {if(metric_choice == "median"){
  #   coord_cartesian(xlim = c(min(cpeddreef$layers[[8]]$data$med_param_reef), max(cpeddreef$layers[[8]]$data$med_param_reef)))
  # } else{
#   coord_cartesian(xlim = c(min(cpeddreef$layers[[8]]$data$mean_param_reef), max(cpeddreef$layers[[8]]$data$mean_param_reef)))
# }} +
labs(x = NULL, y = "Locality", title = NULL) + 
  theme(legend.text = element_text(size = 9))

tavddreef <- tavddreef + 
  {if(metric_choice == "median"){
    scale_x_continuous(breaks = seq(plyr::round_any(min(tavddreef$layers[[9]]$data$med_param_reef), 100, ceiling), 
                                    plyr::round_any(max(tavddreef$layers[[9]]$data$med_param_reef), 100, floor), by = 100))
  } else{
    scale_x_continuous(breaks = seq(plyr::round_any(min(tavddreef$layers[[9]]$data$mean_param_reef), 100, ceiling), 
                                    plyr::round_any(max(tavddreef$layers[[9]]$data$mean_param_reef), 100, floor), by = 100))
  }} +
  labs(x = expression("TAV_IQR"["25-35cm"] * " - TAV_IQR"["15-25cm"] * " (yr)"), y = "Reef", title = NULL) + 
  theme(legend.text = element_text(size = 9))

# layout4 <- "
# A
# B
# B
# "

tavddxspace <- ((tavddsigma - tavddloc) + tavddreef + guide_area()) + plot_layout(guides = "collect", design = layout4) + plot_annotation(tag_levels = "A") & #list(c("A", rep("", length(cpeddsigma$patches$layout$design$t) - 1), 
  #       "B",
  #       "C", rep("", length(cpeddreef$patches$layout$design$t)))), "") &
  theme(legend.text = element_text(size = 9))

ggsave(filename = here::here(paste0("Summer2023/StanOutput/tavddxspace_", metric_choice, ".pdf")),
       plot = tavddxspace,
       width = 6.75,
       height = 6,
       units = "in",
       device = cairo_pdf,
       dpi = 300)

ggsave(filename = here::here(paste0("Summer2023/StanOutput/tavddxspace_", metric_choice, ".png")),
       plot = tavddxspace,
       width = 6.75,
       height = 6,
       units = "in")


medage1525loc <- readRDS(here::here(paste0("Summer2023/StanOutput/HOBS_medagexspace_15to25cm_locintplot_", metric_choice, ".rds")))
medage1525reef <- readRDS(here::here(paste0("Summer2023/StanOutput/HOBS_medagexspace_15to25cm_reefintplot_", metric_choice, ".rds")))
medage1525sigma <- readRDS(here::here(paste0("Summer2023/StanOutput/HOBS_medagexspace_15to25cm_sigmaintplot_", metric_choice, ".rds")))

# #The commented out code is only necessary if using means, in which case the loc plot needs an x-axis break, so is a patchwork
# medage1525loc$patches$layout$design <- c(
#   area(t = 1, l = 1, b = 30, r = 3),
#   area(t = 1, l = 4, b = 30, r = 50),
#   area(t = 27, l = 46, b = 31, r = 50),
#   area(t = 1, l = 50, b = 30, r = 59)
# )
# medage1525loc[[2]] <- medage1525loc[[2]] +
#   {if(metric_choice == "median"){
#     xlim(c(-20, sampsum_15to25cm[, med_param_loc := median(median_age_bp), by = Locality][, max(med_param_loc)] + 50))
#   } else{
#     xlim(c(-20, sampsum_15to25cm[, mean_param_loc := mean(median_age_bp), by = Locality][, max(mean_param_loc)] + 50))
#   }} +
#   {if(metric_choice == "median"){
#     coord_cartesian(xlim = c(-20, sampsum_15to25cm[, med_param_loc := median(median_age_bp), by = Locality][, max(med_param_loc)] + 50))
#   } else{
#     coord_cartesian(xlim = c(-20, sampsum_15to25cm[, mean_param_loc := mean(median_age_bp), by = Locality][, max(mean_param_loc)] + 50))
#   }} +
#   {if(metric_choice == "median"){
#     scale_x_continuous(breaks = seq(0, #plyr::round_any(min(c(plot1$data$ll, Data$mean_param_loc)), 50, floor), 
#                                     plyr::round_any(max(c(medage1525loc[[2]]$data$ll, sampsum_15to25cm[, med_param_loc := median(median_age_bp), by = Locality][, max(med_param_loc)])), 50, ceiling), by = 50))
#   } else{
#     scale_x_continuous(breaks = seq(0, #plyr::round_any(min(c(plot1$data$ll, Data$mean_param_loc)), 50, floor), 
#                                     plyr::round_any(max(c(medage1525loc[[2]]$data$ll, sampsum_15to25cm[, mean_param_loc := mean(median_age_bp), by = Locality][, max(mean_param_loc)])), 50, ceiling), by = 50))
#   }}
#   
# medage1525loc[[4]] <- medage1525loc[[4]] + theme(legend.box.margin = margin(l = 35))

medage1525loc <- medage1525loc + labs(title = NULL, x = NULL)

# medage1525loc$patches$layout$design <- c(
#   area(t = 1, l = 1, b = 31, r = 3),
#   area(t = 1, l = 4, b = 31, r = 50),
#   area(t = 28, l = 47, b = 31, r = 50),
#   area(t = 1, l = 50, b = 31, r = 59)
# )

medage1525reef$patches$layout$design <- c(
  area(t = 1, l = 1, b = 30, r = 3),
  area(t = 1, l = 4, b = 30, r = 50),
  area(t = 29, l = 47, b = 31, r = 50),
  area(t = 1, l = 50, b = 30, r = 59),
  area(t = 31, l = 4, b = 32, r = 59)
)

medage1525sigma <- medage1525sigma + labs(title = NULL, x = NULL) + theme(legend.position = "none")

medage1525sigma$patches$layout$design <- c(
  area(t = 1, l = 1, b = 30, r = 3),
  area(t = 1, l = 4, b = 30, r = 50)
)

medage1525reef[[4]] <- medage1525reef[[4]] + labs(title = NULL, x = NULL) + theme(legend.position = "none")

# layout3 <- "
# AAAAAAAA#
# BBBBBBBBB
# BBBBBBBBB
# "

layout4 <- "
AABB
AABB
CCCD
CCCD
CCCD
CCCD
CCCD
"

medagexspace15to25 <- ((medage1525sigma - medage1525loc) + medage1525reef + guide_area()) + plot_layout(guides = "collect", design = layout4) + plot_annotation(tag_levels = list(c("A", rep("", length(medage1525sigma$patches$layout$design$t) - 1),
                                                                                                                                                                                    "B", #rep("", length(medage1525loc$patches$layout$design$t) - 1),
                                                                                                                                                                                    "C", rep("", length(medage1525reef$patches$layout$design$t) - 1))), "") &
  theme(legend.text = element_text(size = 9))
#medagexspace15to25

ggsave(filename = here::here(paste0("Summer2023/StanOutput/medagexspace15to25_", metric_choice, ".pdf")),
       plot = medagexspace15to25,
       width = 6.75,
       height = 6,
       units = "in",
       device = cairo_pdf,
       dpi = 300)

ggsave(filename = here::here(paste0("Summer2023/StanOutput/medagexspace15to25_", metric_choice, ".png")),
       plot = medagexspace15to25,
       width = 6.75,
       height = 6,
       units = "in")


medage2535loc <- readRDS(here::here(paste0("Summer2023/StanOutput/HOBS_medagexspace_25to35cm_locintplot_", metric_choice, ".rds")))
medage2535reef <- readRDS(here::here(paste0("Summer2023/StanOutput/HOBS_medagexspace_25to35cm_reefintplot_", metric_choice, ".rds")))
medage2535sigma <- readRDS(here::here(paste0("Summer2023/StanOutput/HOBS_medagexspace_25to35cm_sigmaintplot_", metric_choice, ".rds")))

# #The commented out code is only necessary if using means, in which case the loc plot needs an x-axis break, so is a patchwork
# medage2535loc$patches$layout$design <- c(
#   area(t = 1, l = 1, b = 30, r = 3),
#   area(t = 1, l = 4, b = 30, r = 50),
#   area(t = 27, l = 46, b = 31, r = 50),
#   area(t = 1, l = 50, b = 30, r = 59)
# )
# medage2535loc[[2]] <- medage2535loc[[2]] +
#   {if(metric_choice == "median"){
#     xlim(c(-15, sampsum_25to35cm[, med_param_loc := median(median_age_bp), by = Locality][, max(med_param_loc)] + 50))
#   } else{
#     xlim(c(-15, sampsum_25to35cm[, mean_param_loc := mean(median_age_bp), by = Locality][, max(mean_param_loc)] + 50))
#   }} +
#   {if(metric_choice == "median"){
#     coord_cartesian(xlim = c(-15, sampsum_25to35cm[, med_param_loc := median(median_age_bp), by = Locality][, max(med_param_loc)] + 50))
#   } else{
#     coord_cartesian(xlim = c(-15, sampsum_25to35cm[, mean_param_loc := mean(median_age_bp), by = Locality][, max(mean_param_loc)] + 50))
#   }} +
#   {if(metric_choice == "median"){
#     scale_x_continuous(breaks = seq(0, #plyr::round_any(min(c(plot1$data$ll, Data$mean_param_loc)), 50, floor), 
#                                     plyr::round_any(max(c(medage2535loc[[2]]$data$ll, sampsum_25to35cm[, med_param_loc := median(median_age_bp), by = Locality][, max(med_param_loc)])), 50, ceiling), by = 50))
#   } else{
#     scale_x_continuous(breaks = seq(0, #plyr::round_any(min(c(plot1$data$ll, Data$mean_param_loc)), 50, floor), 
#                                     plyr::round_any(max(c(medage2535loc[[2]]$data$ll, sampsum_25to35cm[, mean_param_loc := mean(median_age_bp), by = Locality][, max(mean_param_loc)])), 50, ceiling), by = 50))
#   }}
#   
# medage2535loc[[4]] <- medage2535loc[[4]] + theme(legend.box.margin = margin(l = 20))

medage2535loc <- medage2535loc + labs(title = NULL, x = NULL)

# medage2535loc$patches$layout$design <- c(
#   area(t = 1, l = 1, b = 31, r = 3),
#   area(t = 1, l = 4, b = 31, r = 50),
#   area(t = 28, l = 47, b = 31, r = 50),
#   area(t = 1, l = 50, b = 31, r = 59)
# )

medage2535reef$patches$layout$design <- c(
  area(t = 1, l = 1, b = 30, r = 3),
  area(t = 1, l = 4, b = 30, r = 50),
  area(t = 29, l = 47, b = 31, r = 50),
  area(t = 1, l = 50, b = 30, r = 59),
  area(t = 31, l = 4, b = 32, r = 59)
)

medage2535sigma <- medage2535sigma + labs(title = NULL, x = NULL) + theme(legend.position = "none")

medage2535sigma$patches$layout$design <- c(
  area(t = 1, l = 1, b = 30, r = 3),
  area(t = 1, l = 4, b = 30, r = 50)
)

medage2535reef[[4]] <- medage2535reef[[4]] + labs(title = NULL, x = NULL) + theme(legend.position = "none")

layout4 <- "
AABB
AABB
CCCD
CCCD
CCCD
CCCD
CCCD
"

medagexspace25to35 <- ((medage2535sigma - medage2535loc) + medage2535reef + guide_area()) + plot_layout(guides = "collect", design = layout4) + plot_annotation(tag_levels = list(c("A", rep("", length(medage2535sigma$patches$layout$design$t) - 1),
                                                                                                                                                                                    "B", #rep("", length(medage2535loc$patches$layout$design$t) - 1),
                                                                                                                                                                                    "C", rep("", length(medage2535reef$patches$layout$design$t) - 1))), "") &
  theme(legend.text = element_text(size = 9))
#medagexspace25to35

ggsave(filename = here::here(paste0("Summer2023/StanOutput/medagexspace25to35_", metric_choice, ".pdf")),
       plot = medagexspace25to35,
       width = 6.75,
       height = 6,
       units = "in",
       device = cairo_pdf,
       dpi = 300)

ggsave(filename = here::here(paste0("Summer2023/StanOutput/medagexspace25to35_", metric_choice, ".png")),
       plot = medagexspace25to35,
       width = 6.75,
       height = 6,
       units = "in")


medageddloc <- readRDS(here::here(paste0("Summer2023/StanOutput/HOBS_medageddxspace_locintplot_", metric_choice, ".rds")))
medageddreef <- readRDS(here::here(paste0("Summer2023/StanOutput/HOBS_medageddxspace_reefintplot_", metric_choice, ".rds")))
medageddsigma <- readRDS(here::here(paste0("Summer2023/StanOutput/HOBS_medageddxspace_sigmaintplot_", metric_choice, ".rds")))

medageddsigma <- medageddsigma + labs(title = NULL, x = NULL) + theme(legend.position = "none")

medageddloc <- medageddloc + 
  # #xlim(c(min(medageddreef$layers[[8]]$data$mean_param_reef), max(medageddreef$layers[[8]]$data$mean_param_reef))) +
  # {if(metric_choice == "median"){
  #   scale_x_continuous(breaks = seq(plyr::round_any(min(medageddreef$layers[[8]]$data$med_param_reef), 100, ceiling), 
  #                                   plyr::round_any(max(medageddreef$layers[[8]]$data$med_param_reef), 100, floor), by = 100))
  # } else{
  #   scale_x_continuous(breaks = seq(plyr::round_any(min(medageddreef$layers[[8]]$data$mean_param_reef), 100, ceiling), 
  #                                   plyr::round_any(max(medageddreef$layers[[8]]$data$mean_param_reef), 100, floor), by = 100))
  # }} +
  # {if(metric_choice == "median"){
  #   coord_cartesian(xlim = c(min(medageddreef$layers[[8]]$data$med_param_reef), max(medageddreef$layers[[8]]$data$med_param_reef)))
  # } else{
#   coord_cartesian(xlim = c(min(medageddreef$layers[[8]]$data$mean_param_reef), max(medageddreef$layers[[8]]$data$mean_param_reef)))
# }} +
labs(x = NULL, y = "Locality", title = NULL) + 
  theme(legend.text = element_text(size = 9))

medageddreef <- medageddreef +
  {if(metric_choice == "median"){
    scale_x_continuous(breaks = seq(plyr::round_any(min(medageddreef$layers[[9]]$data$med_param_reef), 100, ceiling), 
                                    plyr::round_any(max(medageddreef$layers[[9]]$data$med_param_reef), 100, floor), by = 100))
  } else{
    scale_x_continuous(breaks = seq(plyr::round_any(min(medageddreef$layers[[9]]$data$mean_param_reef), 100, ceiling), 
                                    plyr::round_any(max(medageddreef$layers[[9]]$data$mean_param_reef), 100, floor), by = 100))
  }} +
  labs(x = expression("Median age"["25-35cm"] * " - Median age"["15-25cm"] * " (years before collection year)"), y = "Reef", title = NULL) + 
  theme(legend.text = element_text(size = 9))

# layout4 <- "
# A
# B
# B
# "

medageddxspace <- ((medageddsigma - medageddloc) + medageddreef + guide_area()) + plot_layout(guides = "collect", design = layout4) + plot_annotation(tag_levels = "A") & #list(c("A", rep("", length(medageddsigma$patches$layout$design$t) - 1), 
  # "B",
  # "C", rep("", length(medageddreef$patches$layout$design$t) - 1))), "") &
  theme(legend.text = element_text(size = 9))

ggsave(filename = here::here(paste0("Summer2023/StanOutput/medageddxspace_", metric_choice, ".pdf")),
       plot = medageddxspace,
       width = 6.75,
       height = 6,
       units = "in",
       device = cairo_pdf,
       dpi = 300)

ggsave(filename = here::here(paste0("Summer2023/StanOutput/medageddxspace_", metric_choice, ".png")),
       plot = medageddxspace,
       width = 6.75,
       height = 6,
       units = "in")


#Results table
files <- c("HOBS_medagexspace_15to25cm.rds", "HOBS_medagexspace_25to35cm.rds", "HOBS_medageddxspace.rds", "HOBS_tavxspace_15to25cm.rds", "HOBS_tavxspace_25to35cm.rds", "HOBS_tavddxspace.rds")
fit <- readRDS(here::here(paste0("StanOutput/", files[1])))
results <- fit$summary()
setDT(results)
results <- results[, .(variable, mean, sd, median, q5, q95)]
results[, locality := c("statewide", "statewide", unique(as.character(locs$locality)), paste0(reefs$locality, " R", reefs$reef_lab), "statewide", "statewide", "statewide")]
setcolorder(results, c("locality", "variable"))
setnames(results, names(results)[c(3:7)], paste0("medage1525_", names(results)[c(3:7)]))

for(file in files[2:length(files)]){
  fit <- readRDS(here::here(paste0("StanOutput/", file)))
  results_i <- fit$summary()
  setDT(results_i)
  results_i <- results_i[, .(mean, sd, median, q5, q95)]
  setnames(results_i, names(results_i), paste0(ifelse(str_detect(file, "medage"), 
                                                      ifelse(str_detect(file, "15"), "medage1525_", 
                                                             ifelse(str_detect(file, "35"), "medage2535_", "medagedd")),
                                                      ifelse(str_detect(file, "15"), "tav1525_", 
                                                             ifelse(str_detect(file, "35"), "tav2535_", "tavdd"))), names(results_i)))
  results <- cbind(results, results_i)
}

saveRDS(results, here::here(paste0("Summer2023/StanOutput/ModelResults_", metric_choice, ".rds")))
openxlsx::write.xlsx(results, here::here(paste0("Summer2023/StanOutput/ModelResults_", metric_choice, ".xlsx")))

sigs <- results[, .(locality, variable, medage1525_median, medage1525_q5, medage1525_q95, medage2535_median, medage2535_q5, medage2535_q95, medageddmedian, medageddq5, medageddq95, 
                    tav1525_median, tav1525_q5, tav1525_q95, tav2535_median, tav2535_q5, tav2535_q95, tavddmedian, tavddq5, tavddq95)]

sigs[, `:=` (medage1525_signif = ifelse(medage1525_q5 < 0 & medage1525_q95 < 0, -1, ifelse(medage1525_q5 > 0 & medage1525_q95 > 0, 1, 0)),
             medage2535_signif = ifelse(medage2535_q5 < 0 & medage2535_q95 < 0, -1, ifelse(medage2535_q5 > 0 & medage2535_q95 > 0, 1, 0)),
             medagedd_signif = ifelse(medageddq5 < 0 & medageddq95 < 0, -1, ifelse(medageddq5 > 0 & medageddq95 > 0, 1, 0)),
             tav1525_signif = ifelse(tav1525_q5 < 0 & tav1525_q95 < 0, -1, ifelse(tav1525_q5 > 0 & tav1525_q95 > 0, 1, 0)),
             tav2535_signif = ifelse(tav2535_q5 < 0 & tav2535_q95 < 0, -1, ifelse(tav2535_q5 > 0 & tav2535_q95 > 0, 1, 0)),
             tavdd_signif = ifelse(tavddq5 < 0 & tavddq95 < 0, -1, ifelse(tavddq5 > 0 & tavddq95 > 0, 1, 0))
             )]

#Table S1----------------------------------------------------------------
madat <- read.xlsx(here::here("Summer2023/EstablishedAcreages_Apr2021_foruse.xlsx"), sheet = "RCPManagedAreas_2019")
madat <- janitor::clean_names(madat)
setDT(madat)

madat_hobs <- madat[str_detect(long_name, "Apalachicola|Seagrasses|Estero|Pierce|Guana"), .(long_name, region, date_of_designation, acres_2019)]
madat_hobs[, `:=` (water_body = fcase(str_detect(long_name, "Apalachicola Bay"), "Apalachicola Bay",
                                      str_detect(long_name, "Apalachicola National"), "Apalachicola River,\nApalachicola Bay",
                                      str_detect(long_name, "Seagrasses"), "Gulf of Mexico\n(Apalachee Bay to\nWaccasassa Bay",
                                      str_detect(long_name, "Estero"), "Estero Bay",
                                      str_detect(long_name, "Pierce"), "Indian River Lagoon",
                                      str_detect(long_name, "Tolomato"), "Guana River, Tolomato\nRiver, Matanzas River,\nPellicer Creek",
                                      str_detect(long_name, "Marsh"), "Guana River, Tolomato\nRiver, Atlantic Ocean\n(Ponte Vedra Beach,\nfrom Sawgrass to\n~1.5km north of Vilano\nBeach)"),
                   study_loc = fcase(str_detect(long_name, "Apalachicola Bay"), "Little St. George Is.",
                                     str_detect(long_name, "Apalachicola National"), "Little St. George Is.,\nGoose Is./East Cove",
                                     str_detect(long_name, "Seagrasses"), "Lone Cabbage",
                                     str_detect(long_name, "Estero"), "Hendry/Mullock Creeks,\nNew Pass,\nBig Hickory",
                                     str_detect(long_name, "Pierce"), "Jack Island",
                                     str_detect(long_name, "Tolomato"), "Guana River,\nMatanzas River,\nPellicer Creek",
                                     str_detect(long_name, "Marsh"), "Guana River"))]

setcolorder(madat_hobs, c("long_name", "water_body", "region", "date_of_designation", "acres_2019", "study_loc"))
madat_hobs[, `:=` (long_name = factor(long_name, levels = c("Apalachicola Bay Aquatic Preserve",
                                                            "Apalachicola National Estuarine Research Reserve",
                                                            "Big Bend Seagrasses Aquatic Preserve",
                                                            "Estero Bay Aquatic Preserve",
                                                            "Indian River-Vero Beach to Ft. Pierce Aquatic Preserve",
                                                            "Guana Tolomato Matanzas National Estuarine Research Reserve",
                                                            "Guana River Marsh Aquatic Preserve")),
                   acres_2019 = round(acres_2019))]
setorder(madat_hobs, long_name)

#format and export table in Excel
wb <- createWorkbook()
addWorksheet(wb, "TableS1")
s1_caption <- createStyle(halign = "center", valign = "top", border = "Bottom", borderStyle = "double", fontName = "Arial", fontSize = 8)
s1_headers_nc1 <- createStyle(halign = "center", valign = "bottom", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 8)
s1_headers_c1 <- createStyle(halign = "left", valign = "bottom", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 8)
s1_body_left <- createStyle(halign = "left", valign = "top", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s1_body_right <- createStyle(halign = "right", valign = "top", numFmt = "#,##0", fontName = "Arial", fontSize = 8)
s1_bottom_left <- createStyle(halign = "left", valign = "top", border = "Bottom", borderStyle = "thin", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s1_bottom_right <- createStyle(halign = "right", valign = "top", numFmt = "#,##0", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 8)
# class(madat_hobs$acres_2019) <- "comma"

title <- "TABLE S1. SAMPLING LOCALITIES AND ORCP MANAGED AREAS"
colnms <- data.table("ORCP Managed Area", "Water Body", "Region", "Date of\nDesignation", "Area\n(acres)", "Study Localities")
writeData(wb, "TableS1", title, startRow = 2, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS1", cols = 2:7, rows = 2)
addStyle(wb, "TableS1", style = s1_caption, rows = 2, cols = 2:7, gridExpand = TRUE)

writeData(wb, "TableS1", colnms, startRow = 3, startCol = 2, colNames = FALSE)
addStyle(wb, "TableS1", style = s1_headers_nc1, rows = 3, cols = 3:7, gridExpand = TRUE)
addStyle(wb, "TableS1", style = s1_headers_c1, rows = 3, cols = 2, gridExpand = TRUE)

writeData(wb, "TableS1", madat_hobs, startRow = 4, startCol = 2, colNames = FALSE)
addStyle(wb, "TableS1", style = s1_body_left, rows = 4:9, cols = c(2:4, 7), gridExpand = TRUE)
addStyle(wb, "TableS1", style = s1_body_right, rows = 4:9, cols = 5:6, gridExpand = TRUE)
addStyle(wb, "TableS1", style = s1_bottom_left, rows = 10, cols = c(2:4, 7), gridExpand = TRUE)
addStyle(wb, "TableS1", style = s1_bottom_right, rows = 10, cols = 5:6, gridExpand = TRUE)

setColWidths(wb, "TableS1", cols = 2:7, widths = c(30, 23, 9, 12, 9, 23))
setRowHeights(wb, "TableS1", rows = 2:10, heights = c(23, 32, 32, 32, 48, 48, 32, 48, 85))


#Tables S2 and S3 prep----------------------------------------------
tp <- fread(here::here("Summer2023/WC_Discrete_TP_Lab_All_MA_Yr_Stats.txt"), sep = "|")
tn <- fread(here::here("Summer2023/WC_Discrete_TN_Lab_All_MA_Yr_Stats.txt"), sep = "|")
temp <- fread(here::here("Summer2023/WC_Discrete_TempW_Field_All_MA_Yr_Stats.txt"), sep = "|")
turb <- fread(here::here("Summer2023/WC_Discrete_Turb_All_All_MA_Yr_Stats.txt"), sep = "|")
sal <- fread(here::here("Summer2023/WC_Discrete_Sal_All_All_MA_Yr_Stats.txt"), sep = "|")
do <- fread(here::here("Summer2023/WC_Discrete_DO_Field_All_MA_Yr_Stats.txt"), sep = "|")

hobsmas <- c("Apalachicola Bay Aquatic Preserve", 
             "Apalachicola National Estuarine Research Reserve", 
             "Big Bend Seagrasses Aquatic Preserve", 
             "Estero Bay Aquatic Preserve", 
             "Indian River-Vero Beach to Ft. Pierce Aquatic Preserve", 
             "Guana Tolomato Matanzas National Estuarine Research Reserve", 
             "Guana River Marsh Aquatic Preserve")

ts2dat <- dplyr::bind_rows(tp[ManagedAreaName %in% hobsmas, ], 
                           tn[ManagedAreaName %in% hobsmas, ],
                           temp[ManagedAreaName %in% hobsmas, ],
                           turb[ManagedAreaName %in% hobsmas, ],
                           sal[ManagedAreaName %in% hobsmas, ],
                           do[ManagedAreaName %in% hobsmas, ])

commonyears <- data.table(ManagedAreaName = character(),
                          years_c = list())

params <- c("Total Phosphorus", "Total Nitrogen", "Water Temperature", "Salinity", "Turbidity", "Dissolved Oxygen")

for(m in unique(ts2dat$ManagedAreaName)){
  years_m <- unique(ts2dat$Year)
  
  for(p in params){
    years_m <- intersect(years_m, ts2dat[ManagedAreaName == m & ParameterName == p, unique(Year)])
    
  }
  commonyears_m <- data.table(ManagedAreaName = m,
                              years_c = list(years_m))
  
  commonyears <- rbind(commonyears, commonyears_m)
  
}

ts2dat_c <- ts2dat[0]

for(t in unique(ts2dat$ManagedAreaName)){
  cyrs <- commonyears[ManagedAreaName == t, years_c[[1]]]
  
  ts2dat_t <- ts2dat[ManagedAreaName == t & ParameterName %in% params & Year %in% cyrs, ]
  
  ts2dat_c <- rbind(ts2dat_c, ts2dat_t)
  
}

ts2dat_c[, `:=` (med_n = median(N_Data, na.rm = TRUE), iqr_nl = quantile(N_Data, probs = 0.25, na.rm = TRUE), iqr_nh = quantile(N_Data, probs = 0.75, na.rm = TRUE),
                 med_min = median(Min, na.rm = TRUE), iqr_minl = quantile(Min, probs = 0.25, na.rm = TRUE), iqr_minh = quantile(Min, probs = 0.75, na.rm = TRUE),
                 med_max = median(Max, na.rm = TRUE), iqr_maxl = quantile(Max, probs = 0.25, na.rm = TRUE), iqr_maxh = quantile(Max, probs = 0.75, na.rm = TRUE),
                 med_med = median(Median, na.rm = TRUE), iqr_medl = quantile(Median, probs = 0.25, na.rm = TRUE), iqr_medh = quantile(Median, probs = 0.75, na.rm = TRUE),
                 med_mn = median(Mean, na.rm = TRUE), iqr_mnl = quantile(Mean, probs = 0.25, na.rm = TRUE), iqr_mnh = quantile(Mean, probs = 0.75, na.rm = TRUE),
                 med_sd = median(StandardDeviation, na.rm = TRUE), iqr_sdl = quantile(StandardDeviation, probs = 0.25, na.rm = TRUE), iqr_sdh = quantile(StandardDeviation, probs = 0.75, na.rm = TRUE),
                 Years = list(c(Year))), by = c("ManagedAreaName", "ParameterName")]

progs <- ts2dat_c[, unlist(strsplit(ProgramIDs, split = ", ")), by = row.names(ts2dat_c)]
progs <- sort(unique(as.integer(progs$V1)))

ts2 <- distinct(ts2dat_c[, .(ManagedAreaName, ParameterName, med_n, iqr_nl, iqr_nh, med_min, iqr_minl, iqr_minh, med_max, iqr_maxl, iqr_maxh, med_med, iqr_medl, iqr_medh, med_mn, iqr_mnl, iqr_mnh, med_sd, iqr_sdl, iqr_sdh, Years)])

# ts2[, `:=` (med_n = round(med_n, 0), iqr_nl = round(med_n, 0), iqr_nh = round(med_n, 0),
#             med_min = round(med_min, 5), iqr_minl = round(iqr_minl, 5), iqr_minh = round(iqr_minh, 5), 
#             med_max = round(med_max, 5), iqr_maxl = round(iqr_maxl, 5), iqr_maxh = round(iqr_maxh, 5), 
#             med_med = round(med_med, 3), iqr_medl = round(iqr_medl, 3), iqr_medh = round(iqr_medh, 3), 
#             med_mn = round(med_mn, 3), iqr_mnl = round(iqr_mnl, 3), iqr_mnh = round(iqr_mnh, 3), 
#             med_sd = round(med_sd, 3), iqr_sdl = round(iqr_sdl, 3), iqr_sdh = round(iqr_sdh, 3))]

ts2a <- melt(ts2, measure.vars = c("med_n", "iqr_nl", "iqr_nh", "med_min", "iqr_minl", "iqr_minh", "med_max", "iqr_maxl", "iqr_maxh", "med_med", "iqr_medl", "iqr_medh", "med_mn", "iqr_mnl", "iqr_mnh", "med_sd", "iqr_sdl", "iqr_sdh"), variable.name = "Statistic (annual)", value.name = "Time series value")

setDT(ts2a)
ts2a[, quantile := ifelse(str_detect(`Statistic (annual)`, "med_"), "50%",
                          ifelse(str_detect(`Statistic (annual)`, "iqr_"), #"25%", "75%")), by = `Statistic (annual)`]
                                 ifelse(str_detect(`Statistic (annual)`, "l"), "25%", "75%"), NA)), by = `Statistic (annual)`]

ts2a[, `Statistic (annual)` := ifelse(str_detect(`Statistic (annual)`, "med_"), 
                                      str_to_title(str_sub(`Statistic (annual)`, 5, -1)), 
                                      str_to_title(str_sub(`Statistic (annual)`, 5, -2))), by = `Statistic (annual)`][, `Statistic (annual)` := fcase(`Statistic (annual)` == "N", "N",
                                                                                                                                                      `Statistic (annual)` == "Min", paste0(`Statistic (annual)`, "imum"),
                                                                                                                                                      `Statistic (annual)` == "Max", paste0(`Statistic (annual)`, "imum"),
                                                                                                                                                      `Statistic (annual)` == "Med", "Median",
                                                                                                                                                      `Statistic (annual)` == "Mn", "Mean",
                                                                                                                                                      `Statistic (annual)` == "Sd", "Standard deviation"), by = `Statistic (annual)`]

ts2b <- pivot_wider(ts2a, id_cols = c("ManagedAreaName", "Years", "Statistic (annual)", "quantile"), names_from = c("ParameterName"), values_from = c("Time series value"))

ts2c <- pivot_wider(ts2b, id_cols = c("ManagedAreaName", "Years", "Statistic (annual)"), names_from = c("quantile"), values_from = c("Total Phosphorus", "Total Nitrogen", "Water Temperature", "Turbidity", "Salinity", "Dissolved Oxygen"))

setorder(ts2c, ManagedAreaName, `Statistic (annual)`)
setcolorder(ts2c, neworder = c("ManagedAreaName", "Years", "Statistic (annual)", "Water Temperature_25%", "Water Temperature_50%", "Water Temperature_75%", "Salinity_25%", "Salinity_50%", "Salinity_75%", "Dissolved Oxygen_25%", "Dissolved Oxygen_50%", "Dissolved Oxygen_75%", "Turbidity_25%", "Turbidity_50%", "Turbidity_75%", "Total Nitrogen_25%", "Total Nitrogen_50%", "Total Nitrogen_75%", "Total Phosphorus_25%", "Total Phosphorus_50%", "Total Phosphorus_75%"))

fwrite(ts2c, here::here("Summer2023/ts2c.csv"))

setDT(ts2c)
ts2c[, `Statistic (annual)` := fcase(`Statistic (annual)` == "N", "N",
                                     `Statistic (annual)` == "Minimum", "Min.",
                                     `Statistic (annual)` == "Maximum", "Max.",
                                     `Statistic (annual)` == "Median", "Median",
                                     `Statistic (annual)` == "Mean", "Mean",
                                     `Statistic (annual)` == "Standard deviation", "St. dev.")]
ts2c[, `:=` (spacer1 = "", spacer2 = "", spacer3 = "", spacer4 = "", spacer5 = "", spacer6 = "")]
setcolorder(ts2c, c("ManagedAreaName", "Years", "Statistic (annual)", 
                    "spacer1", "Water Temperature_25%", "Water Temperature_50%", "Water Temperature_75%", 
                    "spacer2", "Salinity_25%", "Salinity_50%", "Salinity_75%", 
                    "spacer3", "Dissolved Oxygen_25%", "Dissolved Oxygen_50%", "Dissolved Oxygen_75%", 
                    "spacer4", "Turbidity_25%", "Turbidity_50%", "Turbidity_75%", 
                    "spacer5", "Total Nitrogen_25%", "Total Nitrogen_50%", "Total Nitrogen_75%", 
                    "spacer6", "Total Phosphorus_25%", "Total Phosphorus_50%", "Total Phosphorus_75%"))


#Table S2----------------------------------------------------------------
# wb <- createWorkbook()
addWorksheet(wb, "TableS2")
s2_caption <- createStyle(halign = "center", valign = "top", border = "Bottom", borderStyle = "double", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s2_headers_nc1 <- createStyle(halign = "center", valign = "bottom", border = "Bottom", borderStyle = "thin", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s2_headers_c1_1 <- createStyle(halign = "left", valign = "top", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s2_headers_c2_1 <- createStyle(halign = "center", valign = "top", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s2_headers_c1_2 <- createStyle(halign = "left", valign = "bottom", wrapText = TRUE, border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 8)
s2_body_left <- createStyle(halign = "left", valign = "top", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s2_body_right_num <- createStyle(halign = "right", valign = "top", numFmt = "#,##0.00", fontName = "Arial", fontSize = 8)
s2_body_right_int <- createStyle(halign = "right", valign = "top", numFmt = "#,##0", fontName = "Arial", fontSize = 8)
s2_bottom_left <- createStyle(halign = "left", valign = "top", border = "Bottom", borderStyle = "thin", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s2_bottom_right <- createStyle(halign = "right", valign = "top", numFmt = "#,##0.00", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 8)

s2_title_p1 <- "TABLE S2. QUANTILES (25%, 50% AND 75%) OF ANNUAL SUMMARY STATISTICS FOR WATER TEMPERATURE, SALINITY AND DISSOLVED OXYGEN FROM THE SAMPLED ORCP MANAGED AREAS (LISTED IN COUNTER-CLOCKWISE ORDER AROUND FLORIDA, STARTING FROM THE PANHANDLE)"
s2_title_p2 <- "TABLE S2. (CONTINUED)"
s2_colnms1 <- data.table("Managed Area", "Localities", "Years", "Annual Statistic", "", "Water Temperature (\u00B0C)", "", "", "", "Salinity (ppt)", "", "", "", "Dissolved Oxygen (mg/L)", "", "")
s2_colnms2 <- data.table("", "", "", "", "", "25 %", "50 %", "75 %", "", "25 %", "50 %", "75 %", "", "25 %", "50 %", "75%")
s2_foots <- c("AP = Aquatic Preserve; NERR = National Estuarine Research Reserve", "Data source: SEACAR database (https://data.florida-seacar.org/)", paste0("Program IDs: ", capture.output(cat(unlist(ts2dat_c[ParameterName %in% c("Water Temperature", "Salinity", "Dissolved Oxygen"), sort(unique(as.integer(unlist(strsplit(ProgramIDs, split = ", ")), by = row.names(ts2dat_c))))]), sep = ", "))))
s2_mas <- c("Apalachicola Bay AP", "", "", "", "", "", "Apalachicola NERR", "", "", "", "", "", "Big Bend Seagrasses AP", "", "", "", "", "", "Estero Bay AP", "", "", "", "", "", "Indian River-Vero Beach to Ft. Pierce AP", "", "", "", "", "", "Guana Tolomato Matanzas NERR", "", "", "", "", "", "Guana River Marsh AP", "", "", "", "", "")
s2_localities <- c("Little St. George Island", "", "", "", "", "", "Little St. George Island, Goose Island/East Cove", "", "", "", "", "", "Lone Cabbage", "", "", "", "", "", "Hendry/Mullock Creeks, New Pass, Big Hickory", "", "", "", "", "", "Jack Island", "", "", "", "", "", "Pellicer Creek, Matanzas River, Guana River", "", "", "", "", "", "Guana River", "", "", "", "", "")
s2_years <- c("1992, 2000-2008, 2010, 2012-2014, 2018-2022", "", "", "", "", "", "1992, 1998-2023", "", "", "", "", "", "1992-2023", "", "", "", "", "", "1999-2023", "", "", "", "", "", "1997-2023", "", "", "", "", "", "1997-2023", "", "", "", "", "", "1997-2023", "", "", "", "", "")

writeData(wb, "TableS2", s2_title_p1, startRow = 2, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS2", cols = 2:17, rows = 2)
addStyle(wb, "TableS2", style = s2_caption, rows = 2, cols = 2:17, gridExpand = TRUE)

writeData(wb, "TableS2", s2_colnms1, startRow = 3, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS2", cols = 7:9, rows = 3)
mergeCells(wb, "TableS2", cols = 11:13, rows = 3)
mergeCells(wb, "TableS2", cols = 15:17, rows = 3)
addStyle(wb, "TableS2", style = s2_headers_nc1, rows = 3, cols = c(7:9, 11:13, 15:17), gridExpand = TRUE)
addStyle(wb, "TableS2", style = s2_headers_c1_1, rows = 3, cols = 2, gridExpand = TRUE)
addStyle(wb, "TableS2", style = s2_headers_c2_1, rows = 3, cols = 3:5, gridExpand = TRUE)

writeData(wb, "TableS2", s2_colnms2, startRow = 4, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS2", cols = 2, rows = 3:4)
mergeCells(wb, "TableS2", cols = 3, rows = 3:4)
mergeCells(wb, "TableS2", cols = 4, rows = 3:4)
mergeCells(wb, "TableS2", cols = 5, rows = 3:4)
addStyle(wb, "TableS2", style = s2_headers_nc1, rows = 4, cols = 3:17, gridExpand = TRUE)
addStyle(wb, "TableS2", style = s2_headers_c1_2, rows = 4, cols = 2, gridExpand = TRUE)

writeData(wb, "TableS2", ts2c[1:24, 3:15], startRow = 5, startCol = 5, colNames = FALSE)
writeData(wb, "TableS2", s2_mas[1:24], startRow = 5, startCol = 2, colNames = FALSE)
writeData(wb, "TableS2", s2_localities[1:24], startRow = 5, startCol = 3, colNames = FALSE)
writeData(wb, "TableS2", s2_years[1:24], startRow = 5, startCol = 4, colNames = FALSE)
mergeCells(wb, "TableS2", cols = 2, rows = 5:10)
mergeCells(wb, "TableS2", cols = 2, rows = 11:16)
mergeCells(wb, "TableS2", cols = 2, rows = 17:22)
mergeCells(wb, "TableS2", cols = 2, rows = 23:28)
mergeCells(wb, "TableS2", cols = 3, rows = 5:10)
mergeCells(wb, "TableS2", cols = 3, rows = 11:16)
mergeCells(wb, "TableS2", cols = 3, rows = 17:22)
mergeCells(wb, "TableS2", cols = 3, rows = 23:28)
mergeCells(wb, "TableS2", cols = 4, rows = 5:10)
mergeCells(wb, "TableS2", cols = 4, rows = 11:16)
mergeCells(wb, "TableS2", cols = 4, rows = 17:22)
mergeCells(wb, "TableS2", cols = 4, rows = 23:28)
addStyle(wb, "TableS2", style = s2_body_left, rows = 5:28, cols = 2:5, gridExpand = TRUE)
addStyle(wb, "TableS2", style = s2_body_right_int, rows = c(5, 11, 17, 23), cols = 6:17, gridExpand = TRUE)
addStyle(wb, "TableS2", style = s2_body_right_num, rows = c(6:10, 12:16, 18:22, 24:27), cols = 6:17, gridExpand = TRUE)
addStyle(wb, "TableS2", style = s2_bottom_left, rows = 28, cols = 2:5, gridExpand = TRUE)
addStyle(wb, "TableS2", style = s2_bottom_right, rows = 28, cols = 6:17, gridExpand = TRUE)

setColWidths(wb, "TableS2", cols = 2:17, widths = c(13, 13, 10, 7, 1, 6, 6, 6, 1, 6, 6, 6, 1, 6, 6, 6))
setRowHeights(wb, "TableS2", rows = 2:28, heights = c(30, 12, 12, 12, 12, 12, 12, 12, 22, 12, 12, 12, 12, 12, 22, 12, 12, 12, 12, 12, 22, 12, 12, 12, 12, 12, 14))


writeData(wb, "TableS2", s2_title_p2, startRow = 32, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS2", cols = 2:17, rows = 32)
addStyle(wb, "TableS2", style = s2_caption, rows = 32, cols = 2:17, gridExpand = TRUE)

writeData(wb, "TableS2", s2_colnms1, startRow = 33, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS2", cols = 7:9, rows = 33)
mergeCells(wb, "TableS2", cols = 11:13, rows = 33)
mergeCells(wb, "TableS2", cols = 15:17, rows = 33)
addStyle(wb, "TableS2", style = s2_headers_nc1, rows = 33, cols = c(7:9, 11:13, 15:17), gridExpand = TRUE)
addStyle(wb, "TableS2", style = s2_headers_c1_1, rows = 33, cols = 2, gridExpand = TRUE)
addStyle(wb, "TableS2", style = s2_headers_c2_1, rows = 33, cols = 3:5, gridExpand = TRUE)

writeData(wb, "TableS2", s2_colnms2, startRow = 34, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS2", cols = 2, rows = 33:34)
mergeCells(wb, "TableS2", cols = 3, rows = 33:34)
mergeCells(wb, "TableS2", cols = 4, rows = 33:34)
mergeCells(wb, "TableS2", cols = 5, rows = 33:34)
addStyle(wb, "TableS2", style = s2_headers_nc1, rows = 34, cols = 3:17, gridExpand = TRUE)
addStyle(wb, "TableS2", style = s2_headers_c1_2, rows = 34, cols = 2, gridExpand = TRUE)

writeData(wb, "TableS2", ts2c[25:42, 3:15], startRow = 35, startCol = 5, colNames = FALSE)
writeData(wb, "TableS2", s2_mas[25:42], startRow = 35, startCol = 2, colNames = FALSE)
writeData(wb, "TableS2", s2_localities[25:42], startRow = 35, startCol = 3, colNames = FALSE)
writeData(wb, "TableS2", s2_years[25:42], startRow = 35, startCol = 4, colNames = FALSE)
mergeCells(wb, "TableS2", cols = 2, rows = 35:40)
mergeCells(wb, "TableS2", cols = 2, rows = 41:46)
mergeCells(wb, "TableS2", cols = 2, rows = 47:52)
mergeCells(wb, "TableS2", cols = 3, rows = 35:40)
mergeCells(wb, "TableS2", cols = 3, rows = 41:46)
mergeCells(wb, "TableS2", cols = 3, rows = 47:52)
mergeCells(wb, "TableS2", cols = 4, rows = 35:40)
mergeCells(wb, "TableS2", cols = 4, rows = 41:46)
mergeCells(wb, "TableS2", cols = 4, rows = 47:52)
addStyle(wb, "TableS2", style = s2_body_left, rows = 35:51, cols = 2:5, gridExpand = TRUE)
addStyle(wb, "TableS2", style = s2_body_right_int, rows = c(35, 41, 47), cols = 6:17, gridExpand = TRUE)
addStyle(wb, "TableS2", style = s2_body_right_num, rows = c(36:40, 42:46, 48:51), cols = 6:17, gridExpand = TRUE)
addStyle(wb, "TableS2", style = s2_bottom_left, rows = 52, cols = 2:5, gridExpand = TRUE)
addStyle(wb, "TableS2", style = s2_bottom_right, rows = 52, cols = 6:17, gridExpand = TRUE)

writeData(wb, "TableS2", s2_foots, startRow = 53, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS2", cols = 2:17, rows = 53)
mergeCells(wb, "TableS2", cols = 2:17, rows = 54)
mergeCells(wb, "TableS2", cols = 2:17, rows = 55)
addStyle(wb, "TableS2", style = s2_body_left, rows = 53:54, cols = 2, gridExpand = TRUE)
addStyle(wb, "TableS2", style = s2_bottom_left, rows = 55, cols = 2:17, gridExpand = TRUE)

# setColWidths(wb, "TableS2", cols = 2:17, widths = c(10, 10, 10, 7, 1, 5, 5, 5, 1, 5, 5, 5, 1, 5, 5, 5))
setRowHeights(wb, "TableS2", rows = 32:55, heights = c(15, 12, 12, 12, 12, 12, 12, 12, 22, 12, 12, 12, 12, 12, 22, 12, 12, 12, 12, 12, 22, 12, 12, 14))


#Table S3----------------------------------------------------------------
# wb <- createWorkbook()
addWorksheet(wb, "TableS3")
s3_caption <- createStyle(halign = "center", valign = "top", border = "Bottom", borderStyle = "double", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s3_headers_nc1 <- createStyle(halign = "center", valign = "bottom", border = "Bottom", borderStyle = "thin", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s3_headers_c1_1 <- createStyle(halign = "left", valign = "top", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s3_headers_c2_1 <- createStyle(halign = "center", valign = "top", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s3_headers_c1_2 <- createStyle(halign = "left", valign = "bottom", wrapText = TRUE, border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 8)
s3_body_left <- createStyle(halign = "left", valign = "top", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s3_body_right_num <- createStyle(halign = "right", valign = "top", numFmt = "#,##0.00", fontName = "Arial", fontSize = 8)
s3_body_right_int <- createStyle(halign = "right", valign = "top", numFmt = "#,##0", fontName = "Arial", fontSize = 8)
s3_bottom_left <- createStyle(halign = "left", valign = "top", border = "Bottom", borderStyle = "thin", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s3_bottom_right <- createStyle(halign = "right", valign = "top", numFmt = "#,##0.00", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 8)

s3_title_p1 <- "TABLE S3. QUANTILES (25%, 50% AND 75%) OF ANNUAL SUMMARY STATISTICS FOR TURBIDITY, TOTAL NITROGEN AND TOTAL PHOSPHORUS FROM THE SAMPLED ORCP MANAGED AREAS (LISTED IN COUNTER-CLOCKWISE ORDER AROUND FLORIDA, STARTING FROM THE PANHANDLE)"
s3_title_p2 <- "TABLE S3. (CONTINUED)"
s3_colnms1 <- data.table("Managed Area", "Localities", "Years", "Annual Statistic", "", "Turbidity (NTU)", "", "", "", "Total Nitrogen (mg/L)", "", "", "", "Total Phosphorus (mg/L)", "", "")
s3_colnms2 <- data.table("", "", "", "", "", "25 %", "50 %", "75 %", "", "25 %", "50 %", "75 %", "", "25 %", "50 %", "75%")
s3_foots <- c("AP = Aquatic Preserve; NERR = National Estuarine Research Reserve", "Data source: SEACAR database (https://data.florida-seacar.org/)", paste0("Program IDs: ", capture.output(cat(unlist(ts2dat_c[ParameterName %in% c("Turbidity", "Total Nitrogen", "Total Phosphorus"), sort(unique(as.integer(unlist(strsplit(ProgramIDs, split = ", ")), by = row.names(ts2dat_c))))]), sep = ", "))))
s3_mas <- c("Apalachicola Bay AP", "", "", "", "", "", "Apalachicola NERR", "", "", "", "", "", "Big Bend Seagrasses AP", "", "", "", "", "", "Estero Bay AP", "", "", "", "", "", "Indian River-Vero Beach to Ft. Pierce AP", "", "", "", "", "", "Guana Tolomato Matanzas NERR", "", "", "", "", "", "Guana River Marsh AP", "", "", "", "", "")
s3_localities <- c("Little St. George Island", "", "", "", "", "", "Little St. George Island, Goose Island/East Cove", "", "", "", "", "", "Lone Cabbage", "", "", "", "", "", "Hendry/Mullock Creeks, New Pass, Big Hickory", "", "", "", "", "", "Jack Island", "", "", "", "", "", "Pellicer Creek, Matanzas River, Guana River", "", "", "", "", "", "Guana River", "", "", "", "", "")
s3_years <- c("1992, 2000-2008, 2010, 2012-2014, 2018-2022", "", "", "", "", "", "1992, 1998-2023", "", "", "", "", "", "1992-2023", "", "", "", "", "", "1999-2023", "", "", "", "", "", "1997-2023", "", "", "", "", "", "1997-2023", "", "", "", "", "", "1997-2023", "", "", "", "", "")

writeData(wb, "TableS3", s3_title_p1, startRow = 2, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS3", cols = 2:17, rows = 2)
addStyle(wb, "TableS3", style = s3_caption, rows = 2, cols = 2:17, gridExpand = TRUE)

writeData(wb, "TableS3", s3_colnms1, startRow = 3, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS3", cols = 7:9, rows = 3)
mergeCells(wb, "TableS3", cols = 11:13, rows = 3)
mergeCells(wb, "TableS3", cols = 15:17, rows = 3)
addStyle(wb, "TableS3", style = s3_headers_nc1, rows = 3, cols = c(7:9, 11:13, 15:17), gridExpand = TRUE)
addStyle(wb, "TableS3", style = s3_headers_c1_1, rows = 3, cols = 2, gridExpand = TRUE)
addStyle(wb, "TableS3", style = s3_headers_c2_1, rows = 3, cols = 3:5, gridExpand = TRUE)

writeData(wb, "TableS3", s3_colnms2, startRow = 4, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS3", cols = 2, rows = 3:4)
mergeCells(wb, "TableS3", cols = 3, rows = 3:4)
mergeCells(wb, "TableS3", cols = 4, rows = 3:4)
mergeCells(wb, "TableS3", cols = 5, rows = 3:4)
addStyle(wb, "TableS3", style = s3_headers_nc1, rows = 4, cols = 3:17, gridExpand = TRUE)
addStyle(wb, "TableS3", style = s3_headers_c1_2, rows = 4, cols = 2, gridExpand = TRUE)

writeData(wb, "TableS3", ts2c[1:24, c(3, 16:27)], startRow = 5, startCol = 5, colNames = FALSE)
writeData(wb, "TableS3", s3_mas[1:24], startRow = 5, startCol = 2, colNames = FALSE)
writeData(wb, "TableS3", s3_localities[1:24], startRow = 5, startCol = 3, colNames = FALSE)
writeData(wb, "TableS3", s3_years[1:24], startRow = 5, startCol = 4, colNames = FALSE)
mergeCells(wb, "TableS3", cols = 2, rows = 5:10)
mergeCells(wb, "TableS3", cols = 2, rows = 11:16)
mergeCells(wb, "TableS3", cols = 2, rows = 17:22)
mergeCells(wb, "TableS3", cols = 2, rows = 23:28)
mergeCells(wb, "TableS3", cols = 3, rows = 5:10)
mergeCells(wb, "TableS3", cols = 3, rows = 11:16)
mergeCells(wb, "TableS3", cols = 3, rows = 17:22)
mergeCells(wb, "TableS3", cols = 3, rows = 23:28)
mergeCells(wb, "TableS3", cols = 4, rows = 5:10)
mergeCells(wb, "TableS3", cols = 4, rows = 11:16)
mergeCells(wb, "TableS3", cols = 4, rows = 17:22)
mergeCells(wb, "TableS3", cols = 4, rows = 23:28)
addStyle(wb, "TableS3", style = s3_body_left, rows = 5:28, cols = 2:5, gridExpand = TRUE)
addStyle(wb, "TableS3", style = s3_body_right_int, rows = c(5, 11, 17, 23), cols = 6:17, gridExpand = TRUE)
addStyle(wb, "TableS3", style = s3_body_right_num, rows = c(6:10, 12:16, 18:22, 24:27), cols = 6:17, gridExpand = TRUE)
addStyle(wb, "TableS3", style = s3_bottom_left, rows = 28, cols = 2:5, gridExpand = TRUE)
addStyle(wb, "TableS3", style = s3_bottom_right, rows = 28, cols = 6:17, gridExpand = TRUE)

setColWidths(wb, "TableS3", cols = 2:17, widths = c(13, 13, 10, 7, 1, 6, 6, 6, 1, 6, 6, 6, 1, 6, 6, 6))
setRowHeights(wb, "TableS3", rows = 2:28, heights = c(30, 12, 12, 12, 12, 12, 12, 12, 22, 12, 12, 12, 12, 12, 22, 12, 12, 12, 12, 12, 22, 12, 12, 12, 12, 12, 14))

writeData(wb, "TableS3", s3_title_p2, startRow = 32, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS3", cols = 2:17, rows = 32)
addStyle(wb, "TableS3", style = s3_caption, rows = 32, cols = 2:17, gridExpand = TRUE)

writeData(wb, "TableS3", s3_colnms1, startRow = 33, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS3", cols = 7:9, rows = 33)
mergeCells(wb, "TableS3", cols = 11:13, rows = 33)
mergeCells(wb, "TableS3", cols = 15:17, rows = 33)
addStyle(wb, "TableS3", style = s3_headers_nc1, rows = 33, cols = c(7:9, 11:13, 15:17), gridExpand = TRUE)
addStyle(wb, "TableS3", style = s3_headers_c1_1, rows = 33, cols = 2, gridExpand = TRUE)
addStyle(wb, "TableS3", style = s3_headers_c2_1, rows = 33, cols = 3:5, gridExpand = TRUE)

writeData(wb, "TableS3", s3_colnms2, startRow = 34, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS3", cols = 2, rows = 33:34)
mergeCells(wb, "TableS3", cols = 3, rows = 33:34)
mergeCells(wb, "TableS3", cols = 4, rows = 33:34)
mergeCells(wb, "TableS3", cols = 5, rows = 33:34)
addStyle(wb, "TableS3", style = s3_headers_nc1, rows = 34, cols = 3:17, gridExpand = TRUE)
addStyle(wb, "TableS3", style = s3_headers_c1_2, rows = 34, cols = 2, gridExpand = TRUE)

writeData(wb, "TableS3", ts2c[25:42, c(3, 16:27)], startRow = 35, startCol = 5, colNames = FALSE)
writeData(wb, "TableS3", s3_mas[25:42], startRow = 35, startCol = 2, colNames = FALSE)
writeData(wb, "TableS3", s3_localities[25:42], startRow = 35, startCol = 3, colNames = FALSE)
writeData(wb, "TableS3", s3_years[25:42], startRow = 35, startCol = 4, colNames = FALSE)
mergeCells(wb, "TableS3", cols = 2, rows = 35:40)
mergeCells(wb, "TableS3", cols = 2, rows = 41:46)
mergeCells(wb, "TableS3", cols = 2, rows = 47:52)
mergeCells(wb, "TableS3", cols = 3, rows = 35:40)
mergeCells(wb, "TableS3", cols = 3, rows = 41:46)
mergeCells(wb, "TableS3", cols = 3, rows = 47:52)
mergeCells(wb, "TableS3", cols = 4, rows = 35:40)
mergeCells(wb, "TableS3", cols = 4, rows = 41:46)
mergeCells(wb, "TableS3", cols = 4, rows = 47:52)
addStyle(wb, "TableS3", style = s3_body_left, rows = 35:51, cols = 2:5, gridExpand = TRUE)
addStyle(wb, "TableS3", style = s3_body_right_int, rows = c(35, 41, 47), cols = 6:17, gridExpand = TRUE)
addStyle(wb, "TableS3", style = s3_body_right_num, rows = c(36:40, 42:46, 48:51), cols = 6:17, gridExpand = TRUE)
addStyle(wb, "TableS3", style = s3_bottom_left, rows = 52, cols = 2:5, gridExpand = TRUE)
addStyle(wb, "TableS3", style = s3_bottom_right, rows = 52, cols = 6:17, gridExpand = TRUE)

writeData(wb, "TableS3", s3_foots, startRow = 53, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS3", cols = 2:17, rows = 53)
mergeCells(wb, "TableS3", cols = 2:17, rows = 54)
mergeCells(wb, "TableS3", cols = 2:17, rows = 55)
addStyle(wb, "TableS3", style = s3_body_left, rows = 53:54, cols = 2, gridExpand = TRUE)
addStyle(wb, "TableS3", style = s3_bottom_left, rows = 55, cols = 2:17, gridExpand = TRUE)

# setColWidths(wb, "TableS3", cols = 2:17, widths = c(10, 10, 10, 7, 1, 5, 5, 5, 1, 5, 5, 5, 1, 5, 5, 5))
setRowHeights(wb, "TableS3", rows = 32:55, heights = c(15, 12, 12, 12, 12, 12, 12, 12, 22, 12, 12, 12, 12, 12, 22, 12, 12, 12, 12, 12, 22, 12, 12, 14))


#Table S4--------------------------------------------------------
s4_dc <- deadC[yr == 2018, ]
liveyrs <- deadC[yr != 2018 & str_detect(loc, "NA;", negate = TRUE), .(nyr = .N, yrs = capture.output(cat(unlist(unique(yr)), sep = ", "))), by = loc]
s4_dc <- merge(s4_dc, liveyrs[, .(loc, yrs)], by = "loc", all = TRUE)
s4_dc[loc == "Alligator Harbor", loc := "Alligator Harbor*"]
s4_dc[, loc := factor(loc, levels = c("Little St. George Island", "Goose Island/East Cove", "Alligator Harbor*", "Lone Cabbage", "Lemon Bay", "Hendry Creek/Mullock Creek", "New Pass", 
                                      "Big Hickory", "Jack Island", "Pellicer Creek", "Matanzas River", "Guana River"))]
setorder(s4_dc, loc)
s4_dc[, `:=` (genus = c("Crassostrea", "Crassostrea", "Mercenaria", rep("Crassostrea", 9)), yr = NULL)]
setcolorder(s4_dc, c("loc", "genus", "n", "yrs", "f14c_mn", "f14c_sd", "ref_mn", "ref_sd", "deadc", "deadc_sd"))
s4_dc[, `:=` (f14c_mn = plyr::round_any(f14c_mn, 0.0001), f14c_sd = plyr::round_any(f14c_sd, 0.0001),
              ref_mn = plyr::round_any(ref_mn, 0.0001), ref_sd = plyr::round_any(ref_sd, 0.0001),
              deadc = plyr::round_any(deadc, 0.0001), deadc_sd = plyr::round_any(deadc_sd, 0.0001))]

# wb <- createWorkbook()
addWorksheet(wb, "TableS4")
s4_caption <- createStyle(halign = "center", valign = "top", border = "Bottom", borderStyle = "double", fontName = "Arial", fontSize = 8)
s4_headers_nc1 <- createStyle(halign = "center", valign = "bottom", wrapText = TRUE, border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 8)
s4_headers_c1 <- createStyle(halign = "left", valign = "bottom", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 8)
s4_body_left <- createStyle(halign = "left", valign = "top", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s4_body_right_int <- createStyle(halign = "right", valign = "top", numFmt = "#,##0", fontName = "Arial", fontSize = 8)
s4_body_right_num <- createStyle(halign = "right", valign = "top", numFmt = "#,##0.0000", fontName = "Arial", fontSize = 8)
s4_bottom_left <- createStyle(halign = "left", valign = "top", border = "Bottom", borderStyle = "thin", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s4_bottom_right_int <- createStyle(halign = "right", valign = "top", numFmt = "#,##0", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 8)
s4_bottom_right_num <- createStyle(halign = "right", valign = "top", numFmt = "#,##0.0000", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 8)

s4_title <- "TABLE S4. DEAD CARBON ESTIMATES FOR EACH LOCALITY"
s4_colnms <- data.table("Locality\u2020", "Genus", "N", "Collection Year(s)", "Weighted Mean F14Clive", "\u00B1", "Weighted Mean F14Cref", "\u00B1", "Dead C", "\u00B1")
s4_foots <- c("*Location for full-marine salinity clam specimens; no oysters collected.", "\u2020Localities are listed in counter-clockwise geographic order around the state, starting at the panhandle.")

writeData(wb, "TableS4", s4_title, startRow = 2, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS4", cols = 2:11, rows = 2)
addStyle(wb, "TableS4", style = s4_caption, rows = 2, cols = 2:11, gridExpand = TRUE)

writeData(wb, "TableS4", s4_colnms, startRow = 3, startCol = 2, colNames = FALSE)
addStyle(wb, "TableS4", style = s4_headers_nc1, rows = 3, cols = 3:11, gridExpand = TRUE)
addStyle(wb, "TableS4", style = s4_headers_c1, rows = 3, cols = 2, gridExpand = TRUE)

writeData(wb, "TableS4", s4_dc, startRow = 4, startCol = 2, colNames = FALSE)
addStyle(wb, "TableS4", style = s4_body_left, rows = 4:14, cols = c(2:3, 5), gridExpand = TRUE)
addStyle(wb, "TableS4", style = s4_body_right_int, rows = 4:14, cols = 4, gridExpand = TRUE)
addStyle(wb, "TableS4", style = s4_body_right_num, rows = 4:14, cols = 6:11, gridExpand = TRUE)
addStyle(wb, "TableS4", style = s4_bottom_left, rows = 15, cols = c(2:3, 5), gridExpand = TRUE)
addStyle(wb, "TableS4", style = s4_bottom_right_int, rows = 15, cols = 4, gridExpand = TRUE)
addStyle(wb, "TableS4", style = s4_bottom_right_num, rows = 15, cols = 6:11, gridExpand = TRUE)

writeData(wb, "TableS4", s4_foots, startRow = 16, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS4", cols = 2:11, rows = 16)
mergeCells(wb, "TableS4", cols = 2:11, rows = 17)
addStyle(wb, "TableS4", style = s4_body_left, rows = 16, cols = 2:11, gridExpand = TRUE)
addStyle(wb, "TableS4", style = s4_bottom_left, rows = 17, cols = 2:11, gridExpand = TRUE)

setColWidths(wb, "TableS4", cols = 2:11, widths = c(15, 10, 3, 10, 12, 6, 12, 6, 10, 6))
setRowHeights(wb, "TableS4", rows = 2:17, heights = c(20, 25, 24, 24, 15, 15, 15, 24, 15, 15, 15, 15, 15, 15, 15, 15))


# #Table S5 (now Fig. S4)--------------------------------------------------------
# s5dat <- data.table(catno = c("UF 15484", "UF 512435", "UF 512435", "UF 15491"),
#                     colyr = c(1938, 1979, 1979, 1932),
#                     nrloc = c("Little St. George Is.", "Lone Cabbage", "Lone Cabbage", "Big Hickory"),
#                     samp = c("UAL 19523", "UAL 19520", "UAL 19521", "UAL 19522"),
#                     f14c = )
#I decided to present the info from the original Table S5 as a figure instead (Figure S2, code towards the beginning of this script)


#Table S5 (formerly Table S6)--------------------------------------------------------
#note that subscripts in the first row of headings need to be formatted manually in Excel after export
stater <- results[variable %in% c("mu_region", "sigma_hole", "sigma_reef", "sigma_locality"), .(medage1525_median = plyr::round_any(medage1525_median, 0.01), 
                                                                                                medage1525_q5 = plyr::round_any(medage1525_q5, 0.01), 
                                                                                                medage1525_q95 = plyr::round_any(medage1525_q95, 0.01), 
                                                                                                medage2535_median = plyr::round_any(medage2535_median, 0.01), 
                                                                                                medage2535_q5 = plyr::round_any(medage2535_q5, 0.01), 
                                                                                                medage2535_q95 = plyr::round_any(medage2535_q95, 0.01), 
                                                                                                medageddmedian = plyr::round_any(medageddmedian, 0.01), 
                                                                                                medageddq5 = plyr::round_any(medageddq5, 0.01), 
                                                                                                medageddq95 = plyr::round_any(medageddq95, 0.01), 
                                                                                                tav1525_median = plyr::round_any(tav1525_median, 0.01), 
                                                                                                tav1525_q5 = plyr::round_any(tav1525_q5, 0.01), 
                                                                                                tav1525_q95 = plyr::round_any(tav1525_q95, 0.01), 
                                                                                                tav2535_median = plyr::round_any(tav2535_median, 0.01), 
                                                                                                tav2535_q5 = plyr::round_any(tav2535_q5, 0.01), 
                                                                                                tav2535_q95 = plyr::round_any(tav2535_q95, 0.01), 
                                                                                                tavddmedian = plyr::round_any(tavddmedian, 0.01), 
                                                                                                tavddq5 = plyr::round_any(tavddq5, 0.01), 
                                                                                                tavddq95 = plyr::round_any(tavddq95, 0.01)), by = list(locality, variable)]

stater2 <- data.table(variable = rep(c("Median age", "Total age variation"), each = 3),
                      category = rep(c("15-25 cm", "25-35 cm", "Difference\u2020"), 2),
                      mureg_med = unlist(stater[variable == "mu_region", .(medage1525_median, medage2535_median, medageddmedian, tav1525_median, tav2535_median, tavddmedian)]),
                      mureg_q5 = unlist(stater[variable == "mu_region", .(medage1525_q5, medage2535_q5, medageddq5, tav1525_q5, tav2535_q5, tavddq5)]),
                      mureg_q95 = unlist(stater[variable == "mu_region", .(medage1525_q95, medage2535_q95, medageddq95, tav1525_q95, tav2535_q95, tavddq95)]),
                      space1 = c("", "", "", "", "", ""),
                      sigloc_med = unlist(stater[variable == "sigma_locality", .(medage1525_median, medage2535_median, medageddmedian, tav1525_median, tav2535_median, tavddmedian)]),
                      sigloc_q5 = unlist(stater[variable == "sigma_locality", .(medage1525_q5, medage2535_q5, medageddq5, tav1525_q5, tav2535_q5, tavddq5)]),
                      sigloc_q95 = unlist(stater[variable == "sigma_locality", .(medage1525_q95, medage2535_q95, medageddq95, tav1525_q95, tav2535_q95, tavddq95)]),
                      space1 = c("", "", "", "", "", ""),
                      sigreef_med = unlist(stater[variable == "sigma_reef", .(medage1525_median, medage2535_median, medageddmedian, tav1525_median, tav2535_median, tavddmedian)]),
                      sigreef_q5 = unlist(stater[variable == "sigma_reef", .(medage1525_q5, medage2535_q5, medageddq5, tav1525_q5, tav2535_q5, tavddq5)]),
                      sigreef_q95 = unlist(stater[variable == "sigma_reef", .(medage1525_q95, medage2535_q95, medageddq95, tav1525_q95, tav2535_q95, tavddq95)]),
                      space1 = c("", "", "", "", "", ""),
                      sighole_med = unlist(stater[variable == "sigma_hole", .(medage1525_median, medage2535_median, medageddmedian, tav1525_median, tav2535_median, tavddmedian)]),
                      sighole_q5 = unlist(stater[variable == "sigma_hole", .(medage1525_q5, medage2535_q5, medageddq5, tav1525_q5, tav2535_q5, tavddq5)]),
                      sighole_q95 = unlist(stater[variable == "sigma_hole", .(medage1525_q95, medage2535_q95, medageddq95, tav1525_q95, tav2535_q95, tavddq95)]))

# wb <- createWorkbook()
addWorksheet(wb, "TableS5")
s5_caption <- createStyle(halign = "center", valign = "top", border = "Bottom", borderStyle = "double", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s5_headers_nc1_1 <- createStyle(halign = "center", valign = "top", border = "Bottom", borderStyle = "thin", textDecoration = "italic", fontName = "Arial", fontSize = 9)
s5_headers_nc1_2 <- createStyle(halign = "center", valign = "top", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 8)
s5_headers_c1_1 <- createStyle(halign = "left", valign = "top", fontName = "Arial", fontSize = 8)
s5_headers_c2_1 <- createStyle(halign = "center", valign = "top", fontName = "Arial", fontSize = 8)
s5_headers_c1_2 <- createStyle(halign = "left", valign = "top", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 8)
s5_body_left <- createStyle(halign = "left", valign = "top", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s5_body_right <- createStyle(halign = "right", valign = "top", numFmt = "#,##0.00", fontName = "Arial", fontSize = 8)
s5_bottom_left <- createStyle(halign = "left", valign = "top", border = "Bottom", borderStyle = "thin", wrapText = TRUE, fontName = "Arial", fontSize = 8)
s5_bottom_right <- createStyle(halign = "right", valign = "top", numFmt = "#,##0.00", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 8)
# class(madat_hobs$acres_2019) <- "comma"

s5_title <- "TABLE S5. STATEWIDE RESULTS OF HIERARCHICAL MIXED EFFECTS MODELS OF MEDIAN AGE AND TOTAL AGE VARIATION, INCLUDING STATEWIDE MEAN AND STANDARD DEVIATIONS AT LOCALITY, REEF AND SAMPLE HOLE SCALES"
s5_colnms1 <- data.table("Variable", "Category", "e", "", "", "", "\u03c3c", "", "", "", "\u03c3r", "", "", "", "\u03c3h", "", "")
s5_colnms2 <- data.table("", "", "Median", "5 %**", "95 %**", "", "Median", "5 %", "95 %", "", "Median", "5 %", "95 %", "", "Median", "5 %", "95 %")
s5_foots <- c("*CPE = corrected posterior age estimate", "\u2020Differences calculated as (25-35 cm value - 15-25 cm value) for each sample hole", "**5 % and 95 % columns denote the 95 % credible intervals")
writeData(wb, "TableS5", s5_title, startRow = 2, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS5", cols = 2:18, rows = 2)
addStyle(wb, "TableS5", style = s5_caption, rows = 2, cols = 2:18, gridExpand = TRUE)

writeData(wb, "TableS5", s5_colnms1, startRow = 3, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS5", cols = 4:6, rows = 3)
mergeCells(wb, "TableS5", cols = 8:10, rows = 3)
mergeCells(wb, "TableS5", cols = 12:14, rows = 3)
mergeCells(wb, "TableS5", cols = 16:18, rows = 3)
addStyle(wb, "TableS5", style = s5_headers_nc1_1, rows = 3, cols = c(4:6, 8:10, 12:14, 16:18), gridExpand = TRUE)
addStyle(wb, "TableS5", style = s5_headers_c1_1, rows = 3, cols = 2, gridExpand = TRUE)
addStyle(wb, "TableS5", style = s5_headers_c2_1, rows = 3, cols = 3, gridExpand = TRUE)

writeData(wb, "TableS5", s5_colnms2, startRow = 4, startCol = 2, colNames = FALSE)
addStyle(wb, "TableS5", style = s5_headers_nc1_2, rows = 4, cols = 3:18, gridExpand = TRUE)
addStyle(wb, "TableS5", style = s5_headers_c1_2, rows = 4, cols = 2, gridExpand = TRUE)

writeData(wb, "TableS5", stater2, startRow = 5, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS5", cols = 2, rows = 5:7)
mergeCells(wb, "TableS5", cols = 2, rows = 8:10)
addStyle(wb, "TableS5", style = s5_body_left, rows = 5:9, cols = 2:3, gridExpand = TRUE)
addStyle(wb, "TableS5", style = s5_body_right, rows = 5:9, cols = 4:18, gridExpand = TRUE)
addStyle(wb, "TableS5", style = s5_bottom_left, rows = 10, cols = 2:3, gridExpand = TRUE)
addStyle(wb, "TableS5", style = s5_bottom_right, rows = 10, cols = 4:18, gridExpand = TRUE)

writeData(wb, "TableS5", s5_foots, startRow = 11, startCol = 2, colNames = FALSE)
mergeCells(wb, "TableS5", cols = 2:18, rows = 11)
mergeCells(wb, "TableS5", cols = 2:18, rows = 12)
mergeCells(wb, "TableS5", cols = 2:18, rows = 13)
addStyle(wb, "TableS5", style = s5_body_left, rows = 11:12, cols = 2, gridExpand = TRUE)
addStyle(wb, "TableS5", style = s5_bottom_left, rows = 13, cols = 2:18, gridExpand = TRUE)

setColWidths(wb, "TableS5", cols = 2:18, widths = c(10, 10, 5, 5, 5, 1, 5, 5, 5, 1, 5, 5, 5, 1, 5, 5, 5))
setRowHeights(wb, "TableS5", rows = 2:13, heights = c(30, 12, 12, 12, 12, 22, 12, 12, 22, 12, 12, 12))


saveWorkbook(wb, here::here("Summer2023/Durhametal_2023_Tables.xlsx"), overwrite = TRUE)


#Appendix DR1--------------------------------------------------
appendices <- createWorkbook()
addWorksheet(appendices, "DR1")
dr1_title_coldefs <- createStyle(halign = "left", valign = "top", fontName = "Arial", fontSize = 10)
dr1_caption <- createStyle(halign = "left", valign = "top", wrapText = TRUE, fontName = "Arial", fontSize = 10)
dr1_aboveheaders <- createStyle(halign = "center", valign = "top", border = "Bottom", borderStyle = "double", fontName = "Arial", fontSize = 10)
dr1_headers_nc1 <- createStyle(halign = "center", valign = "top", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr1_headers_c1_1 <- createStyle(halign = "left", valign = "top", fontName = "Arial", fontSize = 10)
dr1_headers_c1_2 <- createStyle(halign = "left", valign = "top", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr1_body_left <- createStyle(halign = "left", valign = "top", fontName = "Arial", fontSize = 10)
dr1_body_right_int <- createStyle(halign = "right", valign = "top", numFmt = "###0", fontName = "Arial", fontSize = 10)
dr1_body_right_num <- createStyle(halign = "right", valign = "top", numFmt = "###0.0", fontName = "Arial", fontSize = 10)
dr1_body_right_full <- createStyle(halign = "right", valign = "top", fontName = "Arial", fontSize = 10)
dr1_bottom_left <- createStyle(halign = "left", valign = "top", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr1_bottom_right_int <- createStyle(halign = "right", valign = "top", numFmt = "###0", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr1_bottom_right_num <- createStyle(halign = "right", valign = "top", numFmt = "###0.0", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr1_bottom_right_full <- createStyle(halign = "right", valign = "top", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)

dr1_title <- "Durham et al., 2023, Age variability and time averaging in oyster reef death assemblages"
dr1_cap <- "Appendix DR1. Sample-level geochronology results. All samples are stored in the collections of the Paleontological Research Institution, Ithaca, NY 14850 (www.priweb.org/research-and-collections/research-collection)."
dr1_coldefs <- c("Locality: Name of the study locality where the sample was collected",
                 "Station: Station number assigned by the Paleontological Research Institution",
                 "Lat: Latitude of the sample collection location",
                 "Lon: Longitude of the sample collection location",
                 "IGSN: International GeoSample Number (www.geosamples.org)",
                 "Depth: Burial depth of the sample",
                 "N: Number of specimens dated from each sample",
                 "Median: Median probability calibrated age",
                 "Min: Minimum calibrated age",
                 "Max: Maximum calibrated age",
                 "TAV_IQR: Total age variability")

dr1_colnms1 <- data.table("", "", "", "", "", "", "", "", "Sample_age", "", "")
dr1_colnms2 <- data.table("Locality", "Station", "Lat", "Lon", "IGSN", "Depth", "N", "Median", "Min", "Max", "TAV_IQR")

dr1_dat <- distinct(hobsf14c2[!is.na(tav_iqr), .(locality, station, lat, lon, igsn, depth, n, tav_median = plyr::round_any(tav_median, 0.1), tav_min = plyr::round_any(tav_min, 0.1), tav_max = plyr::round_any(tav_max, 0.1), tav_iqr = plyr::round_any(tav_iqr, 0.01))])

writeData(appendices, "DR1", dr1_title, startRow = 1, startCol = 1, colNames = FALSE)
addStyle(appendices, "DR1", dr1_title_coldefs, rows = 1, cols = 1, gridExpand = TRUE)
writeData(appendices, "DR1", dr1_cap, startRow = 3, startCol = 1, colNames = FALSE)
mergeCells(appendices, "DR1", cols = 1:9, rows = 3:5)
addStyle(appendices, "DR1", dr1_caption, rows = 3:5, cols = 1:9, gridExpand = TRUE)
writeData(appendices, "DR1", dr1_coldefs, startRow = 7, startCol = 1, colNames = FALSE)
addStyle(appendices, "DR1", dr1_title_coldefs, rows = 7:17, cols = 1, gridExpand = TRUE)

addStyle(appendices, "DR1", style = dr1_aboveheaders, rows = 18, cols = 1:11, gridExpand = TRUE)
writeData(appendices, "DR1", dr1_colnms1, startRow = 19, startCol = 1, colNames = FALSE)
addStyle(appendices, "DR1", style = dr1_headers_nc1, rows = 19, cols = 8:10, gridExpand = TRUE)
addStyle(appendices, "DR1", style = dr1_headers_c1_1, rows = 19, cols = 1, gridExpand = TRUE)

writeData(appendices, "DR1", dr1_colnms2, startRow = 20, startCol = 1, colNames = FALSE)
addStyle(appendices, "DR1", style = dr1_headers_nc1, rows = 20, cols = 2:11, gridExpand = TRUE)
addStyle(appendices, "DR1", style = dr1_headers_c1_2, rows = 20, cols = 1, gridExpand = TRUE)

writeData(appendices, "DR1", dr1_dat, startRow = 21, startCol = 1, colNames = FALSE)
addStyle(appendices, "DR1", style = dr1_body_left, rows = 21:144, cols = c(1, 5:6), gridExpand = TRUE)
addStyle(appendices, "DR1", style = dr1_body_right_full, rows = 21:144, cols = 3:4, gridExpand = TRUE)
addStyle(appendices, "DR1", style = dr1_body_right_int, rows = 21:144, cols = c(2, 7), gridExpand = TRUE)
addStyle(appendices, "DR1", style = dr1_body_right_num, rows = 21:144, cols = 8:11, gridExpand = TRUE)
addStyle(appendices, "DR1", style = dr1_bottom_left, rows = 145, cols = c(1, 5:6), gridExpand = TRUE)
addStyle(appendices, "DR1", style = dr1_bottom_right_full, rows = 145, cols = 3:4, gridExpand = TRUE)
addStyle(appendices, "DR1", style = dr1_bottom_right_int, rows = 145, cols = c(2, 7), gridExpand = TRUE)
addStyle(appendices, "DR1", style = dr1_bottom_right_num, rows = 145, cols = 8:11, gridExpand = TRUE)

setColWidths(appendices, "DR1", cols = 1:11, widths = c(25, 8, 12, 13, 10, 8, 4, 8, 8, 8, 8))


#Appendix DR2--------------------------------------------------------------------
# appendices <- createWorkbook()
addWorksheet(appendices, "DR2")
dr2_title_coldefs <- createStyle(halign = "left", valign = "top", fontName = "Arial", fontSize = 10)
dr2_caption <- createStyle(halign = "left", valign = "top", wrapText = TRUE, fontName = "Arial", fontSize = 10)
dr2_aboveheaders <- createStyle(halign = "left", valign = "top", border = "Bottom", borderStyle = "double", fontName = "Arial", fontSize = 10)
dr2_headers_nc1 <- createStyle(halign = "center", valign = "top", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr2_headers_c1 <- createStyle(halign = "left", valign = "top", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr2_body_left <- createStyle(halign = "left", valign = "top", fontName = "Arial", fontSize = 10)
dr2_body_right_int <- createStyle(halign = "right", valign = "top", numFmt = "###0", fontName = "Arial", fontSize = 10)
dr2_body_right_num <- createStyle(halign = "right", valign = "top", numFmt = "###0.0", fontName = "Arial", fontSize = 10)
dr2_body_right_full <- createStyle(halign = "right", valign = "top", fontName = "Arial", fontSize = 10)
dr2_body_right_date <- createStyle(halign = "right", valign = "top", numFmt = "DATE", fontName = "Arial", fontSize = 10)
dr2_bottom_left <- createStyle(halign = "left", valign = "top", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr2_bottom_right_int <- createStyle(halign = "right", valign = "top", numFmt = "###0", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr2_bottom_right_num <- createStyle(halign = "right", valign = "top", numFmt = "###0.0", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr2_bottom_right_full <- createStyle(halign = "right", valign = "top", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr2_bottom_right_date <- createStyle(halign = "right", valign = "top", numFmt = "DATE", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)

dr2_title <- "Durham et al., 2023, Age variability and time averaging in oyster reef death assemblages"
dr2_cap <- "Appendix DR2. High-precision (HP) radiocarbon AMS results (i.e., graphite targets), including both OxCal and empirical posterior distribution calibration results, as applicable. All specimens except those with catalog numbers beginning with \"UF\" are stored in the collections of the Paleontological Research Institution, Ithaca, NY 14850 (www.priweb.org/research-and-collections/research-collection)."
dr2_coldefs <- c("Cat_No: Museum catalog number assigned to the specimen (\"PRI\" and \"UF\" prefixes correspond to the Paleontological Research Institution and the Florida Museum of Natural History, respectively)",
                 "Collection_Date: Date the specimen was collected",
                 "Locality: Name of the study locality where the sample was collected",
                 "Station: Station number assigned by the Paleontological Research Institution",
                 "Lat: Latitude of the specimen collection location",
                 "Lon: Longitude of the specimen collection location",
                 "IGSN: International GeoSample Number (www.geosamples.org)",
                 "Genus: Genus of the radiocarbon-dated specimen",
                 "Live_or_Dead: Categorical variable indicating whether a specimen was alive (or presumed alive) when it was collected",
                 "NAU_ACE_ID: Sample ID assigned by the Arizona Climate and Ecosystems (ACE) Isotope Laboratory at Northern Arizona University",
                 "UCIAMS_ID: Sample ID assigned by the Keck AMS Laboratory at University of California, Irvine",
                 "NAU_ID: Sample ID assigned by Northern Arizona University Quaternary Geochronology Laboratory",
                 "F14C: Fraction modern carbon",
                 "F14C_sd: Standard deviation of the fraction modern carbon value (F14C)",
                 "14C_age: Uncalibrated radiocarbon age",
                 "14C_age_sd: Standard deviation of the uncalibrated radiocarbon age value (14C_age)",
                 "Dead_C: Correction factor used to adjust F14C values for local hardwater or estuarine effects contributing \"old\" carbon 14", 
                 "Dead_C_sd: Standard deviation of the Dead_C value",
                 "F14C_corr: F14C value corrected for \"dead\" carbon content (Dead_C)",
                 "F14C_corr_sd: Standard deviation of the corrected F14C value",
                 "Out_of_range: Boolean variable denoting whether the specimen F14C_corr value was within the range of the calibration curve",
                 "Cal_age_med: Median OxCal likelihood-based calibrated age",
                 "Cal_age_sd: Standard deviation of the OxCal likelihood-based calibrated age",
                 "EPDcal_age_med: Median empirical posterior distribution (EPD) based calibrated age",
                 "EPDcal_age_iqr: Interquartile range (IQR) of the empirical posterior distribution (EPD) based calibrated age",
                 "EPDcal_age_min: Minimum of the empirical posterior distribution (EPD) based calibrated age",
                 "EPDcal_age_max: Maximum of the empirical posterior distribution (EPD) based calibrated age")

dr2_colnms <- data.table("Cat_No", "Collection_Date", "Locality", "Station", "Lat", "Lon", "IGSN", "Genus", "Live_or_Dead", "NAU_ACE_ID", "UCIAMS_ID", 
                         "NAU_ID", "F14C", "F14C_sd", "14C_age", "14C_age_sd", "Dead_C", "Dead_C_sd", "F14C_corr", "F14C_corr_sd", "Out_of_range", "Cal_age_med", 
                         "Cal_age_sd", "EPDcal_age_med", "EPDcal_age_iqr", "EPDcal_age_min", "EPDcal_age_max")

dr2_dat <- distinct(hobsf14c2[analysis == "HP" & !is.na(cat_no), .(cat_no, collection_date, locality, station, lat, lon, igsn, genus, live_or_dead, nau_ace_number, uciams_id, nau_id, 
                                                                   f14c, f14c_sd, x14c_age, x14c_age_sd, f14c_corr, f14c_corr_sd, oor, cal_age_med, cal_age_sd, spec_median, 
                                                                   spec_iqr, spec_min, spec_max)])

dr2_dat[, `:=` (deadc = deadC[yr == year(collection_date) & loc == locality, deadc],
                deadc_sd = deadC[yr == year(collection_date) & loc == locality, deadc_sd]), by = row.names(dr2_dat)]
dr2_dat[oor == TRUE, `:=` (cal_age_med = NA, cal_age_sd = NA, spec_median = NA, spec_iqr = NA, spec_min = NA, spec_max = NA)]

setcolorder(dr2_dat, c("cat_no", "collection_date", "locality", "station", "lat", "lon", "igsn", "genus", "live_or_dead", "nau_ace_number", "uciams_id", "nau_id", 
                       "f14c", "f14c_sd", "x14c_age", "x14c_age_sd", "deadc", "deadc_sd", "f14c_corr", "f14c_corr_sd", "oor", "cal_age_med", "cal_age_sd", "spec_median", 
                       "spec_iqr", "spec_min", "spec_max"))

writeData(appendices, "DR2", dr2_title, startRow = 1, startCol = 1, colNames = FALSE)
writeData(appendices, "DR2", dr2_cap, startRow = 3, startCol = 1, colNames = FALSE)
mergeCells(appendices, "DR2", cols = 1:9, rows = 3:5)
addStyle(appendices, "DR2", dr2_caption, rows = 3:5, cols = 1:9, gridExpand = TRUE)
writeData(appendices, "DR2", dr2_coldefs, startRow = 7, startCol = 1, colNames = FALSE)
addStyle(appendices, "DR2", dr2_title_coldefs, rows = c(1, 7:(length(dr2_colnms) + 6)), cols = 1, gridExpand = TRUE)

writeData(appendices, "DR2", "NA = no data", startRow = (length(dr2_colnms) + 8), startCol = 1, colNames = FALSE)
addStyle(appendices, "DR2", style = dr2_aboveheaders, rows = (length(dr2_colnms) + 8), cols = (1:length(dr2_colnms)), gridExpand = TRUE)
writeData(appendices, "DR2", dr2_colnms, startRow = (length(dr2_colnms) + 9), startCol = 1, colNames = FALSE)
addStyle(appendices, "DR2", style = dr2_headers_nc1, rows = (length(dr2_colnms) + 9), cols = (2:length(dr2_colnms)), gridExpand = TRUE)
addStyle(appendices, "DR2", style = dr2_headers_c1, rows = (length(dr2_colnms) + 9), cols = 1, gridExpand = TRUE)

writeData(appendices, "DR2", dr2_dat, startRow = (length(dr2_colnms) + 10), startCol = 1, colNames = FALSE)
addStyle(appendices, "DR2", style = dr2_body_left, rows = (length(dr2_colnms) + 10):(length(dr2_colnms) + 8 + nrow(dr2_dat)), cols = c(1, 3, 7:9, 12, 19), gridExpand = TRUE)
addStyle(appendices, "DR2", style = dr2_body_right_full, rows = (length(dr2_colnms) + 10):(length(dr2_colnms) + 8 + nrow(dr2_dat)), cols = c(4:6, 10:11, 13:18, 20:length(dr2_colnms)), gridExpand = TRUE)
addStyle(appendices, "DR2", style = dr2_body_right_date, rows = (length(dr2_colnms) + 10):(length(dr2_colnms) + 8 + nrow(dr2_dat)), cols = 2, gridExpand = TRUE)
addStyle(appendices, "DR2", style = dr2_bottom_left, rows = (length(dr2_colnms) + 9 + nrow(dr2_dat)), cols = c(1, 3, 7:9, 12, 19), gridExpand = TRUE)
addStyle(appendices, "DR2", style = dr2_bottom_right_full, rows = (length(dr2_colnms) + 9 + nrow(dr2_dat)), cols = c(4:6, 10:11, 13:18, 20:length(dr2_colnms)), gridExpand = TRUE)
addStyle(appendices, "DR2", style = dr2_bottom_right_date, rows = (length(dr2_colnms) + 9 + nrow(dr2_dat)), cols = 2, gridExpand = TRUE)

setColWidths(appendices, "DR2", cols = 1:length(dr2_colnms), widths = c(10, 18, 30, 8, 12, 12, 10, 12, 15, 15, 15, 15, 12, 12, 15, rep(12, 8), rep(15, 4)))
setRowHeights(appendices, "DR2", rows = 3:5, heights = rep(15, 3))


#Appendix DR3--------------------------------------------------------------------
# appendices <- createWorkbook()
addWorksheet(appendices, "DR3")
dr3_title_coldefs <- createStyle(halign = "left", valign = "top", fontName = "Arial", fontSize = 10)
dr3_caption <- createStyle(halign = "left", valign = "top", wrapText = TRUE, fontName = "Arial", fontSize = 10)
dr3_aboveheaders <- createStyle(halign = "left", valign = "top", border = "Bottom", borderStyle = "double", fontName = "Arial", fontSize = 10)
dr3_headers_nc1 <- createStyle(halign = "center", valign = "top", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr3_headers_c1 <- createStyle(halign = "left", valign = "top", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr3_body_left <- createStyle(halign = "left", valign = "top", fontName = "Arial", fontSize = 10)
dr3_body_right_int <- createStyle(halign = "right", valign = "top", numFmt = "###0", fontName = "Arial", fontSize = 10)
dr3_body_right_num <- createStyle(halign = "right", valign = "top", numFmt = "###0.0", fontName = "Arial", fontSize = 10)
dr3_body_right_full <- createStyle(halign = "right", valign = "top", fontName = "Arial", fontSize = 10)
dr3_body_right_date <- createStyle(halign = "right", valign = "top", numFmt = "DATE", fontName = "Arial", fontSize = 10)
dr3_bottom_left <- createStyle(halign = "left", valign = "top", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr3_bottom_right_int <- createStyle(halign = "right", valign = "top", numFmt = "###0", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr3_bottom_right_num <- createStyle(halign = "right", valign = "top", numFmt = "###0.0", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr3_bottom_right_full <- createStyle(halign = "right", valign = "top", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr3_bottom_right_date <- createStyle(halign = "right", valign = "top", numFmt = "DATE", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)

dr3_title <- "Durham et al., 2023, Age variability and time averaging in oyster reef death assemblages"
dr3_cap <- "Appendix DR3. Low-precision radiocarbon AMS results (i.e., carbonate targets), including both OxCal and empirical posterior distribution calibration results. All specimens are stored in the collections of the Paleontological Research Institution, Ithaca, NY 14850 (www.priweb.org/research-and-collections/research-collection)."
dr3_coldefs <- c("Locality: Name of the locality where the sample was collected.",
                 "Station: Station number assigned by the Paleontological Research Institution",
                 "Lat: Latitude of the specimen collection location",
                 "Lon: Longitude of the specimen collection location",
                 "IGSN: International GeoSample Number (www.geosamples.org)",
                 "Depth: Burial depth of the sample",
                 "Cat_No: Museum catalog number assigned to the specimen",
                 "NAU_ID: Sample ID assigned by Northern Arizona University",
                 "F14C: Fraction modern carbon",
                 "F14C_sd: Standard deviation of the fraction modern carbon value (F14C)",
                 "14C_age: Uncalibrated radiocarbon age",
                 "14C_age_sd: Standard deviation of the uncalibrated radiocarbon age value (14C_age)",
                 "Dead_C: Correction factor to adjust F14C values for local hardwater or estuarine effects contributing \"old\" carbon 14",
                 "Dead_C_sd: Standard deviation of the Dead_C value",
                 "F14C_corr: F14C value corrected for \"dead\" carbon content (Dead_C)",
                 "F14C_corr_sd: Standard deviation of the corrected F14C value",
                 "Out_of_range: Boolean variable denoting whether the specimen F14C_corr value was within the range of the calibration curve",
                 "Cal_age_95CI_low: Lower bound of 95% confidence interval for the OxCal likelihood-based calibrated age",
                 "Cal_age_95CI_high: Upper bound of 95% confidence interval for the OxCal likelihood-based calibrated age",
                 "Cal_age_med: Median OxCal likelihood-based calibrated age",
                 "EPDcal_age_med: Median empirical posterior distribution (EPD) based calibrated age",
                 "EPDcal_age_iqr: Interquartile range (IQR) of the empirical posterior distribution (EPD) based calibrated age",
                 "EPDcal_age_min: Minimum of the empirical posterior distribution (EPD) based calibrated age",
                 "EPDcal_age_max: Maximum of the empirical posterior distribution (EPD) based calibrated age")

dr3_colnms <- data.table("Locality", "Station", "Lat", "Lon", "IGSN", "Depth", "Cat_No", "NAU_ACE_ID", "NAU_ID", "F14C", "F14C_sd", "14C_age", 
                         "14C_age_sd", "Dead_C", "Dead_C_sd", "F14C_corr", "F14C_corr_sd", "Out_of_range", "Cal_age_med", "Cal_age_sd", 
                         "EPDcal_age_med", "EPDcal_age_iqr", "EPDcal_age_min", "EPDcal_age_max")

dr3_dat <- distinct(hobsf14c2[analysis == "LP" & !is.na(cat_no) & !is.na(oor), .(locality, station, lat, lon, igsn, depth, cat_no, nau_ace_number, nau_id, f14c, 
                                                                                 f14c_sd, x14c_age, x14c_age_sd, f14c_corr, f14c_corr_sd, oor, cal_age_med, 
                                                                                 cal_age_sd, spec_median, spec_iqr, spec_min, spec_max)])

dr3_dat[, `:=` (deadc = deadC[yr == 2018 & loc == locality, deadc],
                deadc_sd = deadC[yr == 2018 & loc == locality, deadc_sd]), by = row.names(dr3_dat)]
dr3_dat[oor == TRUE, `:=` (cal_age_med = NA, cal_age_sd = NA, spec_median = NA, spec_iqr = NA, spec_min = NA, spec_max = NA)]

setcolorder(dr3_dat, c("locality", "station", "lat", "lon", "igsn", "depth", "cat_no", "nau_ace_number", "nau_id", "f14c", "f14c_sd", "x14c_age", 
                       "x14c_age_sd", "deadc", "deadc_sd", "f14c_corr", "f14c_corr_sd", "oor", "cal_age_med", "cal_age_sd", "spec_median", "spec_iqr", 
                       "spec_min", "spec_max"))

writeData(appendices, "DR3", dr3_title, startRow = 1, startCol = 1, colNames = FALSE)
writeData(appendices, "DR3", dr3_cap, startRow = 3, startCol = 1, colNames = FALSE)
mergeCells(appendices, "DR3", cols = 1:9, rows = 3:5)
addStyle(appendices, "DR3", dr3_caption, rows = 3:5, cols = 1:9, gridExpand = TRUE)
writeData(appendices, "DR3", dr3_coldefs, startRow = 7, startCol = 1, colNames = FALSE)
addStyle(appendices, "DR3", dr3_title_coldefs, rows = c(1, 7:(length(dr3_colnms) + 6)), cols = 1, gridExpand = TRUE)

writeData(appendices, "DR3", "NA = no data", startRow = (length(dr3_colnms) + 8), startCol = 1, colNames = FALSE)
addStyle(appendices, "DR3", style = dr3_aboveheaders, rows = (length(dr3_colnms) + 8), cols = (1:length(dr3_colnms)), gridExpand = TRUE)
writeData(appendices, "DR3", dr3_colnms, startRow = (length(dr3_colnms) + 9), startCol = 1, colNames = FALSE)
addStyle(appendices, "DR3", style = dr3_headers_nc1, rows = (length(dr3_colnms) + 9), cols = (2:length(dr3_colnms)), gridExpand = TRUE)
addStyle(appendices, "DR3", style = dr3_headers_c1, rows = (length(dr3_colnms) + 9), cols = 1, gridExpand = TRUE)

writeData(appendices, "DR3", dr3_dat, startRow = (length(dr3_colnms) + 10), startCol = 1, colNames = FALSE)
addStyle(appendices, "DR3", style = dr3_body_left, rows = (length(dr3_colnms) + 10):(length(dr3_colnms) + 8 + nrow(dr3_dat)), cols = c(1, 5:8, 17), gridExpand = TRUE)
addStyle(appendices, "DR3", style = dr3_body_right_full, rows = (length(dr3_colnms) + 10):(length(dr3_colnms) + 8 + nrow(dr3_dat)), cols = c(2:4, 9:16, 18:length(dr3_colnms)), gridExpand = TRUE)
addStyle(appendices, "DR3", style = dr3_bottom_left, rows = (length(dr3_colnms) + 9 + nrow(dr3_dat)), cols = c(1, 5:8, 17), gridExpand = TRUE)
addStyle(appendices, "DR3", style = dr3_bottom_right_full, rows = (length(dr3_colnms) + 9 + nrow(dr3_dat)), cols = c(2:4, 9:16, 18:length(dr3_colnms)), gridExpand = TRUE)

setColWidths(appendices, "DR3", cols = 1:length(dr3_colnms), widths = c(20, 8, 15, 15, 15, 10, 12, 12, rep(10, 3), 18, 15, rep(10, 3), rep(15, 8)))
setRowHeights(appendices, "DR3", rows = 3:5, heights = rep(15, 3))


#Appendix DR4--------------------------------------------------------------------
# appendices <- createWorkbook()
addWorksheet(appendices, "DR4")
dr4_title_coldefs <- createStyle(halign = "left", valign = "top", fontName = "Arial", fontSize = 10)
dr4_caption <- createStyle(halign = "left", valign = "top", wrapText = TRUE, fontName = "Arial", fontSize = 10)
dr4_aboveheaders <- createStyle(halign = "left", valign = "top", border = "Bottom", borderStyle = "double", fontName = "Arial", fontSize = 10)
dr4_headers_nc1 <- createStyle(halign = "center", valign = "top", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr4_headers_c1 <- createStyle(halign = "left", valign = "top", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr4_body_left <- createStyle(halign = "left", valign = "top", fontName = "Arial", fontSize = 10)
dr4_body_right_int <- createStyle(halign = "right", valign = "top", numFmt = "###0", fontName = "Arial", fontSize = 10)
dr4_body_right_num <- createStyle(halign = "right", valign = "top", numFmt = "###0.0", fontName = "Arial", fontSize = 10)
dr4_body_right_full <- createStyle(halign = "right", valign = "top", fontName = "Arial", fontSize = 10)
dr4_body_right_date <- createStyle(halign = "right", valign = "top", numFmt = "DATE", fontName = "Arial", fontSize = 10)
dr4_bottom_left <- createStyle(halign = "left", valign = "top", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr4_bottom_right_int <- createStyle(halign = "right", valign = "top", numFmt = "###0", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr4_bottom_right_num <- createStyle(halign = "right", valign = "top", numFmt = "###0.0", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr4_bottom_right_full <- createStyle(halign = "right", valign = "top", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)
dr4_bottom_right_date <- createStyle(halign = "right", valign = "top", numFmt = "DATE", border = "Bottom", borderStyle = "thin", fontName = "Arial", fontSize = 10)

dr4_title <- "Durham et al., 2023, Age variability and time averaging in oyster reef death assemblages"
dr4_cap <- "Appendix DR4. The calibration curve values used to calibrate the radiocarbon ages in this study. The values from 1950 and older are from the Marine20 curve using a constant regional marine reservoir correction, ΔR = −134 ± 26 yr. The values from 1950.5 to 2022 were generated by fitting a Bayesian generalized additive model to 14C measurements of marine carbonates from across the Gulf of Mexico and western Caribbean region, then calculating the mean and standard deviation of 1000 draws from the posterior distribution at each half-year interval. See the main text for additional details."
dr4_coldefs <- c("Year: Year values of the calibration curve at which the calibration curve model posterior distribution was sampled",
                 "F14C: Average of 1000 draws from the posterior distribution of fraction modern carbon (F14C) at the corresponding year value",
                 "F14C_sd: Standard deviation of 1000 draws from the posterior distribution of fraction modern carbon (F14C) at the corresponding year value")
                 
dr4_colnms <- data.table("Year", "F14C", "F14C_sd")

writeData(appendices, "DR4", dr4_title, startRow = 1, startCol = 1, colNames = FALSE)
writeData(appendices, "DR4", dr4_cap, startRow = 3, startCol = 1, colNames = FALSE)
mergeCells(appendices, "DR4", cols = 1:9, rows = 3:6)
addStyle(appendices, "DR4", dr4_caption, rows = 3:6, cols = 1:9, gridExpand = TRUE)
writeData(appendices, "DR4", dr4_coldefs, startRow = 8, startCol = 1, colNames = FALSE)
addStyle(appendices, "DR4", dr4_title_coldefs, rows = c(1, 8:(length(dr4_colnms) + 7)), cols = 1, gridExpand = TRUE)

addStyle(appendices, "DR4", style = dr4_aboveheaders, rows = (length(dr4_colnms) + 8), cols = (1:length(dr4_colnms)), gridExpand = TRUE)
writeData(appendices, "DR4", dr4_colnms, startRow = (length(dr4_colnms) + 9), startCol = 1, colNames = FALSE)
addStyle(appendices, "DR4", style = dr4_headers_nc1, rows = (length(dr4_colnms) + 9), cols = (2:length(dr4_colnms)), gridExpand = TRUE)
addStyle(appendices, "DR4", style = dr4_headers_c1, rows = (length(dr4_colnms) + 9), cols = 1, gridExpand = TRUE)

writeData(appendices, "DR4", predictedCurve, startRow = (length(dr4_colnms) + 10), startCol = 1, colNames = FALSE)
addStyle(appendices, "DR4", style = dr4_body_right_full, rows = (length(dr4_colnms) + 10):(length(dr4_colnms) + 8 + nrow(predictedCurve)), cols = 1:3, gridExpand = TRUE)
addStyle(appendices, "DR4", style = dr4_bottom_right_full, rows = (length(dr4_colnms) + 9 + nrow(predictedCurve)), cols = 1:3, gridExpand = TRUE)

setColWidths(appendices, "DR4", cols = 1:length(dr4_colnms), widths = c(8, 12, 12))
setRowHeights(appendices, "DR4", rows = 3:6, heights = rep(15, 4))


saveWorkbook(appendices, here::here("Summer2023/Durhametal_2023_Appendices.xlsx"), overwrite = TRUE)



#Figure S1----------------------------------------------------------------------
#Note, uses some of the same base layers loaded for Figure 1 above

currents %>% 
  mutate(x_nudge = case_when(OBJECTID == 96 ~ -4,
                             OBJECTID == 98 ~ 6,
                             OBJECTID == 99 ~ 1,
                             OBJECTID == 138 ~ -2.5,
                             TRUE ~ 0),
         y_nudge = case_when(OBJECTID == 96 ~ 2,
                             OBJECTID == 98 ~ -1,
                             OBJECTID == 99 ~ -1.75,
                             OBJECTID == 138 ~ 0.5,
                             TRUE ~ 0)
  ) ->
  currents

currents$NAME[which(currents$NAME == "South Equatorial")] <- "South\nEquatorial"
currents$NAME[which(currents$NAME == "North Equatorial")] <- "North\nEquatorial"
currents$NAME[which(currents$NAME == "Gulf Stream")] <- "Gulf\nStream"

outbox <- data.table(lat = c(6.8, 6.8, 33.2, 33.2), lon = c(-100.15, -52.9, -52.9, -100.15))
outbox2 <- st_as_sf(outbox, coords = c("lon", "lat"), crs = st_crs(4326))
outbox3 <- st_combine(outbox2)
outbox4 <- st_cast(outbox3, "POLYGON")

FigS2 <- ggplot() +
  geom_sf(data = subset(world_map_data, world_map_data$continent %in% c("North America", "South America")), lwd = 0.25, fill = "grey75", color = "grey10") +
  geom_sf(data = subset(currents, currents$OBJECTID %in% c(96, 98, 99, 138)), fill = "lightblue1", color = "lightblue3", alpha = 0.5, lwd = 0.25, inherit.aes = FALSE) +
  geom_sf_text(data = subset(currents, currents$OBJECTID %in% c(96, 98, 99, 138)), aes(label = NAME), nudge_x = subset(currents, currents$OBJECTID %in% c(96, 98, 99, 138))$x_nudge, nudge_y = subset(currents, currents$OBJECTID %in% c(96, 98, 99, 138))$y_nudge, size = 6*0.35, family = "Arial", inherit.aes = FALSE) +
  geom_sf(data = comb14c2_sf, aes(fill = study, shape = source), color = "black", alpha = 0.5, size = 2.75, stroke = 0.5, inherit.aes = FALSE) +
  geom_sf(data = outbox4, lwd = 0.3, color = "black", fill = "transparent") +
  coord_sf(xlim = c(-98, -55), ylim = c(8, 32)) +
  theme_minimal() +
  scale_shape_manual(values = c("coral" = 21, "otolith" = 22, "shell" = 23, "surface water" = 24)) +
  labs(x = "Lon.", y = "Lat.", fill = "Study", shape = "Sample type") +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0, nrow = 7, byrow = TRUE, override.aes = list(shape = 22, color = "transparent")),
         shape = guide_legend(title.position = "top", title.hjust = 0, nrow = 7)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.background = element_rect(color = "black", size = 0.5, fill = "white"),
        plot.background = element_rect(color = "transparent", size = 0.5, fill = "white"),
        axis.title = element_text(size = 7, colour = "black", family = "Arial"), 
        axis.text = element_text(size = 7, colour = "black", family = "Arial"), 
        legend.text = element_text(size = 6, colour = "black", family = "Arial"), 
        legend.title = element_text(size = 7, colour = "black", family = "Arial", margin = margin(c(0, 0, 1, 0), unit = "mm")),
        legend.margin=margin(c(0, 1, 0, 0), unit = "cm"),
        legend.spacing.y = unit(0, 'mm'),
        legend.position = "bottom",
        legend.justification = "left",
        legend.spacing.x = unit(1, "mm"))

FigS2 <- FigS2 + 
  ggspatial::annotation_scale(
    location = "bl",
    text_cex = .pt/5.5,
    pad_x = unit(0.05, "in"), pad_y = unit(0.05, "in"),
    line_width = 0.5,
    height = unit(0.05, "in"),
    bar_cols = c("grey60", "white")) +
  ggspatial::annotation_north_arrow(
    location = "bl", which_north = "true",
    pad_x = unit(0.1, "in"), pad_y = unit(0.15, "in"),
    height = unit(0.3, "in"), width = unit(0.3, "in"),
    style = ggspatial::north_arrow_fancy_orienteering(
      fill = c("grey40", "white"),
      line_col = "grey20"))


# ggsave(here::here("Summer2023/Figures/FigS2.pdf"),
#        FigS2,
#        width = 43*(1/6), # 1 pica = ~1/6 inch; 3.5, #column width for 2-column page layout GSA Pub spec
#        height = 5,
#        units = "in",
#        device = grDevices::cairo_pdf)

ggsave(here::here("Summer2023/Figures/FigS2.png"),
       FigS2,
       width = 43*(1/6), # 1 pica = ~1/6 inch; 3.5, #column width for 2-column page layout GSA Pub spec
       height = 6,
       units = "in")


#Figure S2---------------------------------------------------------------------
ggplot(data = comb14c2_sf_mod[inIQR_1yr == TRUE, ], aes(x = year, y = F14C_calc, color = study, shape = source)) +
  geom_ribbon(data = predictedCurve, aes(x = year, ymin = f14c_mn - f14c_sd, ymax = f14c_mn + f14c_sd), alpha = 0.2, inherit.aes = FALSE) +
  geom_ribbon(data = newcurvepreds, aes(x = year, ymin = f14c_mn - f14c_sd, ymax = f14c_mn + f14c_sd), fill = "pink", alpha = 0.75, inherit.aes = FALSE) +
  geom_ribbon(data = newcurvepreds2, aes(x = year, ymin = f14c_mn - f14c_sd, ymax = f14c_mn + f14c_sd), fill = "purple1", alpha = 0.5, inherit.aes = FALSE) +
  geom_point() +
  geom_line(data = predictedCurve, aes(x = year, y = f14c_mn), lwd = 0.75, inherit.aes = FALSE) +
  geom_line(data = newcurvepreds, aes(x = year, y = f14c_mn), color = "red", lwd = 0.75, inherit.aes = FALSE) +
  geom_line(data = newcurvepreds2, aes(x = year, y = f14c_mn), color = "purple", lwd = 0.75, inherit.aes = FALSE) +
  theme_bw() +
  labs(x = "Known-age year", y = expression(paste(F^14, "C")), color = "Study", shape = "Source") +
  coord_cartesian(xlim = c(1750, 2023), ylim = c(0.75, 1.2)) +
  theme(axis.title = element_text(size = 7, colour = "black", family = "Arial"), 
        axis.text = element_text(size = 7, colour = "black", family = "Arial"), 
        legend.text = element_text(size = 7, colour = "black", family = "Arial"), 
        legend.title = element_text(size = 7, colour = "black", family = "Arial"),
        text = element_text(size = 7, colour = "black", family = "Arial"))


modcompare <- ggplot(data = comb14c2_sf_mod[inIQR_3yr == TRUE, ], aes(x = year, y = F14C_calc, color = study, shape = source)) +
  geom_ribbon(data = firstnewcurve, aes(x = year, ymin = estimate - esterror, ymax = estimate + esterror), fill = "green", alpha = 0.75, inherit.aes = FALSE) +
  geom_ribbon(data = newcurvepreds, aes(x = year, ymin = f14c_mn - f14c_sd, ymax = f14c_mn + f14c_sd), fill = "pink", alpha = 0.75, inherit.aes = FALSE) +
  geom_ribbon(data = newcurvepreds2, aes(x = year, ymin = f14c_mn - f14c_sd, ymax = f14c_mn + f14c_sd), fill = "purple1", alpha = 0.5, inherit.aes = FALSE) +
  geom_point() +
  geom_line(data = firstnewcurve, aes(x = year, y = estimate), color = "green4", lwd = 0.75, inherit.aes = FALSE) +
  geom_line(data = newcurvepreds, aes(x = year, y = f14c_mn), color = "red", lwd = 0.75, inherit.aes = FALSE) +
  geom_line(data = newcurvepreds2, aes(x = year, y = f14c_mn), color = "purple", lwd = 0.75, inherit.aes = FALSE) +
  theme_bw() +
  labs(x = "Known-age year", y = expression(paste(F^14, "C")), color = "Study", shape = "Source") +
  coord_cartesian(xlim = c(1750, 2023), ylim = c(0.75, 1.2)) +
  theme(axis.title = element_text(size = 7, colour = "black", family = "Arial"), 
        axis.text = element_text(size = 7, colour = "black", family = "Arial"), 
        legend.text = element_text(size = 7, colour = "black", family = "Arial"), 
        legend.title = element_text(size = 7, colour = "black", family = "Arial"),
        text = element_text(size = 7, colour = "black", family = "Arial"))

FigS2 <- ggplot(data = comb14c2_sf, aes(x = year, y = F14C_calc, fill = study, shape = source)) +
  geom_ribbon(data = newcurvepreds3, aes(x = year, ymin = f14c_mn - f14c_sd, ymax = f14c_mn + f14c_sd), fill = "pink", alpha = 0.4, inherit.aes = FALSE) +
  # geom_ribbon(data = predictedCurve_may, aes(x = `#year`, ymin = F14C - F14Csd, ymax = F14C + F14Csd), fill = "pink", alpha = 0.5, inherit.aes = FALSE) +
  geom_ribbon(data = newcurvepreds3_oto, aes(x = year, ymin = f14c_mn - f14c_sd, ymax = f14c_mn + f14c_sd), fill = "green", alpha = 0.4, inherit.aes = FALSE) +
  geom_ribbon(data = newcurvepreds3_cor, aes(x = year, ymin = f14c_mn - f14c_sd, ymax = f14c_mn + f14c_sd), fill = "lightblue", alpha = 0.75, inherit.aes = FALSE) +
  geom_point(color = "black", alpha = 0.5, size = 2.25, stroke = 0.5) +
  geom_line(data = newcurvepreds3, aes(x = year, y = f14c_mn, color = "Shell"), lwd = 0.75, inherit.aes = FALSE) +
  # geom_line(data = predictedCurve_may, aes(x = `#year`, y = F14C), color = "red", lwd = 0.75, inherit.aes = FALSE) +
  geom_line(data = newcurvepreds3_oto, aes(x = year, y = f14c_mn, color = "Otolith"), lwd = 0.75, inherit.aes = FALSE) +
  geom_line(data = newcurvepreds3_cor, aes(x = year, y = f14c_mn, color = "Coral"), lwd = 0.5, inherit.aes = FALSE) +
  geom_ribbon(data = predictedCurve[year <= 1950, ], aes(x = year, ymin = f14c_mn - f14c_sd, ymax = f14c_mn + f14c_sd), fill = "grey20", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = predictedCurve[year <= 1950, ], aes(x = year, y = f14c_mn, color = "Marine20"), lwd = 0.5, inherit.aes = FALSE) +
  theme_bw() +
  labs(x = "Known-age year", y = expression(paste(F^14, "C")), fill = "Study", shape = "Sample\ntype", color = "Prediction\ntype") +
  coord_cartesian(xlim = c(1750, 2023), ylim = c(0.88, 1.18)) +
  scale_shape_manual(values = c("coral" = 21, "otolith" = 22, "shell" = 23, "surface water" = 24)) +
  scale_color_manual(values = c("Coral" = "dodgerblue", "Otolith" = "green4", "Shell" = "firebrick", "Marine20" = "black")) +
  guides(fill = guide_legend(order = 3, title.position = "top", title.hjust = 0, nrow = 7, byrow = TRUE, override.aes = list(shape = 22, color = "transparent")),
         shape = guide_legend(order = 2, title.position = "top", title.hjust = 0, nrow = 7),
         color = guide_legend(order = 1, title.position = "top", title.hjust = 0, nrow = 7, override.aes = list(linetype = rep("solid", 4), shape = rep(NA, 4)))) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.background = element_rect(color = "black", size = 0.5, fill = "white"),
        plot.background = element_rect(color = "transparent", size = 0.5, fill = "white"),
        axis.title = element_text(size = 7, colour = "black", family = "Arial"), 
        axis.text = element_text(size = 7, colour = "black", family = "Arial"), 
        legend.text = element_text(size = 6, colour = "black", family = "Arial"), 
        legend.title = element_text(size = 7, colour = "black", family = "Arial", margin = margin(c(0, 0, 1, 0), unit = "mm")),
        legend.margin=margin(c(0, 0.25, 0, 0), unit = "cm"),
        legend.spacing.y = unit(0, 'mm'),
        legend.position = "bottom",
        legend.justification = "left",
        legend.spacing.x = unit(1, "mm"))

ggsave(filename = here::here("Summer2023/Figures/FigS2.png"),
       plot = FigS2,
       width = 7,
       height = 6,
       units = "in",
       dpi = 300)


#Sample size analysis------------------------------------------------------------
#Code modified from John Handley's code (individual reef .Rmd files) in "CZ325_Sampling_task_2" directory
pilotsamps <- unique(distinct(hobsf14c2[n_tot > 20 & oor == FALSE, ])$igsn)
ss_dat <- distinct(hobsf14c2[igsn %in% pilotsamps & oor == FALSE, .(loc, reef_num, igsn, analysis, spec_field_id, analysis_id, locality, station, lat, lon, collection_date, depth, cat_no, genus, live_or_dead, ace_file, nau_ace_number, uciams_id, nau_id, f14c, f14c_sd, remove, comment, file_name, depth_int, f14c_corr, f14c_corr_sd)])
ss_post <- all_posteriors[name %in% unique(ss_dat$spec_field_id) & analysis == "LP" & deadcref == "meanref" & probability > 0, ]

ss_results = data.table()
ss_results_sum = data.table()
set.seed(1234)
#WARNING - this next step takes a long time to run!
for(s in pilotsamps){
  ss_post_s <- ss_post[igsn == s, ]
  for( n in 2:length(unique(ss_post_s$name))){
    for( ii in 1:100){
      random_sample = sample(unique(ss_post_s$name), size = n, replace = FALSE)
      # tmp = ss_post_s %>% filter(name %in% random_sample)
      # ss_results_ii <- data.table(igsn = s,
      #                             n = n, 
      #                             cpe = CPE(tmp)$cpe,
      #                             tav_iqr = CPE(tmp)$tav_iqr,
      #                             Median = Median(TAV2(tmp)))
      ss_results_ii <- ss_post_s[name %in% random_sample, .(igsn = unique(igsn),
                                                            n = length(unique(name)),
                                                            cpe = CPE(ss_post_s[name %in% random_sample, ])$cpe,
                                                            tav_iqr = CPE(ss_post_s[name %in% random_sample, ])$tav_iqr,
                                                            Median = Median(TAV2(ss_post_s[name %in% random_sample, ])))]
      ss_results = rbind(ss_results, ss_results_ii)
    }
  }
  # pop_values = last(ss_results_ii)
  # median.tmp = ss_results_ii %>% group_by(n)
  # 
  # cpe.tmp = ss_results_ii %>% group_by(n) %>% dplyr::summarize( cpe_Bias = Bias(cpe, pop_values$cpe), cpe_SD = StandardDeviation(cpe, pop_values$cpe), cpe_CV =  StandardDeviation(cpe, pop_values$cpe)/pop_values$cpe)
  # 
  # tav_iqr.tmp = ss_results_ii %>% group_by(n) %>% dplyr::summarize( tav_iqr_Bias = Bias(tav_iqr, pop_values$tav_iqr), tav_iqr_SD = StandardDeviation(tav_iqr, pop_values$tav_iqr), tav_iqr_CV =  StandardDeviation(tav_iqr, pop_values$tav_iqr)/pop_values$tav_iqr) %>% select(-n)
  # 
  # ss_results_sum_s = round(cbind(cpe.tmp,tav_iqr.tmp), 3)
  # ss_results_sum_s$igsn <- s
  # ss_results_sum <- rbind(ss_results_sum, ss_results_sum_s)
  
  print(paste0("done with ", s))
}

# maxn <- ss_results[, .(nmax = max(n)), by = igsn]
# ss_results_sum <- data.table()
# for(sn in maxn$igsn){
#   nmax <- maxn[igsn == sn, nmax]
#   popvals <- distinct(ss_results[n == nmax & igsn == sn, .(n, cpe, tav_iqr, Median)])
#   ss_results_sum_sn <- ss_results[igsn == sn & !is.na(cpe) & !is.na(tav_iqr) & !is.na(Median), .(igsn = unique(igsn),
#                                                                                                  cpe_Bias = Bias(cpe, popvals$cpe), 
#                                                                                                  cpe_SD = StandardDeviation(cpe, popvals$cpe), 
#                                                                                                  cpe_CV =  StandardDeviation(cpe, popvals$cpe)/popvals$cpe,
#                                                                                                  tav_iqr_Bias = Bias(tav_iqr, popvals$tav_iqr), 
#                                                                                                  tav_iqr_SD = StandardDeviation(tav_iqr, popvals$tav_iqr), 
#                                                                                                  tav_iqr_CV =  StandardDeviation(tav_iqr, popvals$tav_iqr)/popvals$tav_iqr,
#                                                                                                  median_Bias = Bias(Median, popvals$Median),
#                                                                                                  median_SD = StandardDeviation(Median, popvals$Median),
#                                                                                                  median_CV = StandardDeviation(Median, popvals$Median)/popvals$Median), by = n]
#   
#   ss_results_sum <- rbind(ss_results_sum, ss_results_sum_sn)
# }

ss_orig <- data.table()
for(sn in pilotsamps){
  ss_orig_sn <- ss_post[igsn == sn & name %in% ss_dat[is.na(ace_file), unique(spec_field_id)], .(igsn = sn,
                                                                                                 n = length(unique(name)), 
                                                                                                 cpe = CPE(ss_post[igsn == sn & name %in% ss_dat[is.na(ace_file), unique(spec_field_id)], ])$cpe,
                                                                                                 tav_iqr = CPE(ss_post[igsn == sn & name %in% ss_dat[is.na(ace_file), unique(spec_field_id)], ])$tav_iqr,
                                                                                                 Median = Median(TAV2(ss_post[igsn == sn & name %in% ss_dat[is.na(ace_file), unique(spec_field_id)], ])))]
  ss_orig <- rbind(ss_orig, ss_orig_sn)
}

ss_results_meds <- ss_results[, .(Median = median(Median), tav_iqr = median(tav_iqr), cpe = median(cpe)), by = list(igsn, n)]

ss_results2 <- melt(ss_results, id.vars = c("igsn", "n"), measure.vars = c("cpe", "tav_iqr", "Median"), variable.name = "metric", value.name = "value")
ss_results_meds2 <- melt(ss_results_meds, id.vars = c("igsn", "n"), measure.vars = c("cpe", "tav_iqr", "Median"), variable.name = "metric", value.name = "value")

ggplot(data = ss_results, aes(x = n, y = tav_iqr)) +
  geom_jitter(width = 0.2) +
  geom_smooth(color = "orange") +
  theme_bw() +
  facet_wrap(~igsn)

ggplot(data = ss_results, aes(x = n, y = cpe)) +
  geom_jitter(width = 0.2) +
  geom_smooth(color = "orange") +
  theme_bw() +
  facet_wrap(~igsn)

ggplot(data = ss_post[igsn == "IEPRI005Q" & name %in% c("GI-EC_R4H3S1_101", "GI-EC_R4H3S1_44", "GI-EC_R4H3S1_94"), ], aes(x = year, y = probability, color = name)) +
  geom_point() +
  # coord_cartesian(xlim = c(1750, 2023)) +
  facet_wrap(~igsn, scales = "free")


#Figure S11----------------------------------------------------------------------

medplot <- ggplot(data = ss_results, aes(x = n, y = Median)) +
  geom_jitter(width = 0.1, fill = "grey85", color = "grey75", alpha = 0.4, shape = 21, size = 1.5) +
  geom_point(data = ss_results_meds, aes(x = n, y = Median), fill = "dodgerblue", color = "grey20", shape = 23, size = 2) +
  geom_vline(xintercept = 7, color = "black", lwd = 0.5) +
  # geom_smooth(color = "orange") +
  geom_point(data = ss_orig, aes(x = n, y = Median), fill = "firebrick", color = "black", shape = 24, size = 2) +
  labs(x = "Number of specimens", y = "Median sample age (year)") +
  theme_bw() +
  facet_wrap(~igsn, ncol = 1, strip.position = "top", scales = "free_y", labeller = as_labeller(c("IEPRI005E" = "GI-EC Reef 2, IEPRI005E", 
                                                                                                  "IEPRI005Q" = "GI-EC Reef 3, IEPRI005Q", 
                                                                                                  "IEPRI006M" = "PC Reef 3, IEPRI006M", 
                                                                                                  "IEPRI00BE" = "GR Reef 3, IEPRI00BE"))) +
  theme(axis.title = element_text(size = 7, colour = "black", family = "Arial"), 
        axis.text = element_text(size = 7, colour = "black", family = "Arial"), 
        legend.text = element_text(size = 7, colour = "black", family = "Arial"), 
        legend.title = element_text(size = 7, colour = "black", family = "Arial"),
        text = element_text(size = 7, colour = "black", family = "Arial"),
        strip.background = element_rect(color = "transparent", fill = "white"),
        strip.text = element_text(color = "black", hjust = 0, family = "Arial", size = 7))

tavplot <- ggplot(data = ss_results, aes(x = n, y = tav_iqr)) +
  geom_jitter(width = 0.1, fill = "grey85", color = "grey75", alpha = 0.4, shape = 21, size = 1.5) +
  geom_point(data = ss_results_meds, aes(x = n, y = tav_iqr), fill = "dodgerblue", color = "grey20", shape = 23, size = 2) +
  geom_vline(xintercept = 7, color = "black", lwd = 0.5) +
  # geom_smooth(color = "orange") +
  geom_point(data = ss_orig, aes(x = n, y = tav_iqr), fill = "firebrick", color = "black", shape = 24, size = 2) +
  labs(x = "Number of specimens", y = "Sample total age variation (years)") +
  theme_bw() +
  facet_wrap(~igsn, ncol = 1, strip.position = "top", scales = "free_y", labeller = as_labeller(c("IEPRI005E" = "", #GI-EC Reef 2, IEPRI005E
                                                                                                  "IEPRI005Q" = "", #GI-EC Reef 3, IEPRI005Q
                                                                                                  "IEPRI006M" = "", #PC Reef 3, IEPRI006M
                                                                                                  "IEPRI00BE" = ""))) + #GR Reef 3, IEPRI00BE
  theme(axis.title = element_text(size = 7, colour = "black", family = "Arial"), 
        axis.text = element_text(size = 7, colour = "black", family = "Arial"), 
        legend.text = element_text(size = 7, colour = "black", family = "Arial"), 
        legend.title = element_text(size = 7, colour = "black", family = "Arial"),
        text = element_text(size = 7, colour = "black", family = "Arial"),
        strip.background = element_rect(color = "transparent", fill = "white"),
        strip.text = element_text(color = "black", hjust = 0, family = "Arial", size = 7))

# cpeplot <- ggplot(data = ss_results, aes(x = n, y = cpe)) +
#   geom_jitter(width = 0.1, fill = "grey85", color = "grey75", shape = 21, size = 1.5) +
#   geom_point(data = ss_results_meds, aes(x = n, y = Median, group = igsn), fill = "dodgerblue", color = "grey20", shape = 23, size = 2) +
#   geom_vline(xintercept = 7, color = "black", lwd = 0.5) +
#   geom_smooth(color = "orange") +
#   labs(x = "Number of specimens", y = "Corrected posterior age estimate (years)") +
#   theme_bw() +
#   facet_wrap(~igsn, ncol = 1, strip.position = "right")


FigS6 <- medplot + tavplot + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")

ggsave(filename = here::here("Summer2023/Figures/FigS6.png"),
       plot = FigS6,
       width = 6,
       height = 8,
       units = "in",
       dpi = 300)

