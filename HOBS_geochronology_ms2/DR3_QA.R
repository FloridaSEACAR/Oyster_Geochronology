
#This script shows the discrepancy in corrected F14C and SD values from an error in an Excel formula I made in 
#the formulas caluclating these values for Guana River, Jack Island, Goose Island East Cove, and Lone Cabbage localities 
#in '14C results & Posteriors_LP_Marine20_Feb 2022.xlsx' ('All regions' tab). The values are re-calculated and graphed
#for comparison with the 'All regions' values and then I show that the recalculated values match the values in the 
#source tab for the Guana River data (i.e., the 'NE regions' tab), demonstrating that the problem was isolated to the 
#aggregation of regional data (meaning the correct values were used in Quan's calibrations and this issue is isolated
#to my work composing the data appendices).
#
#Stephen Durham, 03/10/2023

library(tidyverse)
library(data.table)
library(openxlsx)
library(sf)
library(metaDigitise)
library(rintcal)
library(ggpmisc)
library(brms)
# library(plotKML)
library(mapview)
mapviewOptions(fgb = FALSE)


dr3 <- read.xlsx(here::here("HOBS_geochronology_ms2/HOBSGeo_Apps_202208_final.xlsx"), sheet = "DR3", startRow = 32)
setDT(dr3)

ggplot(data = dr3, aes(x = Locality, y = F14C_corr_sd, color = Locality)) +
  geom_jitter()

ggplot(data = dr3, aes(x = Locality, y = F14C_sd, color = Locality)) +
  geom_jitter()

dr3[, corr2 := F14C/(1 - Dead_C)]
dr3[, corrsd2 := sqrt(((F14C_sd/(1 - Dead_C))^2) + (Dead_C_sd * (F14C/(1 - Dead_C)^2))^2)]

ggplot(data = dr3, aes(x = Locality, y = corrsd2, color = Locality)) +
  geom_jitter() +
  coord_cartesian(ylim = c(0, max(dr3$F14C_corr_sd))) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

dr3[, `:=` (correq = ifelse(F14C_corr == corr2, 1, 0), corrsdeq = ifelse(F14C_corr_sd == corrsd2, 1, 0)), by = row.names(dr3)]

dr3[correq == 0, unique(Locality)]
dr3[correq == 1, unique(Locality)]

dr3[corrsdeq == 0, unique(Locality)]
dr3[corrsdeq == 1, unique(Locality)]


neregtest <- read.xlsx("C:/Users/durham_s/OneDrive - Florida Department of Environmental Protection/Documents/Historical oyster body size project/Geochronology PO/Data/14C results & Posteriors_LP_Marine20_Feb 2022.xlsx", sheet = "NE Region", startRow = 2)
neregtest <- neregtest[2:nrow(neregtest),]
setDT(neregtest)

gr_dr3 <- dr3[Locality == "Guana River", ]
gr_nrt <- neregtest[Locality == "GR", c(1:10, 13:17)]
gr_nrt[, X14 := as.numeric(X14)]

grtest <- merge(gr_dr3[, .(Locality, NAU_ID, F14C, F14C_sd, F14C_corr, F14C_corr_sd, corr2, corrsd2)],
                gr_nrt[, .(Sample.name, F14C, `±`, `Dead.C-corrected.F14C`, X14)], by.x = c("NAU_ID"), by.y = c("Sample.name"))

grtest[, `:=` (F14C_eq = F14C.x == F14C.y,
               F14C_sd_eq = F14C_sd == `±`,
               F14C_corr_eq = F14C_corr == `Dead.C-corrected.F14C`,
               F14C_corr_sd_eq = F14C_corr_sd == X14,
               corr2_eq = corr2 == `Dead.C-corrected.F14C`,
               corrsd2_eq = corrsd2 == X14)]

unique(grtest$F14C_eq) #all TRUE
unique(grtest$F14C_sd_eq) #all TRUE
unique(grtest$F14C_corr_eq) #all FALSE
unique(grtest$F14C_corr_sd_eq) #all FALSE
unique(grtest$corr2_eq) #all TRUE
unique(grtest$corrsd2_eq) #Some TRUE, some FALSE, but all FALSE are the same through at least 9 decimal places


#Which localities have out-of-range F14C and/or F14C_corr?
oor <- dr3[F14C > 1.1588 | corr2 > 1.1588, ]
table(dr3[F14C > 1.1588, Locality])
table(dr3[corr2 > 1.1588, Locality])
oor_sum <- data.table(Locality = unique(oor$Locality),
                      F14C = table(dr3[F14C > 1.1588, Locality]),
                      F14C_corr = table(dr3[corr2 > 1.1588, Locality]))
oor_sum[, `:=` (F14C.V1 = NULL, F14C_corr.V1 = NULL)]
setnames(oor_sum, c("F14C.N", "F14C_corr.N"), c("N_F14C", "N_F14C_corr"))
fwrite(oor_sum, here::here(paste0("OutOfRangeSpecSum_", Sys.Date(), ".csv")))

#Closely look at which reefs/samples had specimens with out-of-range results
dr3[, `:=` (oor_F14C = ifelse(F14C > 1.1588, 1, 0), oor_F14C_corr = ifelse(corr2 > 1.1588, 1, 0))]
xwalk <- fread(here::here("HOBS_geochronology_ms2/SpecNames_Xwalk.csv"))
dr3n <- merge(dr3, xwalk[, .(Sample_ID, Cat_No)])
dr3n[, reef := str_sub(Sample_ID, 1, str_locate(Sample_ID, "_")[1] + 2), by = Sample_ID]
dr3n[, reef2 := reef]

reefsum <- distinct(dr3n[, `:=` (n_tot = .N,
                                 n_oor_F14C = sum(oor_F14C),
                                 n_oor_F14C_corr = sum(oor_F14C_corr)), by = reef][, .(Locality, reef, n_tot, n_oor_F14C, n_oor_F14C_corr)][, `:=` (prop_oor_F14C = n_oor_F14C/n_tot, prop_oor_F14C_corr = n_oor_F14C_corr/n_tot)])
setorder(reefsum, Locality, reef)

dr3n[, reefn := which(sort(unique(dr3n$reef)) == reef2), by = reef2]

ggplot(dr3n) +
  geom_vline(xintercept = 1, color = "black", lwd = 0.75) +
  geom_vline(xintercept = 1.1588, color = "firebrick", lwd = 0.75) +
  geom_point(aes(x = F14C, y = reefn + 0.15, color = Locality, shape = Depth), alpha = 0.8, size = 1.5) +
  geom_point(aes(x = corr2, y = reefn - 0.15, color = Locality, shape = Depth), alpha = 0.4, size = 1.5) +
  theme_bw() +
  scale_y_continuous(breaks = unique(dr3n$reefn), labels = unique(dr3n$reef2), trans = "reverse") +
  labs(y = "Reef", x = "F14C (alpha = 0.8)\nF14C_corr (alpha = 0.4)", shape = "Burial depth")


#Need to identify potential lots that could help us flesh out the 1965-1975 portion of the bomb-pulse curve and/or our Dead C corrections

pri <- fread("C:/Users/durham_s/Downloads/PRI_Class_Bival.csv")
colnames(pri)
pri2 <- pri[, .(PRICatalogNumber, PRIStationNumber, IGSNNumber, Count, Order, Family, Genus, Species, Subspecies, County, Latitude, Longitude)]
setnames(pri2, c("PRICatalogNumber", "PRIStationNumber", "Count"), c("CatNum", "Station", "Preparation"))
pri2[, Institution := "PRI"]
pri2[, `:=` (Station = as.character(Station), Preparation = as.character(Preparation))]

bmsm <- fread("C:/Users/durham_s/Downloads/BMSM_Class_Bival.csv")
colnames(bmsm)
bmsm2 <- bmsm[, .(BMSMNo, Order, Family, Genus, Subgenus, Species, Subspecies, County, Site, Latitude1, Longitude1, preparations, StartDate, DepthMIN)]
setnames(bmsm2, c("BMSMNo", "Site", "preparations", "Latitude1", "Longitude1", "DepthMIN"), c("CatNum", "Station", "Preparation", "Latitude", "Longitude", "Depth"))
bmsm2[, Institution := "BMSM"]
bmsm2[, `:=` (Station = as.character(Station), Preparation = as.character(Preparation))]

flmnh <- fread("C:/Users/durham_s/Downloads/FLMNH_Class_Bival.csv")
colnames(flmnh)
flmnh2 <- flmnh[, .(UFID, Order, Family, Genus, Species, Subspecies, Preparations, County, Locality, Latitude, Longitude, StartDate, Depth1, DepthUnit)]
setnames(flmnh2, c("UFID", "Locality", "Preparations", "Depth1"), c("CatNum", "Station", "Preparation", "Depth"))
flmnh2[, Institution := "FLMNH"]
flmnh2[, `:=` (Station = as.character(Station), Preparation = as.character(Preparation))]

fwc <- fread("C:/Users/durham_s/Downloads/FWC_Class_Bival.csv")
colnames(fwc)
fwc2 <- fwc[, .(CatalogFSBCI, Preparations, Order, Family, Genus, Species, Subspecies, Date, Verbatimlocality, Latitude, Longitude, Minimumdepth)]
setnames(fwc2, c("CatalogFSBCI", "Verbatimlocality", "Preparations", "Minimumdepth", "Date"), c("CatNum", "Station", "Preparation", "Depth", "StartDate"))
fwc2[, Institution := "FWC"]
fwc2[, `:=` (Station = as.character(Station), Preparation = as.character(Preparation))]

allm <- bind_rows(pri2, bmsm2, flmnh2, fwc2)
setcolorder(allm, c("Institution", "IGSNNumber", "CatNum", "County", "StartDate", "Preparation", "Depth", "DepthUnit", "Station", "Order", "Family", "Genus", "Subgenus", "Species", "Subspecies", "Latitude", "Longitude"))
allm2 <- distinct(allm[!is.na(Longitude) & !is.na(Latitude), ])
# allm3 <- st_as_sf(allm2, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

rm(pri, pri2, bmsm, bmsm2, flmnh, flmnh2, fwc, fwc2)

cv <- allm2[Genus == "Crassostrea" & Species == "virginica", ]
vens <- allm2[Family == "Veneridae" & Genus %in% c("Mercenaria", "Chione"), ]
old <- allm2[year(StartDate) %in% c(1965:1975), ]
new <- allm2[year(StartDate) %in% c(2018:2023), ]

hobs_sta <- distinct(dr3[, .(Station, Lat, Lon)])
hobs_sta <- st_as_sf(hobs_sta, coords = c("Lon", "Lat"), crs = 4326)
buffs <- st_buffer(hobs_sta, 2000)

cv2 <- distinct(cv[!is.na(Longitude) & !is.na(Latitude), ])
cv3 <- st_as_sf(cv2, coords = c("Longitude", "Latitude"), crs = 4326)
for(i in seq(1, nrow(cv2))){
  cv2[i, nearHOBS := ifelse(TRUE %in% st_intersects(buffs, cv3$geometry[i], sparse = FALSE), 1, 0)]
}
# cv2_hobs <- st_intersects(buffs, cv2)
# cv2_hobs <- cv2[which(st_contains(buffs, cv2, sparse = FALSE)),]
# cv2[buffs, , op = st_intersects]

View(subset(cv2_hobs, cv2_hobs$CatNum %in% cv2_hobs$CatNum[which(cv2_hobs$CatNum %in% setdiff(unique(cv2_hobs$CatNum), unique(str_sub(dr3$Cat_No, 4, -1))))]))
mapview::mapview(subset(cv2_hobs, cv2_hobs$CatNum %in% cv2_hobs$CatNum[which(cv2_hobs$CatNum %in% setdiff(unique(cv2_hobs$CatNum), unique(str_sub(dr3$Cat_No, 4, -1))))]))

vens2 <- distinct(vens[!is.na(Longitude) & !is.na(Latitude), ])
vens3 <- st_as_sf(vens2, coords = c("Longitude", "Latitude"), crs = 4326)
for(i in seq(1, nrow(vens2))){
  vens2[i, nearHOBS := ifelse(TRUE %in% st_intersects(buffs, vens3$geometry[i], sparse = FALSE), 1, 0)]
}
unique(vens2$nearHOBS)

mapview::mapview(subset(vens3, vens3$CatNum %in% vens2[nearHOBS == 1, unique(CatNum)]), col.regions = c("Mercenaria" = "blue", "Chione" = "red"))
View(subset(vens3, vens3$CatNum %in% vens2[nearHOBS == 1, unique(CatNum)]))


old2 <- distinct(old[str_detect(Preparation, "et", negate = TRUE) & Family != "", ])
old3 <- st_as_sf(old2, coords = c("Longitude", "Latitude"), crs = 4326)
for(i in seq(1, nrow(old2))){
  old2[i, nearHOBS := ifelse(TRUE %in% st_intersects(buffs, old3$geometry[i], sparse = FALSE), 1, 0)]
}
unique(old2$nearHOBS)

#which Families have at least 50 lots?
old2[, nFam := .N, by = Family]
o50_old <- old2[nFam > 50, unique(Family)]
o50_old_use <- setdiff(o50_old, c("Lucinidae", "Nuculanidae", "Propeamussiidae", "Unionidae"))
old2[!is.na(Depth) & DepthUnit != "", Depth_m := fcase(DepthUnit == "fathoms", Depth * 1.8288,
                                                       DepthUnit == "feet", Depth / 3.28084,
                                                       DepthUnit == "meters", Depth)]
old2b <- old2[Family %in% o50_old_use & (is.na(Depth) | Depth_m < 50), ]
old2b[!is.na(Genus) & !is.na(Species), fullName := paste0(Genus, " ", Species)]
table(old2b$fullName)
old2c <- old2b[str_detect(fullName, "Undet.", negate = TRUE) & Species != "" & Species != "sp.", ]
oldspo5 <- as.data.table(table(old2c$fullName)[which(table(old2c$fullName) > 5)])
old2c <- old2c[fullName %in% oldspo5$V1, ]

#Generate filtered lists for each museum

uflist <- old2c[Institution == "FLMNH", ]
setorder(uflist, fullName, StartDate)
uflist[, `:=` (IGSNNumber = NULL, Subgenus = NULL, nFam = NULL)]
setcolorder(uflist, c("Institution", "fullName", "CatNum", "StartDate", "County", "Preparation", "Depth", "DepthUnit", "Depth_m", "Station", "Order", "Family",
                      "Genus", "Species", "Subspecies", "Latitude", "Longitude"))
write.xlsx(uflist, here::here(paste0("CalCurve/DurhamSR_FLBivalves_oldUF_", Sys.Date(), ".xlsx")), sheetName = "FLMNH", colNames = TRUE, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)

bmlist <- old2c[Institution == "BMSM", ]
setorder(bmlist, fullName, StartDate)
bmlist[, `:=` (IGSNNumber = NULL, Subgenus = NULL, nFam = NULL)]
setcolorder(bmlist, c("Institution", "fullName", "CatNum", "StartDate", "County", "Preparation", "Depth", "DepthUnit", "Depth_m", "Station", "Order", "Family",
                      "Genus", "Species", "Subspecies", "Latitude", "Longitude"))
write.xlsx(bmlist, here::here(paste0("CalCurve/DurhamSR_FLBivalves_oldBM_", Sys.Date(), ".xlsx")), sheetName = "BMSM", colNames = TRUE, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)

fwclist <- old2c[Institution == "FWC", ]
setorder(fwclist, fullName, StartDate)
fwclist[, `:=` (IGSNNumber = NULL, Subgenus = NULL, nFam = NULL)]
setcolorder(fwclist, c("Institution", "fullName", "CatNum", "StartDate", "County", "Preparation", "Depth", "DepthUnit", "Depth_m", "Station", "Order", "Family",
                      "Genus", "Species", "Subspecies", "Latitude", "Longitude"))
write.xlsx(fwclist, here::here(paste0("CalCurve/DurhamSR_FLBivalves_oldFWC_", Sys.Date(), ".xlsx")), sheetName = "FWC", colNames = TRUE, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)

mapview::mapview(subset(old3, old3$Family %in% o50_old & old3$CatNum %in% old2[nearHOBS == 1, unique(CatNum)]))
View(subset(old3, old3$Family %in% o50_old & old3$CatNum %in% old2[nearHOBS == 1, unique(CatNum)]))



allm2[!is.na(Depth) & DepthUnit != "", Depth_m := fcase(DepthUnit %in% c("fathoms", "fathom"), Depth * 1.8288,
                                                        DepthUnit %in% c("feet", "feet ", "ft", "foot"), Depth / 3.28084,
                                                        DepthUnit %in% c("meters", "meter"), Depth,
                                                        DepthUnit == "inches", Depth * 39.3701)]
allm2[, nFam := .N, by = Family]
o50_allm2 <- allm2[nFam > 50, unique(Family)]
allm3 <- distinct(allm2[str_detect(Preparation, "et", negate = TRUE) & !is.na(Longitude) & !is.na(Latitude) & (is.na(Depth) | Depth_m < 50) & Family %in% o50_allm2 & Family != "" & !is.na(StartDate), ])
allm3 <- allm3[Family %in% setdiff(unique(allm3$Family), c("Lucinidae", "Nuculanidae", "Propeamussiidae", "Unionidae")), ]
allm3[!is.na(Genus) & !is.na(Species), fullName := paste0(Genus, " ", Species)]
table(allm3$fullName)
allm3 <- allm3[str_detect(fullName, "Undet.", negate = TRUE) & Species != "" & Species != "sp.", ]
allm3spo5 <- as.data.table(table(allm3$fullName)[which(table(allm3$fullName) > 5)])
allm3 <- allm3[fullName %in% allm3spo5$V1, ]
setorder(allm3, fullName, StartDate)
allm3[, `:=` (Subgenus = NULL)]
setcolorder(allm3, c("Institution", "fullName", "IGSNNumber", "CatNum", "StartDate", "County", "Preparation", "Depth", "DepthUnit", "Depth_m", "Station", "Order", "Family",
                     "Genus", "Species", "Subspecies", "Latitude", "Longitude"))
new <- distinct(allm3[year(StartDate) %in% c(2018:2023), ])



#Are any of the other lots from similar locations to the "new" lots of same Family?
#Use only families that show up in both the older and newest lots
# sharedFams <- intersect(unique(old$Family), unique(new$Family))
# new_sta <- distinct(new2[, .(CatNum, Family, Station, Latitude, Longitude)])
new_sta <- distinct(new[, .(Latitude, Longitude)])
new_sta <- st_as_sf(new_sta, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)
buffs_new <- st_buffer(new_sta, 1500)
buffs_new$ind <- seq(1, nrow(buffs_new))
allm4 <- st_as_sf(allm3, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

# setnames(buffs_new, c("CatNum", "Family"), c("CatNo", "Fam"))

# mapview::mapview(buffs_new) +
#   mapview::mapview(subset(allm3, allm3$CatNum %in% allm2[Family %in% sharedFams, unique(CatNum)]))

# setnames(allm3, c("CatNum", "Family"), c("CatNo", "Fam"))

# allm2[, `:=` (nearnew = ifelse(TRUE %in% st_intersects(buffs_new$geometry[which(buffs_new$Fam == Family)], allm3$geometry[which(allm3$CatNo == CatNum)], sparse = FALSE), 1, 0),
#               nearcats = list(list(subset(buffs_new, buffs_new$Fam == Family)$CatNo[which(st_intersects(buffs_new$geometry[which(buffs_new$Fam == Family)], 
#                                                                                                         allm3$geometry[which(allm3$CatNo == CatNum)], sparse = FALSE) == TRUE)]))), 
#                    by = row.names(allm2)]

newlocs <- allm3[0]
newlocs[, `:=` (indval = integer(), nearnew = NULL, nearcats = NULL)]
for(x in seq(1, nrow(buffs_new))){
  # new_x <- allm3[buffs_new[x,], -c(which(colnames(allm3) %in% c("nearnew", "nearcats"))), op = st_intersects]
  new_x <- allm4[buffs_new[x,], , op = st_intersects]
  new_x <- st_drop_geometry(new_x)
  # setnames(new_x, c("Fam", "CatNo"), c("Family", "CatNum"))
  # setnames(new_x, "CatNo", "CatNum")
  # new_x$Latitude <- newlocs$Latitude[x]
  # new_x$Longitude <- newlocs$Longitude[x]
  new_x$indval <- x
  # new_x$nearcats <- newlocs$nearcats[x]
  # new_x$nearnew <- newlocs$nearnew[x]
  newlocs <- rbind(newlocs, new_x)
}

newlocs_dry <- newlocs[str_detect(Preparation, "ry") | Preparation == "", ]
newlocs_dry[, Year := year(StartDate)]
saveRDS(newlocs_dry, here::here(paste0("CalCurve/DryPrep_NewLocs_wOldColls_", Sys.Date(), ".rds")))

newlocs_dry_fam <- distinct(newlocs_dry[, .(.N), by = c("indval", "Family", "Year")])
setorder(newlocs_dry_fam, indval, Family, Year)
newlocs_dry_fam_sum <- distinct(newlocs_dry_fam[!is.na(Year), .(minYear = min(Year), maxYear = max(Year), N_yr = .N, N_spec = sum(N)), by = c("indval", "Family")])
newlocs_dry_fam_sum2 <- newlocs_dry_fam_sum[minYear <= 1980 & maxYear >= 2018 & N_yr >= 3, ]

newlocs_dry_gen <- distinct(newlocs_dry[, .(.N), by = c("indval", "Family", "Genus", "Year")])
setorder(newlocs_dry_gen, indval, Family, Genus, Year)
newlocs_dry_gen_sum <- distinct(newlocs_dry_gen[!is.na(Year), .(minYear = min(Year), maxYear = max(Year), N_yr = .N, N_spec = sum(N)), by = c("indval", "Family", "Genus")])
newlocs_dry_gen_sum2 <- newlocs_dry_gen_sum[minYear <= 1980 & maxYear >= 2018 & N_yr >= 3, ]

newlocs_dry_sp <- distinct(newlocs_dry[, .(.N), by = c("indval", "Family", "Genus", "Species", "fullName", "Year")])
setorder(newlocs_dry_sp, fullName, indval, Family, Genus, Species, Year)
newlocs_dry_sp_sum <- distinct(newlocs_dry_sp[!is.na(Year), .(minYear = min(Year), maxYear = max(Year), N_yr = .N, N_spec = sum(N)), by = c("indval", "Family", "Genus", "Species", "fullName")])
newlocs_dry_sp_sum2 <- newlocs_dry_sp_sum[minYear <= 1980 & maxYear >= 2018 & N_yr >= 3, ]

length(unique(newlocs_dry_sp_sum2$indval))
fwrite(newlocs_dry_sp_sum2, here::here(paste0("newlocs_dry_sp_sum2_", Sys.Date(), ".csv")))

mapview::mapview(subset(buffs_new, buffs_new$ind %in% unique(newlocs_dry_fam_sum2$indval)), col.regions = "green", layer.name = "Family") +
  mapview::mapview(subset(buffs_new, buffs_new$ind %in% unique(newlocs_dry_gen_sum2$indval)), col.regions = "firebrick", layer.name = "Genus") +
    mapview::mapview(subset(buffs_new, buffs_new$ind %in% unique(newlocs_dry_sp_sum2$indval)), col.regions = "dodgerblue", layer.name = "Species")


#Generate filtered lists for each museum

uflist2 <- newlocs_dry[Institution == "FLMNH" & fullName %in% unique(newlocs_dry_sp_sum2$fullName) & indval %in% unique(newlocs_dry_sp_sum2$indval), ]
setorder(uflist2, fullName, indval, StartDate)
uflist2[, `:=` (IGSNNumber = NULL, nFam = NULL)]
setcolorder(uflist2, c("Institution", "fullName", "CatNum", "indval", "Year", "StartDate", "County", "Preparation", "Depth", "DepthUnit", "Depth_m", "Station", "Order", "Family",
                       "Genus", "Species", "Subspecies", "Latitude", "Longitude"))
write.xlsx(uflist2, here::here(paste0("CalCurve/DurhamSR_FLBivalves_UF_reps_1500m_", Sys.Date(), ".xlsx")), sheetName = "FLMNH", colNames = TRUE, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)

bmlist2 <- newlocs_dry[Institution == "BMSM" & fullName %in% unique(newlocs_dry_sp_sum2$fullName) & indval %in% unique(newlocs_dry_sp_sum2$indval), ]
setorder(bmlist2, fullName, indval, StartDate)
bmlist2[, `:=` (IGSNNumber = NULL, Subspecies = NULL, nFam = NULL)]
setcolorder(bmlist2, c("Institution", "fullName", "CatNum", "indval", "Year", "StartDate", "County", "Preparation", "Depth", "DepthUnit", "Depth_m", "Station", "Order", "Family",
                       "Genus", "Species", "Latitude", "Longitude"))
write.xlsx(bmlist2, here::here(paste0("CalCurve/DurhamSR_FLBivalves_BM_reps_1500m_", Sys.Date(), ".xlsx")), sheetName = "BMSM", colNames = TRUE, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)

fwclist2 <- newlocs_dry[Institution == "FWC" & fullName %in% unique(newlocs_dry_sp_sum2$fullName) & indval %in% unique(newlocs_dry_sp_sum2$indval), ]
setorder(fwclist2, fullName, indval, StartDate)
fwclist2[, `:=` (IGSNNumber = NULL, Subspecies = NULL, nFam = NULL)]
setcolorder(fwclist2, c("Institution", "fullName", "CatNum", "indval", "Year", "StartDate", "County", "Preparation", "Depth", "DepthUnit", "Depth_m", "Station", "Order", "Family",
                       "Genus", "Species", "Latitude", "Longitude"))
write.xlsx(fwclist2, here::here(paste0("CalCurve/DurhamSR_FLBivalves_FWC_reps_1500m_", Sys.Date(), ".xlsx")), sheetName = "FWC", colNames = TRUE, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)



#Extract delta14C info from figures---------------------------

data <- metaDigitise(dir = here::here("CalCurve/")) 
data2 <- metaDigitise(dir = here::here("CalCurve/"), summary = FALSE) 
druffel1980 <- setDT(data2$scatterplot$Druffel1980_Fig2.png)
MoyerandGrottoli2011 <- setDT(data2$scatterplot$MoyerandGrottoli2011_Fig3.png)

druffel1980[, nearestYear := round(x)]
druffel1980[c(18, 19, 21), nearestYear := c(1955, 1956, 1958)] #Correct years where it looks like the rounding is wrong, based on looking at the figure.
fwrite(druffel1980, here::here(paste0("CalCurve/druffel1980_Fig2_metaDigitize_", Sys.Date(), ".csv")))

MoyerandGrottoli2011[, nearestYear := round(x)]
fwrite(MoyerandGrottoli2011, here::here(paste0("MoyerandGrottoli2011_Fig3_metaDigitize_", Sys.Date(), ".csv")))

comb14c <- fread(here::here("CalCurve/14Cdata_FL_Gulf_Caribbean_20230324.csv"))

ggplot(comb14c[!is.na(delta14C), ]) +
  geom_point(aes(x = year, y = delta14C, color = study)) +
  geom_smooth(aes(x = year, y = delta14C), color = "black", method = "loess") +
  theme_bw()


#Compiled calibration curve data--------------------------

comb14c <- fread(here::here("CalCurve/14Cdata_FL_Gulf_Caribbean_20230329.csv"))
comb14c[lon > 1, lon := lon * -1]
comb14c_sf <- st_as_sf(comb14c[!is.na(lon) & !is.na(lat), ], coords = c("lon", "lat"), crs = 4326)
mapview::mapview(comb14c_sf)

colnames(comb14c)

ggplot(comb14c[!is.na(delta14C) & !is.na(year), ]) +
  geom_point(aes(x = year, y = delta14C, color = study, shape = source)) +
  theme_bw()

comb14c[, ind := seq(1, nrow(comb14c))]

comb14c[!is.na(`14Cage`), `:=` (F14C_calc = age.F14C(`14Cage`, sdev = `14Cage_sd`)[,1],
                                F14C_calc_sd = age.F14C(`14Cage`, sdev = `14Cage_sd`)[,2])]

comb14c[!is.na(pMC), `:=` (C14age_calc = pMC.age(pMC, sdev = pMC_sd)[[1]][1], 
                           C14age_sd_calc = pMC.age(pMC, sdev = pMC_sd)[[1]][2]), by = ind] #comb14c[!is.na(pMC), ind]

comb14c[, year_bp := 2023 - year]

comb14c[!is.na(delta14C) & !is.na(year_bp) & is.na(F14C_calc), F14C_calc := D14C.F14C(delta14C, year_bp)]

comb14c[!is.na(C14age_calc) & is.na(F14C_calc), `:=` (F14C_calc = age.F14C(C14age_calc, sdev = C14age_sd_calc)[,1],
                                                      F14C_calc_sd = age.F14C(C14age_calc, sdev = C14age_sd_calc)[,2])]
comb14c[!is.na(year) & !is.na(F14C) & is.na(F14C_calc), `:=` (F14C_calc = F14C, F14C_calc_sd = F14C_sd)]

comb14c[!is.na(year) & !is.na(F14C) & !is.na(F14C_calc), delta14C_calc := F14C.D14C(F14C_calc, year_bp)]

ggplot(comb14c[!is.na(year) & !is.na(F14C_calc), ]) +
  geom_point(aes(x = year, y = F14C_calc, color = study, shape = source)) +
  geom_point(data = comb14c[!is.na(year) & !is.na(F14C_calc) & str_detect(species, "Mercenaria sp"), ], aes(x = year, y = F14C_calc), shape = 1, color = "black", stroke = 1, size = 1.5) +
  theme_bw()

# ggplot(comb14c[!is.na(year) & !is.na(F14C) & !is.na(F14C_calc), ]) +
#   geom_point(aes(x = delta14C, y = delta14C_calc, color = study)) +
#   geom_smooth(aes(x = delta14C, y = delta14C_calc), method = "lm", color = "black", lwd = 0.5) +
#   theme_bw() +
#   coord_cartesian(ylim = c(-55, -30))

ggplot(comb14c[!is.na(year) & !is.na(F14C) & !is.na(F14C_calc), ], aes(x = delta14C, y = delta14C_calc, color = study)) +
  geom_point() +
  stat_poly_line(color = "black", lwd = 0.5) +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")), color = "black") +
  theme_bw() +
  coord_cartesian(ylim = c(-55, -30))

ggplot(comb14c[!is.na(year) & !is.na(F14C) & !is.na(F14C_calc), ], aes(x = F14C, y = F14C_calc, color = study)) +
  geom_point() +
  stat_poly_line(color = "black", lwd = 0.5) +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")), color = "black") +
  theme_bw() #+
  # coord_cartesian(ylim = c(-55, -30))


gamtest <- brm(bf(F14C_calc ~ s(year)), data = comb14c[!is.na(year), ], family = gaussian(), 
               cores = 4, seed = 456, chains = 4, iter = 5000, warmup = 1500, thin = 3, 
               control = list(adapt_delta = 0.99, max_treedepth = 15), backend = "cmdstanr", threads = threading(2))
gamtest_ms <- conditional_smooths(gamtest)
plot(gamtest_ms)

ggplot(comb14c[!is.na(year), ]) +
  geom_ribbon(data = gamtest_ms$mu, aes(x = year, ymin = lower__ + 1, ymax = upper__ + 1), alpha = 0.2) +
  geom_point(aes(x = year, y = F14C_calc, color = study, shape = source)) +
  geom_point(data = comb14c[!is.na(year) & !is.na(F14C_calc) & str_detect(species, "Mercenaria sp"), ], aes(x = year, y = F14C_calc), shape = 1, color = "black", stroke = 1, size = 1.5) +
  geom_line(data = gamtest_ms$mu, aes(x = year, y = estimate__ + 1)) +
  geom_point(data = comb14c[!is.na(year) & study == "Andrews et al. 2020" & F14C_calc <= 0.94, ], aes(x = year, y = F14C_calc), shape = 2, size = 1) +
  theme_bw()


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

# comb14c[, decade := fcase(year < 1900, 1,
#                           year >= 1900 & year < 1910, 3,
#                           year >= 1910 & year < 1920, 5,
#                           year >= 1920 & year < 1930, 7,
#                           year >= 1930 & year < 1940, 9,
#                           year >= 1940 & year < 1950, 11,
#                           year >= 1950 & year < 1960, 13,
#                           year >= 1960 & year < 1970, 15,
#                           year >= 1970 & year < 1980, 17,
#                           year >= 1980 & year < 1990, 19,
#                           year >= 1990 & year < 2000, 21,
#                           year >= 2000 & year < 2010, 23,
#                           year >= 2010 & year < 2020, 25,
#                           year >= 2020 & year < 2030, 27)]

comb14c_sf2 <- st_as_sf(comb14c[!is.na(year) & !is.na(lon) & !is.na(lat), ], coords = c("lon", "lat"), crs = 4326)

currents <- st_read(here::here("CalCurve/MajorOceanCurrents/Major_Ocean_Currents.shp"))
currents <- st_make_valid(currents)
currents <- st_transform(currents, crs = 4326)
# currents2 <- st_read("C:/Users/durham_s/Downloads/export25.kml")
# currents2 <- kml_open("C:/Users/durham_s/Downloads/export25.kml")

mapview::mapview(currents) +
  mapview::mapview(subset(comb14c_sf2, !is.na(comb14c_sf2$year)), zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), label = "sta", layer.name = "F14C_calc", alpha = 0, alpha.regions = 0, legend = TRUE, popup = FALSE, highlight = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year < 1900), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "pre-1900", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1900 & year < 1910), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1900s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1910 & year < 1920), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1910s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1920 & year < 1930), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1920s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1930 & year < 1940), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1930s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1940 & year < 1950), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1940s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1950 & year < 1960), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1950s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1960 & year < 1970), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1960s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1970 & year < 1980), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1970s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1980 & year < 1990), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1980s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 1990 & year < 2000), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "1990s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 2000 & year < 2010), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "2000s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 2010 & year < 2020), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "2010s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c_sf2, !is.na(comb14c_sf2$year) & year >= 2020 & year < 2030), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(0.75, 1.35, 0.1), layer.name = "2020s", legend = FALSE)


#Remove points north of the southern Bahamas and east of the Caribbean Sea from the data, based on the ocean currents
comb14c2 <- comb14c[!is.na(year) & lon < -77.87 | lat < 26.43 & lon < -57, ]
comb14c2_sf <- st_as_sf(comb14c2[!is.na(year) & !is.na(lon) & !is.na(lat), ], coords = c("lon", "lat"), crs = 4326)
minF14C <- plyr::round_any(range(comb14c2_sf$F14C_calc)[1], 0.05, f = floor)
maxF14C <- plyr::round_any(range(comb14c2_sf$F14C_calc)[2], 0.05, f = ceiling)

mapview::mapview(currents) +
  mapview::mapview(subset(comb14c2_sf, !is.na(comb14c2_sf$year)), zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), label = "sta", layer.name = "F14C_calc", alpha = 0, alpha.regions = 0, legend = TRUE, popup = FALSE, highlight = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year < 1900), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "pre-1900", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1900 & year < 1910), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1900s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1910 & year < 1920), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1910s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1920 & year < 1930), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1920s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1930 & year < 1940), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1930s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1940 & year < 1950), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1940s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1950 & year < 1960), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1950s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1960 & year < 1970), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1960s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1970 & year < 1980), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1970s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1980 & year < 1990), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1980s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 1990 & year < 2000), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "1990s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 2000 & year < 2010), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "2000s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 2010 & year < 2020), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "2010s", legend = FALSE) +
  mapview::mapview(st_jitter(subset(comb14c2_sf, !is.na(comb14c2_sf$year) & year >= 2020 & year < 2030), factor = 0.001), label = "year", zcol = "F14C_calc", at = seq(minF14C, maxF14C, 0.1), layer.name = "2020s", legend = FALSE)

comb14c2_sf$current <- st_nearest_feature(comb14c2_sf, currents)
comb14c2_sf$current2 <- lapply(comb14c2_sf$current, function(x) currents$NAME[x])
comb14c2_sf$current2 <- unlist(comb14c2_sf$current2)
comb14c2_sf$current2[which(is.na(comb14c2_sf$current2))] <- "Gulf Stream"

gamtest7 <- brm(bf(F14C_calc ~ s(year, k = 20) + current + source), data = subset(comb14c2_sf, !is.na(comb14c2_sf$year)), family = gaussian(), 
                cores = 4, seed = 456, chains = 4, iter = 5000, warmup = 1500, thin = 3, 
                control = list(adapt_delta = 0.99, max_treedepth = 15), backend = "cmdstanr", threads = threading(2))
gamtest3_ms <- conditional_smooths(gamtest3)
plot(gamtest3_ms)

ggplot(comb14c[!is.na(year), ]) +
  geom_ribbon(data = gamtest5_ms$mu, aes(x = year, ymin = lower__ + 1, ymax = upper__ + 1), alpha = 0.2) +
  geom_point(aes(x = year, y = F14C_calc, color = study, shape = source)) +
  geom_point(data = comb14c[!is.na(year) & !is.na(F14C_calc) & str_detect(species, "Mercenaria sp"), ], aes(x = year, y = F14C_calc), shape = 1, color = "black", stroke = 1, size = 1.5) +
  geom_line(data = gamtest5_ms$mu, aes(x = year, y = estimate__ + 1)) +
  geom_point(data = comb14c[!is.na(year) & study == "Andrews et al. 2020" & F14C_calc <= 0.94, ], aes(x = year, y = F14C_calc), shape = 2, size = 1) +
  theme_bw()


