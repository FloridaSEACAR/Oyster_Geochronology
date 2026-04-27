
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
# mapviewOptions(fgb = FALSE)
library(patchwork)
library(colorspace)
library(ggpubr)
library(plotly)
library(broom.mixed)


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

pri <- fread(here::here("HOBS_geochronology_ms2/PRI_Class_Bival.csv"))
colnames(pri)
pri2 <- pri[, .(PRICatalogNumber, PRIStationNumber, IGSNNumber, Count, Order, Family, Genus, Species, Subspecies, County, Latitude, Longitude)]
setnames(pri2, c("PRICatalogNumber", "PRIStationNumber", "Count"), c("CatNum", "Station", "Preparation"))
pri2[, Institution := "PRI"]
pri2[, `:=` (Station = as.character(Station), Preparation = as.character(Preparation), locDetails = NA, Notes = NA, Depth = NA, CatNum = as.character(CatNum), StartDate = NA)]

bmsm <- fread(here::here("HOBS_geochronology_ms2/BMSM_Class_Bival.csv"))
colnames(bmsm)
bmsm2 <- bmsm[, .(BMSMNo, Order, Family, Genus, Species, Subspecies, County, Site, Latitude1, Longitude1, preparations, StartDate, DepthMIN, Remarks)]
setnames(bmsm2, c("BMSMNo", "Site", "preparations", "Latitude1", "Longitude1", "DepthMIN", "Remarks"), c("CatNum", "Station", "Preparation", "Latitude", "Longitude", "Depth", "Notes"))
bmsm2[, Institution := "BMSM"]
bmsm2[, `:=` (Station = as.character(Station), Preparation = as.character(Preparation), locDetails = NA, Depth = as.character(Depth), CatNum = as.character(CatNum), 
              StartDate = as.Date(StartDate), IGSNNumber = NA)]

flmnh <- fread(here::here("HOBS_geochronology_ms2/FLMNH_Class_Bival.csv"))
colnames(flmnh)
flmnh2 <- flmnh[, .(UFID, Order, Family, Genus, Species, Subspecies, Preparations, County, Locality, Latitude, Longitude, StartDate, Depth1, DepthUnit)]
setnames(flmnh2, c("UFID", "Locality", "Preparations", "Depth1"), c("CatNum", "Station", "Preparation", "Depth"))
flmnh2[, Institution := "FLMNH"]
flmnh2[, `:=` (Station = as.character(Station), Preparation = as.character(Preparation), locDetails = NA, Notes = NA, Depth = paste0(Depth, " ", DepthUnit), CatNum = as.character(CatNum), 
               StartDate = as.Date(StartDate), IGSNNumber = NA)]
flmnh2[, DepthUnit := NULL]

fwc <- fread(here::here("HOBS_geochronology_ms2/FWC_Class_Bival.csv"))
colnames(fwc)
fwc2 <- fwc[, .(CatalogFSBCI, Preparations, Order, Family, Genus, Species, Subspecies, Date, Verbatimlocality, Latitude, Longitude, Minimumdepth)]
setnames(fwc2, c("CatalogFSBCI", "Verbatimlocality", "Preparations", "Minimumdepth", "Date"), c("CatNum", "Station", "Preparation", "Depth", "StartDate"))
fwc2[, Institution := "FWC"]
fwc2[, `:=` (Station = as.character(Station), County = NA, Preparation = as.character(Preparation), locDetails = NA, Notes = NA, Depth = as.character(Depth), CatNum = as.character(CatNum), 
             StartDate = as.Date(StartDate), IGSNNumber = NA)]

nmnh1 <- fread(here::here("HOBS_geochronology_ms2/NMNH_Class_Bival.csv"))
nmnh1[, `:=` (Cruise = as.character(Cruise), `Station Number` = as.character(`Station Number`))]
nmnh2 <- fread(here::here("HOBS_geochronology_ms2/nmnhsearch-20230404211837.csv"))
nmnh2[, `:=` (Cruise = as.character(Cruise), `Station Number` = as.character(`Station Number`))]
nmnh3 <- fread(here::here("HOBS_geochronology_ms2/nmnhsearch-20230404211806.csv"))
nmnh3[, `:=` (Cruise = as.character(Cruise), `Station Number` = as.character(`Station Number`))]
nmnh4 <- fread(here::here("HOBS_geochronology_ms2/nmnhsearch-20230404211724.csv"))
nmnh4[, `:=` (Cruise = as.character(Cruise), `Station Number` = as.character(`Station Number`))]
nmnh5 <- fread(here::here("HOBS_geochronology_ms2/nmnhsearch-20230404211627.csv"))
nmnh5[, `:=` (Cruise = as.character(Cruise), `Station Number` = as.character(`Station Number`))]
nmnh6 <- fread(here::here("HOBS_geochronology_ms2/nmnhsearch-20230404211613.csv"))
nmnh6[, `:=` (Cruise = as.character(Cruise), `Station Number` = as.character(`Station Number`))]
nmnh7 <- fread(here::here("HOBS_geochronology_ms2/nmnhsearch-20230404211440.csv"))
nmnh7[, `:=` (Cruise = as.character(Cruise), `Station Number` = as.character(`Station Number`))]
nmnh8 <- fread(here::here("HOBS_geochronology_ms2/nmnhsearch-20230404211408.csv"))
nmnh8[, `:=` (Cruise = as.character(Cruise), `Station Number` = as.character(`Station Number`))]
nmnh9 <- fread(here::here("HOBS_geochronology_ms2/nmnhsearch-20230404211251.csv"))
nmnh9[, `:=` (Cruise = as.character(Cruise), `Station Number` = as.character(`Station Number`))]
nmnh10 <- fread(here::here("HOBS_geochronology_ms2/nmnhsearch-20230404211229.csv"))
nmnh10[, `:=` (Cruise = as.character(Cruise), `Station Number` = as.character(`Station Number`))]
nmnh11 <- fread(here::here("HOBS_geochronology_ms2/nmnhsearch-20230404211102.csv"))
nmnh11[, `:=` (Cruise = as.character(Cruise), `Station Number` = as.character(`Station Number`))]
nmnh <- bind_rows(nmnh1, nmnh2, nmnh3, nmnh4, nmnh5, nmnh6, nmnh7, nmnh8, nmnh9, nmnh10, nmnh11)
rm(nmnh1, nmnh2, nmnh3, nmnh4, nmnh5, nmnh6, nmnh7, nmnh8, nmnh9, nmnh10, nmnh11)
nmnh <- janitor::clean_names(nmnh)
colnames(nmnh)
nmnh[, `:=` (Genus = as.character(ifelse(scientific_name == "", "NA", str_sub(scientific_name, 1, str_locate(scientific_name, " ")[1] - 1))), 
             Species = as.character(ifelse(scientific_name == "", "NA", ifelse(str_detect(scientific_name, ", "), 
                                                                               str_sub(scientific_name, str_locate_all(scientific_name, " [:lower:]")[[1]][length(str_locate_all(scientific_name, " [:lower:]")[[1]])], 
                                                                                       str_locate_all(scientific_name, " ")[[1]][which(str_locate_all(scientific_name, " ")[[1]] == str_locate_all(scientific_name, ", ")[[1]][2])[1] - 1] - 1), 
                                                                               str_sub(scientific_name, str_locate(scientific_name, " [:lower:]")[2], -1))))),
     by = scientific_name]

nmnh2 <- nmnh[, .(catalog_number, specimen_count, preparation_details_preparation_remarks, order, family, Genus, Species, station_number, date_collected, district_county, city_town, precise_locality, centroid_latitude, centroid_longitude, depth_m, depth_notes, notes)]
nmnh2[, `:=` (Preparation = paste0(specimen_count, "; ", preparation_details_preparation_remarks), 
              Station = paste0("city: ", city_town, "; depth (m): ", depth_m, " (", depth_notes, "); ", precise_locality),
              Institution = "NMNH",
              Subspecies = NA)]
setnames(nmnh2, c("catalog_number", "depth_m", "date_collected", "district_county", "station_number", "order", "family", "centroid_latitude", "centroid_longitude", "notes"), 
                c("CatNum", "Depth", "StartDate", "County", "locDetails", "Order", "Family", "Latitude", "Longitude", "Notes"))
nmnh2[, `:=` (Station = as.character(Station), Preparation = as.character(Preparation), CatNum = as.character(CatNum), StartDate = as.Date(StartDate, tryFormats = c("%e %b %Y", "%d %b %Y", "%Y")), IGSNNumber = NA)]
nmnh2 <- distinct(nmnh2[, .(CatNum, County, Depth, Family, Genus, Institution, Latitude, locDetails, Longitude, Notes, Order, Preparation, Species, StartDate, Station, Subspecies)])

allm <- bind_rows(pri2, bmsm2, flmnh2, fwc2, nmnh2)
setcolorder(allm, c("Institution", "IGSNNumber", "CatNum", "County", "StartDate", "Preparation", "Depth", "Station", "Order", "Family", "Genus", "Species", "Subspecies", "Latitude", "Longitude", "locDetails"))
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

comb14c <- fread(here::here("CalCurve/14Cdata_FL_Gulf_Caribbean_20230329.csv"))

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

# #make sure text encoding is UTF-8 so mapview can work (following: https://github.com/r-spatial/mapview/issues/194#issuecomment-471172491)
# char_utf8 <- function(mapdat){
#   # delete sticky geometry list column (will re-add later)
#   # tst3 = st_set_geometry(comb14c_sf2, NULL)
#   tst3 = st_set_geometry(mapdat, NULL)
#   
#   # function to check for 
#   f = function(x) x == "character"
#   
#   # which columns are of class character?
#   charcols = sapply(sapply(tst3, class), f)
#   
#   # convert those columns from latin1 to utf8 encoding
#   # note that it is actually necessary to know the original encoding (here guessed to be latin1)
#   tst3[, charcols] = sapply(tst3[, charcols], iconv, from = "latin1", to = "UTF-8", sub = "")
#   tst3[, (charcols) := iconv(charcols)]
#   for(i in which(charcols)){
#     tst3[, (i) := iconv(tst3[, i, with = F], from = "latin1", to = "UTF-8", sub = "")]
#   }
#   
#   # re-add geometry column and convert back to sf
#   tst3$geometry = st_geometry(mapdat)
#   tst3 = st_as_sf(tst3)
#   return(tst3)
# }

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
                control = list(adapt_delta = 0.99, max_treedepth = 15), backend = "cmdstanr", threads = threading(2), file = here::here("gamtest7b.rds")) #Not sure what happened to "gamtest7.rds" but it started plotting weird so I re-ran the model and saved it as "gamtest7b.rds" and the plotting issues were fixed
gamtest7_ms <- conditional_smooths(gamtest7, int_conditions = list(year = seq(range(subset(comb14c2_sf, !is.na(comb14c2_sf$year))$year)[1], range(subset(comb14c2_sf, !is.na(comb14c2_sf$year))$year)[2])))
gamtest7_ms$`mu: s(year,k=20)`$estimate2 <- gamtest7_ms$`mu: s(year,k=20)`$estimate__ + 1.022722
setnames(gamtest7_ms$`mu: s(year,k=20)`, "year", "yr")
plot(gamtest7_ms)

ggplot(comb14c2[!is.na(year), ]) +
  geom_ribbon(data = gamtest7_ms$`mu: s(year,k=20)`, aes(x = yr, ymin = lower__ + 1.022722, ymax = upper__ + 1.022722), alpha = 0.2) +
  geom_point(aes(x = year, y = F14C_calc, color = study, shape = source)) +
  geom_point(data = comb14c2[!is.na(year) & !is.na(F14C_calc) & str_detect(species, "Mercenaria sp"), ], aes(x = year, y = F14C_calc), shape = 1, color = "black", stroke = 1, size = 1.5) +
  geom_line(data = gamtest7_ms$`mu: s(year,k=20)`, aes(x = yr, y = estimate2)) +
  geom_point(data = comb14c2[!is.na(year) & study == "Andrews et al. 2020" & F14C_calc <= 0.94, ], aes(x = year, y = F14C_calc), shape = 2, size = 1) +
  theme_bw()

inits <- data.table(year = seq(range(subset(comb14c2_sf, !is.na(comb14c2_sf$year))$year)[1], range(subset(comb14c2_sf, !is.na(comb14c2_sf$year))$year)[2], by = 0.5))
inits[, `:=` (current = rep(65, nrow(inits)), source = rep("otolith", nrow(inits)))]
preds <- predict(gamtest7, newdata = inits)
preds2 <- data.table(estimate = preds[,1], esterror = preds[,2], q2_5 = preds[,3], q97_5 = preds[,4])
preds2[, year := inits$year]
setcolorder(preds2, "year")

ggplot(comb14c2[!is.na(year), ]) +
  geom_ribbon(data = preds2, aes(x = year, ymin = estimate - esterror, ymax = estimate + esterror), alpha = 0.2) +
  geom_point(aes(x = year, y = F14C_calc, color = study, shape = source)) +
  geom_point(data = comb14c2[!is.na(year) & !is.na(F14C_calc) & str_detect(species, "Mercenaria sp"), ], aes(x = year, y = F14C_calc), shape = 1, color = "black", stroke = 1, size = 1.5) +
  geom_line(data = preds2, aes(x = year, y = estimate)) +
  geom_point(data = comb14c2[!is.na(year) & study == "Andrews et al. 2020" & F14C_calc <= 0.94, ], aes(x = year, y = F14C_calc), shape = 2, size = 1) +
  theme_bw()


max(gamtest7_ms$`mu: s(year,k=20)`$estimate2) #1.141562, which is lower than the peak of Quan's curve (1.1588)...

f14c_quan <- fread(here::here("CalCurve/QuanBombCurve_20230112.csv"), sep = "|")

ggplot(comb14c2[!is.na(year), ]) +
  geom_ribbon(data = gamtest7_ms$`mu: s(year,k=20)`, aes(x = year, ymin = lower__ + 1.022722, ymax = upper__ + 1.022722), alpha = 0.2) +
  geom_point(aes(x = year, y = F14C_calc, color = study, shape = source)) +
  geom_point(data = comb14c2[!is.na(year) & !is.na(F14C_calc) & str_detect(species, "Mercenaria sp"), ], aes(x = year, y = F14C_calc), shape = 1, color = "black", stroke = 1, size = 1.5) +
  geom_line(data = gamtest7_ms$`mu: s(year,k=20)`, aes(x = year, y = estimate__ + 1.022722)) +
  geom_point(data = comb14c2[!is.na(year) & study == "Andrews et al. 2020" & F14C_calc <= 0.94, ], aes(x = year, y = F14C_calc), shape = 2, size = 1) +
  geom_line(data = f14c_quan, aes(x = Year, y = F14C), lty = "dashed", color = "firebrick", lwd = 0.75) +
  theme_bw()

studycolslist <- qualitative_hcl(28)
studycols <- setNames(studycolslist, sort(unique(c(comb14c2$study, f14c_quan$study))))
f14c_quan[, study := factor(study, levels = names(studycols))]
comb14c2[, study := factor(study, levels = names(studycols))]
f14c_quan[, source := factor(source, levels = c("coral", "otolith", "shell", "surface water"))]
comb14c2[, source := factor(source, levels = c("coral", "otolith", "shell", "surface water"))]

qcurve <- ggplot(data = f14c_quan, aes(x = Year, y = F14C, color = study, shape = source)) +
  geom_line(aes(x = Year, y = F14C), lty = "solid", color = "firebrick", lwd = 0.75, inherit.aes = FALSE) +
  geom_point() +
  theme_bw() +
  labs(x = NULL, y = "F14C", color = "Study", shape = "Source") +
  coord_cartesian(xlim = c(1950, 2023), ylim = c(0.75, 1.2))

scurve <- ggplot(data = comb14c2[!is.na(year), ], aes(x = year, y = F14C_calc, color = study, shape = source)) +
  geom_ribbon(data = gamtest7_ms$`mu: s(year,k=20)`, aes(x = year, ymin = lower__ + 1.022722, ymax = upper__ + 1.022722), alpha = 0.2, inherit.aes = FALSE) +
  geom_point() +
  geom_point(data = comb14c2[!is.na(year) & !is.na(F14C_calc) & str_detect(species, "Mercenaria sp"), ], aes(x = year, y = F14C_calc), shape = 1, color = "black", stroke = 1, size = 1.5, inherit.aes = FALSE) +
  geom_line(data = gamtest7_ms$`mu: s(year,k=20)`, aes(x = year, y = estimate2), lwd = 0.75, inherit.aes = FALSE) +
  theme_bw() +
  labs(x = NULL, y = "F14C_calc", color = "Study", shape = "Source") +
  coord_cartesian(xlim = c(1750, 2023), ylim = c(0.75, 1.2))

qscurves <- ggplot() +
  geom_ribbon(data = gamtest7_ms$`mu: s(year,k=20)`, aes(x = year, ymin = lower__ + 1.022722, ymax = upper__ + 1.022722), alpha = 0.2) +
  geom_line(data = gamtest7_ms$`mu: s(year,k=20)`, aes(x = year, y = estimate2), lwd = 0.75) +
  geom_line(data = f14c_quan, aes(x = Year, y = F14C), lty = "solid", color = "firebrick", lwd = 0.75) +
  theme_bw() +
  labs(x = "year", y = "F14C") +
  coord_cartesian(xlim = c(1950, 2023), ylim = c(0.75, 1.2)) +
  theme(legend.position = "none")

curvecompare <- qcurve / scurve / qscurves + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & 
  scale_colour_manual(values = studycols, drop = FALSE) & theme(legend.text = element_text(size = 6),
                                                                legend.title = element_text(size = 8),
                                                                axis.text = element_text(size = 8),
                                                                axis.title = element_text(size = 8),
                                                                plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))

ggsave(filename = here::here(paste0("CalCurve/CurveComparePlot_", Sys.Date(), ".png")),
       plot = curvecompare,
       width = 7,
       height = 6,
       units = "in",
       dpi = 300)



#Compare Apalachicola Bay oyster results from Hadden & Durham work----------------------------------
apaCompare <- comb14c2[study %in% c("Durham et al. 2023", "Hadden et al. 2018", "Hadden & Cherkinsky 2017"), ]
apaCompare2 <- apaCompare[study == "Durham et al. 2023" & str_detect(sta, "George|Goose|8728"), ]
apaCompare3 <- apaCompare[(study == "Durham et al. 2023" & sta %in% unique(apaCompare2$sta)) | str_detect(study, "Hadden"), ]

apaCompare3[, sta2 := fcase(sta %in% c("Indian Pass", "Indian Lagoon", "NA (Little St. George Is.)", "St. Vincent Sound"), "West Apalachicola Bay",
                            sta %in% c("10578 (Little St. George Is.)", "10577 (Little St. George Is.)"), "Mid Apalachicola Bay",
                            sta %in% c("Cat Point", "10576 (Goose Island/East Cove)", "10570 (Goose Island/East Cove)"), "East Apalachicola Bay",
                            sta == "Dog Island", "Dog Island",
                            sta %in% c("Alligator Harbor", "8728 (NA)"), "Alligator Harbor")]
apaCompare3[, sta2 := factor(sta2, levels = c("West Apalachicola Bay", "Mid Apalachicola Bay", "East Apalachicola Bay", "Dog Island", "Alligator Harbor"))]

apaCompPlot <- ggplot(data = apaCompare3, aes(x = sta2, y = F14C_calc, fill = as.factor(year), shape = study)) +
  geom_point(size = 2, color = "black") +
  ylim(0.75, 1.2) +
  theme_bw() +
  labs(x = "Locality", fill = "Year", shape = "Study") +
  scale_fill_discrete_sequential(palette = "inferno", nmax = 10, order = c(3:7), breaks = c("1938", "1941", "2016", "2017", "2020")) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  scale_shape_manual(values = c("Durham et al. 2023" = 21, "Hadden & Cherkinsky 2017" = 22, "Hadden et al. 2018" = 23)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = here::here(paste0("CalCurve/ApalachicolaBayComparePlot_", Sys.Date(), ".png")),
       plot = apaCompPlot,
       width = 6,
       height = 4,
       units = "in",
       dpi = 300)

#Test calculations of Dead C-------------------------------------------

comb14c2[, `:=` (w_i = 1/(F14C_sd^2), w_i_calc = 1/(F14C_calc_sd^2))]
# comb14c2[year]
# comb14c2[, `:=` (deadC_r = , deadC_o = )]

apaCompare3[is.na(F14C_calc_sd), F14C_calc_sd := apaCompare3[!is.na(F14C_calc_sd), max(F14C_calc_sd)]]
apaCompare3[, `:=` (w_i = 1/(F14C_sd^2), 
                    w_i_calc = 1/(F14C_calc_sd^2), 
                    refval = gamtest7_ms$`mu: s(year,k=20)`$estimate2[which(gamtest7_ms$`mu: s(year,k=20)`$yr == year)],
                    relyear = year - (min(apaCompare3$year) - 1)), by = year]
apaCompare3[year == 1938, mean_w := sum(w_i_calc * F14C_calc)/sum(w_i_calc), by = sta2]
apaCompare3[year == 1941, mean_w := sum(w_i_calc * F14C_calc)/sum(w_i_calc), by = sta2]
apaCompare3[year > 2015, mean_w := sum(w_i_calc * F14C_calc)/sum(w_i_calc), by = sta2]

apaCompare3[, `:=` (deadC = (1 - mean_w)/refval)]
apaCompare3[, `:=` (F14C_calc_corr = F14C_calc/(1 - deadC))]

apaCompPlot_corr <- ggplot(data = apaCompare3, aes(x = sta2, y = F14C_calc_corr, fill = as.factor(year), shape = study)) +
  geom_point(size = 2, color = "black") +
  ylim(0.75, 1.2) +
  theme_bw() +
  labs(x = "Locality", fill = "Year", shape = "Study") +
  scale_fill_discrete_sequential(palette = "inferno", nmax = 10, order = c(3:7), breaks = c("1938", "1941", "2016", "2017", "2020")) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  scale_shape_manual(values = c("Durham et al. 2023" = 21, "Hadden & Cherkinsky 2017" = 22, "Hadden et al. 2018" = 23)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

apaCompPlot_corr2 <- ggplot(data = apaCompare3, aes(x = relyear, y = F14C_calc_corr, fill = sta2, shape = study)) +
  geom_point(size = 2, color = "black") +
  ylim(0.75, 1.2) +
  theme_bw() +
  labs(x = "Year", fill = "Locality", shape = "Study") +
  scale_fill_discrete_sequential(palette = "inferno", nmax = 10, order = c(3:7), breaks = levels(apaCompare3$sta2)) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  scale_shape_manual(values = c("Durham et al. 2023" = 21, "Hadden & Cherkinsky 2017" = 22, "Hadden et al. 2018" = 23)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

apaCompPlot_deadc <- ggplot(data = apaCompare3, aes(x = relyear, y = deadC, fill = sta2, shape = study)) +
  geom_point(size = 2, color = "black") +
  # ylim(0.75, 1.2) +
  theme_bw() +
  labs(x = "Year", fill = "Locality", shape = "Study") +
  scale_fill_discrete_sequential(palette = "inferno", nmax = 10, order = c(3:7), breaks = levels(apaCompare3$sta2)) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  scale_shape_manual(values = c("Durham et al. 2023" = 21, "Hadden & Cherkinsky 2017" = 22, "Hadden et al. 2018" = 23)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

apadr3 <- dr3[str_detect(Locality, "George|Goose"), ]
apadr3[, estyr := (F14C - 0.93)/0.00098]

apaCompPlot_f14c <- ggplot(data = apaCompare3, aes(x = relyear, y = F14C_calc)) +
  geom_point(aes(x = relyear, y = F14C_calc, fill = sta2, shape = study), size = 2, color = "black") +
  geom_smooth(method = "lm", color = "black") +
  ggpubr::stat_regline_equation() +
  # geom_point(data = dr3[])
  ylim(0.75, 1.2) +
  theme_bw() +
  labs(x = "Year", fill = "Locality", shape = "Study") +
  scale_fill_discrete_sequential(palette = "inferno", nmax = 10, order = c(3:7), breaks = levels(apaCompare3$sta2)) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  scale_shape_manual(values = c("Durham et al. 2023" = 21, "Hadden & Cherkinsky 2017" = 22, "Hadden et al. 2018" = 23)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

apaComPlots <- apaCompPlot_f14c / apaCompPlot_deadc / apaCompPlot_corr2 + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & 
  coord_cartesian(xlim = c(78, 85)) &
  theme(legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))


#Retrieve Gov Apalachicola data for temp, sal, carbon & flow pulled using the "ds-pipelines-targets-example-wqp.Rproj" pipeline.
# apaw <- fread(here::here("ds-pipelines-targets-example-wqp-main/SavedData/ApalachicolaBay_TempSalFlowC_2023-04-07.txt"), sep = "|")
# ahw <- fread(here::here("ds-pipelines-targets-example-wqp-main/SavedData/AlligatorHarbortoAPA_TempSalFlowC_2023-04-27.txt"), sep = "|")
# apaw <- rbind(apaw, ahw)
# apaw_sites <- fread(here::here("ds-pipelines-targets-example-wqp-main/SavedData/ApalachicolaBay_TempSalFlowC_SiteInfo_2023-04-07.txt"), sep = "|")
# ahw_sites <- fread(here::here("ds-pipelines-targets-example-wqp-main/SavedData/AlligatorHarbortoAPA_TempSalFlowC_SiteInfo_2023-04-27.txt"), sep = "|")
# apaw_sites <- rbind(apaw_sites, ahw_sites)
apaw <- fread(here::here("ds-pipelines-targets-example-wqp-main/SavedData/AlligatorHarbortoAPA_TempSalFlowC_2023-04-27.txt"), sep = "|")
apaw_sites <- fread(here::here("ds-pipelines-targets-example-wqp-main/SavedData/AlligatorHarbortoAPA_TempSalFlowC_SiteInfo_2023-04-27.txt"), sep = "|")
apaw_red <- distinct(apaw[, .(MonitoringLocationIdentifier, parameter, ActivityStartDate, ActivityStartDateTime, CharacteristicName, ResultMeasureValue, ResultMeasure.MeasureUnitCode, ResultCommentText)])
apaw2 <- merge(apaw_red, distinct(apaw_sites[, .(site_no, MonitoringLocationTypeName, dec_lon_va, dec_lat_va)]), by.x = "MonitoringLocationIdentifier", by.y = "site_no", all.x = TRUE)
apaw2 <- apaw2[str_detect(CharacteristicName, "stage", negate = TRUE) & 
               str_detect(CharacteristicName, "Flow, severity (choice list)", negate = TRUE) & 
               ResultMeasure.MeasureUnitCode != "ft/sec" & 
               ResultCommentText == "" &
               !is.na(ResultMeasureValue), ]

apaw_sites_sf <- st_as_sf(apaw_sites, coords = c("dec_lon_va", "dec_lat_va"), crs = 4326)
apaw_sites2 <- merge(apaw_sites_sf, distinct(apaw_red[, .(MonitoringLocationIdentifier, parameter)]), by.x = "site_no", by.y = "MonitoringLocationIdentifier", all.x = TRUE)
mapview(subset(apaw_sites2, apaw_sites2$parameter == "carbon"), layer.name = "carbon", col.regions = "dodgerblue") +
  mapview(subset(apaw_sites2, apaw_sites2$parameter == "flow"), layer.name = "flow", col.regions = "firebrick") +
  mapview(subset(apaw_sites2, apaw_sites2$parameter == "salinity"), layer.name = "salinity", col.regions = "orange") +
  mapview(subset(apaw_sites2, apaw_sites2$parameter == "temperature"), layer.name = "temperature", col.regions = "forestgreen")

# apaw2[!is.na(ActivityStartDateTime), asdate := date(ActivityStartDateTime)]
# apaw2[!is.na(ActivityStartDate) & is.na(ActivityStartDateTime), asdate := date(ActivityStartDate)]
apaw2[, asdate := as_date(ActivityStartDate)]
apaw2[, `:=` (MonitoringLocationIdentifier2 = MonitoringLocationIdentifier, parameter2 = parameter)]
sumdat <- distinct(apaw2[, .(n = .N, 
                             begin = min(asdate), #ifelse(is.na(ActivityStartDateTime), min(as_date(ActivityStartDate)), min(as_date(ActivityStartDateTime))), 
                             end = max(asdate)), #ifelse(is.na(ActivityStartDateTime), max(as_date(ActivityStartDate)), max(as_date(ActivityStartDateTime))), 
                             # units = apaw2[parameter2 == parameter & MonitoringLocationIdentifier2 == MonitoringLocationIdentifier, unique(ResultMeasure.MeasureUnitCode)]),
                         by = list(parameter, MonitoringLocationIdentifier)])
for(i in seq(1, nrow(sumdat))){
  sumdat[i, units := list(list(apaw2[parameter2 == sumdat[i, parameter] & MonitoringLocationIdentifier2 == sumdat[i, MonitoringLocationIdentifier], unique(ResultMeasure.MeasureUnitCode)]))]
}

# sumdat[, `:=` (begin = as_date(begin), end = as_date(end))]
sumdat[, yr_rng := year(end) - year(begin)]

# locs <- sumdat[n >= 50 & yr_rng >= 10, .(parameter, MonitoringLocationIdentifier)]
locs <- distinct(sumdat[, .(parameter, MonitoringLocationIdentifier)])

mapview(subset(apaw_sites2, apaw_sites2$parameter == "carbon" & apaw_sites2$MonitoringLocationIdentifier %in% locs$MonitoringLocationIdentifier), layer.name = "carbon", col.regions = "dodgerblue") +
  mapview(subset(apaw_sites2, apaw_sites2$parameter == "flow" & apaw_sites2$MonitoringLocationIdentifier %in% locs$MonitoringLocationIdentifier), layer.name = "flow", col.regions = "firebrick") +
  mapview(subset(apaw_sites2, apaw_sites2$parameter == "salinity" & apaw_sites2$MonitoringLocationIdentifier %in% locs$MonitoringLocationIdentifier), layer.name = "salinity", col.regions = "orange") +
  mapview(subset(apaw_sites2, apaw_sites2$parameter == "temperature" & apaw_sites2$MonitoringLocationIdentifier %in% locs$MonitoringLocationIdentifier), layer.name = "temperature", col.regions = "forestgreen")

apaw_use <- apaw2[MonitoringLocationIdentifier %in% locs$MonitoringLocationIdentifier & !is.na(ResultMeasureValue), ]

ggplot(apaw_use, aes(x = asdate, y = ResultMeasureValue, color = ResultMeasure.MeasureUnitCode)) +
  geom_line() +
  geom_point(shape = 1) +
  theme_bw() +
  facet_wrap(~ parameter, ncol = 1, scales = "free_y")

disc <- fread("C:/Users/steph/Dropbox/SEACAR/WQ_Summaries/Combined_WQ_WC_NUT_row-2021-Nov-30.csv")
apadisc <- disc[ManagedAreaName %in% c("Apalachicola National Estuarine Research Reserve", "Alligator Harbor Aquatic Preserve") & ParameterName %in% c("Salinity", "Water Temperature", "Turbidity", "Total Suspended Solids, TSS"), ]
setnames(apadisc, c("ProgramID", "ProgramName", "ProgramLocationID", "Activity_Start_Date_Time", "ActivityDepth", 
                    "ActivityDepthUnit", "ParameterName", "ResultValue", "Units", "Latitude", "Longitude"), 
         c("OrganizationIdentifier", "OrganizationFormalName", "MonitoringLocationIdentifier", "ActivityStartDateTime", 
           "ResultDepthHeightMeasure.MeasureValue", "ResultDepthHeightMeasure.MeasureUnitCode", "CharacteristicName", 
           "ResultMeasureValue", "ResultMeasure.MeasureUnitCode", "dec_lat_va", "dec_lon_va"))
apadisc[, `:=` (parameter = str_to_lower(CharacteristicName), 
                OrganizationIdentifier = as.character(OrganizationIdentifier),
                ResultMeasureValue = as.numeric(ResultMeasureValue), 
                ResultDepthHeightMeasure.MeasureValue = as.numeric(ResultDepthHeightMeasure.MeasureValue), 
                ActivityStartDateTime = as_datetime(ActivityStartDateTime))]


apadisc[!is.na(Relative_Depth), MonitoringLocationIdentifier := paste0(MonitoringLocationIdentifier, "_", Relative_Depth)]

# apadisc_apaw <- distinct(apadisc[, .(OrganizationIdentifier, OrganizationFormalName, MonitoringLocationIdentifier, ActivityStartDateTime, 
#                                      ResultDepthHeightMeasure.MeasureValue, ResultDepthHeightMeasure.MeasureUnitCode, CharacteristicName, 
#                                      ResultMeasureValue, ResultMeasure.MeasureUnitCode, parameter, dec_lat_va, dec_lon_va)])
apadisc_apaw <- distinct(apadisc[, .(MonitoringLocationIdentifier, ActivityStartDateTime, CharacteristicName, 
                                     ResultMeasureValue, ResultMeasure.MeasureUnitCode, parameter, dec_lat_va, dec_lon_va)])
apadisc_apaw[, asdate := date(ActivityStartDateTime)]

apadisc_apaw[, RMV := mean(ResultMeasureValue), by = list(MonitoringLocationIdentifier, asdate, parameter)][, ResultMeasureValue := RMV][, RMV := NULL]
apadisc_apaw <- distinct(apadisc_apaw)

apadisc_sf <- st_as_sf(apadisc, coords = c("dec_lon_va", "dec_lat_va"), crs = 4326)
apadisc_sf2 <- distinct(apadisc_sf[, c("MonitoringLocationIdentifier", "parameter", "geometry")])
apadisc_apaw[, MonitoringLocationTypeName := NA]

apaw_use2 <- bind_rows(apaw_use, apadisc_apaw)
apaw_use2 <- distinct(apaw_use2[!is.na(ResultMeasureValue) & ResultMeasure.MeasureUnitCode != "g/kg" & ResultMeasure.MeasureUnitCode != "ft/sec", ])
apaw_sites3 <- bind_rows(apaw_sites2, apadisc_sf2)

apaw_use2[parameter == "salinity", ResultMeasureValue := fcase(ResultMeasure.MeasureUnitCode %in% c("ppt", "ppth", "PSU"), ResultMeasureValue,
                                                               ResultMeasure.MeasureUnitCode == "g/L", ResultMeasureValue * 1000,
                                                               ResultMeasure.MeasureUnitCode == "ppb", ResultMeasureValue / 1000 )][parameter == "salinity", ResultMeasure.MeasureUnitCode := "ppt"]
apaw_use2[parameter == "flow", ResultMeasureValue := fcase(ResultMeasure.MeasureUnitCode == "m3/sec", ResultMeasureValue/0.028316832,
                                                           ResultMeasure.MeasureUnitCode %in% c("ft3/sec", "ft3/s"), ResultMeasureValue,
                                                           ResultMeasure.MeasureUnitCode == "Mgd", ResultMeasureValue * 1.5472286365101,
                                                           ResultMeasure.MeasureUnitCode == "gal/min", ResultMeasureValue * 0.0022280092365745)][parameter == "flow", ResultMeasure.MeasureUnitCode := "cfs"]
apaw_use2[parameter == "carbon", ResultMeasure.MeasureUnitCode := fcase(ResultMeasure.MeasureUnitCode %in% c("mg/l", "per mil"), "mg/l",
                                                                        ResultMeasure.MeasureUnitCode == "pct modern", "pct modern")]

apaw_use2 <- distinct(apaw_use2)
apaw_use2[parameter == "water temperature", `:=` (parameter = "temperature", ResultMeasure.MeasureUnitCode = "deg C")]
apaw_use2[ResultMeasure.MeasureUnitCode == "mg/L", ResultMeasure.MeasureUnitCode := "mg/l"]

ggplot(apaw_use2, aes(x = asdate, y = ResultMeasureValue, color = ResultMeasure.MeasureUnitCode)) +
  geom_line() +
  geom_point(shape = 1) +
  theme_bw() +
  facet_wrap(~ parameter, ncol = 1, scales = "free_y")


shore <- st_read(here::here("FLCounties/Counties_-_Detailed_Shoreline.shp"))
shore <- st_make_valid(shore)
shore <- st_transform(shore, crs = 4326)
bath <- st_read(here::here("Bathymetry/bathym.shp"))
bath <- st_make_valid(bath)
bath <- st_transform(bath, crs = 4326)
wbid <- st_read(here::here("IWR_WBIDs_Run61/WBID_Run61.shp"))
wbid <- st_make_valid(wbid)
wbid <- st_transform(wbid, crs = 4326)
apaw_use2_sf <- distinct(apaw_use2[, .(MonitoringLocationIdentifier, dec_lon_va, dec_lat_va)])
apaw_use2_sf <- st_as_sf(apaw_use2_sf, coords = c("dec_lon_va", "dec_lat_va"), crs = 4326)
wbid_apa <- wbid[apaw_use2_sf$geometry, , op = st_intersects]
apaw_use2[, MonitoringLocationIdentifier := as.character(MonitoringLocationIdentifier)]
wbid_apa$centroid <- lapply(wbid_apa$geometry, function(x) st_centroid(x))
wbid_apa$cent_lon <- lapply(wbid_apa$centroid, function(x) x[[1]][1])
wbid_apa$cent_lon <- as.numeric(wbid_apa$cent_lon)

for(w in unique(wbid_apa$WBID)){
  apaw_use2[MonitoringLocationIdentifier %in% unique(apaw_use2_sf[subset(wbid_apa, wbid_apa$WBID == w)$geometry, , op = st_intersects]$MonitoringLocationIdentifier), `:=` (WBID = w, MonitoringLocationTypeName = subset(wbid_apa, wbid_apa$WBID == w)$BODY), by = MonitoringLocationIdentifier]
  
}

#Use the hobs_sta object to identify water quality data from within 1km of the HOBS localities
apaCompare3_sf <- st_as_sf(apaCompare3, coords = c("lon", "lat"), crs = 4326)
wq_buffs <- st_buffer(apaCompare3_sf, 1000)
apaw_use2[, `:=` (keepforHOBS = fcase(dec_lon_va > -85.28201 & dec_lat_va > 29.67753 & dec_lon_va < -85.1791 & dec_lat_va < 29.69915, 1,
                                      dec_lon_va > -85.03282 & dec_lat_va > 29.59517 & dec_lon_va < -84.98524 & dec_lat_va < 29.62679, 1,
                                      dec_lon_va > -84.9004 & dec_lat_va > 29.71014 & dec_lon_va < -84.87791 & dec_lat_va < 29.72519, 1,
                                      dec_lon_va > -84.78893 & dec_lat_va > 29.69102 & dec_lon_va < -84.76181 & dec_lat_va < 29.7095, 1,
                                      dec_lon_va > -84.71675 & dec_lat_va > 29.77010 & dec_lon_va < -84.55830 & dec_lat_va < 29.88735, 1,
                                      dec_lon_va > -84.44519 & dec_lat_va > 29.90154 & dec_lon_va < -84.39438 & dec_lat_va < 29.93367, 1,
                                      default = 0),
                  area = fcase(dec_lon_va > -85.28201 & dec_lat_va > 29.67753 & dec_lon_va < -85.1791 & dec_lat_va < 29.69915, "West Apalachicola Bay",
                               dec_lon_va > -85.03282 & dec_lat_va > 29.59517 & dec_lon_va < -84.98524 & dec_lat_va < 29.62679, "Mid Apalachicola Bay",
                               dec_lon_va > -84.9004 & dec_lat_va > 29.71014 & dec_lon_va < -84.87791 & dec_lat_va < 29.72519, "East Apalachicola Bay 2",
                               dec_lon_va > -84.78893 & dec_lat_va > 29.69102 & dec_lon_va < -84.76181 & dec_lat_va < 29.7095, "East Apalachicola Bay 1",
                               dec_lon_va > -84.71675 & dec_lat_va > 29.77010 & dec_lon_va < -84.55830 & dec_lat_va < 29.88735, "Dog Island",
                               dec_lon_va > -84.44519 & dec_lat_va > 29.90154 & dec_lon_va < -84.39438 & dec_lat_va < 29.93367, "Alligator Harbor",
                               default = NA))]
apaw_use2[keepforHOBS == 1, keepInd := seq(1, nrow(apaw_use2[keepforHOBS == 1, ]))]
apaw_use2_sf <- st_as_sf(apaw_use2, coords = c("dec_lon_va", "dec_lat_va"), crs = 4326)

for(i in unique(wq_buffs$samp)){
  apaw_use2[MonitoringLocationIdentifier %in% unique(apaw_use2_sf[subset(wq_buffs, wq_buffs$samp == i)$geometry, , op = st_intersects]$MonitoringLocationIdentifier), nearHOBS := 1, by = MonitoringLocationIdentifier]
}
apaw_use2[is.na(nearHOBS), nearHOBS := 0]

apaw_stas <- distinct(apaw_use2_sf[, c("MonitoringLocationIdentifier", "geometry")])
apaw_stas$geo <- paste0(apaw_stas$geometry)
st_geometry(apaw_stas) <- NULL
apaw_umli <- data.table(geo = unique(apaw_stas$geo),
                        umli = seq(1, length(unique(apaw_stas$geo))))
apaw_stas2 <- merge(apaw_stas, apaw_umli, by = "geo", all.x = TRUE)
apaw_use2 <- merge(apaw_use2, apaw_stas2, by = "MonitoringLocationIdentifier", all.x = TRUE)
# apaw_use2[, month := date(floor_date(ActivityStartDateTime, unit = "months"))]
# apaw_use2[is.na(month), month := date(floor_date(ActivityStartDate, unit = "months"))]
apaw_use2[, month := date(floor_date(asdate, unit = "months"))]



#Look at long Chattahoochee R. flow time series from USGS 02358000 and relationship to salinity---------------------

chatflow <- fread(here::here("usgs02358000.txt"))
names(chatflow) <- as.character(chatflow[1, ])
chatflow <- chatflow[seq(3, nrow(chatflow)), ]
setnames(chatflow, c("26908_00060_00003", "26908_00060_00003_cd", "26909_00065_00003", "26909_00065_00003_cd"), c("discharge_cfs", "discharge_cd", "gageHeight_ft", "gageHeight_cd"))
chatflow[, datetime2 := as_date(datetime)]
apaw_use2b <- merge(apaw_use2, chatflow[, .(datetime2, discharge_cfs, discharge_cd, gageHeight_ft, gageHeight_cd)], by.x = "asdate", by.y = "datetime2", all.x = TRUE, all.y = TRUE)

ggplot(apaw_use2b[parameter == "salinity" & !is.na(ResultMeasureValue) & keepforHOBS == 1, ], aes(x = discharge_cfs, y = ResultMeasureValue)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~area, ncol = 2)

# apaw_use2b[, `:=` (asdate_2d = asdate + 2, 
#                    asdate_4d = asdate + 4, 
#                    asdate_6d = asdate + 6, 
#                    asdate_8d = asdate + 8,
#                    asdate_10d = asdate + 10,
#                    asdate_12d = asdate + 12,
#                    asdate_14d = asdate + 14,
#                    asdate_21d = asdate + 21,
#                    asdate_30d = asdate + 30)]
# 
# apaw_use2c <- copy(apaw_use2b)
# apaw_use2c[, `:=` (asdate_2d = NULL, asdate_4d = NULL, asdate_6d = NULL, asdate_8d = NULL, asdate_10d = NULL, asdate_12d = NULL, asdate_14d = NULL, asdate_21d = NULL, asdate_30d = NULL)]
# for(n in seq(27, 35)){
#   shift_n <- distinct(apaw_use2b[, c(23, ..n)]) #parameter == "salinity" & !is.na(ResultMeasureValue) & keepforHOBS == 1
#   setnames(shift_n, "discharge_cfs", paste0("discharge", str_sub(colnames(apaw_use2b)[n], str_locate(colnames(apaw_use2b)[n], "_")[1], -1)))
#   apaw_use2c <- merge(apaw_use2c, shift_n, by.x = "asdate", by.y = colnames(apaw_use2b)[n], all.x = TRUE, all.y = TRUE)
#   
# }
# 
# apaw_use2c[, `:=` (discharge_cfs = as.numeric(discharge_cfs), discharge_2d = as.numeric(discharge_2d), discharge_4d = as.numeric(discharge_4d), discharge_6d = as.numeric(discharge_6d),
#                    discharge_8d = as.numeric(discharge_8d), discharge_10d = as.numeric(discharge_10d), discharge_12d = as.numeric(discharge_12d), 
#                    discharge_14d = as.numeric(discharge_14d), discharge_21d = as.numeric(discharge_21d), discharge_30d = as.numeric(discharge_30d))]
# apaw_use2c[!is.na(ResultMeasureValue), RMV_mn := mean(ResultMeasureValue), by = list(parameter, keepforHOBS, area, asdate)]
# apaw_use2c[!is.na(ResultMeasureValue), RMV_mn_m := mean(ResultMeasureValue), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_2d), discharge_2dmn_m := mean(discharge_2d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_4d), discharge_4dmn_m := mean(discharge_4d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_6d), discharge_6dmn_m := mean(discharge_6d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_8d), discharge_8dmn_m := mean(discharge_8d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_10d), discharge_10dmn_m := mean(discharge_10d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_12d), discharge_12dmn_m := mean(discharge_12d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_14d), discharge_14dmn_m := mean(discharge_14d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_21d), discharge_21dmn_m := mean(discharge_21d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_30d), discharge_30dmn_m := mean(discharge_30d), by = list(parameter, keepforHOBS, area, month)]
# #medians
# apaw_use2c[!is.na(ResultMeasureValue), RMV_md_m := median(ResultMeasureValue), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_2d), discharge_2dmd_m := median(discharge_2d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_4d), discharge_4dmd_m := median(discharge_4d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_6d), discharge_6dmd_m := median(discharge_6d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_8d), discharge_8dmd_m := median(discharge_8d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_10d), discharge_10dmd_m := median(discharge_10d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_12d), discharge_12dmd_m := median(discharge_12d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_14d), discharge_14dmd_m := median(discharge_14d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_21d), discharge_21dmd_m := median(discharge_21d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_30d), discharge_30dmd_m := median(discharge_30d), by = list(parameter, keepforHOBS, area, month)]
# #minimum sal & maximum flow
# apaw_use2c[!is.na(ResultMeasureValue), RMV_min_m := min(ResultMeasureValue), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_2d), discharge_2dmax_m := max(discharge_2d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_4d), discharge_4dmax_m := max(discharge_4d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_6d), discharge_6dmax_m := max(discharge_6d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_8d), discharge_8dmax_m := max(discharge_8d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_10d), discharge_10dmax_m := max(discharge_10d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_12d), discharge_12dmax_m := max(discharge_12d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_14d), discharge_14dmax_m := max(discharge_14d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_21d), discharge_21dmax_m := max(discharge_21d), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_30d), discharge_30dmax_m := max(discharge_30d), by = list(parameter, keepforHOBS, area, month)]
# #5% quantile sal & 95% quantile flow
# apaw_use2c[!is.na(ResultMeasureValue), RMV_q05_m := quantile(ResultMeasureValue, probs = 0.05), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_cfs), discharge_q95_m := quantile(discharge_cfs, probs = 0.95), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_2d), discharge_2dq95_m := quantile(discharge_2d, probs = 0.95), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_4d), discharge_4dq95_m := quantile(discharge_4d, probs = 0.95), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_6d), discharge_6dq95_m := quantile(discharge_6d, probs = 0.95), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_8d), discharge_8dq95_m := quantile(discharge_8d, probs = 0.95), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_10d), discharge_10dq95_m := quantile(discharge_10d, probs = 0.95), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_12d), discharge_12dq95_m := quantile(discharge_12d, probs = 0.95), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_14d), discharge_14dq95_m := quantile(discharge_14d, probs = 0.95), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_21d), discharge_21dq95_m := quantile(discharge_21d, probs = 0.95), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[!is.na(discharge_30d), discharge_30dq95_m := quantile(discharge_30d, probs = 0.95), by = list(parameter, keepforHOBS, area, month)]
# apaw_use2c[, Month := factor(month(asdate), levels = seq(1, 12))]
# 
# ggplot(apaw_use2c[parameter == "salinity" & !is.na(ResultMeasureValue) & keepforHOBS == 1, ], aes(x = discharge_2d, y = RMV_mn)) +
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~area, ncol = 2)
# 
# ggplot(distinct(apaw_use2c[parameter == "salinity" & !is.na(ResultMeasureValue) & keepforHOBS == 1, .(discharge_q95_m, RMV_q05_m, Month, area)]), 
#        aes(x = discharge_q95_m, y = RMV_q05_m)) +
#   geom_point(aes(color = Month)) +
#   geom_smooth(method = "lm") +
#   stat_regline_equation(aes(label = ..eq.label..), label.x.npc = 0.5, label.y.npc = 0.95) +
#   stat_regline_equation(aes(label = ..rr.label..), label.x.npc = 0.5, label.y.npc = 0.85) + 
#   theme_bw() +
#   scale_x_continuous(trans = "log") +
#   facet_wrap(~area, ncol = 2)
# 
# ggplot(distinct(apaw_use2c[parameter == "salinity" & !is.na(ResultMeasureValue) & keepforHOBS == 1, .(discharge_2dq95_m, RMV_q05_m, Month, area)]), 
#        aes(x = discharge_2dq95_m, y = RMV_q05_m)) +
#   geom_point(aes(color = Month)) +
#   geom_smooth(method = "lm") +
#   stat_regline_equation(aes(label = ..eq.label..), label.x.npc = 0.5, label.y.npc = 0.95) +
#   stat_regline_equation(aes(label = ..rr.label..), label.x.npc = 0.5, label.y.npc = 0.85) +
#   theme_bw() +
#   scale_x_continuous(trans = "log") +
#   facet_wrap(~area, ncol = 2)
# 
# ggplot(distinct(apaw_use2c[parameter == "salinity" & !is.na(ResultMeasureValue) & keepforHOBS == 1, .(discharge_4dq95_m, RMV_q05_m, Month, area)]), 
#        aes(x = discharge_4dq95_m, y = RMV_q05_m)) +
#   geom_point(aes(color = Month)) +
#   geom_smooth(method = "lm") +
#   stat_regline_equation(aes(label = ..eq.label..), label.x.npc = 0.5, label.y.npc = 0.95) +
#   stat_regline_equation(aes(label = ..rr.label..), label.x.npc = 0.5, label.y.npc = 0.85) +
#   theme_bw() +
#   scale_x_continuous(trans = "log") +
#   facet_wrap(~area, ncol = 2)
# 
# ggplot(distinct(apaw_use2c[parameter == "salinity" & !is.na(ResultMeasureValue) & keepforHOBS == 1, .(discharge_6dq95_m, RMV_q05_m, Month, area)]), 
#        aes(x = discharge_6dq95_m, y = RMV_q05_m)) +
#   geom_point(aes(color = Month)) +
#   geom_smooth(method = "lm") +
#   stat_regline_equation(aes(label = ..eq.label..), label.x.npc = 0.5, label.y.npc = 0.95) +
#   stat_regline_equation(aes(label = ..rr.label..), label.x.npc = 0.5, label.y.npc = 0.85) +
#   theme_bw() +
#   scale_x_continuous(trans = "log") +
#   facet_wrap(~area, ncol = 2)
# 
# ggplot(distinct(apaw_use2c[parameter == "salinity" & !is.na(ResultMeasureValue) & keepforHOBS == 1, .(discharge_8dq95_m, RMV_q05_m, Month, area)]), 
#        aes(x = discharge_8dq95_m, y = RMV_q05_m)) +
#   geom_point(aes(color = Month)) +
#   geom_smooth(method = "lm") +
#   stat_regline_equation(aes(label = ..eq.label..), label.x.npc = 0.5, label.y.npc = 0.95) +
#   stat_regline_equation(aes(label = ..rr.label..), label.x.npc = 0.5, label.y.npc = 0.85) +
#   theme_bw() +
#   scale_x_continuous(trans = "log") +
#   facet_wrap(~area, ncol = 2)
# 
# ggplot(distinct(apaw_use2c[parameter == "salinity" & !is.na(ResultMeasureValue) & keepforHOBS == 1, .(discharge_10dq95_m, RMV_q05_m, Month, area)]), 
#        aes(x = discharge_10dq95_m, y = RMV_q05_m)) +
#   geom_point(aes(color = Month)) +
#   geom_smooth(method = "lm") +
#   stat_regline_equation(aes(label = ..eq.label..), label.x.npc = 0.5, label.y.npc = 0.95) +
#   stat_regline_equation(aes(label = ..rr.label..), label.x.npc = 0.5, label.y.npc = 0.85) +
#   theme_bw() +
#   scale_x_continuous(trans = "log") +
#   facet_wrap(~area, ncol = 2)
# 
# ggplot(distinct(apaw_use2c[parameter == "salinity" & !is.na(ResultMeasureValue) & keepforHOBS == 1, .(discharge_12dq95_m, RMV_q05_m, Month, area)]), 
#        aes(x = discharge_12dq95_m, y = RMV_q05_m)) +
#   geom_point(aes(color = Month)) +
#   geom_smooth(method = "lm") +
#   stat_regline_equation(aes(label = ..eq.label..), label.x.npc = 0.5, label.y.npc = 0.95) +
#   stat_regline_equation(aes(label = ..rr.label..), label.x.npc = 0.5, label.y.npc = 0.85) +
#   theme_bw() +
#   scale_x_continuous(trans = "log") +
#   facet_wrap(~area, ncol = 2)
# 
# ggplot(distinct(apaw_use2c[parameter == "salinity" & !is.na(ResultMeasureValue) & keepforHOBS == 1, .(discharge_14dq95_m, RMV_q05_m, Month, area)]), 
#        aes(x = discharge_14dq95_m, y = RMV_q05_m)) +
#   geom_point(aes(color = Month)) +
#   geom_smooth(method = "lm") +
#   stat_regline_equation(aes(label = ..eq.label..), label.x.npc = 0.5, label.y.npc = 0.95) +
#   stat_regline_equation(aes(label = ..rr.label..), label.x.npc = 0.5, label.y.npc = 0.85) +
#   theme_bw() +
#   scale_x_continuous(trans = "log") +
#   facet_wrap(~area, ncol = 2)

#What is the best flow offset for each area?
apaw_use2b[!is.na(ResultMeasureValue), RMV_q05_m := quantile(ResultMeasureValue, probs = 0.05), by = list(parameter, keepforHOBS, area, month)]
saldisregs <- data.table(flowOffset = integer(),
                         area = character(),
                         effect = character(),
                         component = character(),
                         group = character(),
                         level = character(),
                         term = character(),
                         estimate = numeric(),
                         std.error = numeric(),
                         conf.low = numeric(),
                         conf.high = numeric())

#Model salinity v. flow relationship with offsets of +0 to +14 days
for(x in seq(0, 25)){
  saldisregs_x <- data.table(asdate = seq(min(apaw_use2b$asdate), max(apaw_use2b$asdate), by = 1))
  saldisregs_x <- merge(saldisregs_x, distinct(apaw_use2b[!is.na(discharge_cfs), .(asdate, discharge_cfs)]), by = "asdate", all = TRUE)
  saldisregs_x[, asdate := asdate + x]
  saldisregs_x <- merge(saldisregs_x, distinct(apaw_use2b[!is.na(discharge_cfs) & keepforHOBS == 1, .(asdate, parameter, area, RMV_q05_m)]), by = "asdate", all = TRUE)
  saldisregs_x[, `:=` (discharge_cfs = as.numeric(discharge_cfs), area = factor(area, levels = c("West Apalachicola Bay", "Mid Apalachicola Bay", "East Apalachicola Bay 2", "East Apalachicola Bay 1", "Dog Island", "Alligator Harbor"), ordered = TRUE))]
  if(x == 0){
    for(a in c(paste0(saldisregs_x[parameter == "salinity" & !is.na(RMV_q05_m), unique(area)]))){
      if(a == "All"){
        mod_x <- brm(discharge_cfs ~ 1 + RMV_q05_m + (1 + RMV_q05_m | area), family = "lognormal", data = saldisregs_x[parameter == "salinity" & !is.na(RMV_q05_m), ],
                     cores = 4, threads = threading(2), iter = 5000, warmup = 1000, thin = 2, control = c(adapt_delta = 0.9, max_treedepth = 18), backend = "cmdstanr", 
                     seed = 453, file = here::here(paste0("salinity_v_discharge_", formatC(x, width = 2, format = "d", flag = "0"), "off_All.rds")))
        
        # mod_x2 <- brm(discharge_cfs | trunc(l = 0) ~ 1 + RMV_q05_m + (1 + RMV_q05_m | area), family = "gaussian", data = saldisregs_x[parameter == "salinity" & !is.na(RMV_q05_m), ],
        #              cores = 4, threads = 4, iter = 5000, warmup = 1000, thin = 2, control = c(adapt_delta = 0.8, max_treedepth = 10), backend = "cmdstanr")
        #              
        # mod_x3 <- brm(discharge_cfs ~ 1 + RMV_q05_m + (1 + RMV_q05_m | area), family = "lognormal", data = saldisregs_x[parameter == "salinity" & !is.na(RMV_q05_m), ],
        #               cores = 4, threads = 4, iter = 5000, warmup = 1000, thin = 2, control = c(adapt_delta = 0.8, max_treedepth = 10), backend = "cmdstanr")
      } else{
        # mod_x <- brm(discharge_cfs ~ 1 + RMV_q05_m + (1 + RMV_q05_m | area), family = "lognormal", data = saldisregs_x[parameter == "salinity" & !is.na(RMV_q05_m) & area == a, ],
        #              cores = 4, threads = threading(2), iter = 5000, warmup = 1000, thin = 2, control = c(adapt_delta = 0.9, max_treedepth = 12), backend = "cmdstanr", seed = 453,
        #              file = here::here(paste0("salinity_v_discharge_", x, "_off_", str_replace_all(a, " ", ""), ".rds")))
        mod_x <- brm(formula = RMV_q05_m ~ 1 + discharge_cfs, family = "lognormal", data = saldisregs_x[parameter == "salinity" & !is.na(RMV_q05_m) & area == a, ], 
                     cores = 4, threads = threading(2), iter = 5000, warmup = 1000, thin = 2, control = c(adapt_delta = 0.9, max_treedepth = 18), backend = "cmdstanr", 
                     seed = 453, file = here::here(paste0("discharge_v_salinity_", formatC(x, width = 2, format = "d", flag = "0"), "off_", str_replace_all(a, " ", ""), ".rds")))
      }
      
      mod_x_res <- tidy(mod_x, effects = c("fixed", "ran_vals", "ran_pars"))
      setDT(mod_x_res)
      mod_x_res[, `:=` (flowOffset = x, area = a)]
      
      saldisregs <- rbind(saldisregs, mod_x_res)
    }
    
  } else{
    for(a in c(paste0(saldisregs_x[parameter == "salinity" & !is.na(RMV_q05_m), unique(area)]))){
      if(a == "All"){
        mod_x <- update(mod_x, formula = discharge_cfs ~ 1 + RMV_q05_m + (1 + RMV_q05_m | area), family = "lognormal", newdata = saldisregs_x[parameter == "salinity" & !is.na(RMV_q05_m), ], 
                        cores = 4, threads = threading(2), iter = 5000, warmup = 1000, thin = 2, control = c(adapt_delta = 0.9, max_treedepth = 18), backend = "cmdstanr", 
                        seed = 453, file = here::here(paste0("salinity_v_discharge_", formatC(x, width = 2, format = "d", flag = "0"), "off_All.rds")))
      } else{
        mod_x <- update(mod_x, formula = RMV_q05_m ~ 1 + discharge_cfs, family = "lognormal", newdata = saldisregs_x[parameter == "salinity" & !is.na(RMV_q05_m) & area == a, ], 
                        cores = 4, threads = threading(2), iter = 5000, warmup = 1000, thin = 2, control = c(adapt_delta = 0.9, max_treedepth = 18), backend = "cmdstanr", 
                        seed = 453, file = here::here(paste0("discharge_v_salinity_", formatC(x, width = 2, format = "d", flag = "0"), "off_", str_replace_all(a, " ", ""), ".rds")))
      }
  
      mod_x_res <- tidy(mod_x, effects = c("fixed", "ran_vals", "ran_pars"))
      setDT(mod_x_res)
      mod_x_res[, `:=` (flowOffset = x, area = a)]
      
      saldisregs <- rbind(saldisregs, mod_x_res)
    }
  }
}

saldisregs[, sig := ifelse(conf.low < 0 & conf.high > 0, 0, 1)]
saldisregs_sum <- saldisregs[!is.na(level), .(minSE = min(std.error)), by = list(term, level)]
saldisregs_sum[, `:=` (flowOffset = saldisregs[!is.na(level) & std.error == minSE, flowOffset], sig = saldisregs[!is.na(level) & std.error == minSE, sig]), by = minSE]

#All of the lowest SD's for area effects occurred with either model 11 or model 12; does one of these models fit best overall based on the LOO criterion?
bestmods <- data.table(area = character(), bestmod = list())

for(a in c("All", paste0(saldisregs_x[parameter == "salinity" & !is.na(RMV_q05_m), unique(area)]))){
  aacro <- fcase(a == "All", "All",
                 a == "West Apalachicola Bay", "WAB",
                 a == "East Apalachicola Bay 2", "EAB2",
                 a == "Alligator Harbor", "AH",
                 a == "Mid Apalachicola Bay", "MAB",
                 a == "East Apalachicola Bay 1", "EAB1",
                 a == "Dog Island", "DI")
  mods <- list.files(here::here())[which(str_detect(list.files(here::here()), "salinity_v_discharge_") & str_detect(list.files(here::here()), ifelse(a == "All", "All.rds", str_replace_all(a, " ", ""))))]
  modobjs <- c()
  for(mod in mods){
    assign(paste0("off", formatC(which(mods == mod) - 1, width = 2, format = "d", flag = "0")), readRDS(here::here(paste0(mod))))
    modobjs <- append(modobjs, as.name(paste0("off", formatC(which(mods == mod) - 1, width = 2, format = "d", flag = "0")))) 
  }
  
  assign(paste0(aacro, "_sums"), lapply(modobjs, function(x) summary(eval(x)))) #check rhats
  hasloo <- lapply(modobjs, function(x) length(eval(x)$criteria))
  lapply(modobjs[which(hasloo == 0)], function(x) assign(paste0(x), add_criterion(eval(x), criterion = "loo", model_name = paste0(x), overwrite = TRUE, file = str_sub(mods[which(str_detect(mods, paste0("_", str_sub(paste0(x), 4, -1), "off")))], 1, -5), force_save = TRUE)))
  
  # #Potentially faster to use foreach, but I didn't finish this code
  # library(doFuture)
  # plan(multisession, workers = 14)
  # foreach(m = modobjs) %dofuture% {
  #   assign(paste0(m), add_criterion(eval(m), criterion = "loo", model_name = paste0(m), overwrite = TRUE, file = str_sub(modobjs[which(str_detect(modobjs, paste0("_", str_sub(paste0(m), 4, -1), "_")))], 1, -5), force_save = TRUE))
  # }
  
  for(mod in mods){
    assign(paste0("off", formatC(which(mods == mod) - 1, width = 2, format = "d", flag = "0")), readRDS(here::here(paste0(mod))))
  }
  assign(paste0(aacro, "_loo"), loo_compare(off00, off01, off02, off03, off04, off05, off06, off07, off08, off09, off10, off11, off12, off13, off14, off15, off16, off17, off18, off19, off20, off21, off22, off23, off24, off25))
  bestmod_a <- data.table(area = a,
                          bestmod = list(list(names(which(eval(as.name(paste0(aacro, "_loo")))[, 1] == 0 & eval(as.name(paste0(aacro, "_loo")))[, 2] == 0)))))
  bestmods <- rbind(bestmods, bestmod_a)
  rm(off00, off01, off02, off03, off04, off05, off06, off07, off08, off09, off10, off11, off12, off13, off14, off15, off16, off17, off18, off19, off20, off21, off22, off23, off24, off25)
  print(paste0("Done with ", a))
}

View(bestmods)
offsets <- setorder(distinct((apaw_use2b[!is.na(area), .(avg_lon = mean(dec_lon_va)), by = area])), avg_lon)
offsets <- merge(offsets, bestmods, by = "area")
offsets[, `:=` (area2 = area, 
                offset1 = fcase(area == "Alligator Harbor", 3, 
                                area == "Dog Island", 20, 
                                area == "East Apalachicola Bay 1", 4, 
                                area == "East Apalachicola Bay 2", 14, 
                                area == "Mid Apalachicola Bay", 8, 
                                area == "West Apalachicola Bay", 21), 
                offset2 = fcase(area == "Alligator Harbor", 3, 
                                area == "Dog Island", 20, 
                                area == "East Apalachicola Bay 1", 13, 
                                area == "East Apalachicola Bay 2", 14, 
                                area == "Mid Apalachicola Bay", 8, 
                                area == "West Apalachicola Bay", 21))]


saldis_woff <- data.table(area = rep(unique(offsets$area), each = length(seq(min(apaw_use2b$asdate), max(apaw_use2b$asdate), by = 1))), 
                          asdate = rep(seq(min(apaw_use2b$asdate), max(apaw_use2b$asdate), by = 1), length(unique(offsets$area))))
saldis_woff <- merge(saldis_woff, distinct(apaw_use2b[!is.na(discharge_cfs), .(asdate, discharge_cfs)]), by = c("asdate"), all = TRUE)
saldis_woff[, asdate_woff := asdate + offsets$offset1, by = asdate]
saldis_woff[, asdate := NULL]
setnames(saldis_woff, "asdate_woff", "asdate")
saldis_woff <- merge(saldis_woff, distinct(apaw_use2b[parameter == "salinity" & keepforHOBS == 1 & !is.na(ResultMeasureValue) & !is.na(area), .(asdate, parameter, area, ResultMeasureValue, RMV_q05_m, WBID, umli, month, dec_lon_va, dec_lat_va, geo)]), by = c("asdate", "area"), all.x = TRUE)
setorder(saldis_woff, asdate, area)


apaCompare3[, date2 := fcase(date == "Feb-38", as_date("1938-02-01"),
                             date == "Jul-41", as_date("1941-07-01"),
                             date == "Aug-41", as_date("1941-08-01"),
                             date == "16-Dec", as_date("2016-12-01"),
                             date == "17-Apr", as_date("2017-04-01"),
                             date == "9/3/2020", as_date("2020-09-01"),
                             date == "8/11/2020", as_date("2020-08-01"),
                             date == "7/30/2019", as_date("2019-08-01"),
                             date == "1/31/1938", as_date("1938-02-01")), by = date]

saldis_woff[, `:=` (area2 = area, year = year(asdate), discharge_cfs = as.numeric(discharge_cfs))]
apaCompare3[, yr := year(date2)]
apaCompare3[sta2 == "East Apalachicola Bay", sta2 := ifelse(sta == "Cat Point", "East Apalachicola Bay 2", "East Apalachicola Bay 1")]
apaCompare3[, sta2 := factor(sta2, levels = c("West Apalachicola Bay", "Mid Apalachicola Bay", "East Apalachicola Bay 2", "East Apalachicola Bay 1", "Dog Island", "Alligator Harbor"), ordered = TRUE)]

surftest <- distinct(saldis_woff[!is.na(discharge_cfs) & !is.na(year), .(lon_mn = mean(dec_lon_va, na.rm = TRUE), flow_mn = mean(discharge_cfs), q05_mn = mean(RMV_q05_m, na.rm = TRUE)), by = list(area, year)])
surftest[, lon_mn := plyr::round_any(lon_mn, 0.1)]
apaC3 <- distinct(apaCompare3[, .(lon_mn = plyr::round_any(lon2, 0.1), deadC_mn = mean(deadC)), by = list(sta2, yr)])
setnames(apaC3, c("sta2", "yr"), c("area", "year"))
surftest2 <- merge(surftest, apaC3, by = c("area", "year"), all = TRUE)
surftest2[is.nan(lon_mn.x), lon_mn.x := NA]
surftest2[is.nan(q05_mn), q05_mn := NA]
surftest2[is.na(lon_mn.x), lon_mn.x := lon_mn.y]
surftest2[, lon_mn.y := NULL]
setnames(surftest2, "lon_mn.x", "lon_mn")
surftest2[, relyear := year - min(year)] # min(year) is 1922

deadcmod2 <- brm(formula = deadC_mn ~ 1 + lon_mn * flow_mn + (1 | relyear), family = gaussian, data = surftest2[!is.na(deadC_mn), ],
                 cores = 4, chains = 4, warmup = 1000, iter = 5000, thin = 2, control = c(adapt_delta = 0.99, max_treedepth = 18), 
                 backend = "cmdstanr", threads = threading(2), seed = 453, file = here::here("deadcmod_v2.rds"))

for(i in seq(1, length(unique(apaCompare3$date2)))){
  date_i <- unique(apaCompare3$date2)[i]
  if(year(date_i) %in% unique(saldis_woff$year)){
    stas <- apaCompare3[date2 == date_i, unique(sta2)]
    saldis_woff[year == year(date_i) & area %in% stas, use_exact := 1]
    saldis_woff[year == year(date_i) & area %in% stas, `:=` (use_prioryr = 1,
                                                             deadC_mn = apaCompare3[yr == year(date_i) & sta2 == area, mean(deadC)],
                                                             deadC_sd = apaCompare3[yr == year(date_i) & sta2 == area, sd(deadC)],
                                                             flow_mn = saldis_woff[year <= year(date_i) & year >= year(date_i) - 1 & !is.na(discharge_cfs) & area2 == area, mean(discharge_cfs)],
                                                             flow_sd = saldis_woff[year <= year(date_i) & year >= year(date_i) - 1 & !is.na(discharge_cfs) & area2 == area, sd(discharge_cfs)],
                                                             sal_mn = saldis_woff[year <= year(date_i) & year >= year(date_i) - 1 & !is.na(ResultMeasureValue) & area2 == area, mean(ResultMeasureValue)],
                                                             sal_sd = saldis_woff[year <= year(date_i) & year >= year(date_i) - 1 & !is.na(ResultMeasureValue) & area2 == area, sd(ResultMeasureValue)],
                                                             q05sal_mn = saldis_woff[year <= year(date_i) & year >= year(date_i) - 1 & !is.na(RMV_q05_m) & area2 == area, mean(RMV_q05_m)],
                                                             q05sal_sd = saldis_woff[year <= year(date_i) & year >= year(date_i) - 1 & !is.na(RMV_q05_m) & area2 == area, sd(RMV_q05_m)]),
                by = area]
  }
}

View(saldis_woff[!is.na(q05sal_mn), ])




# # model Chat. R. flow & salinity v. deadC-----------------------------------------
# salmod <- brm(formula = deadC_mn ~ 1 + me(minsal_mn, minsal_sd) + (1 | year:area), family = gaussian, data = distinct(saldis_woff[use_prioryr == 1 & !is.na(salinity), .(year, sal_mn, sal_sd, minsal_mn, minsal_sd, deadC_mn, area)]),
#               cores = 4, chains = 4, iter = 5000, control = c(adapt_delta = 0.99, max_treedepth = 20), backend = "cmdstanr", threads = threading(2), seed = 453, file = here::here("salmod.rds"))
# summary(salmod)
# 
# salmod2 <- brm(formula = deadC_mn ~ 1 + minsal_mn + (1 | year:area), family = gaussian, data = distinct(saldis_woff[use_prioryr == 1 & !is.na(salinity), .(year, sal_mn, sal_sd, minsal_mn, minsal_sd, deadC_mn, area)]),
#                cores = 4, chains = 4, iter = 5000, control = c(adapt_delta = 0.99, max_treedepth = 20), backend = "cmdstanr", threads = threading(2), seed = 453, file = here::here("salmod2.rds"))
# summary(salmod2)
# 
# flowmod1 <- brm(formula = deadC_mn ~ 1 + flow_mn + (1 | year:area), family = gaussian, data = distinct(saldis_woff[use_prioryr == 1 & !is.na(flow_mn), .(year, flow_mn, deadC_mn, area)]),
#                 cores = 4, chains = 4, iter = 5000, control = c(adapt_delta = 0.99, max_treedepth = 20), backend = "cmdstanr", threads = threading(2), seed = 453, file = here::here("flowmod1.rds"))
# summary(flowmod1)
# 
# # plausibility check
# testdat <- data.table(minsal_mn = seq(0, 30), year = 2016, area = "East Apalachicola Bay 2")
# testdat2 <- cbind(testdat, predict(salmod2, testdat)) #sal when deadC is ~0.078 is ~5

# model longitude v. deadC-----------------------------------------

ggplot(data = distinct(apaCompare3[, c("sta2", "deadC", "year", "study")]), aes(x = sta2, y = deadC, fill = as.factor(year), shape = study)) +
  geom_point(size = 2, color = "black") +
  geom_label(aes(label = year), position = position_jitter(width = 0.55, height = 0.007), size = 2) +
  ylim(-0.02, 0.08) +
  theme_bw() +
  labs(x = "Locality", fill = "Year", shape = "Study") +
  scale_fill_discrete_sequential(palette = "inferno", nmax = 10, order = c(2:7), breaks = c("1938", "1941", "2016", "2017", "2019", "2020")) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  scale_shape_manual(values = c("Durham et al. 2023" = 21, "Hadden & Cherkinsky 2017" = 22, "Hadden et al. 2018" = 23)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = distinct(apaCompare3[, c("lon", "deadC", "year", "study")]), aes(x = lon, y = deadC, fill = as.factor(year), shape = study)) +
  geom_point(size = 2, color = "black") +
  # geom_label(aes(label = year), position = position_jitter(width = 0.55, height = 0.007), size = 2) +
  ylim(-0.02, 0.08) +
  theme_bw() +
  labs(x = "Longitude", fill = "Year", shape = "Study") +
  scale_fill_discrete_sequential(palette = "inferno", nmax = 10, order = c(2:7), breaks = c("1938", "1941", "2016", "2017", "2019", "2020")) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  scale_shape_manual(values = c("Durham et al. 2023" = 21, "Hadden & Cherkinsky 2017" = 22, "Hadden et al. 2018" = 23)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

apaCompare3[, postbomb := ifelse(decade %in% c("2010s", "2020s"), 1, 0)]
apaCompare3[, postbomb := as.logical(postbomb)]
apaCompare3[, lon2 := fcase(sta2 == "Alligator Harbor", -84.41000,
                            sta2 == "Dog Island", -84.61000,
                            sta2 == "East Apalachicola Bay 1", -84.77293,
                            sta2 == "East Apalachicola Bay 2", -84.89000,
                            sta2 == "Mid Apalachicola Bay", -85.00554,
                            sta2 == "West Apalachicola Bay", -85.22000)]

# lonmod <- brm(formula = deadC ~ 1 + lon + (1 | postbomb), family = gaussian, data = apaCompare3,
#               cores = 4, chains = 4, warmup = 1000, iter = 5000, thin = 2, control = c(adapt_delta = 0.99, max_treedepth = 18), backend = "cmdstanr", threads = threading(2), seed = 453, file = here::here("lonmod.rds"))
# summary(lonmod)
# 
# lonmod2 <- brm(formula = deadC ~ 1 + lon2 + (1 | postbomb), family = gaussian, data = apaCompare3,
#                cores = 4, chains = 4, warmup = 1000, iter = 5000, thin = 2, control = c(adapt_delta = 0.99, max_treedepth = 18), backend = "cmdstanr", threads = threading(2), seed = 453, file = here::here("lonmod2.rds"))
# summary(lonmod2)
# 
# lonmod3 <- brm(formula = deadC ~ 1 + lon:postbomb, family = gaussian, data = apaCompare3,
#                cores = 4, chains = 4, warmup = 1000, iter = 5000, thin = 2, control = c(adapt_delta = 0.99, max_treedepth = 18), backend = "cmdstanr", threads = threading(2), seed = 453, file = here::here("lonmod3.rds"))
# summary(lonmod3)

surftest <- distinct(saldis_woff[!is.na(dec_lon_va) & !is.na(discharge_cfs) & !is.na(deadC_mn) & !is.na(year) & !is.na(RMV_q05_m), .(area, dec_lon_va, discharge_cfs, RMV_q05_m, deadC_mn, year)])

deadcmod <- brm(formula = deadC_mn ~ 1 + dec_lon_va * minsal_mn + (1 | year), family = gaussian, data = surftest,
                cores = 4, chains = 4, warmup = 1000, iter = 5000, thin = 2, control = c(adapt_delta = 0.99, max_treedepth = 18), 
                backend = "cmdstanr", threads = threading(2), seed = 453, file = here::here("deadcmod.rds"))

surftest[, relyear := year - min(surftest$year)]

surfdat <- data.table(dec_lon_va = rep(seq(-85.222, -84.401, 0.05), each = 16), minsal_mn = rep(seq(0, 30, 2), 17), year = 2016)
surfdat2 <- cbind(surfdat, predict(deadcmod, surfdat))
Dead_C <- matrix(surfdat2$Estimate, nrow = 17, ncol = 16, byrow = TRUE,
                 dimnames = list(c(unique(surfdat2$dec_lon_va)),
                                 c(unique(surfdat2$minsal_mn))))

surfplot <- plot_ly(x = unique(surfdat2$minsal_mn), 
                    y = unique(surfdat2$dec_lon_va), 
                    z = ~Dead_C,
                    contours = list(z = list(show = TRUE,
                                             start = round(min(Dead_C), 2),
                                             end = round(max(Dead_C), 2),
                                             size = 0.01,
                                             usecolormap=TRUE,
                                             highlightcolor="#ff0000",
                                             project = list(z = TRUE)))) %>% 
  add_surface(contours = list(z = list(show = TRUE,
                                       start = round(min(Dead_C), 2),
                                       end = round(max(Dead_C), 2),
                                       size = 0.01,
                                       usecolormap=TRUE,
                                       highlightcolor="#ff0000",
                                       project = list(z = TRUE)))) %>% 
  layout(scene = list(xaxis = list(title = "Salinity"),
                      yaxis = list(title = "Longitude"),
                      zaxis = list(title = "Dead C"))) #only works if "~" is removed from z argument in plot_ly()


surfdat3 <- data.table(dec_lon_va = rep(seq(-85.222, -84.401, 0.05), each = 16*9), minsal_mn = rep(seq(0, 30, 2), 17*9), year = rep(seq(1940, 2020, 10), each = 16))
surfdat4 <- cbind(surfdat3, predict(deadcmod, surfdat3, allow_new_levels = TRUE))

yr <- 1940
Dead_C <- matrix(surfdat4[year == yr, Estimate], nrow = 17, ncol = 16, byrow = TRUE,
                 dimnames = list(c(surfdat4[year == yr, unique(dec_lon_va)]),
                                 c(surfdat4[year == yr, unique(minsal_mn)])))

surfplot2 <- plot_ly(x = surfdat4[year == yr, unique(minsal_mn)], 
                     y = surfdat4[year == yr, unique(dec_lon_va)], 
                     z = ~Dead_C,
                     contours = list(z = list(show = TRUE,
                                              start = round(min(Dead_C), 2),
                                              end = round(max(Dead_C), 2),
                                              size = 0.01,
                                              usecolormap=TRUE,
                                              highlightcolor="#ff0000",
                                              project = list(z = TRUE)))) %>% 
  add_surface(contours = list(z = list(show = TRUE,
                                       start = round(min(Dead_C), 2),
                                       end = round(max(Dead_C), 2),
                                       size = 0.01,
                                       usecolormap=TRUE,
                                       highlightcolor="#ff0000",
                                       project = list(z = TRUE)))) %>% 
  layout(scene = list(title = yr,
                      xaxis = list(title = "Salinity"),
                      yaxis = list(title = "Longitude"),
                      zaxis = list(title = "Dead C"))) #only works if "~" is removed from z argument in plot_ly()


surfplot3 <- surfplot %>% add_surface(x = surfdat4[year == yr, unique(minsal_mn)], 
                                      y = surfdat4[year == yr, unique(dec_lon_va)], 
                                      z = Dead_C + 0.3)























#split data components by parameter and merge----------------------

apaw_flow <- distinct(apaw_use2[parameter == "flow", .(MonitoringLocationIdentifier, umli, parameter, area, keepforHOBS, keepInd, nearHOBS, month, dec_lon_va,
                                                       dec_lat_va, ResultMeasureValue, ResultMeasure.MeasureUnitCode, WBID, MonitoringLocationTypeName)])
apaw_flow[, RMV := mean(ResultMeasureValue), by = list(umli, month, parameter)][, ResultMeasureValue := RMV][, RMV := NULL]
apaw_flow <- distinct(apaw_flow[, .(umli, month, area, keepforHOBS, keepInd, nearHOBS, ResultMeasureValue, dec_lat_va, dec_lon_va, MonitoringLocationTypeName, WBID)])
setnames(apaw_flow, c("ResultMeasureValue"), c("flow"))
# apaw_flow[, `:=` (lat_rd = plyr::round_any(dec_lat_va, 0.01), lon_rd = plyr::round_any(dec_lon_va, 0.01), date = date(ActivityStartDateTime))]
# apaw_flow[, RMV := mean(flow), by = month][, flow := RMV][, RMV := NULL]

apaw_carbon <- distinct(apaw_use2[parameter == "carbon", .(MonitoringLocationIdentifier, umli, parameter, area, keepforHOBS, keepInd, nearHOBS, month, dec_lon_va,
                                                           dec_lat_va, ResultMeasureValue, ResultMeasure.MeasureUnitCode, WBID, MonitoringLocationTypeName)])
apaw_carbon[, RMV := mean(ResultMeasureValue), by = list(umli, month, parameter)][, ResultMeasureValue := RMV][, RMV := NULL]
apaw_carbon <- distinct(apaw_carbon[, .(umli, month, area, keepforHOBS, keepInd, nearHOBS, ResultMeasureValue, dec_lat_va, dec_lon_va, MonitoringLocationTypeName, WBID)])
setnames(apaw_carbon, c("ResultMeasureValue"), c("carbon"))
# apaw_carbon[, `:=` (lat_rd = plyr::round_any(dec_lat_va, 0.01), lon_rd = plyr::round_any(dec_lon_va, 0.01), date = date(ActivityStartDateTime))]
# apaw_carbon[, RMV := mean(carbon), by = date][, carbon := RMV][, RMV := NULL]

apaw_salinity <- distinct(apaw_use2[parameter == "salinity", .(MonitoringLocationIdentifier, umli, parameter, area, keepforHOBS, keepInd, nearHOBS, month, dec_lon_va,
                                                               dec_lat_va, ResultMeasureValue, ResultMeasure.MeasureUnitCode, WBID, MonitoringLocationTypeName)])
apaw_salinity[, RMV := mean(ResultMeasureValue), by = list(umli, month, parameter)][, ResultMeasureValue := RMV][, RMV := NULL]
apaw_salinity <- distinct(apaw_salinity[, .(umli, month, area, keepforHOBS, keepInd, nearHOBS, ResultMeasureValue, dec_lat_va, dec_lon_va, MonitoringLocationTypeName, WBID)])
setnames(apaw_salinity, c("ResultMeasureValue"), c("salinity"))
# apaw_salinity[, `:=` (lat_rd = plyr::round_any(dec_lat_va, 0.01), lon_rd = plyr::round_any(dec_lon_va, 0.01), date = date(ActivityStartDateTime))]

apaw_temperature <- distinct(apaw_use2[parameter == "temperature", .(MonitoringLocationIdentifier, umli, parameter, area, keepforHOBS, keepInd, nearHOBS, month, dec_lon_va,
                                                                     dec_lat_va, ResultMeasureValue, ResultMeasure.MeasureUnitCode, WBID, MonitoringLocationTypeName)])
apaw_temperature[, RMV := mean(ResultMeasureValue), by = list(umli, month, parameter)][, ResultMeasureValue := RMV][, RMV := NULL]
apaw_temperature <- distinct(apaw_temperature[, .(umli, month, area, keepforHOBS, keepInd, nearHOBS, ResultMeasureValue, dec_lat_va, dec_lon_va, MonitoringLocationTypeName, WBID)])
setnames(apaw_temperature, c("ResultMeasureValue"), c("temperature"))
# apaw_temperature[, `:=` (lat_rd = plyr::round_any(dec_lat_va, 0.01), lon_rd = plyr::round_any(dec_lon_va, 0.01), date = date(ActivityStartDateTime))]

apaw_tss <- distinct(apaw_use2[parameter == "total suspended solids, tss", .(MonitoringLocationIdentifier, umli, parameter, area, keepforHOBS, keepInd, nearHOBS, month, dec_lon_va,
                                                                     dec_lat_va, ResultMeasureValue, ResultMeasure.MeasureUnitCode, WBID, MonitoringLocationTypeName)])
apaw_tss[, RMV := mean(ResultMeasureValue), by = list(umli, month, parameter)][, ResultMeasureValue := RMV][, RMV := NULL]
apaw_tss <- distinct(apaw_tss[, .(umli, month, area, keepforHOBS, keepInd, nearHOBS, ResultMeasureValue, dec_lat_va, dec_lon_va, MonitoringLocationTypeName, WBID)])
setnames(apaw_tss, c("ResultMeasureValue"), c("tss"))

apaw_turbidity <- distinct(apaw_use2[parameter == "turbidity", .(MonitoringLocationIdentifier, umli, parameter, area, keepforHOBS, keepInd, nearHOBS, month, dec_lon_va,
                                                                             dec_lat_va, ResultMeasureValue, ResultMeasure.MeasureUnitCode, WBID, MonitoringLocationTypeName)])
apaw_turbidity[, RMV := mean(ResultMeasureValue), by = list(umli, month, parameter)][, ResultMeasureValue := RMV][, RMV := NULL]
apaw_turbidity <- distinct(apaw_turbidity[, .(umli, month, area, keepforHOBS, keepInd, nearHOBS, ResultMeasureValue, dec_lat_va, dec_lon_va, MonitoringLocationTypeName, WBID)])
setnames(apaw_turbidity, c("ResultMeasureValue"), c("turbidity"))

apaw_use3 <- merge(apaw_temperature, apaw_salinity, 
                   by = c("umli", "month", "area", "keepforHOBS", "keepInd", "nearHOBS", "dec_lat_va", "dec_lon_va", "MonitoringLocationTypeName", "WBID"), all.x = TRUE, all.y = TRUE)
apaw_use3 <- merge(apaw_use3, apaw_carbon, 
                   by = c("umli", "month", "area", "keepforHOBS", "keepInd", "nearHOBS", "dec_lat_va", "dec_lon_va", "MonitoringLocationTypeName", "WBID"), all.x = TRUE, all.y = TRUE)
apaw_use3 <- merge(apaw_use3, apaw_flow, 
                   by = c("umli", "month", "area", "keepforHOBS", "keepInd", "nearHOBS", "dec_lat_va", "dec_lon_va", "MonitoringLocationTypeName", "WBID"), all.x = TRUE, all.y = TRUE)
apaw_use3 <- merge(apaw_use3, apaw_tss, 
                   by = c("umli", "month", "area", "keepforHOBS", "keepInd", "nearHOBS", "dec_lat_va", "dec_lon_va", "MonitoringLocationTypeName", "WBID"), all.x = TRUE, all.y = TRUE)
apaw_use3 <- merge(apaw_use3, apaw_turbidity, 
                   by = c("umli", "month", "area", "keepforHOBS", "keepInd", "nearHOBS", "dec_lat_va", "dec_lon_va", "MonitoringLocationTypeName", "WBID"), all.x = TRUE, all.y = TRUE)

# apaw_use3[, `:=` (lat_rd2 = ordered(lat_rd), lon_rd2 = ordered(lon_rd))]
# apaw_use3 <- apaw_use3[!is.na(lon_rd) & !is.na(lat_rd), ]
apaw_use3 <- apaw_use3[!is.na(dec_lon_va) & !is.na(dec_lat_va), ]
apaw_use3 <- merge(apaw_use3, wbid_apa[, c("WBID", "cent_lon")], by = "WBID", all.x = TRUE)
wbidorder <- wbid_apa[, c("WBID", "cent_lon")]
setorder(wbidorder, cent_lon)
apaw_use3[, WBID2 := factor(WBID, levels = wbidorder$WBID)]
apaw_use3_sf <- st_as_sf(apaw_use3, coords = c("dec_lon_va", "dec_lat_va"), crs = 4326)

#are turbidity and salinity and flow correlated?
tsf <- apaw_use3[!is.na(turbidity) & !is.na(salinity), ]
ggplot(tsf, aes(x = salinity, y = turbidity)) +
  geom_point() +
  theme_bw()

ggplot(tsf, aes(x = salinity, y = turbidity, color = as.factor(umli))) +
  geom_point() +
  guides(color = "none") +
  theme_bw()

ggplot(tsf[umli %in% names(sort(table(tsf$umli)))[(length(names(sort(table(tsf$umli)))) - 5):length(names(sort(table(tsf$umli))))], ], aes(x = salinity, y = turbidity, color = as.factor(umli))) +
  geom_point() +
  guides(color = "none") +
  theme_bw() +
  facet_wrap(~umli)


tsf2 <- apaw_use3[!is.na(tss) & !is.na(salinity), ]
ggplot(tsf2, aes(x = salinity, y = tss)) +
  geom_point() +
  theme_bw()

apaw_stas2_acf <- subset(apaw_sites3, str_detect(apaw_sites3$MonitoringLocationDescriptionText, ".hattahoochee") | 
                                      str_detect(apaw_sites3$MonitoringLocationDescriptionText, ".lint") | 
                                      str_detect(apaw_sites3$MonitoringLocationDescriptionText, ".palachicola") | 
                                      str_detect(apaw_sites3$MonitoringLocationDescriptionText, ".hipola") |
                                      str_detect(apaw_sites3$MonitoringLocationDescriptionText, ".arrabelle") |
                                      str_detect(apaw_sites3$MonitoringLocationDescriptionText, "New R") |
                                      str_detect(apaw_sites3$MonitoringLocationDescriptionText, ".rooked") |
                                      str_detect(apaw_sites3$MonitoringLocationDescriptionText, "St. George Sound") |
                                      str_detect(apaw_sites3$station_nm, ".hattahoochee") | 
                                      str_detect(apaw_sites3$station_nm, ".lint") | 
                                      str_detect(apaw_sites3$station_nm, ".palachicola") | 
                                      str_detect(apaw_sites3$station_nm, ".hipola") |
                                      str_detect(apaw_sites3$station_nm, ".arrabelle") |
                                      str_detect(apaw_sites3$station_nm, "New R") |
                                      str_detect(apaw_sites3$station_nm, ".rooked") |
                                      str_detect(apaw_sites3$station_nm, "St. George Sound") |
                                      str_detect(apaw_sites3$MonitoringLocationDescriptionText, ".HATTAHOOCHEE") | 
                                      str_detect(apaw_sites3$MonitoringLocationDescriptionText, ".LINT") | 
                                      str_detect(apaw_sites3$MonitoringLocationDescriptionText, ".PALACHICOLA") | 
                                      str_detect(apaw_sites3$MonitoringLocationDescriptionText, ".HIPOLA") |
                                      str_detect(apaw_sites3$MonitoringLocationDescriptionText, ".ARRABELLE") |
                                      str_detect(apaw_sites3$MonitoringLocationDescriptionText, "NEW R") |
                                      str_detect(apaw_sites3$MonitoringLocationDescriptionText, ".ROOKED") |
                                      str_detect(apaw_sites3$MonitoringLocationDescriptionText, "ST. GEORGE SOUND") |
                                      str_detect(apaw_sites3$station_nm, ".HATTAHOOCHEE") | 
                                      str_detect(apaw_sites3$station_nm, ".LINT") | 
                                      str_detect(apaw_sites3$station_nm, ".PALACHICOLA") |
                                      str_detect(apaw_sites3$station_nm, ".HIPOLA") |
                                      str_detect(apaw_sites3$station_nm, ".ARRABELLE") |
                                      str_detect(apaw_sites3$station_nm, "NEW R") |
                                      str_detect(apaw_sites3$station_nm, ".ROOKED") |
                                      str_detect(apaw_sites3$station_nm, "ST. GEORGE SOUND") |
                                      apaw_sites3$MonitoringLocationIdentifier %in% unique(subset(apaw_use2, !is.na(apaw_use2$area))$MonitoringLocationIdentifier) |
                                      apaw_sites3$MonitoringLocationIdentifier %in% unique(apadisc_apaw$MonitoringLocationIdentifier))[, c("site_no", "station_nm", "MonitoringLocationIdentifier", "MonitoringLocationName", "MonitoringLocationTypeName", "MonitoringLocationDescriptionText")]
apaw_stas2_acf <- merge(apaw_stas2_acf, apaw_stas2, by = "MonitoringLocationIdentifier", all.x = TRUE)
apaw_stas2_acf <- distinct(apaw_stas2_acf)

mapview(apaw_use3_sf, col.regions = "dodgerblue", layer.name = "all") +
  mapview(subset(apaw_use3_sf, apaw_use3_sf$umli %in% unique(apaw_stas2_acf$umli)), col.regions = "firebrick", layer.name = "acf") #+
  #mapview(subset(apaw_use3_sf, apaw_use3_sf$umli == 3684), col.regions = "green", layer.name = "1907")

apaw_stas2_acf2 <- copy(apaw_stas2_acf)
st_geometry(apaw_stas2_acf2) <- NULL
setDT(apaw_stas2_acf2)
apaw_stas2_acf2[, watershed := fcase(str_detect(MonitoringLocationDescriptionText, ".hattahoochee"), "Chattahoochee River",
                                     str_detect(MonitoringLocationDescriptionText, ".lint"), "Flint River",
                                     str_detect(MonitoringLocationDescriptionText, ".palachicola"), "Apalachicola River",
                                     str_detect(MonitoringLocationDescriptionText, ".hipola"), "Chipola River",
                                     str_detect(MonitoringLocationDescriptionText, ".arrabelle"), "Carrabelle River",
                                     str_detect(MonitoringLocationDescriptionText, "New R"), "New River",
                                     str_detect(MonitoringLocationDescriptionText, ".rooked"), "Crooked River",
                                     str_detect(MonitoringLocationDescriptionText, "St. George Sound"), "St. George Sound",
                                     str_detect(station_nm, ".hattahoochee"), "Chattahoochee River",
                                     str_detect(station_nm, ".lint"), "Flint River",
                                     str_detect(station_nm, ".palachicola"), "Apalachicola River",
                                     str_detect(station_nm, ".hipola"), "Chipola River",
                                     str_detect(station_nm, ".arrabelle"), "Carrabelle River",
                                     str_detect(station_nm, "New R"), "New River",
                                     str_detect(station_nm, ".rooked"), "Crooked River",
                                     str_detect(station_nm, "St. George Sound"), "St. George Sound",
                                     str_detect(MonitoringLocationDescriptionText, ".HATTAHOOCHEE"), "Chattahoochee River",
                                     str_detect(MonitoringLocationDescriptionText, ".LINT"), "Flint River",
                                     str_detect(MonitoringLocationDescriptionText, ".PALACHICOLA"), "Apalachicola River",
                                     str_detect(MonitoringLocationDescriptionText, ".HIPOLA"), "Chipola River",
                                     str_detect(MonitoringLocationDescriptionText, ".ARRABELLE"), "Carrabelle River",
                                     str_detect(MonitoringLocationDescriptionText, "NEW R"), "New River",
                                     str_detect(MonitoringLocationDescriptionText, ".ROOKED"), "Crooked River",
                                     str_detect(MonitoringLocationDescriptionText, "ST. GEORGE SOUND"), "St. George Sound",
                                     str_detect(station_nm, ".HATTAHOOCHEE"), "Chattahoochee River",
                                     str_detect(station_nm, ".LINT"), "Flint River",
                                     str_detect(station_nm, ".PALACHICOLA"), "Apalachicola River",
                                     str_detect(station_nm, ".HIPOLA"), "Chipola River",
                                     str_detect(station_nm, ".ARRABELLE"), "Carrabelle River",
                                     str_detect(station_nm, "NEW R"), "New River",
                                     str_detect(station_nm, ".ROOKED"), "Crooked River",
                                     str_detect(station_nm, "ST. GEORGE SOUND"), "St. George Sound",
                                     default = "estuary"), by = row.names(apaw_stas2_acf2)]

apaw_use3b <- merge(apaw_use3, distinct(apaw_stas2_acf2[, .(umli, watershed)]), by = "umli", all.x = TRUE)
apaw_use4 <- apaw_use3b[umli %in% unique(apaw_stas2_acf$umli), ]

apaw_use4_sum <- distinct(apaw_use4[, .(wshed = unique(watershed),
                                        n = .N,
                                        begin = min(month),
                                        end = max(month),
                                        temp = ifelse(FALSE %in% is.na(temperature), "Y", "N"),
                                        sal = ifelse(FALSE %in% is.na(salinity), "Y", "N"),
                                        flow = ifelse(FALSE %in% is.na(flow), "Y", "N")),
                                    by = umli])
setorder(apaw_use4_sum, wshed, begin, n)

# sumdat <- distinct(apaw2[, .(n = .N, 
#                              begin = min(asdate), #ifelse(is.na(ActivityStartDateTime), min(as_date(ActivityStartDate)), min(as_date(ActivityStartDateTime))), 
#                              end = max(asdate)), #ifelse(is.na(ActivityStartDateTime), max(as_date(ActivityStartDate)), max(as_date(ActivityStartDateTime))), 
#                          # units = apaw2[parameter2 == parameter & MonitoringLocationIdentifier2 == MonitoringLocationIdentifier, unique(ResultMeasure.MeasureUnitCode)]),
#                          by = list(parameter, MonitoringLocationIdentifier)])

ggplot(apaw_use4[watershed == "Apalachicola River" & !is.na(flow), ], aes(x = month, y = flow, color = dec_lat_va)) +
  geom_point() +
  theme_bw()

ggplot(apaw_use4[!is.na(flow), ], aes(x = month, y = flow, color = dec_lat_va)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~watershed, ncol = 1, scales = "free_y")

ggplot(apaw_use4[!is.na(flow), ], aes(x = dec_lat_va, y = flow, group = dec_lat_va, color = watershed)) +
  geom_boxplot() +
  geom_point() +
  theme_bw() +
  coord_flip()
  # facet_wrap(~watershed, ncol = 1, scales = "free_y")


apaw_use3 <- merge(apaw_use3, distinct(apaw_carbon[, .(month, carbon)]), by = "month", all.x = TRUE, all.y = TRUE)
apaw_use3 <- merge(apaw_use3, distinct(apaw_flow[, .(month, flow)]), by = "month", all.x = TRUE, all.y = TRUE)

ggplot(apaw_use3[!is.na(flow) & !is.na(salinity), ], aes(x = flow, y = salinity, color = lon_rd)) +
  geom_point() +
  theme_bw() +
  scale_color_continuous_sequential(palette = "Inferno") +
  coord_cartesian(xlim = c(0, 5000)) +
  facet_wrap(~ lat_rd2)

ggplot(apaw_use3[!is.na(carbon) & !is.na(salinity), ], aes(x = carbon, y = salinity, color = lon_rd)) +
  geom_point() +
  theme_bw() +
  scale_color_continuous_sequential(palette = "Inferno") +
  coord_cartesian(xlim = c(0, 1)) +
  facet_wrap(~ lat_rd2)

ggplot(apaw_use3[!is.na(salinity), ], aes(x = date, y = salinity, color = WBID, shape = MonitoringLocationTypeName)) +
  geom_point() +
  geom_smooth(aes(x = date, y = salinity), method = "lm", inherit.aes = FALSE) +
  theme_bw() +
  scale_color_discrete_sequential(palette = "Inferno") +
  facet_wrap(~ lat_rd2)

# apaw_use3[!is.na(salinity), month := floor_date(date, unit = "month")]
# apaw_use3[!is.na(salinity), `:=` (minsal = min(salinity), maxsal = max(salinity)), by = list(month, lat_rd, lon_rd)]
apaw_use3[!is.na(salinity), `:=` (minsal = min(salinity), maxsal = max(salinity)), by = list(keepInd, month)]
apaw_use3_sal_sf <- st_as_sf(apaw_use3[!is.na(salinity), ], coords = c("dec_lon_va", "dec_lat_va"), crs = 4326)

shore_apa <- st_crop(shore, st_bbox(wbid_apa))
saltest <- subset(apaw_use3_sal_sf, year(apaw_use3_sal_sf$month) > 2005 & year(apaw_use3_sal_sf$month) < 2015)
saltest <- distinct(saltest[, c("month", "minsal", "maxsal", "geometry")])

ggplot(data = subset(saltest, year(saltest$month) == 2010)) +
  geom_sf(data = shore_apa) +
  geom_sf(aes(fill = minsal), color = "black", shape = 21) +
  theme_void() +
  scale_fill_continuous_sequential(palette = "Inferno") +
  facet_wrap(~ month, nrow = 4)

ggplot(data = saltest) +
  geom_sf(data = shore_apa, fill = "grey90", color = "grey70") +
  geom_sf(aes(fill = minsal), color = "black", shape = 21) +
  theme_void() +
  scale_fill_continuous_sequential(palette = "Inferno") +
  facet_wrap(~ month, ncol = 12)

# apaw_use3[!is.na(date), month := floor_date(date, unit = "month")]
# apaw_use3[!is.na(flow), `:=` (minflow = min(flow), maxflow = max(flow)), by = list(month, lat_rd, lon_rd)]
# apaw_use3[!is.na(carbon), `:=` (mincarb = min(carbon), maxcarb = max(carbon)), by = list(month, lat_rd, lon_rd)]
# apaw_use3[!is.na(temperature), `:=` (mintemp = min(temperature), maxtemp = max(temperature)), by = list(month, lat_rd, lon_rd)]

apaw_use3[!is.na(flow), `:=` (minsal = min(flow), maxsal = max(flow)), by = month]
apaw_use3[!is.na(carbon), `:=` (minsal = min(carbon), maxsal = max(carbon)), by = month]
apaw_use3[!is.na(temperature), `:=` (minsal = min(temperature), maxsal = max(temperature)), by = list(keepInd, month)]


# #Use the hobs_sta object to identify water quality data from within 500m of the HOBS localities
# apaCompare3_sf <- st_as_sf(apaCompare3, coords = c("lon", "lat"), crs = 4326)
# wq_buffs <- st_buffer(apaCompare3_sf, 500)
# apaw_use3_sf <- st_as_sf(apaw_use3[!is.na(lon_rd) & !is.na(lat_rd), ], coords = c("lon_rd", "lat_rd"), crs = 4326)
# # for(i in seq(1, nrow(apaw_use3))){
# #   apaw_use3[i, nearHOBS := ifelse(TRUE %in% st_intersects(wq_buffs, apaw_use3_sf$geometry[i], sparse = FALSE), 1, 0)]
# # }
# # saveRDS(apaw_use3, here::here(paste0("apaw_use3_", Sys.Date(), ".rds")))
# apaw_use3_nearHOBS <- readRDS(here::here("apaw_use3_2023-04-15.rds"))
# # st_geometry(apaw_use3_nearHOBS) <- NULL
# apaw_use3_sf <- st_as_sf(apaw_use3, coords = c("lon_rd", "lat_rd"), crs = 4326)
# apaw_use3_sf2 <- merge(apaw_use3_sf, distinct(apaw_use3_nearHOBS[, .(month, lat_rd2, lon_rd2, MonitoringLocationTypeName, nearHOBS)]), all.x = TRUE)
# 
# mapview(subset(apaw_use3_sf2, apaw_use3_sf2$nearHOBS == 1)) +
#   mapview(apaCompare3_sf, col.regions = "red")
# 
# # apaw_use4 <- apaw_use3_sf2
# # st_geometry(apaw_use4) <- NULL
# # setDT(apaw_use4)

apaCompare3[, date2 := fcase(date == "Feb-38", as_date("1938-02-01"),
                             date == "Jul-41", as_date("1941-07-01"),
                             date == "Aug-41", as_date("1941-08-01"),
                             date == "16-Dec", as_date("2016-12-01"),
                             date == "17-Apr", as_date("2017-04-01"),
                             date == "9/3/2020", as_date("2020-09-01"),
                             date == "8/11/2020", as_date("2020-08-01"),
                             date == "7/30/2019", as_date("2019-08-01"),
                             date == "1/31/1938", as_date("1938-02-01")), by = date]

# apaw_use4_sf <- st_as_sf(apaw_use4, coords = c("lon_rd2", "lat_rd2"), crs = 4326)
# apaw_use3_sf2$closest <- "NA" #https://stackoverflow.com/questions/49200458/find-nearest-features-using-sf-in-r
# for(pt in seq_len(nrow(apaw_use3_sf2))){
#   apaw_use3_sf2$closest[pt] <- apaCompare3_sf[which.min(
#     st_distance(apaCompare3_sf, apaw_use3_sf2[pt,])),]$sta2
# }
# apaw_use5 <- copy(apaw_use3_sf2)
# st_geometry(apaw_use5) <- NULL
# setDT(apaw_use5)
# apaw_use5[, closest2 := levels(apaCompare3$sta2)[as.numeric(closest)], by = closest]
# apaw_use5[, sta2 := closest2]

# for(i in seq(1, length(unique(apaCompare3$date2)))){
#   date_i <- unique(apaCompare3$date2)[i]
#   if(date_i %in% unique(apaw_use5$month)){
#     stas <- apaCompare3[date2 == date_i, unique(sta2)]
#     apaw_use5[month == date_i & closest2 %in% stas, use_exact := 1]
#     apaw_use5[month <= date_i & month > floor_date(date_i - 365, unit = "month") & closest2 %in% stas, `:=` (use_prioryr = 1,
#                                                                                                              deadC_mn = apaCompare3[date2 == date_i & sta2 == closest2, mean(deadC)],
#                                                                                                              deadC_sd = apaCompare3[date2 == date_i & sta2 == closest2, sd(deadC)],
#                                                                                                              # deadC = apaCompare3[date2 == date_i & sta2 == closest2, deadC],
#                                                                                                              flow_mn = apaw_use5[month <= date_i & month > floor_date(date_i - 365, unit = "month") & !is.na(flow) & sta2 == closest2, mean(flow)],
#                                                                                                              flow_sd = apaw_use5[month <= date_i & month > floor_date(date_i - 365, unit = "month") & !is.na(flow) & sta2 == closest2, sd(flow)],
#                                                                                                              carbon_mn = apaw_use5[month <= date_i & month > floor_date(date_i - 365, unit = "month") & !is.na(carbon) & sta2 == closest2, mean(carbon)],
#                                                                                                              carbon_sd = apaw_use5[month <= date_i & month > floor_date(date_i - 365, unit = "month") & !is.na(carbon) & sta2 == closest2, sd(carbon)],
#                                                                                                              sal_mn = apaw_use5[month <= date_i & month > floor_date(date_i - 365, unit = "month") & !is.na(salinity) & sta2 == closest2, mean(salinity)],
#                                                                                                              sal_sd = apaw_use5[month <= date_i & month > floor_date(date_i - 365, unit = "month") & !is.na(salinity) & sta2 == closest2, sd(salinity)],
#                                                                                                              minsal_mn = apaw_use5[month <= date_i & month > floor_date(date_i - 365, unit = "month") & !is.na(minsal) & sta2 == closest2, mean(salinity)],
#                                                                                                              minsal_sd = apaw_use5[month <= date_i & month > floor_date(date_i - 365, unit = "month") & !is.na(minsal) & sta2 == closest2, sd(salinity)],
#                                                                                                              temp_mn = apaw_use5[month <= date_i & month > floor_date(date_i - 365, unit = "month") & !is.na(temperature) & sta2 == closest2, mean(temperature)],
#                                                                                                              temp_sd = apaw_use5[month <= date_i & month > floor_date(date_i - 365, unit = "month") & !is.na(temperature) & sta2 == closest2, sd(temperature)]),
#               by = closest2]
#   }
# }

apaw_use3[, `:=` (area2 = area, year = year(month))]
apaCompare3[, yr := year(date2)]
apaCompare3[sta2 == "East Apalachicola Bay", sta2 := ifelse(sta == "Cat Point", "East Apalachicola Bay 2", "East Apalachicola Bay 1")]
apaCompare3[, sta2 := factor(sta2, levels = c("West Apalachicola Bay", "Mid Apalachicola Bay", "East Apalachicola Bay 2", "East Apalachicola Bay 1", "Dog Island", "Alligator Harbor"), ordered = TRUE)]

for(i in seq(1, length(unique(apaCompare3$date2)))){
  date_i <- unique(apaCompare3$date2)[i]
  if(year(date_i) %in% unique(apaw_use3$year)){
    stas <- apaCompare3[date2 == date_i, unique(sta2)]
    apaw_use3[year == year(date_i) & area %in% stas, use_exact := 1]
    apaw_use3[year == year(date_i) & area %in% stas, `:=` (use_prioryr = 1,
                                                           deadC_mn = apaCompare3[yr == year(date_i) & sta2 == area, mean(deadC)],
                                                           deadC_sd = apaCompare3[yr == year(date_i) & sta2 == area, sd(deadC)],
                                                           # deadC = apaCompare3[date2 == date_i & sta2 == closest2, deadC],
                                                           flow_mn = apaw_use3[year <= year(date_i) & year >= year(date_i) - 1 & !is.na(flow) & area2 == area, mean(flow)],
                                                           flow_sd = apaw_use3[year <= year(date_i) & year >= year(date_i) - 1 & !is.na(flow) & area2 == area, sd(flow)],
                                                           carbon_mn = apaw_use3[year <= year(date_i) & year >= year(date_i) - 1 & !is.na(carbon) & area2 == area, mean(carbon)],
                                                           carbon_sd = apaw_use3[year <= year(date_i) & year >= year(date_i) - 1 & !is.na(carbon) & area2 == area, sd(carbon)],
                                                           sal_mn = apaw_use3[year <= year(date_i) & year >= year(date_i) - 1 & !is.na(salinity) & area2 == area, mean(salinity)],
                                                           sal_sd = apaw_use3[year <= year(date_i) & year >= year(date_i) - 1 & !is.na(salinity) & area2 == area, sd(salinity)],
                                                           minsal_mn = apaw_use3[year <= year(date_i) & year >= year(date_i) - 1 & !is.na(minsal) & area2 == area, mean(salinity)],
                                                           minsal_sd = apaw_use3[year <= year(date_i) & year >= year(date_i) - 1 & !is.na(minsal) & area2 == area, sd(salinity)],
                                                           temp_mn = apaw_use3[year <= year(date_i) & year >= year(date_i) - 1 & !is.na(temperature) & area2 == area, mean(temperature)],
                                                           temp_sd = apaw_use3[year <= year(date_i) & year >= year(date_i) - 1 & !is.na(temperature) & area2 == area, sd(temperature)]),
              by = area]
  }
}


# ggplot(apaw_use5[use_exact == 1 & !is.na(flow), ], aes(x = flow_mn, y = deadC_mn, color = month)) +
#   geom_point() +
#   theme_bw()
# 
# ggplot(apaw_use5[use_exact == 1 & !is.na(salinity), ], aes(x = sal_mn, y = deadC_mn, color = month)) +
#   geom_point() +
#   theme_bw()
# 
# ggplot(distinct(apaw_use5[use_exact == 1 & !is.na(salinity), .(month, sal_mn, deadC_mn, closest2)]), aes(x = sal_mn, y = deadC_mn, fill = closest2)) +
#   geom_jitter(width = 0.1, height = 0.0001, shape = 21, size = 2.5, color = "black") +
#   theme_bw()
# 
# ggplot(distinct(apaw_use5[use_exact == 1 & !is.na(salinity), .(month, sal_mn, minsal_mn, deadC_mn, closest2)]), aes(x = minsal_mn, y = deadC_mn)) + #, fill = closest2
#   geom_point(shape = 21, size = 2.5, color = "black") +
#   geom_smooth(aes(x = minsal_mn, y = deadC_mn), method = "lm", inherit.aes = FALSE) +
#   # geom_point(aes(x = sal_mn, y = deadC_mn, fill = closest2), shape = 22, size = 2.5, color = "black", inherit.aes = FALSE) +
#   stat_regline_equation(label.y = 0.004, label.x = 18, aes(label = ..eq.label..)) +
#   stat_regline_equation(label.y = 0.0025, label.x = 18, aes(label = ..rr.label..)) +
#   theme_bw()
# 
# ggplot(data = distinct(apaCompare3_sf[, c("sta2", "deadC", "year", "study")]), aes(x = sta2, y = deadC, fill = as.factor(year), shape = study)) +
#        geom_point(size = 2, color = "black") +
#        geom_label(aes(label = year), position = position_jitter(width = 0.55, height = 0.007), size = 2) +
#        ylim(-0.02, 0.08) +
#        theme_bw() +
#        labs(x = "Locality", fill = "Year", shape = "Study") +
#        scale_fill_discrete_sequential(palette = "inferno", nmax = 10, order = c(2:7), breaks = c("1938", "1941", "2016", "2017", "2019", "2020")) +
#        guides(fill = guide_legend(override.aes = list(shape=21))) +
#        scale_shape_manual(values = c("Durham et al. 2023" = 21, "Hadden & Cherkinsky 2017" = 22, "Hadden et al. 2018" = 23)) +
#        theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(distinct(apaw_use3[use_exact == 1 & !is.na(salinity), .(month, sal_mn, deadC_mn, area)]), aes(x = sal_mn, y = deadC_mn, fill = area)) +
  geom_jitter(width = 0.1, height = 0.0001, shape = 21, size = 2.5, color = "black") +
  theme_bw()

ggplot(distinct(apaw_use3[use_prioryr == 1 & !is.na(salinity), .(year, sal_mn, deadC_mn, area)]), aes(x = sal_mn, y = deadC_mn, fill = area)) +
  geom_jitter(width = 0.1, height = 0.0001, shape = 21, size = 2.5, color = "black") +
  theme_bw()

ggplot(distinct(apaw_use3[use_prioryr == 1 & !is.na(salinity), .(year, sal_mn, minsal_mn, deadC_mn, area)]), aes(x = minsal_mn, y = deadC_mn)) + #, fill = closest2
  geom_point(shape = 21, size = 2.5, color = "black") +
  geom_smooth(aes(x = minsal_mn, y = deadC_mn), method = "lm", inherit.aes = FALSE) +
  # geom_point(aes(x = sal_mn, y = deadC_mn, fill = closest2), shape = 22, size = 2.5, color = "black", inherit.aes = FALSE) +
  stat_regline_equation(label.y = 0.004, label.x = 18, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 0.0025, label.x = 18, aes(label = ..rr.label..)) +
  theme_bw()

# model salinity v. deadC-----------------------------------------
salmod <- brm(formula = deadC_mn ~ 1 + me(minsal_mn, minsal_sd) + (1 | year:area), family = gaussian, data = distinct(apaw_use3[use_prioryr == 1 & !is.na(salinity), .(year, sal_mn, sal_sd, minsal_mn, minsal_sd, deadC_mn, area)]),
              cores = 4, chains = 4, iter = 5000, control = c(adapt_delta = 0.99, max_treedepth = 20), backend = "cmdstanr", threads = threading(2), seed = 453, file = here::here("salmod.rds"))
summary(salmod)

salmod2 <- brm(formula = deadC_mn ~ 1 + minsal_mn + (1 | year:area), family = gaussian, data = distinct(apaw_use3[use_prioryr == 1 & !is.na(salinity), .(year, sal_mn, sal_sd, minsal_mn, minsal_sd, deadC_mn, area)]),
              cores = 4, chains = 4, iter = 5000, control = c(adapt_delta = 0.99, max_treedepth = 20), backend = "cmdstanr", threads = threading(2), seed = 453, file = here::here("salmod2.rds"))
summary(salmod2)

# plausibility check
testdat <- data.table(minsal_mn = seq(0, 30), year = 2016, area = "East Apalachicola Bay 2")
testdat2 <- cbind(testdat, predict(salmod2, testdat)) #sal when deadC is ~0.078 is ~5

# model longitude v. deadC-----------------------------------------

ggplot(data = distinct(apaCompare3[, c("sta2", "deadC", "year", "study")]), aes(x = sta2, y = deadC, fill = as.factor(year), shape = study)) +
  geom_point(size = 2, color = "black") +
  geom_label(aes(label = year), position = position_jitter(width = 0.55, height = 0.007), size = 2) +
  ylim(-0.02, 0.08) +
  theme_bw() +
  labs(x = "Locality", fill = "Year", shape = "Study") +
  scale_fill_discrete_sequential(palette = "inferno", nmax = 10, order = c(2:7), breaks = c("1938", "1941", "2016", "2017", "2019", "2020")) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  scale_shape_manual(values = c("Durham et al. 2023" = 21, "Hadden & Cherkinsky 2017" = 22, "Hadden et al. 2018" = 23)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = distinct(apaCompare3[, c("lon", "deadC", "year", "study")]), aes(x = lon, y = deadC, fill = as.factor(year), shape = study)) +
  geom_point(size = 2, color = "black") +
  # geom_label(aes(label = year), position = position_jitter(width = 0.55, height = 0.007), size = 2) +
  ylim(-0.02, 0.08) +
  theme_bw() +
  labs(x = "Longitude", fill = "Year", shape = "Study") +
  scale_fill_discrete_sequential(palette = "inferno", nmax = 10, order = c(2:7), breaks = c("1938", "1941", "2016", "2017", "2019", "2020")) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  scale_shape_manual(values = c("Durham et al. 2023" = 21, "Hadden & Cherkinsky 2017" = 22, "Hadden et al. 2018" = 23)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

apaCompare3[, postbomb := ifelse(decade %in% c("2010s", "2020s"), 1, 0)]
apaCompare3[, postbomb := as.logical(postbomb)]
apaCompare3[, lon2 := fcase(sta2 == "Alligator Harbor", -84.41000,
                            sta2 == "Dog Island", -84.61000,
                            sta2 == "East Apalachicola Bay 1", -84.77293,
                            sta2 == "East Apalachicola Bay 2", -84.89000,
                            sta2 == "Mid Apalachicola Bay", -85.00554,
                            sta2 == "West Apalachicola Bay", -85.22000)]

lonmod <- brm(formula = deadC ~ 1 + lon + (1 | postbomb), family = gaussian, data = apaCompare3,
              cores = 4, chains = 4, warmup = 1000, iter = 5000, thin = 2, control = c(adapt_delta = 0.99, max_treedepth = 18), backend = "cmdstanr", threads = threading(2), seed = 453, file = here::here("lonmod.rds"))
summary(lonmod)

lonmod2 <- brm(formula = deadC ~ 1 + lon2 + (1 | postbomb), family = gaussian, data = apaCompare3,
              cores = 4, chains = 4, warmup = 1000, iter = 5000, thin = 2, control = c(adapt_delta = 0.99, max_treedepth = 18), backend = "cmdstanr", threads = threading(2), seed = 453, file = here::here("lonmod2.rds"))
summary(lonmod2)

lonmod3 <- brm(formula = deadC ~ 1 + lon:postbomb, family = gaussian, data = apaCompare3,
               cores = 4, chains = 4, warmup = 1000, iter = 5000, thin = 2, control = c(adapt_delta = 0.99, max_treedepth = 18), backend = "cmdstanr", threads = threading(2), seed = 453, file = here::here("lonmod3.rds"))
summary(lonmod3)

surftest <- distinct(apaw_use3[!is.na(dec_lon_va) & !is.na(minsal_mn) & !is.na(deadC_mn) & !is.na(year), .(dec_lon_va, minsal_mn, deadC_mn, year)])

deadcmod <- brm(formula = deadC_mn ~ 1 + dec_lon_va * minsal_mn + (1 | year), family = gaussian, data = surftest,
                cores = 4, chains = 4, warmup = 1000, iter = 5000, thin = 2, control = c(adapt_delta = 0.99, max_treedepth = 18), 
                backend = "cmdstanr", threads = threading(2), seed = 453, file = here::here("deadcmod.rds"))

surftest[, relyear := year - min(surftest$year)]

surfdat <- data.table(dec_lon_va = rep(seq(-85.222, -84.401, 0.05), each = 16), minsal_mn = rep(seq(0, 30, 2), 17), year = 2016)
surfdat2 <- cbind(surfdat, predict(deadcmod, surfdat))
Dead_C <- matrix(surfdat2$Estimate, nrow = 17, ncol = 16, byrow = TRUE,
                  dimnames = list(c(unique(surfdat2$dec_lon_va)),
                                  c(unique(surfdat2$minsal_mn))))

surfplot <- plot_ly(x = unique(surfdat2$minsal_mn), 
                    y = unique(surfdat2$dec_lon_va), 
                    z = ~Dead_C,
                    contours = list(z = list(show = TRUE,
                                             start = round(min(Dead_C), 2),
                                             end = round(max(Dead_C), 2),
                                             size = 0.01,
                                             usecolormap=TRUE,
                                             highlightcolor="#ff0000",
                                             project = list(z = TRUE)))) %>% 
  add_surface(contours = list(z = list(show = TRUE,
                                       start = round(min(Dead_C), 2),
                                       end = round(max(Dead_C), 2),
                                       size = 0.01,
                                       usecolormap=TRUE,
                                       highlightcolor="#ff0000",
                                       project = list(z = TRUE)))) %>% 
  layout(scene = list(xaxis = list(title = "Salinity"),
                      yaxis = list(title = "Longitude"),
                      zaxis = list(title = "Dead C"))) #only works if "~" is removed from z argument in plot_ly()


surfdat3 <- data.table(dec_lon_va = rep(seq(-85.222, -84.401, 0.05), each = 16*9), minsal_mn = rep(seq(0, 30, 2), 17*9), year = rep(seq(1940, 2020, 10), each = 16))
surfdat4 <- cbind(surfdat3, predict(deadcmod, surfdat3, allow_new_levels = TRUE))

yr <- 1940
Dead_C <- matrix(surfdat4[year == yr, Estimate], nrow = 17, ncol = 16, byrow = TRUE,
                 dimnames = list(c(surfdat4[year == yr, unique(dec_lon_va)]),
                                 c(surfdat4[year == yr, unique(minsal_mn)])))

surfplot2 <- plot_ly(x = surfdat4[year == yr, unique(minsal_mn)], 
                     y = surfdat4[year == yr, unique(dec_lon_va)], 
                     z = ~Dead_C,
                     contours = list(z = list(show = TRUE,
                                              start = round(min(Dead_C), 2),
                                              end = round(max(Dead_C), 2),
                                              size = 0.01,
                                              usecolormap=TRUE,
                                              highlightcolor="#ff0000",
                                              project = list(z = TRUE)))) %>% 
  add_surface(contours = list(z = list(show = TRUE,
                                       start = round(min(Dead_C), 2),
                                       end = round(max(Dead_C), 2),
                                       size = 0.01,
                                       usecolormap=TRUE,
                                       highlightcolor="#ff0000",
                                       project = list(z = TRUE)))) %>% 
  layout(scene = list(title = yr,
                      xaxis = list(title = "Salinity"),
                      yaxis = list(title = "Longitude"),
                      zaxis = list(title = "Dead C"))) #only works if "~" is removed from z argument in plot_ly()


surfplot3 <- surfplot %>% add_surface(x = surfdat4[year == yr, unique(minsal_mn)], 
                                      y = surfdat4[year == yr, unique(dec_lon_va)], 
                                      z = Dead_C + 0.3)

