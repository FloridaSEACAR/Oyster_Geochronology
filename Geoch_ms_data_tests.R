
library(data.table)
library(tidyverse)
#install.packages("ggThemeAssist")

# sampsum <- data.table::fread(here::here('sampleSummary_7-10-2020.csv'))
# sampsum2  <- data.table::fread(here::here('sampleSummary_7-28-2020.csv'))
# specdat <- data.table::fread(here::here('HOBS_PhaseI_Posteriors.csv'))
# specsum <- data.table::fread(here::here('specimenSummary_7-10-2020.csv'))
specdat <- fread(here::here("HOBS_Posteriors_all.csv"))
sampsum <- fread(here::here("sampleSummary_2021-10-20.csv"))
sampsum[, depth := str_sub(RandL, start = -8, end = -1)]
specsum <- fread(here::here("specimenSummary_2021-10-20.csv"))
MKetal_specdat <- data.table::fread(here::here('AppendixDR2_PosteriorDistributions_dataonly.csv'))
MKetal_specsum <- data.table::fread(here::here('AppendixDR1_dataSummaryTable_dataonly.csv'))

names(specsum)<-str_replace_all(names(specsum), c(" " = "_"))
specsum$median_age_bp <- abs(specsum$median_age - 2019)
specsum$loc <- substring(specsum$name, 1, 2)
locs_abbrev <- unique(specsum$loc)
locs_abbrev_corr <- c("BH", "GI-EC", "GR", "HC-MC", "JI-WC", "LB", "LC", "LSG", "MR", "NP", "PC")
locs_full <- c("Big Hickory Island", "Goose Island & East Cove", "Guana River", "Hendry & Mullock Creeks", 
               "Jack Island", "Lemon Bay", "Lone Cabbage Island", "Little St. George Island", "Matanzas River",
               "New Pass", "Pellicer Creek")
specsum$loc_full <- specsum$loc
for(i in seq_along(locs_abbrev)){
  specsum <- specsum %>% dplyr::mutate_at(c("loc_full"), stringr::str_replace_all, pattern = locs_abbrev[i], replacement = locs_full[i])
}
for(i in seq_along(locs_abbrev)){
  specsum <- specsum %>% dplyr::mutate_at(c("loc"), stringr::str_replace_all, pattern = locs_abbrev[i], replacement = locs_abbrev_corr[i])
}
specsum$reefID <- paste0(specsum$loc, " R", specsum$reef)


specdat$value_bp <- abs(specdat$value - 2019)
specdat$loc <- substring(specdat$name, 1, 2)
specdat$loc_full <- specdat$loc
for(i in seq_along(locs_abbrev)){
  specdat <- specdat %>% dplyr::mutate_at(c("loc_full"), stringr::str_replace_all, pattern = locs_abbrev[i], replacement = locs_full[i])
}
for(i in seq_along(locs_abbrev)){
  specdat <- specdat %>% dplyr::mutate_at(c("loc"), stringr::str_replace_all, pattern = locs_abbrev[i], replacement = locs_abbrev_corr[i])
}
specdat$reefID <- paste0(specdat$loc, " R", specdat$reef)

## From script written by Matt Kosnik for our pilot LP14C data
## clear out the really 0 probability ages
specdat <- specdat[(specdat$probability > 0),]
MKetal_specdat <- MKetal_specdat[(MKetal_specdat$Probability > 0),]

## check to see what the probabilities sum to 1 (they do not)
( cSum <- aggregate(probability ~ name, data=specdat, FUN=sum) )
( cSum2 <- aggregate(Probability ~ Specimen, data=MKetal_specdat, FUN=sum) )

## make the probabilities for each specimen sum to 1
## all these are half year probabilities so can all get uniformly scaled
specdat2 <- merge(specdat,cSum, by='name')
specdat2$aProb <- specdat2$probability.x/specdat2$probability.y

MKetal_specdat2 <- merge(MKetal_specdat, cSum2, by='Specimen')
MKetal_specdat2$aProb <- MKetal_specdat2$Probability.x/MKetal_specdat2$Probability.y

## check to see what the probabilies sum to 1
( cSum <- aggregate(aProb ~ name, data=specdat2, FUN=sum) )
( cSum2 <- aggregate(aProb ~ Specimen, data=MKetal_specdat2, FUN=sum) )


# sampsum2$locality <- substring(sampsum2$sample, 1, 2)
# sampsum$reef <- substring(sampsum$sample, 4, 5)
# sampsum$depth <- substring(sampsum$sample, 7, 13)
sampsum$median_age_bp <- abs(sampsum$Sample_median_age - 2019)
# names(sampsum)<-str_replace_all(names(sampsum), c(" " = "_"))

vlines <- data.frame(xint = c(median(subset(sampsum, sampsum$depth == "15-25cm")$median_age_bp), 
                              median(subset(sampsum, sampsum$depth == "25-35cm")$median_age_bp)), 
                    depth = c('15-25cm', '25-35cm'))

# vlines2 <- data.frame(xint = c(median(subset(sampsum2, sampsum2$depth == "15-25cm")$median_age_bp), 
#                               median(subset(sampsum2, sampsum2$depth == "25-35cm")$median_age_bp)), 
#                      depth = c('15-25cm', '25-35cm'))


sampsum_nolc <- subset(sampsum, sampsum$locality != "LC")
vlines_nolc <- data.frame(xint = c(median(subset(sampsum_nolc, sampsum_nolc$depth == "15-25cm")$median_age_bp), 
                                   median(subset(sampsum_nolc, sampsum_nolc$depth == "25-35cm")$median_age_bp)), 
                          depth = c('15-25cm', '25-35cm'))


#Time-averaging by median sample age (RandL version)

Ketal <- data.frame(median_age_bp = 1058, corrected_posterior_age_estimate = 1830)
label1 <- c(expression(atop(NA, atop(textstyle('Kowalewski et al. 2018'), 
                                    textstyle(paste('(', italic('Tucetona pectinata'), ')'))))))

Detal <- data.frame(median_age_bp = c(123, 283, 146, 148, 160, 2043), 
                    corrected_posterior_age_estimate = c(68, 2670, 16, 27, 16, 1937))
label2 <- c(expression(atop(NA, atop(textstyle('Dominguez et al. 2016'), 
                                     textstyle(paste('(', italic('Fulvia tenuicostata'), ')'))))))

ggplot(sampsum, aes(x = median_age_bp, y = corrected_posterior_age_estimate)) +
  theme_classic() +
  scale_fill_brewer(palette = 'Dark2') +
  scale_color_brewer(palette = 'Dark2') +
  geom_hline(yintercept = c(1, 10, 100, 1000)) +
  geom_vline(data = vlines, aes(xintercept = xint, color = depth), show.legend = FALSE) +
  geom_point(aes(fill = depth), color = "black", shape = 21, size = 3) +
  geom_point(data = Ketal, fill = 'dark grey', color = 'black', shape = 21, size = 3) +
  annotate('text', x = Ketal$median_age_bp[1] - 250, y = Ketal$corrected_posterior_age_estimate[1] + 600, 
           label = label1) +
  geom_point(data = Detal, fill = 'dark grey', color = 'black', shape = 22, size = 3) +
  annotate('text', x = Detal$median_age_bp[6] - 250, y = Detal$corrected_posterior_age_estimate[6] + 600, 
           label = label2) +
  scale_y_continuous(trans='log10') +
  scale_x_reverse() +
  labs(x = 'Median age (yr)',
       y = 'Corrected posterior age estimate \n (scale of time-averaging)',
       fill = 'Burial depth') + 
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12, colour = "black"), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))
  

#Time-averaging by median sample age (True sample version)

Ketal <- data.frame(median_age_bp = 1058, Sample_corrected_posterior_age_estimate = 1830)
label1 <- c(expression(atop(NA, atop(textstyle('Kowalewski et al. 2018'), 
                                     textstyle(paste('(', italic('Tucetona pectinata'), ')'))))))

Detal <- data.frame(median_age_bp = c(123, 283, 146, 148, 160, 2043), 
                    Sample_corrected_posterior_age_estimate = c(68, 2670, 16, 27, 16, 1937))
label2 <- c(expression(atop(NA, atop(textstyle('Dominguez et al. 2016'), 
                                     textstyle(paste('(', italic('Fulvia tenuicostata'), ')'))))))

ggplot(sampsum2, aes(x = median_age_bp, y = Sample_corrected_posterior_age_estimate)) +
  theme_classic() +
  scale_fill_brewer(palette = 'Dark2') +
  scale_color_brewer(palette = 'Dark2') +
  geom_hline(yintercept = c(1, 10, 100, 1000)) +
  geom_vline(data = vlines2, aes(xintercept = xint, color = depth), show.legend = FALSE) +
  geom_point(aes(fill = depth), color = "black", shape = 21, size = 3) +
  geom_point(data = Ketal, fill = 'dark grey', color = 'black', shape = 21, size = 3) +
  annotate('text', x = Ketal$median_age_bp[1] - 250, y = Ketal$Sample_corrected_posterior_age_estimate[1] + 600, 
           label = label1) +
  geom_point(data = Detal, fill = 'dark grey', color = 'black', shape = 22, size = 3) +
  annotate('text', x = Detal$median_age_bp[6] - 250, y = Detal$Sample_corrected_posterior_age_estimate[6] + 600, 
           label = label2) +
  scale_y_continuous(trans='log10') +
  scale_x_reverse() +
  labs(x = 'Median age (yr)',
       y = 'Corrected posterior age estimate \n (scale of time-averaging)',
       fill = 'Burial depth') + 
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12, colour = "black"), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))

#Boxplots of median sample ages and time-averaging
sampsum2$reefID <- paste0(sampsum2$locality, "_R", sampsum2$reef)
sampsum2$holeID <- stringr::str_extract(sampsum2$sample, ".._R.H.")
sampsum2$holeN <- stringr::str_extract(sampsum2$holeID, "R.H.")
sampsum2$holeN <- stringr::str_extract(sampsum2$holeN, "H.")

ggplot(sampsum2, aes(x = depth, y = median_age_bp, color = holeN)) +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2') +
  geom_errorbar(aes(ymin = median_age_bp - Sample_corrected_posterior_age_estimate/2, 
                    ymax = median_age_bp + Sample_corrected_posterior_age_estimate/2),
                position = position_dodge(width = 0.5), width = 0.3) +
  geom_point(position = position_dodge(width = 0.5)) +
  ylim(0, 60) +
  facet_wrap(~reefID)


#Posterior distributions by Locality
L1 <- subset(specdat2, specdat2$depth == "15-25cm")
L2 <- subset(specdat2, specdat2$depth == "25-35cm")
colnames(MKetal_specdat2)[1] <- "name"
colnames(MKetal_specdat2)[3] <- "value_bp"
MKetal_specdat2$loc <- "Kowalewski et al. 2018"
MKetal_specdat2_tuc <- subset(MKetal_specdat2, MKetal_specdat2$Taxon == "Tucetona")
colnames(MKetal_specsum)[3] <- "name"
MKetal_specsum$loc <- "Kowalewski et al. 2018"
MKetal_specsum_tuc <- subset(MKetal_specsum, MKetal_specsum$taxon == "Tucetona")

ggplot(NULL, aes(x = value_bp, y = aProb, group = name)) +
  geom_hline(yintercept = 0, color = 'grey', size = 0.5) +
  geom_area(data = L1, aes(y = aProb, fill = depth), position = 'identity', alpha = 0.45, color = 'black', size = 0.2) + 
  geom_area(data = L2, aes(y = aProb, fill = depth), position = 'identity', alpha = 0.45, color = 'black', size = 0.2) +
  geom_area(data = MKetal_specdat2_tuc, aes(y = aProb), fill = 'dark grey', position = 'identity', alpha = 0.45, color = 'black', size = 0.2) +
  geom_point(data = specsum, aes(x = median_age_bp, y = -0.05, fill = depth), color = 'black', stroke = 0.3, size = 1.5, shape = 24) +
  geom_point(data = MKetal_specsum_tuc, aes(x = medianAge, y = -0.05), fill = 'dark grey', color = 'black', stroke = 0.3, size = 1.5, shape = 24) +
  scale_fill_brewer(palette = 'Dark2') +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw() +
  facet_wrap(~factor(loc, levels = c('LSG', 'GI-EC', 'LC', 'LB', 'HC-MC', 'NP', 'BH', 'JI-WC', 'PC', 'MR', 'GR', 'Kowalewski et al. 2018'),
                     ordered = TRUE)) +
  scale_x_reverse() +
  xlim(60, 0) +
  labs(x = 'Age (yr)',
       y = 'Annual probability',
       fill = 'Burial depth') + 
  theme(panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank")) 
#        axis.title = element_text(face = "bold"), axis.text = element_text(size = 10, colour = "black"))
  
#Posterior distributions by reefID
ggplot(NULL, aes(x = value_bp, y = aProb, group = name)) +
  geom_hline(yintercept = 0, color = 'grey', size = 0.5) +
  geom_area(data = L1, aes(y = aProb, fill = depth), position = 'identity', alpha = 0.45, color = "black", size = 0.2) + 
  geom_area(data = L2, aes(y = aProb, fill = depth), position = 'identity', alpha = 0.45, color = "black", size = 0.2) +
  geom_point(data = specsum, aes(x = median_age_bp, y = -0.05, fill = depth), color = 'black', stroke = 0.3, size = 1.5, shape = 24) +
  scale_fill_brewer(palette = 'Dark2') +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw() +
  facet_wrap(~factor(reefID, levels = c('LSG R1', 'LSG R2', 'LSG R3', 
                                        'GI-EC R1', 'GI-EC R2', 'GI-EC R4',
                                        'LC R1', 'LC R2', 'LC R3', 'LC R5', 
                                        'LB R1', 'LB R2', 'LB R3',
                                        'HC-MC R1', 'HC-MC R2',
                                        'NP R2', 'NP R3',
                                        'BH R1', 'BH R2', 'BH R3',
                                        'JI-WC R1', 'JI-WC R2',
                                        'PC R1', 'PC R2', 'PC R3',
                                        'MR R2', 'MR R3',
                                        'GR R1', 'GR R2', 'GR R3'),
                     ordered = TRUE)) +
  #facet_grid(rows = vars(loc)) +
  scale_x_reverse() +
  xlim(60, 0) +
  labs(x = 'Age (yr)',
       y = 'Annual probability',
       fill = 'Burial depth') + 
  theme(panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"))



