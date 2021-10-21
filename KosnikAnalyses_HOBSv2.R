## read data
pData <- read.csv(here::here('HOBS_PhaseI_Posteriors.csv'), as.is=TRUE)
sData <- read.csv(here::here('HOBS_PhaseI_SpecInfo.csv'), as.is=TRUE)

#load June 2020 updated museum specimen posteriors
pData <- read.csv(here::here("MuseumSpecPosteriors_updated062021.csv"), as.is = TRUE)
pData$name <- str_replace(pData$name, "UF ", "UF")

#Pilot/test data

pData <- read.csv(here::here('PILOT_Posteriors-Table 1.csv'), as.is=TRUE)
sData <- read.csv(here::here('PILOT_HOBS_14CPilot_SpecInfo.csv'), as.is=TRUE)

#Final HOBS Phase I and II data
library(tidyverse)
library(data.table)
pData <- fread(here::here("HOBS_Posteriors_all.csv"))
sData <- fread(here::here("HOBS_Geochronology_DAspecimensVWW.csv"))
sData[, Spec_ID := str_replace(Spec_ID, "S1_", "S1-")]
sData[, Spec_ID := str_replace(Spec_ID, "S2_", "S2-")]

## COLUMN / PAGE SIZES FOR FIGURES
pageWidthOne <- 2.33
pageWidthTwo <- 3.5
pageWidthFull <- 7.125
pageHeight <- 9.5

summary(pData)
summary(sData)

library(stringr)
library(Hmisc)

## make matchable name between sData and pData
#sData$name <- paste(sData$SpecimenID, sData$Sample_depth, sep='  ')
sData$name <- paste(sData$Spec_ID)

setdiff(unique(sData$name),unique(pData$name))
setdiff(unique(pData$name),unique(sData$name))

nrow(sData)
nrow(pData)

## clear out the really 0 probability ages
pData <- pData[(pData$probability > 0),]

## check to see what the probabilities sum to 1 (they do not)
( cSum <- aggregate(probability ~ name, data=pData, FUN=sum) )

## make the probabilities for each specimen sum to 1
## all these are half year probabilities so can all get uniformly scaled
pData2 <- merge(pData, cSum, by='name')
pData2$aProb <- pData2$probability.x/pData2$probability.y
pData2 <- pData2[,c('name','value','aProb','depth','reef')]
pData2 <- pData2[,c('name','value','aProb')] #version to run when analyzing the museum specimens

## check to see what the probabilies sum to 1
( cSum <- aggregate(aProb ~ name, data=pData2, FUN=sum) )

## make site, excavation depth and sample variables
pData2$site <- substring(pData2$name,0,2)
#pData2$depth <- substring(pData2$name,(nchar(pData2$name)-6),nchar(pData2$name))
pData2$RandL <- paste0(pData2$site, " R", pData2$reef, " ", pData2$depth)
pData2$sample <- paste0(pData2$site, "_", stringr::str_extract(pData2$name, "R.H.S."))

sData$site <- substring(sData$name,0,2)
sData$RandL <- paste0(sData$site, " R", sData$Reef_num, " ", sData$Sample_depth)
sData$sample <- paste0(sData$site, "_", stringr::str_extract(sData$name, "R.H.S."))

## list of unique specimen names
SPECIMENS <- unique(pData2$name)



## plot specimen probabilities
pdf(here::here('plot-Specimens.pdf'), width=pageWidthFull, height=pageHeight/2)
plot(aProb ~ value, data= pData2, type='n', xlab='years (AD)', ylab='posterior probability')
mtext('Specimen age probability', line=-1)
for (s in 1:length(SPECIMENS)) {
	sad <- pData2[(pData2$name == SPECIMENS[s]),]
	polygon(c(min(sad$value),sad$value,max(sad$value)), c(0,sad$aProb,0), col= rgb(0,0,0,0.1))
}

dev.off()

head(pData2)

## make sample level age posterior probabilities
sAgeDist <- aggregate(aProb ~ sample + value, data=pData2, FUN=sum)
SAMPLES <- unique(sAgeDist$sample)
nsamp <- sData %>% dplyr::count(sample)
sAgeDist_nsamp <- dplyr::left_join(sAgeDist, nsamp, by = "sample")
## adjust for the fact that there are x specimens per sample.
sAgeDist_nsamp$aProb <- sAgeDist_nsamp$aProb/sAgeDist_nsamp$n



## PLOT

# colors to use...
pCol <- c(rgb(0,0.1,0.4,0.3), rgb(0.4,0.1,0,0.3), rgb(0,0.4,0.1,0.3), rgb(0.4,0,0.1,0.3))

## plot sample age probabilities
pdf(here::here('plot-Samples.pdf'), width=pageWidthFull, height=pageHeight/2)
plot(aProb ~ value, data= sAgeDist_nsamp, type='n', xlab='years (AD)', ylab='posterior probability')
mtext('Sample age probability', line=-1)
for (s in 1:length(SAMPLES)) {
	sad <- sAgeDist_nsamp[(sAgeDist_nsamp$sample == SAMPLES[s]),]
	polygon(c(min(sad$value),sad$value,max(sad$value)), c(0,sad$aProb,0), col= pCol[s])
	peakProb <- sad[(sad$aProb == max(sad$aProb)),]
	
	text(peakProb$value+1, peakProb$aProb+0.001, peakProb$sample)
}
dev.off()


## SPECIMEN SUMMARY TABLE
cols <- c('name', 'median age', 'age range', 'minimum age', 'maximum age', 'specimen age variability (IQR)')
specimenSummary <- data.frame(matrix(ncol=length(cols), nrow=length(SPECIMENS)))
colnames(specimenSummary) <- cols
specimenSummary$name <- SPECIMENS

for (s in 1:length(SPECIMENS)) {

	subSpecDist <- pData2[(pData2$name == SPECIMENS[s]),]
	specimenQuantiles <- wtd.quantile(subSpecDist$value, subSpecDist$aProb*10000)
	specimenSummary[s,'median age'] <- specimenQuantiles[3]
	specimenSummary[s,'minimum age'] <- specimenQuantiles[1]
	specimenSummary[s,'maximum age'] <- specimenQuantiles[5]
	specimenSummary[s,'age range'] <- specimenQuantiles[5] - specimenQuantiles[1]
	specimenSummary[s,'specimen age variability (IQR)'] <- specimenQuantiles[4] - specimenQuantiles[2]
	
}

depths <- pData2[, c("name", "sample", "depth", "reef", "RandL")]
depths2 <- dplyr::distinct(depths)
specimenSummary <- dplyr::left_join(specimenSummary, depths2, by = "name")
site <- substring(specimenSummary$name,0,2)
specimenSummary <- cbind(site,specimenSummary)

write.csv(specimenSummary, file=here::here('specimenSummary_2021-10-20.csv'))

## SAMPLE SUMMARY TABLE
#stats <- c('RandL','RandL_n', 'RandL_mean_age', 'RandL_median_age', 'RandL_age_range', 'RandL_minimum_age', 'RandL_maximum_age', 'RandL_median_specimen_age_variability_(IQR)', 'RandL_total_age_variability_(IQR)', 'RandL_median_specimen_age-estimation_error', 'RandL_corrected_posterior_age_estimate') 
stats <- c('sample', 'Sample_n', 'Sample_mean_age', 'Sample_median_age', 'Sample_age_range', 'Sample_minimum_age', 
           'Sample_maximum_age', 'Sample_median_specimen_age_variability_(IQR)', 'Sample_total_age_variability_(IQR)',
           'Sample_median_specimen_age-estimation_error', 'Sample_corrected_posterior_age_estimate')

sampleSummary <- data.frame(matrix(nrow=length(SAMPLES), ncol=length(stats)))
colnames(sampleSummary) <- stats
sampleSummary$sample <- SAMPLES

for (s in 1:length(SAMPLES)) {
	
	specSum <- specimenSummary[(specimenSummary$sample == SAMPLES[s]),]
	
	sampleSummary[s,'Sample_n'] <- nrow(specSum)
	sampleSummary[s,'Sample_median_specimen_age-estimation_error'] <- median(specSum[,'specimen age variability (IQR)'])
	sampleSummary[s,'Sample_median_specimen_age_variability_(IQR)'] <- IQR(specSum[,'median age'])

	subSampDist <- sAgeDist[(sAgeDist$sample == SAMPLES[s]),]

	sampleSummary[s,'Sample_mean_age'] <- wtd.mean(subSampDist$value, subSampDist$aProb*10000)
	
	sampleQuantiles <- wtd.quantile(subSampDist$value, subSampDist$aProb*10000)
	sampleSummary[s,'Sample_minimum_age'] <- sampleQuantiles[1]
	sampleSummary[s,'Sample_median_age'] <- sampleQuantiles[3]
	sampleSummary[s,'Sample_maximum_age'] <- sampleQuantiles[5]
	sampleSummary[s,'Sample_age_range'] <- sampleQuantiles[5] - sampleQuantiles[1]
	sampleSummary[s,'Sample_total_age_variability_(IQR)'] <- sampleQuantiles[4] - sampleQuantiles[2]
}

	sampleSummary[,'Sample_corrected_posterior_age_estimate'] <- sampleSummary[,'Sample_total_age_variability_(IQR)'] - sampleSummary[,'Sample_median_specimen_age-estimation_error']
  depths3 <- unique(depths2[,2:5])
	sampleSummary <- dplyr::left_join(sampleSummary, depths3, by = "sample")
	
write.csv(sampleSummary, file=here::here('sampleSummary_2021-10-20.csv'))

sampleSummary
