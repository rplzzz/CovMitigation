## Run from package top level

library('dplyr')
library('tidyr')
library('here')

datadir <- here('data-raw')

#### Import, aggregate and clean up census data :  Result - vaMyAgeBands
#Import, clean & aggregate demographic data 2018 Census Bureau Estimates-- 130 counties and cities in virginia
#import census data from Weldon Cooper Center
vaDemo <- read.csv(file.path(datadir,"Census_2018_AgeSexEstimates_forVA_2019-08.csv"),
                   sep = ",", header = TRUE, stringsAsFactors = FALSE)

#Collapse Age bands into a smaller set for analysis
vaDemo$Age0to5 <- vaDemo$MUnder5 + vaDemo$FUnder5
vaDemo$Age5to19 <- vaDemo$M5to9    + vaDemo$M10to14 + vaDemo$M15to19 +
  vaDemo$F5to9    + vaDemo$F10to14 + vaDemo$F15to19 
vaDemo$Age20to39 <- vaDemo$M20to24  + vaDemo$M25to29 + vaDemo$M30to34 + vaDemo$M35to39  + 
  vaDemo$F20to24  + vaDemo$F25to29 + vaDemo$F30to34 + vaDemo$F35to39  
vaDemo$Age40to59 <- vaDemo$M40to44 + vaDemo$M45to49 + vaDemo$M50to54  + vaDemo$M55to59 + 
  vaDemo$F40to44 + vaDemo$F45to49 + vaDemo$F50to54  + vaDemo$F55to59 
vaDemo$Age60to74 <- vaDemo$M60to64 + vaDemo$M65to69  + vaDemo$M70to74  + 
  vaDemo$F60to64 + vaDemo$F65to69  + vaDemo$F70to74
vaDemo$Age75plus <- vaDemo$M75to79  + vaDemo$M80to84  + vaDemo$M85andover +
  vaDemo$F75to79  + vaDemo$F80to84  + vaDemo$F85andover 

# select the age bands, current program does not use them but might be useful in a age banded simulation (NEXT VERSION)
vaMyAgeBands <- vaDemo %>%
  select(Locality, Total, Age0to5, Age5to19, Age20to39, Age40to59, Age60to74, Age75plus)

#### join the agebanded census population data with the RESPIRATORY DRG market share data into ageMarketShare (% values) and ageMarketFraction (fraction values)
mktShareRespDRG <-  read.csv(file.path(datadir,"marketShareRespDRG_UVA_VCU.csv"),
                             sep = ",", header = TRUE, stringsAsFactors = FALSE)
# join the dataframes- this is for market share as a percentage, just for respiratory DRGS
ageMarketShareResp <- full_join(mktShareRespDRG, vaMyAgeBands, by = "Locality")
# drop age band columns and hospitalizations from patients living in other states
ageMarketShareResp <- ageMarketShareResp %>%
  select(Locality, shareUVA, shareVCU) %>%
  filter(Locality != "Virginia")

# recreate the same dataframe with market share as a fraction
ageMarketFractionResp <- ageMarketShareResp
ageMarketFractionResp$shareUVA <- ageMarketFractionResp$shareUVA/100
ageMarketFractionResp$shareVCU <- ageMarketFractionResp$shareVCU/100
names(ageMarketFractionResp) <-c("Locality", "fractionRespUVA", "fractionRespVCU")

#### Repeat for the ALL DRG market share values
# join the agebanded population data with the ALL DRR market share data into ageMarketShare (% values) and ageMarketFraction (fraction values)
mktShareAllDRG <-  read.csv(file.path(datadir,"marketShareAllDRG_UVA_VCU.csv"),
                            sep = ",", header = TRUE, stringsAsFactors = FALSE)

# join the dataframes
ageMarketShareAll<- full_join(mktShareAllDRG, vaMyAgeBands, by = "Locality")

# select only the columns needed
ageMarketShareAll <- ageMarketShareAll %>%
  select(Locality, shareUVA, shareVCU) %>%
  filter(Locality != "Virginia")

# recreate the same dataframe with market share as a fraction
ageMarketFractionAll <- ageMarketShareAll
ageMarketFractionAll$shareUVA <- ageMarketFractionAll$shareUVA/100
ageMarketFractionAll$shareVCU <- ageMarketFractionAll$shareVCU/100
names(ageMarketFractionAll) <-c("Locality", "fractionAllUVA", "fractionAllVCU")

marketFractionFinal <- full_join(ageMarketFractionAll, ageMarketFractionResp, by = "Locality")
marketFractionFinal <- gather(marketFractionFinal, key = TypeAMCDRG, value = marketShare, fractionAllUVA:fractionRespVCU)

#### Read in and join up the county and map location data
vaFIPScodes <- read.csv(file.path(datadir,"VA_county_FIPScodes.csv"), 
                        sep = ",", header = TRUE, stringsAsFactors = FALSE)     #read in the county FIPS 
sampleCounties <- inner_join(vaFIPScodes, countySelection, vaFIPScodes, by = "Locality")                   #join fips codes with counties selected
names(sampleCounties) <- c("fips", "Locality", "TypeAMCDRG", "marketShare")                                #label correctly so plot_map() can use


usethis::use_data(vaMyAgeBands, marketFractionFinal, sampleCounties)
