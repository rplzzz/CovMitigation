# COVID-19 SIRS projections

###########################################################################
# PART ONE PART ONE PART ONE ----------------------------------------------


# Install packages --------------------------------------------------------
require(tidyverse)
require(deSolve)
options(scipen = 999)
library(usmap)
library(CovMitigation)

# Simulation Principles ---------------------------------------------------
#The goal is to simulate the spread of the COVID infection in each virginia county and then
#use that information plus other assummptions about health care utilization to estmate demand

#the simulation runs in two loops(external loop and internal loop)
#the extermal loop contains variables for key values like  R0, S(susceptibleFraction), M factor for mitigation, D for duration

# the internal loop runs the SIR model for x days for each county in virginia and stores the output
#important variables to calculate include beta

#The E loop has two components, a dataframe containing baseline parameter values and a vector containing one
#of the baseline values that will be varied.

#Exmple: we store baseline values of R0, D, gamma, Sfraction, T (doubling time), beta) these are the fixed values for an E run
# But to control the number and sensitivity analysis of an E run, we supply a vector of length n  and we use
# the values of this vector to run the E loop n times using the values from teh vector



# Start Simulation Run Here -----------------------------------------------

#####    Step One                                   ################
#####    Select VCU or UVA and select market Share  ################
typeAMCDRG <- "fractionAllUVA"   #or "fractionAllVCU"  "fractionRespUVA" "fractionRespVCU" "fractionAllUVA"
selectMarketShare  <-  0.2    #select a fracton from 0.0 to 0.99

#filter the counties who meet the market share threshold
countySelection <- marketFractionFinal %>%
  filter(marketFractionFinal$TypeAMCDRG ==typeAMCDRG & marketFractionFinal$marketShare >= selectMarketShare)


#create the county level map- what counties refer to UVA- a visual display tool

#make the map
referralMap <- plot_usmap(data = sampleCounties, include = "VA", values = "marketShare", color = "white", size = 0.2) + 
  scale_fill_continuous(name = "AMC Market Share", label = scales::comma) + 
  theme(legend.position = "right") + theme(legend.text=element_text(size=16)) + theme(legend.title=element_text(size=20)) +
  ggtitle(paste("Counties with AMC Market Share Exceeding", selectMarketShare*100,"Percent")) + theme(plot.title = element_text(size = 24, face = "bold"))
print(referralMap)


#####    Step Two                                  ################
#####    Verify the Fixed Variable Values          ################
#These are the fixed variable settings that will be used unless overridden
R0 <- 2.5                               #R0----used to back calculate beta
durationInfection <- 7                  #duration of infectivity
gamma <- 1/durationInfection            #recovery time
S <- 1                                #fraction of total population (N) susceptible to infection
T <- 7                                  #doubling time (not currently used)
M <- 1                                  #mitigation fraction
symptoFraction <- 0.50                  #fraction of cases that are infected

#This is a datafrome of fixed variable label names
ErunVariableLabels <- as.data.frame(c("RO", "durationInfection", "gamma", "fractionS", "doubleTime", "mitigation fraction"))
names(ErunVariableLabels) ="Labels"
ErunVariableLabels$Labels <-as.character(ErunVariableLabels$Labels)

# make some assumptions about hospitalization and ICU and mechanical ventilation
hospFraction <- 0.03   # MMWR(20.7-31.4%)  range I think likely is half that 10.4 - 15.7% with a mean of 13%
# but even this is likely an overestimate so I'll start with half of that 6.5%

ratioICUtoAcute     <-   0.238          # from the same MMWR 121/508 cases admitted went to the ICU- that is 23.8%  of all admissions
# 26.4% from the 201 person case series Risk Factors Associated with ARDS, C Wu, JAMA 3/13/20
fractionMV          <-   0.328         # from the 201 person case series Risk Factors Associated with ARDS, C Wu, JAMA 3/13/20
hospFractionNIPPV   <-   0.2605         # from the 201 person case series Risk Factors Associated with ARDS, C Wu, JAMA 3/13/20
hospFractionIMV     <-   0.0675         # from the 201 person case series Risk Factors Associated with ARDS, C Wu, JAMA 3/13/20

# time lags
timeToSympto      <- 3              #days to symptomatic infection
timeToAcuteHosp   <- 5              #days to symptomatic infection
timeToICU         <- 6              #days to symptomatic infection
timeToMV          <- 6              #days to symptomatic infection
timeToNIPPV       <- 6             #days to symptomatic infection
timeToIMV         <- 6              #days to symptomatic infection


#####    Step Three                                 ################
#####    Select Sensitivity Variable and Values     ################
#####    Options: R0, fractionS, Duration           ################
#Turn on the sensitivity vector for External run using one of these three statement groups
#have not built for gamma (inverse of duration) or doubling time
#You must comment out the other two statement groups

# scenarioVector <- as.data.frame(c(1.75, 2.25,  2.75))
# names(scenarioVector) ="value_R0"
# runCountELoop <- nrow(scenarioVector)
# scenarioLabel <- ErunVariableLabels[1,1]

# scenarioVector <- as.data.frame(c(4,5,7,10))
# names(scenarioVector) ="value_D"
# runCountELoop <- nrow(scenarioVector)
# sscenarioLabel <- ErunVariableLabels[2,1]
# 
# scenarioVector <- as.data.frame(c(0.3, 0.4, 0.5, 0.6))
# names(scenarioVector) ="value_S"
# runCountELoop <- nrow(scenarioVector)
# scenarioLabel <- ErunVariableLabels[4,1]

scenarioVector <- as.data.frame(c(5, 7, 10))
names(scenarioVector) ="doubleTime"
runCountELoop <- nrow(scenarioVector)
scenarioLabel <- ErunVariableLabels[5,1]


# scenarioVector <- as.data.frame(c(0.05, 0.1, 0.2, 0.3))
# names(scenarioVector) ="value_S"
# runCountELoop <- nrow(scenarioVector)
# scenarioLabel <- ErunVariableLabels[6,1]



#####    Step Four                                  ################
#####    Scenario Calculations Now Run              ################

# #Set the market share selections as well- decide whether to use UVA or VCU data and whether to use market fraction for resp or all drgs
# # countySelection <- marketFractionFinal %>%
# #   filter(marketFractionFinal$TypeAMCDRG =="fractionAllUVA" & marketFractionFinal$marketShare >= selectMarketShare)
# 
#  countySelection <- marketFractionFinal %>%
#    filter(marketFractionFinal$TypeAMCDRG =="" & marketFractionFinal$marketShare >= selectMarketShare)


# create a null dataframe to contain consecutive county level outpts
# # setting the initial values for the elements we want to track
sir_values_collect  <-data.frame(
  time=double(),
  S = double(),
  I = double(),
  R = double(),
  newCases  = double(),
  scenarioLabel    = character(),    #this will be scenarioLabel
  scenarioValue    = character(),   #this will be scenarioValue
  AMCDRG           = character(),         #this will be typeAMCDRG
  marketFraction   = double(),     #this will be selectMarketShare
  locality         = character(),        #this will be the county location
  population       = double(),          #thsi will be the baseline population within a county
  symptoInfection  = double(), 
  daysToSympto     = double(), 
  acuteHosp        = double(), 
  daysToAcuteHosp  = double(), 
  icuHosp          = double(), 
  daysToicuHosp    = double(), 
  anyMV            = double(), 
  daysToAnyMV      = double(), 
  NIPPVHosp        = double(), 
  daysToNIPPVHosp  = double(), 
  IMVHosp          = double(), 
  daysToIMVHosp    = double() 
)


##### this is the beginning of the extermal loop
for (EloopRun in 1:runCountELoop){
  scenarioValue <- (scenarioVector[EloopRun,1])
  #R0 <- scenarioValue
  T <- scenarioValue
  
  #write the basic SIR equations and create a null dataframe t0 hold the county-by-county SIR calculations
  
  #This is the start of the internal loop
  #The internal loop runs through the counties
  #the county selection dataframe contains the counties we need to assess
  #so first we need to know how long the loop should be by counting the rows in the countySelection dataframe
  #call this variable countyCount
  #use it to set the length of the internal loop
  
  lengthCounty <- nrow(countySelection)
  for(countyRuns in 1:lengthCounty){
    
    #determine the locality and population variables for this county run. we will store them later
    locality <- as.character(countySelection[countyRuns,1])                              # this sets the label for location 
    population <- as.numeric(vaMyAgeBands[vaMyAgeBands$Locality == locality,2])          # use the corresponding population for this location
    
    #begin the SIR equation now
    #set the initial values and labels for each scenario
    
    N <- population * S
    parameters_values <- c(
      #beta <- R0/(N*durationInfection) * M,      # infectious contact rate (/person/day) * mitigation factor
      beta =  (1 + (log(2)*durationInfection)/T)/(N*durationInfection) * M,
      gamma = 1/durationInfection                # recovery rate (/day)
    )
    
    
    initial_values <- c(
      #S = 10000,                  # number of susceptibles at time = 0
      S = population * S,          # number of susceptibles at time = 0
      I =   1,                     # number of infectious at time = 0
      R =   0                      # number of recovered (and immune) at time = 0
    )
    
    #set the duration
    time_values <- seq(0, 180)      # days
    
    #check the parameter values
    parameters_values
    initial_values
    time_values
    
    sir_values_1 <- ode(
      y = initial_values,
      times = time_values,
      func = sir_equations,
      parms = parameters_values 
    )
    
    #calculate the new Cases as the daily drop in susceptible cases
    sir_values_1 <- as.data.frame(sir_values_1)
    newCase <- as.data.frame(diff(sir_values_1$S)*-1)          #difference the S cases to get the drop but order of diff is wrong so * negative 1
    zeroDay <- as.data.frame(0)                                #create 0 zero value for Day zero - no new cases that day
    names(newCase) <- "newCases"                               #name them for binding purposes
    names(zeroDay) <- "newCases"                               #name them for binding purposes
    newCase <- rbind(zeroDay, newCase)                         #bind them
    
    sir_values_1 <- cbind(sir_values_1, newCase)                       #join back to the sir_values_1 to report out
    
    
    
    
    #store the output and attach labels for the locations and scenarios
    #sir_values_1 <- as.data.frame(sir_values_1)
    sir_values_1$scenarioLabel    <- scenarioLabel
    sir_values_1$scenarioValue    <- scenarioValue
    sir_values_1$AMCDRG           <- typeAMCDRG
    sir_values_1$marketFraction   <- selectMarketShare
    sir_values_1$locality         <- locality
    sir_values_1$population       <- population
    
    # calculate hospitalization numbers and days to hospitalization etc
    sir_values_1$symptoInfection <- sir_values_1$newCases * symptoFraction
    sir_values_1$daysToSympto   <- sir_values_1$time + timeToSympto
    sir_values_1$acuteHosp <- sir_values_1$newCases * symptoFraction * hospFraction * (1-ratioICUtoAcute)
    sir_values_1$daysToAcuteHosp   <- sir_values_1$time + timeToAcuteHosp
    sir_values_1$icuHosp   <- sir_values_1$newCases * symptoFraction * hospFraction * ratioICUtoAcute
    sir_values_1$daysToicuHosp   <- sir_values_1$time + timeToICU
    sir_values_1$anyMV     <- sir_values_1$newCases * symptoFraction * hospFraction* fractionMV
    sir_values_1$daysToAnyMV   <- sir_values_1$time + timeToMV
    sir_values_1$NIPPVHosp <- sir_values_1$newCases * symptoFraction *  hospFraction * hospFractionNIPPV
    sir_values_1$daysToNIPPVHosp   <- sir_values_1$time + timeToNIPPV
    sir_values_1$IMVHosp   <- sir_values_1$newCases * symptoFraction  *  hospFraction * hospFractionIMV
    sir_values_1$daysToIMVHosp   <- sir_values_1$time + timeToIMV
    
    sir_values_collect <- rbind(sir_values_collect, sir_values_1)
    
  } 
}



# Part 2: Calculating hospitalization and ICU admission rates -------------
# take the SIR model output and make a copy for graphing

inpatientEstimates <- sir_values_collect


# epidemic curve of new cases in VA that will eventually show up at UVA
newCasesBoundToUVA <- inpatientEstimates %>%
  mutate(scenarioValue = factor(scenarioValue)) %>%
  select(time, scenarioLabel, scenarioValue, marketFraction, locality, newCases) %>%
  group_by(scenarioLabel, scenarioValue, time) %>%
  summarize(newCasesPerDay = sum(newCases))
plotNewCasesBoundToUVA <-ggplot(newCasesBoundToUVA, aes(x = time, y = newCasesPerDay, group = scenarioValue, colour = scenarioValue)) + geom_line(size=1.25) +
  ggtitle("New Community Infections in UVA Referral Base Epidemiologic Curve") + 
  scale_x_continuous(breaks=seq(0, 180, by= 5 )) +geom_vline(xintercept = c(0, 30, 60, 90, 120, 150, 180)) 
print(plotNewCasesBoundToUVA)


newSymptoCasesBoundToUVA <- inpatientEstimates %>%
  mutate(scenarioValue = factor(scenarioValue)) %>%
  select(daysToSympto, scenarioLabel, scenarioValue, marketFraction, locality, symptoInfection) %>%
  group_by(scenarioLabel, scenarioValue, daysToSympto) %>%
  summarize(newCasesPerDay = sum(symptoInfection))
plotNewSymptoCasesBoundToUVA <-ggplot(newSymptoCasesBoundToUVA, aes(x = daysToSympto, y = newCasesPerDay, group = scenarioValue, colour = scenarioValue)) + geom_line(size=1.25) +
  ggtitle("New Symptomatic Community Infections in UVA Referral Base Epidemiologic Curve") + ylab("New Symptomatic Infections in UVA Referral Area") +
  scale_x_continuous(breaks=seq(0, 180, by= 5 )) +geom_vline(xintercept = c(0, 30, 60, 90, 120, 150, 180)) 
print(plotNewSymptoCasesBoundToUVA)


newAcuteHospBoundToUVA <- inpatientEstimates %>%
  mutate(scenarioValue = factor(scenarioValue)) %>%
  select(daysToAcuteHosp, scenarioLabel, scenarioValue, marketFraction, locality, acuteHosp) %>%
  group_by(scenarioLabel, scenarioValue, daysToAcuteHosp) %>%
  summarize(newAcutePerDay = sum(acuteHosp))
plotNewAcuteHospBoundToUVA <-ggplot(newAcuteHospBoundToUVA, aes(x = daysToAcuteHosp, y = newAcutePerDay, group = scenarioValue, colour = scenarioValue)) + geom_line(size=1.25)+
  ggtitle("New Acute Care Admissions for UVA Epidemiologic Curve") + 
  scale_x_continuous(breaks=seq(0, 180, by= 5 ))  + geom_vline(xintercept = c(0, 30, 60, 90, 120, 150, 180))
print(plotNewAcuteHospBoundToUVA)

newICUHospBoundToUVA <- inpatientEstimates %>%
  mutate(scenarioValue = factor(scenarioValue)) %>%
  select(daysToicuHosp, scenarioLabel, scenarioValue, marketFraction, locality, icuHosp) %>%
  group_by(scenarioLabel, scenarioValue, daysToicuHosp) %>%
  summarize(newICUPerDay = sum(icuHosp))
plotNewICUHospBoundToUVA <-ggplot(newICUHospBoundToUVA, aes(x = daysToicuHosp, y = newICUPerDay, group = scenarioValue, colour = scenarioValue)) + geom_line(size=1.25)+
  ggtitle("New ICU Admissions for UVA Epidemiologic Curve") + 
  scale_x_continuous(breaks=seq(0, 180, by= 5 ))  + geom_vline(xintercept = c(0, 30, 60, 90, 120, 150, 180))
print(plotNewICUHospBoundToUVA)


newIMVHospBoundToUVA <- inpatientEstimates %>%
  mutate(scenarioValue = factor(scenarioValue)) %>%
  select(daysToIMVHosp, scenarioLabel, scenarioValue, marketFraction, locality, IMVHosp) %>%
  group_by(scenarioLabel, scenarioValue, daysToIMVHosp) %>%
  summarize(newIMVPerDay = sum(IMVHosp))
plotNewIMVHospBoundToUVA <-ggplot(newIMVHospBoundToUVA, aes(x = daysToIMVHosp, y = newIMVPerDay, group = scenarioValue, colour = scenarioValue)) + geom_line(size=1.25)+
  ggtitle("New Invasive Mechanical Ventilation Admissions for UVA Epidemiologic Curve") + 
  xlab("time (days") + ylab("New Daily Admissions for IMV") +
  scale_y_continuous(breaks=1:40) +
  scale_x_continuous(breaks=seq(0, 180, by= 5 ))  + geom_vline(xintercept = c(0, 30, 60, 90, 120, 150, 180))
print(plotNewIMVHospBoundToUVA)

