######################################
##    Functions for DGM R version   ##
##      Brooke Davis                ##
##    started March 10, 2017        ##
##   Brooke.davis@dfo-mpo.gc.ca     ##
######################################

# Misc useful function -- adaptation of dplyr %in% function
`%notin%` <- function(x,y) !(x %in% y) 

# Read in all csv data files
read.data <- function(){
  # Load Data
  # Base period Exploitation Rates
  BPER <<- read.csv("DataIn/BasePeriod_ER.csv")
  
  #Regional Distributions coefficients multipool, only for 1979
  RDistMP <<- read.csv("DataIn/CohortRegionDist_MultiPool.csv")
  
  # can extract max age and number of periods, regions from this file
  # extract ages
  Ages <<- unique(RDistMP$Age)
  # Max Age
  MaxAge <<- max(Ages)
  # only ages 2-maxAge are modelled
  NAges <<- length(Ages)
  #Number of stocks
  Stocks <<- unique(RDistMP$StockID)
  NS <<- length(unique(Stocks))
  # Number of periods
  NP <<- length(unique(RDistMP$Period))
  # number of Regions
  NR <<- length(unique(RDistMP$FromReg))
  
  #Regional Distribution coefficients, single pool
  RDistSP <<- read.csv("DataIn/CohortRegionDist_SinglePool.csv")
  
  # Base year cohort abundance by age and markings (for hatchery fish)
  BaseCohort <<- read.csv("DataIn/CohortSize_Age.csv")
  # save marking types
  AllNames <<- names(BaseCohort)
  MarkTypes <<- AllNames[AllNames %notin% c("StockID", "Age", "Abundance")]
 
  
  # Read in fisheries info
  FisheryInfo <<- read.csv("DataIn/CTCFisheries.csv")
  # Get number of gears -- in future might be different number per region
  Gears <<- unique(FisheryInfo$Sector)
  NG <<- length(Gears)
  #Number of fisheries
  NF <<- length(unique(FisheryInfo$FisheryID))
  
  #Read in CTC stock info
  StocksInfo <<- read.csv("DataIn/CTCStocks.csv")
  
  # Fishery period distributions - "Fish allocation over time periods" for each year
  FisheryDist <<- read.csv("DataIn/Fishery_Dist_Period.csv")
  
  # Fishery harvest Rate scalars
  HRS <<- read.csv("DataIn/Fishery_HrScalar.csv")
  
  # Marked retention -- ***need to come back to this******
  #MarkedRet <<- read.csv("DataIn/Fishery_MarkedRetention.csv")
  
  # Fishery release mortality, for 1979
  ReleaseMort <<- read.csv("DataIn/Fishery_Mortality.csv")
  
  # Fishery size limits -- for 1979, for one fishery for 1979, 1989, 1991
  SizeLims <<- read.csv("DataIn/Fishery_SizeLimit.csv")
  
  # Group cohort survival
  Grp_Surv <<- read.csv("DataIn/GroupSurvival.csv")
  
  # Regional Sector Allocations -- proportion of fish "allocated" to each sector (net, sport etc)
  ## *** come back to this **** strange coding
  SectorAlloc <<- read.csv("DataIn/Region_Sector_Alloc.csv")
  
  # Sampling Rates  ** what is this used for? -- for sampling moduule
  #SampingRates <<- read.csv("DataIn/SamplingRates.csv")
  
  # Size at Age by stock
  Size <<- read.csv("DataIn/Stock_Age_Size.csv")
  
  # Maturation by stock, age -- period in which they mature is in StocksInfo
  Maturation <<- read.csv("DataIn/Stock_MatRate.csv")
  
  # Stock PSM -- represented as survival rate (Surv_Rate)
  PreSpawnSurv <<- read.csv("DataIn/Stock_PrespawnMort.csv")
  
  # Terminal harvest rate by stock, fishery, age
  TermHR <<- read.csv("DataIn/Stock_TermHR.csv")
  
  # Terminal harvest rate scalar *** what is this for? why is it only for base year??*******
  TermHRS <<- read.csv("DataIn/Stock_TermHRScalar.csv")
  
  # Proportion of each class vulnerable to a fishery
  PV <<- read.csv("DataIn/Prop_Age_Vulnerable.csv")
  
  # Terminal fishery mort rates
  TermRelMort <<- read.csv("DataIn/Stock_TermMortRetention.csv") 
} # end read.data function



#########################################################
# Not sure if thsi is right, but this is what is in the VB code
AABM.TAC <- function(AI){
  # create vector to store total TAC, length NR
  Total_TAC <- NULL
  # TAC in array NR * NSec
  TAC <- array(dim=c(NR, NG))
  # REgion 1
    if (AI[1] < 1.005){
      Total_TAC[1] = 17000 + 110500 * AI[1]
      TAC[1,1] = 110500 * AI[1] * 0.8 #troll TAC
    }
    if (AI[1] >= 1.005 & AI[1] < 1.2) {
      Total_TAC[1] = -114750 + 242250 * AI[1]
      TAC[1,1] = (-131750 + 242250 * AI[1]) * 0.8 #troll TAC
    }
    if (AI[1] >= 1.2 & AI[1] < 1.5) {
      Total_TAC[1] = 17000 + 151721 * AI[1]
      TAC[1,1] = 151721 * AI[1] * 0.8 #troll TAC
    }
    if (AI[1] >= 1.5) {
      Total_TAC[1] = 17000 + 164364 * AI[1]
      TAC[1,1] = 164364 * AI[1] * 0.8 #troll TAC
    }
    
    TAC[1,2] = TAC[1,1] * 0.25 #sport TAC
    TAC[1,3] = 17000 # net TAC
  
 # Region 2
    if (AI[2] < 1.205) {
      Total_TAC[2] = 130000 * AI[2]
      TAC[2,1] = Total_TAC[2] * 0.8 #troll TAC
    }
    if (AI[2] >= 1.205 & AI[2] < 1.5) {
      Total_TAC[2] = -20000 + 146667 * AI[2]
      TAC[2,1] = Total_TAC[2] * 0.8 #troll TAC
    }
    if (AI[2] >= 1.5) {
      Total_TAC[2] = 145892 * AI[2]
      TAC[2,1] = Total_TAC[2] * 0.8 #troll TAC
    }
    
    TAC[2,2] = TAC[2,1] * 0.25 #sport TAC
    TAC[2,3] = 0 # net fishery is ISBM

# Region 3
    if (AI[3] < 0.5) {
      Total_TAC[3] = 128347 * AI[2]
      TAC[3,1] = Total_TAC[3] * 0.8 #troll TAC
    }
    if (AI[3] >= 0.5 & AI[3] < 1.0) {
      Total_TAC[3] = 149739 * AI[3]
      TAC[3,1] = Total_TAC[3] * 0.8 #troll TAC
    }
    if (AI[3] >= 1.0) {
      Total_TAC[3] = 171130 * AI[3]
      TAC[3,1] = Total_TAC[3] * 0.8 #troll TAC
    }
    
    TAC[3,2] = TAC[3,1] * 0.25 #sport TAC
    TAC[3,3] = 0 # net fishery is ISBM
    
    # just return total over regions
    out <- Total_TAC
    
    
    #### Want to return a TAC by r,y,f according to sector and region

} # end AABM.TAC function

#####################################################################3
## DGM TAC function takes "fishery+catch" (N*BPER*PV -- numercator of AI calc) 
## and AI to get TAC according to HRI, rather than just a number
 
AABM.TAC.HRI <- function(numer, AI) {
  # store TAC, HRI in vectors length 3 (each region)
  HRI <- NULL
  TAC <- NULL
  
  ## Region 1
  if(ai[1] < 1.005){
  HRI[1] = 0.371
  } else if (ai[1] >= 1.005 & ai[1] < 1.2) {
  HRI[1] = 0.3795 * ai[1] - 0.0104
  } else if (ai[1] >= 1.2 & ai[1] <= 1.5) {
  HRI[1] = 0.51
  } else if (ai[1] > 1.5) {
  HRI[1] = 0.5525
  }
  
  # Region 2
  if (ai[2] < 1.205) {
  HRI[2] = 0.757
  } else if (ai[2] >= 1.205 & ai[2] <= 1.5) {
  HRI[2] = 0.3153 * ai[2] + 0.3771
  } else if (ai[2] > 1.5) {
  HRI[2] = 0.85
  }
  
  # Region 3
  if (ai[3] < 0.5) {
  HRI[3] = 0.21
  } else if (ai[3] >= 0.5 & ai[3] <= 1.0) {
  HRI[3] = 0.245
  } else if (ai[3] > 1.0) {
  HRI[3] = 0.28
  }
  

  TAC = 1.25 * numer * HRI
  
  # region one Net is AABM, needs additional TAC
  TAC[1] <- TAC[1]+ 17000 
  
  TAC
} # end AABM.TAC.HRI function


###################################################
## Proportion of each stock age of legal size
## Can be done outside of simulation
# Get proportion of fish from each stock and age that will be sub-legal
GetPRel <- function(NY, NP, NR, NG, NS, NAges){
  PRel <- array(dim=c(NY, NP, NR, NG, NS, NAges))
  # Preterm ff's
  Preterm_ff <- FisheryInfo$FisheryID[which(FisheryInfo$Type !="TERM")]
  for(ff in Preterm_ff){  
    region <- FisheryInfo$Region[which(FisheryInfo$FisheryID==ff)]
    gear <- which(Gears == FisheryInfo$Sector[which(FisheryInfo$FisheryID==ff)])
    for(yy in 1:NY){
      # get size limits for most recent year
      DatYrs <- SizeLims$Year[which(SizeLims$FisheryID==ff)]
      Smin <- SizeLims$Min[which(SizeLims$FisheryID==ff & SizeLims$Year==max(DatYrs[DatYrs<=Years[yy]]))]
      Smax <- SizeLims$Max[which(SizeLims$FisheryID==ff & SizeLims$Year==max(DatYrs[DatYrs<=Years[yy]]))]
        for(pp in 1:NP){
          for(ss in 1:NS){
            for(aa in 1:NAges){
              # Get Size Distribution parameters for this stock, age, period
              mu <- Size$Mean[which(Size$StockID==Stocks[ss] & Size$Age==Ages[aa] & Size$Period==pp)]
              sigma <- Size$SD[which(Size$StockID==Stocks[ss] & Size$Age==Ages[aa] & Size$Period==pp)]
              PRel[yy, pp,  region, gear, ss, aa] <- pnorm(Smin, mu, sigma) + 1 - pnorm(Smax, mu, sigma)
            } # end age loop
          } # end stock loop
        } # end period loop
      } # end year loop
    } # end fisheries loop
PRel
} # end GetPRel function

#########################################
## Same for terminal fisheries
#########################################
# for terminal fisheries will just be by year, stock, age
GetPRel_Term <- function( NY, NS, NAges){
  PRel <- array(dim=c( NY, NS, NAges))
  Term_ff <- FisheryInfo$FisheryID[which(FisheryInfo$Type == "TERM")]
  for(yy in 1:NY){
    for(ff in Term_ff){  
      stk <- TermHR$StockID[which(TermHR$FisheryID==ff)]
      # Use limits for most recent year
      DatYrs <- TermRelMort$Year[which(TermRelMort$FisheryID==ff)]
      Smin <- TermRelMort$SizeLimit_Lower[which(TermRelMort$FisheryID==ff & TermRelMort$Year==max(DatYrs[DatYrs<=Years[yy]]))]
      Smax <- TermRelMort$SizeLimit_Upper[which(TermRelMort$FisheryID==ff & TermRelMort$Year==max(DatYrs[DatYrs<=Years[yy]]))]
      if(is.na(Smax)){
        Smax <- 9999
      }
      for(aa in 1:NAges){
        # Get Size Distribution parameters for this stock, age, assuming terminal fishing occurs in period 3
        mu <- Size$Mean[which(Size$StockID==stk & Size$Age==Ages[aa] & Size$Period==3)]
        sigma <- Size$SD[which(Size$StockID==stk & Size$Age==Ages[aa] & Size$Period==3)]
        PRel[yy, stk, aa] <- pnorm(Smin, mu, sigma) + 1 - pnorm(Smax, mu, sigma)
      } # end age loop
    } # end fisheries loop
  } # end year loop
  PRel
} # end GetPRel function

#############################
##  Get Age 1 to Adult Surv
#############################

# Given age 1-5 survival and age 2-5 maturation, get age 1 to
# adult survival rate

GetAge1Surv <- function(S=c(0.5, 0.6, 0.7, 0.8, 0.9), M=c(NA, 0.0737, 0.0573, 0.7087, 1.0000)){
  # Make flexible in case have more age classes in future
  N <- length(S)
  SEQ <- NULL
  SEQ[1] <- NA
  SEQ[2] <- S[1]*S[2]
  for(i in 3:5){
    SEQ[i] <- SEQ[i-1]*(1-M[i-1])*S[i]
  }
  sum(SEQ*M, na.rm=T)
}
