####################################################
###    Chinook DGM R Test version                 ##
###       Brooke Davis                            ##
###      Started March 10, 2017                   ##
###     brooke.davis@dfo-mpo.gc.ca                ##
### Re-start 5/12 trying to re-create DGM results ##
####################################################

# Set Working Directory
setwd("C:/Users/DavisBr/Documents/Chinook MSE/R Test")

# Load packages
library(dplyr)

# source functions
FPath <- "R/Functions.r"
source(FPath)

#read in data
read.data()

#####################################
# Set years for which are simulating
Years <- 1979 : 1985
NY <- length(Years)
######################################
# Set to Single pool or multipool
#PoolType <- "MP"
PoolType <- "SP"
######################################
#** maybe read in from control file lateR?


###################################
#Set up arrays to store values
###################################

# Stick to indexing order
# Year | period | region | fishery | Sector | stock | age

#Cohort -- starting abundance before distribution, 
#represents "new" fish to model that need to be distributed
# note Cohort is indexed by age (1-5), where N is indexed by adult age (2-MaxAge)
Cohort <- array(dim=c(NY, NS, NAges+1)) 
#abundance N -- indexes by stock, year, period, region, adult ages age (ages 2-MaxAge indexed as 1-(MaxAge-1))
N <- array(dim=c(NY, NP, NR, NS, NAges)) 
# have seperated for each step to make checks easier
# after pre-terminal fishing
N1 <- array(dim=c(NY, NP, NR, NS, NAges))
# after release mortality
N2 <- array(dim=c(NY, NP, NR, NS, NAges)) 
# after dropoff mortality
N3 <- array(dim=c(NY, NP, NR, NS, NAges)) 
# after maturity (only changes in period where stock matures)
N4 <- array(dim=c(NY, NP, NR, NS, NAges)) 
# Harvest rate array -- indexed by stock, year, period, age, gear
h <- array(dim=c( NY, NP, NR, NG, NS, NAges)) 
# TAC by year, period, region, gear, stock, ages -- only for AABM
TAC <- array(dim=c( NY, NP, NR, NG))
# Catch C -- indexed by stock, year, period, region, age, sector
Catch <- array(dim=c( NY, NP, NR, NG, NS, NAges))
# Keep terminal catch separate
Catch_Term <- array(dim=c(NY, NS, NAges))
# Abundance Index
AI <- array(dim=c(NY, NR))
# Proportion of sub-legal size (and therefore released), for non-terminal fisheries
# Since can be done outside of simulation, just calculate using function
PRel <- GetPRel(NY,NP, NR, NG, NS, NAges)
# For terminal fisheries only vary by year, stock, and age
PRel_Term <- GetPRel_Term(NY, NS, NAges)
# Releases
Release <- array(dim=c(NY, NP, NR, NG, NS, NAges))
Release_Term <- array(dim=c(NY, NS, NAges))
# Release Mortality
RelMort <- array(dim=c(NY, NP, NR, NG, NS, NAges))
RelMort_Term <- array(dim=c(NY, NS, NAges))
# DropOff Mortality
Dropoff <- array(dim=c(NY, NP, NR, NG, NS, NAges))
Dropoff_Term <- array(dim=c(NY, NS, NAges))
# Mature Fish
Mature <- array(dim=c(NY, NP, NR, NS, NAges))
# Escapement (after terminal fishing)
Escape <- array(dim=c(NY, NS, NAges))
# Spawners (after pre-spawn mort)
Spawners <- array(dim=c(NY, NS))
# Recruits
Recruits <- array(dim=c(NY, NS))

############################
## Start Simualtion here  ##
############################
# Cycle over years
for(yy in 1:length(Years)){
  ################################
  ## Set inital Cohort values   ##
  ################################
  ## If first year set first year cohort abundance, if not first year need to add new year 2's to pop and distribute
  if(yy==1){
    # Put values from BaseCohort into Cohort array
    for(ss in 1:NS){
      for(aa in 2:MaxAge){
        # Ages in array will start from 2
        Cohort[yy,ss,aa] <- BaseCohort$Abundance[which(BaseCohort$StockID==Stocks[ss] & BaseCohort$Age==(aa))]
      } # end age loop
      # need to get age one abundance using Age 2 survival
      # Need to know survival group
      grp <- StocksInfo$Surv_Grp[which(StocksInfo$StockID==Stocks[ss])]
      # Get survival from age 1 to 2
      Surv_1_2 <- Grp_Surv$SR_1_2[which(Grp_Surv$Grp==as.character(grp) & Grp_Surv$Year==Years[yy])]
      Cohort[yy, ss, 1] <-  Cohort[yy, ss, 2] / Surv_1_2
    } # end stock loop
  } else if(yy==2) { # end yy==1
    # Need age two's for second year, assume same as first year
     Cohort[2,,2] <- Cohort[1,,2]
  } else if(yy>2) { # end yy==2
    ####################
    ## Year 2+ Cohort  #
    ####################
    # when year > 2 need to fill in age 2's using last year's 1's applying survival
    for(ss in 1:NS){
      # get age 1-2 survival for appropriate survival group
      grp <- StocksInfo$Surv_Grp[which(StocksInfo$StockID==Stocks[ss])]
      S <- as.numeric(Grp_Surv[which(Grp_Surv$Grp==as.character(grp) & Grp_Surv$Year==Years[yy]), "SR_1_2"])
      Cohort[yy, ss, 2] <- Cohort[yy-1, ss, 1] * S
    } # end stock loop
  } # end yy >2 if


  
  ##########################
  ## Dist across regions  ##
  ##########################
  for(pp in 1:NP){
    # Current data only for 1979, will use for all years
    ##############
    ## Period 1 ##
    ##############
    # distribute cohort
    if(pp==1){ 
      #distribute cohort (new fish to model) using SP coeffs
      ## also need to include adults from last generation
      for(rr in 1:NR){
        for(ss in 1:NS){
          # In year 1 "cohort" includes adults, since all age classes "new" to model
            if(yy==1){
              for(aa in 1:NAges){
                # in year one N comes only from cohort, all ages are in cohort
                N[ yy, 1, rr, ss, aa] <- RDistSP[which(RDistSP$StockID==Stocks[ss] & RDistSP$Age==Ages[aa] & RDistSP$Period==1), 
                                               which(names(RDistSP)==paste("P", rr, sep="_"))]  *  Cohort[yy,ss,aa+1]
              } # end age loop
            } else { # end yy == 1
                # in subsequent only distribute age 2's ( which are entry 1 in N array) in this way
                N[yy, 1, rr, ss, 1] <-  RDistSP[which(RDistSP$StockID==Stocks[ss] & RDistSP$Age==2 & RDistSP$Period==1), 
                                                         which(names(RDistSP)==paste("P", rr, sep="_"))]  *  Cohort[yy,ss,2]
                # distribute surviving adults from last year
                # get appropriate group survival to get adults from last year
                grp <- StocksInfo$Surv_Grp[which(StocksInfo$StockID==Stocks[ss])]
                # for this to work survival needs to be order by column from young to old
                S <- as.numeric(Grp_Surv[which(Grp_Surv$Grp==as.character(grp) & Grp_Surv$Year==Years[yy]), 
                                         which(names(Grp_Surv) %in% paste("SR",(Ages-1),Ages, sep="_"))])
                
                ######################################################
                ## Distribute adults from last year's final period  ##
                ######################################################
                if(PoolType=="SP"){
                  for(aa in 2:NAges){
                    # for single pool need to sum all of age and stock across regions
                    N[ yy, 1, rr, ss, aa] <- RDistSP[which(RDistSP$StockID==Stocks[ss] & RDistSP$Age==Ages[aa] & RDistSP$Period==pp), 
                                                        which(names(RDistSP)==paste("P", rr, sep="_"))]  *  sum(N4[ yy-1, NP, , ss,aa-1 ]) * S[aa]
                  } # end age loop
                  
                } else if(PoolType=="MP"){ # end single pool if
                  # check that regions are in right order or else will mess up multiplication
                  Check <- RDistMP$FromReg[which(RDistMP$StockID==Stocks[ss] & RDistMP$Age==Ages[aa] & RDistMP$Period==pp) ]
                  if(F %in% (Check == c(1:NR))){ print("STOP -- MP DIST Data in Wrong Order")}
                  for(aa in 2:NAges){
                    N[ yy, 1, rr, ss, aa] <- sum( RDistMP[which(RDistMP$StockID==Stocks[ss] & RDistMP$Age==Ages[aa] & RDistMP$Period==pp),
                                                         which(names(RDistMP)==paste("P", rr, sep="_"))]  *  N4[ yy-1, NP, ,ss, aa-1] * S[aa] )
                  } # end age loop
 
                } # end multipool loop
            } # end yy>1
        } # end stock loop
      } # end region loop
    } else { # pp>1
      ##############
      ## Period 2+ #
      ##############
      # redistribute according to adult pop at end of last period
      for(rr in 1:NR){
        for(ss in 1:NS){
          for(aa in 1:NAges){
            # Single or multi pool?
            if(PoolType=="SP"){ 
              # Now get N according to distribution and Survival
              N[ yy, pp, rr, ss, aa] <- RDistSP[which(RDistSP$StockID==Stocks[ss] & RDistSP$Age==Ages[aa] & RDistSP$Period==pp), 
                                                    which(names(RDistSP)==paste("P", rr, sep="_"))]  *  sum(N4[ yy, pp-1, , ss, aa]) 
            } else{ #multipool dist sum over all regions could come from
              # check that regions are in right order or else will mess up multiplication
              Check <- RDistMP$FromReg[which(RDistMP$StockID==Stocks[ss] & RDistMP$Age==Ages[aa] & RDistMP$Period==pp) ]
              if(F %in% (Check == c(1:NR))){ print("STOP -- MP DIST Data in Wrong Order")}
              N[ yy, pp, rr, ss, aa] <- sum( RDistMP[which(RDistMP$StockID==Stocks[ss] & RDistMP$Age==Ages[aa] & RDistMP$Period==pp),
                                                    which(names(RDistMP)==paste("P", rr, sep="_"))]  *  N4[ yy, pp-1, ,ss, aa])
            } # end multipool dist
          } # end age loop
        } # end stock loop
      } # end region loop
      
    } # end pp > 1
    
    ##############################################
    ##  Setting AABM/ISBM TAC/Catch in Period 1 ##
    ##############################################
    if(pp==1){
      #based on period one N, set AABM TAC
        # If it's the first year AI=1
        if(yy==1){
          AI[yy,] <- 1
          # set up denominator now for future AI calcs
          # vector of length regions
          denom <- NULL
          numer <- NULL
          for(rr in 1:NR){
            # which is driver fishery for that region?
            DriverFishery <- FisheryInfo$FisheryID[which(FisheryInfo$Region==rr & FisheryInfo$IsALDriver=="Y")]
            # Vector of age proportions vulnerabel by age
            V <- as.numeric(PV[which(PV$FisheryID==DriverFishery & PV$Year==Years[1]), which(names(PV) %in% paste("PV", Ages, sep="_"))])
            # test if multiplying by BPER makes a difference 
            BPER_P1 <- as.matrix(BPER[which(BPER$Period==1 & BPER$FisheryID==DriverFishery),
                 which(names(BPER) %in% paste("ER_Age", Ages, sep=""))])
            # still need numerator for TAC equation, same as denominator
            denom[rr] <- numer[rr] <-  sum( (BPER_P1*N[1,1,rr,,]) %*% V )
            # denom[rr] <- sum( N[1,1,rr,,] %*% V )
          } # end region
        } else { # if year isn't==1
          for(rr in 1:NR){
            # which is driver fishery for that region?
            DriverFishery <- FisheryInfo$FisheryID[which(FisheryInfo$Region==rr & FisheryInfo$IsALDriver=="Y")]
            # pick out proportion for that year and driver fishery
            V <- as.numeric(PV[which(PV$FisheryID==DriverFishery ), which(names(PV) %in% paste("PV", Ages, sep="_"))])
            BPER_P1 <- as.matrix(BPER[which(BPER$Period==1 & BPER$FisheryID==DriverFishery),
                                        which(names(BPER) %in% paste("ER_Age", Ages, sep=""))])
            numer[rr] <- sum( (BPER_P1*N[yy,1,rr,,]) %*% V )
            AI[yy,rr] <- numer[rr]/denom[rr]
            #AI[yy,rr] <- sum( N[yy,1,rr,,] %*% V )/denom[rr]
          } # end region loop
        } #end year==1
      
        # Get TAC from function
        #TAC_R <- AABM.TAC(AI[yy,])
        TAC_R <- AABM.TAC.HRI(numer=numer, ai=AI[yy,])
        
        # Allocate regional TAC by period and gear distributions
        for(rr in 1:NR){
          for(gg in 1:NG){
            # get fishery ID from FisheryInfo
            FID <- FisheryInfo$FisheryID[which(FisheryInfo$Region==rr & FisheryInfo$Sector==Gears[gg])]
            # check that it is AABM
            if(FisheryInfo$Type[which(FisheryInfo$FisheryID==FID)]=="AABM"){
            # get most recent data year
            DatYrs <- FisheryDist$Year[which(FisheryDist$FisheryID==FID)]
            # first Allocate by period
            PeriodDist <-  as.numeric(FisheryDist[which(FisheryDist$FisheryID==FID & FisheryDist$Year==max(DatYrs[DatYrs<=Years[yy]])), 
                                                                    which(names(FisheryDist) %in% paste("Prop_Period_", 1:NP, sep=""))])
            # Now apply sector allocs
            DatYrs <-  SectorAlloc$Year[which(SectorAlloc$Region==rr)]
            GearDist <- as.numeric( SectorAlloc[which(SectorAlloc$Region==rr  & SectorAlloc$Year==max(DatYrs[DatYrs<=Years[yy]])), 
                                                         paste("Alloc", gg, sep="_")] )
            TAC[yy, , rr, gg] <- PeriodDist * GearDist * TAC_R[rr]
            # TAC will be empty for rr=2/g=3 and r=3/g=3 because these are ISBM
            } # end AABM if
          } # end gear loop
        } # end region loop
     
    } # end pp==1    
   
      # Now ISBM  -- Needs to be calculated each period!
      # which fisheries are ISBM?
      ff_ISBM <- FisheryInfo$FisheryID[which(FisheryInfo$Type=="ISBM")]
      for(ff in ff_ISBM){
        # extract fishery HRS
        # need to find closest lower year with dat
        DatYrs <-  HRS$Year[which(HRS$FisheryID==ff)]
        # choose largest of the years from same or previous years
        delta <- HRS$HRS[which(HRS$FisheryID==ff & HRS$Year==max(DatYrs[DatYrs<=Years[yy]]))] 
        # which region is the fishery in?
        region <- FisheryInfo$Region[which(FisheryInfo$FisheryID==ff)]
        gear <- which(Gears == FisheryInfo$Sector[which(FisheryInfo$FisheryID==ff)])
        for(ss in 1:NS){
            # harvest rate is BPER * HRS (delta), matrix NP*NA
            h[yy,,region,gear,ss, ] <- as.matrix(BPER[which(BPER$StockID==Stocks[ss] & BPER$FisheryID==ff),
                                    which(names(BPER) %in% paste("ER_Age", Ages, sep=""))]  *  delta)
            # Set TAC for period x Age based on period 1 abundance
            # Fishery Distributions across periods (vector lenth NP)
            # check which years have data
            #DatYrs <- FisheryDist$Year[which(FisheryDist$FisheryID==ff)]
            # Choose larger year of those that are smaller
            # Tao is vector length NP May 31/2017 don't need period allocations
            #Tao <- as.numeric(FisheryDist[which(FisheryDist$FisheryID==ff & FisheryDist$Year==max(DatYrs[DatYrs<=Years[yy]])), 
                              # which(names(FisheryDist) %in% paste("Prop_Period_", 1:NP, sep=""))])
            # Proportion vulnerable by age (Vector length NAges)
            # will need to alter if get more years of data
            V <- as.numeric(PV[which(PV$FisheryID==ff & PV$Year==Years[yy]), which(names(PV) %in% paste("PV", Ages, sep="_"))])
            # TAC here will be an NP X NAges matrix -- removed vulnerability
            # Tao(NPx1) X N(1xNAges) * h (NP X NAges)
            # this is where can add variability, no need to calculate TAC, since based on harvest rates
            #Catch[yy,,region, gear, ss, ] <- (cbind(Tao) %*% N[yy,pp,region,ss, ]) * h[yy,,region,gear,ss, ]
            Catch[yy,pp,region, gear, ss, ] <- h[yy,pp,region,gear,ss, ] * V * N[yy,pp,region,ss, ]  #* (1-PRel[yy,pp,region,gear,ss,])
            # values are tiny!
            # release and dropoff mortality
            # get release mortality, dropoff mortality
            # Get release mortality and dropoff mortality rates
            DatYrs <- TermRelMort$Year[which(TermRelMort$StockID==Stocks[ss])]
            eta <- TermRelMort$RelMort_subl[which(TermRelMort$StockID==Stocks[ss] & TermRelMort$Year==max(DatYrs[DatYrs<=Years[yy]]))]
            theta <- TermRelMort$Dropoff_Mort[which(TermRelMort$StockID==Stocks[ss] & TermRelMort$Year==max(DatYrs[DatYrs<=Years[yy]]))]
            Release[yy,,region, gear, ss, ] <- Catch[yy, , region, gear, ss, ] * PRel[yy, , region, gear, ss,  ] /(1-PRel[yy, , region, gear, ss, ])
            RelMort[yy,,region, gear, ss, ] <- Release[yy, , region, gear, ss, ] * eta
            # drop off mort applies to total catch
            Dropoff[yy,,region, gear, ss, ] <- (Catch[yy, , region, gear, ss, ] + Release[yy,,region, gear, ss, ]) * theta
          } # End Stock Loop
          
        } # end ISBM ff loop
     
   
    
    ##################################
    ####   Pre-terminal Fishing   ####
    ##################################
    
    ## Account for Catch according to ISBM catches and AABM TAC's and current stock, ages distributions
    
    # Turn AABM TAC into catch by stock, age
    # Loop over AABM fisheries
    ff_AABM <- FisheryInfo$FisheryID[which(FisheryInfo$Type=="AABM")]
    for(ff in ff_AABM){
      # which region is the fishery in?
      region <- FisheryInfo$Region[which(FisheryInfo$FisheryID==ff)]
      gear <- which(Gears == FisheryInfo$Sector[which(FisheryInfo$FisheryID==ff)])
      # extract release and dropoff mort rates
      #sublegal release mort, eta
      eta <- ReleaseMort[which(ReleaseMort$FisheryID==ff), "Rel_MortRate_Sub"]
      # Dropoff
      theta <-  ReleaseMort[which(ReleaseMort$FisheryID==ff), "Dropoff_MortRate"]
      for(ss in 1:NS){
        for(aa in 1:NAges){
          # get vulnerabilities
          #V <- as.numeric(PV[which(PV$FisheryID==ff & PV$Year==Years[yy]), which(names(PV) %in% paste("PV", Ages, sep="_"))])
          # could add variability here
          # Catch by stock age, depends on proportions of that stock and age
          # Catch[yy, pp, region, gear, ss, aa] <- V[aa] * N[ yy, pp, region, ss, aa]*(1-PRel[yy, pp, region, gear, ss, aa ]) / 
          #   sum( ((1-PRel[yy, pp, region, gear, ,  ]) * N[ yy, pp, region, , ]) %*% V) * TAC[yy, pp, region, gear]
          # Release[yy,pp,region, gear, ss, aa] <- Catch[yy, pp, region, gear, ss, aa] * PRel[yy, pp, region, gear, ss, aa ] /(1-PRel[yy, pp, region, gear, ss, aa ])
          # RelMort[yy,pp,region, gear, ss, aa] <- Release[yy, pp, region, gear, ss, aa] * eta
          # # drop off mort applies to total catch
          # Dropoff[yy,pp,region, gear, ss, aa] <- (Catch[yy, pp, region, gear, ss, aa] + Release[yy,pp,region, gear, ss, aa]) * theta
          
          # New Catch eqns 5/31/2017 to match DGM MOdel
          Catch[yy, pp, region, gear, ss, aa] <-  N[ yy, pp, region, ss, aa]*(1-PRel[yy, pp, region, gear, ss, aa ]) / 
               sum( ((1-PRel[yy, pp, region, gear, ,  ]) * N[ yy, pp, region, , ])) * TAC[yy, pp, region, gear]
          ## Catch * Encounter Rate -- prop(sublegal)/prop(Legal)
          Release[yy,pp,region, gear, ss, aa] <- Catch[yy, pp, region, gear, ss, aa] * PRel[yy, pp, region, gear, ss, aa ] /(1-PRel[yy, pp, region, gear, ss, aa ])
          RelMort[yy,pp,region, gear, ss, aa] <- Release[yy, pp, region, gear, ss, aa] * eta
          Dropoff[yy,pp,region, gear, ss, aa] <- (Catch[yy, pp, region, gear, ss, aa] + Release[yy,pp,region, gear, ss, aa]) * theta
        } # end Age loop
      } # end stock loop
    } # end fishery loop
    

  # Go through regions, stocks ages, and remove fishing, incidental mort.
    for(rr in 1:NR){
      for(ss in 1:NS){
        for(aa in 1:NAges){
          # Remove fishing over all gears
          N1[ yy, pp, rr, ss, aa] <- N[yy, pp, rr, ss, aa] - sum(Catch[yy, pp, rr, , ss, aa], na.rm=T)
          # Remove release Mort over all gears
          N2[ yy, pp, rr, ss, aa] <- N1[yy, pp, rr, ss, aa] #- sum(RelMort[yy, pp, rr, , ss, aa], na.rm=T)
          # Remove dropoff Mort
          N3[ yy, pp, rr, ss, aa] <- N2[yy, pp, rr, ss, aa] #- sum(Dropoff[yy, pp, rr, , ss, aa], na.rm=T)
        } # end Age loop
      } # end stock loop
    } # end region loop
    
    #################
    ## Maturation  ##
    #################
    # remove mature fish  -- need to change so that it is by region (required for multi-pool)
    # Check which stocks mature in this period, if any
    ss_mature <- which(Stocks == StocksInfo$StockID[which(StocksInfo$Mat_Period==pp)])
    if(length(ss_mature) >= 1){
      for(ss in ss_mature){
        # get maturation rates as vector, across ages (length(NAges))
        # don't care which region these mature fish came from for escapement
        MatRates <- sapply( Ages, function(x) Maturation$Maturation[which(Maturation$StockID==Stocks[ss] & Maturation$Age==x)] )
        # for multipool need to extract these mature fish from each region
        for(rr in 1:NR){
          Mature[yy, pp, rr, ss, ] <-  N3[yy, pp, rr ,ss, ] * MatRates
          # Now remove these from population
          N4[yy, pp, rr, ss, ] <- N3[yy, pp, rr, ss,  ] -  Mature[yy, pp, rr, ss, ]
        } # end region loop
      } # end stock loop
    } else {  # end mature if, else if no stocks maturing this period
      # population after maturity same as before
      N4 <- N3
    } # end else non mature this period
    
  } # End Period loop

  ############################################
  ##  Terminal Fishing, PSM, recruitment    ##
  ############################################
  # Cycle over stocks
  for(ss in 1:NS){
    # check that tehre is a fishery associated with this stock
    if(Stocks[ss] %in% TermHR$StockID){
      # extract HR for age classes as vector
      HR <- sapply ( Ages, function(x) TermHR[which(TermHR$StockID==Stocks[ss]), paste("HR_Age", x, sep="")])
      # Get HRS for most recent year
      ff <- TermHR$FisheryID[which(TermHR$StockID==Stocks[ss])]
      DatYrs <- HRS$Year[which(HRS$FisheryID==ff)]
      HRScaler <- HRS$HRS[which(HRS$FisheryID==ff & HRS$Year == max(DatYrs[DatYrs<=Years[yy]]))] 
      # apply HRS and HRScaler
      # need to sum N over all regions, for period 3
      Catch_Term[yy,ss, ] <- apply(Mature[yy, 3, ,ss,], 2, sum) * HR * HRScaler
      # Min size is so small, rarely have releses
      Release_Term[yy, ss, ] <- Catch_Term[yy, ss, ] * PRel_Term[yy, ss, ] /(1-PRel_Term[yy, ss,  ])
      # get release mortality, dropoff mortality
      DatYrs <- TermRelMort$Year[which(TermRelMort$StockID==Stocks[ss])]
      eta <- TermRelMort$RelMort_subl[which(TermRelMort$StockID==Stocks[ss] & TermRelMort$Year==max(DatYrs[DatYrs<=Years[yy]]))]
      theta <- TermRelMort$Dropoff_Mort[which(TermRelMort$StockID==Stocks[ss] & TermRelMort$Year==max(DatYrs[DatYrs<=Years[yy]]))]
      #calculate release mortality
      RelMort_Term[yy, ss, ] <- Release_Term[yy, ss, ] * eta
      # drop off mort applies to total catch
      Dropoff_Term[yy, ss,] <- (Catch_Term[yy,  ss, ] + Release_Term[yy, ss, ]) * theta
    # just in case there isn't a terminal fishery for that given stocks
    } else { # else if not terminal fishery
      Catch_Term[yy, ss, ] <- 0
    } # end terminal fishery if
    Escape[yy, ss, ] <- apply(Mature[yy, 3, ,ss,], 2, sum) - Catch_Term[yy,ss, ] -  RelMort_Term[yy, ss, ] - Dropoff_Term[yy, ss,]
    # apply PSM, in future will likely be by stock, so put inside stock loop
    PS_Surv <- PreSpawnSurv$Surv_Rate[which(PreSpawnSurv$Year==Years[yy])]
    Spawners[yy, ss] <- sum(Escape[yy, ss, ] * PS_Surv)
    #Recruitment
    # Is it a hatchery or natural stock?
    if(as.character(StocksInfo$Hatchery[which(StocksInfo$StockID==Stocks[ss])])=="N"){
      # Get Ricker parameters for this stock
      a <- StocksInfo$Ricker_A[which(StocksInfo$StockID==Stocks[ss])]
      B <- StocksInfo$Ricker_B[which(StocksInfo$StockID==Stocks[ss])]
      RickerRecruits <-  Spawners[yy, ss] * exp( a*(1- Spawners[yy, ss]/B))
      Recruits[yy, ss] <- ifelse(RickerRecruits > B , B, RickerRecruits)
      # translate adult recruits to age 1's for next year
      # get age 1 to adult survival
      # need vectors of S and M for this stock, year
      grp <- StocksInfo$Surv_Grp[which(StocksInfo$StockID==Stocks[ss])]
      S <- as.numeric(Grp_Surv[which(Grp_Surv$Grp==as.character(grp) & Grp_Surv$Year==Years[yy]), paste("SR", 0:(MaxAge-1), 1:MaxAge, sep="_")])
      # for now Maturation rates are only for one year
      # first entry is NA becasuse age 1's dont' mature at all
      M <- c(NA,  Maturation$Maturation[which(Maturation$StockID==Stocks[ss] & Maturation$Age==2:MaxAge)])
      Age1Surv <- GetAge1Surv(S=S, M=M)
      if (yy < NY){
        Cohort[yy+1, ss, 1] <- Recruits[yy, ss] / Age1Surv
      } # end if not final year
    } else { # else if hatchery stock
      # for now just do what DGM did, multiply by alpha from "associated stock"
      a <- StocksInfo$Ricker_A[which(StocksInfo$Assoc_Stock==as.character(Stocks[ss]))]
      # for hatchery these go directly into age 1's for next year
      if (yy < NY){
        Cohort[yy+1, ss, 1] <-  Spawners[yy, ss] * exp(a)
      } # end if not final year
    } # end hatchery if
    
   } # end stock loop
    

} # End Year Cycle

########################
##   End Simulation   ##
########################
 
# Export Values to csv files
# Abundance

# Cohort abundance (beginnng of period 1)
Cohort_Tab <- data.frame(Year=numeric(), Stock=character(), Age=numeric(), Abund=numeric() )
for(yy in 1:NY){
  for(ss in 1:NS){
    for(aa in 1:MaxAge)
      if(aa <= 2){
        NewRow <- data.frame(Year==Years[yy], Stock=Stocks[ss], Age=aa, Abund=Cohort[yy,ss,aa])
      } else {
        NewRow <- data.frame(Year==Years[yy], Stock=Stocks[ss], Age=aa, Abund=Cohort[yy,ss,aa])
      }
  }
}
