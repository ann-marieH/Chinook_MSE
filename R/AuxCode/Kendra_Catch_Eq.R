# Chinook Fishery Catch Equations
# Kendra Holt, May 16


# Specify variables & parameter values =====================================================
# Create cohort data for 5 stock of ages 2 to 5 (rows=stocks, columns = ages; Cohort[s,a])
Cohort<-matrix(c(2000,1000,500,200,6000,4000,2000,1000,950,500,250,150,3000,1000,500,300,4000,3000,2000,1000),5,4,byrow=T)
# Create mean size and sd size by stock and age (rows=stocks, columns = ages; muSize[s,a])
muSize<-matrix(c(557,728,829,889,460,660,797,890,476,635,737,802,476,635,737,802,562,781,886,936),5,4,byrow=T)
sigSize<-matrix(c(63,82,94,101,rep(86,4),rep(89,4),67,74,77,78,rep(86,4)),5,4,byrow=T)
# Set min size vulnerable to gear and minimum size limit for retention

dRate<-0.05 # dropoff rate
relMortRate<-0.3

# AABM Fisheries =================================================================

# Specify fishery-specific TAC
TAC <- 500
minVuln<-500 # assume knife-edge recruitment to fishing gear at 350 mm
minSizeLim<-600

# Note: size limits can be adjusted for AABM fisheries by changing the minSizeLim value


# Calculate proprtion of each stock-age cohort that are legal
propLegal<- 1-pnorm(minSizeLim,mean=muSize,sd=sigSize)
#propLegal <- 1-(1 / (1+exp(-1.7*((minSizeLim- muSize) / sigSize)))) # approximation of cummulative norm. probablity (may be faster for simulation???)

# Calculate the proportion of each stock that are sublegal (i.e., above min vulnerability to gear but below size limit)
# This is just 1-propNotLegal
propNotLegal<-pnorm(minSizeLim,mean=muSize,sd=sigSize)
#propNotLegal<- 1 / (1+exp(-1.7*((minSizeLim- muSize) / sigSize))) 


#Proportion too small to be caught in gear
propNotVuln<-pnorm(minVuln,mean=muSize,sd=sigSize)
#propNotVuln<- 1 / (1+exp(-1.7*((minVuln- muSize) / sigSize)))
propVuln <- 1-propNotVuln

#proportion of fish vulnerable, but also sublegal
propSubLegal<- propNotLegal - propNotVuln

A1<- propNotVuln
A2 <- propSubLegal
A3 <- propLegal
  
# Catch brough on Board from each pop will be proportional to proportion that are vulnerable

COB <- TAC * ((1+A2/(A2+A3)) * ((A2+A3) * Cohort) / sum((A2+A3) * Cohort))
Landed <- COB * (A3/(A2+A3))


# Calculate landed catch (i.e., allocate TAC to stocks and ages)
Catch <- TAC * ((propLegal*Cohort) / sum(propLegal*Cohort))


# Calculate sub-legal mortality
   # (for each age- and stock-, the ectRate is an expansion factor to get from Catch to Sub-legal releases based on assumtpion 
    # of equal vulnerability of all fish above minVuln to gear)
ectRate<-propSubLegal / propLegal  
Releases_sublegal<-ectRate * Catch
RelMortalities<-Releases_sublegal*relMortRate

COB <- TAC * (1+ectRate) * ((A2+A3) * Cohort) / sum((A2+A3) * Cohort)
rel_Alt <- COB - Catch


# calculate drop-off mortality
Dropoff_legal <- Catch * dRate
Dropoff_sublegal<-Releases_sublegal * dRate

(Catch+Releases_sublegal) * dRate

# Calculate total incidental mortality (release + drop-off for legal and sublegal)
IncMort<- RelMortalities + Dropoff_legal + Dropoff_sublegal

browser()

# ISBM Fisheries ===============================================================

# Specify fishery-specific base period exploitation rate, by stock and age (assuming total mortality, incl. release & dropoff)
BPER<-matrix(c(0.00024,0.00162,0.0007,0.00074,0,0.00981,0.02111,0.02126,0,0.00014,0.00147,0,0,0.01923,0.02072,0.02645,0.0001,0.014,0.031,0.035),5,4,byrow=T)

BP_minSizeLim<-600
scalar<-0.8

# Note: size limits for IIBM fisheries must be kept at base period size limit at this point; 
    # I have not yet figured out a way to dientangle the vulnerability and size limit effects within BPER

# Calculate legal mortalities based on scalar relative to base period
LegalMortalities<-Cohort*propLegal*BPER*scalar

# in DGm they use PV instead of propLegal, so don't need propLegal to try and match it
# this will be called catch for now to match DGM, but we think it is double accounting for dropoff 

# Calculate portion of legal mortalities that are catch
Catch<-LegalMortalities/(1+dRate)

#Calculate portion of legal mortalities that are drop-off
Dropoff_legal<-LegalMortalities - Catch

# Calculate sublegal mortalities (release and drop-off combined)
SubLegalMortalities<-Cohort*propNotLegal*BPER*scalar

# Calculate total incidental mortality (release + drop-off for legal and sublegal)
IncMort<-SubLegalMortalities+Dropoff_legal

