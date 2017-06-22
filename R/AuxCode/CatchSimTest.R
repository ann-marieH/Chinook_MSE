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

# Note: size limits can be adjusted for AABM fisheries by changing the minSizeLim valu# Calculate proprtion of each stock-age cohort that are legal
propLegal<- A3 <-  1-pnorm(minSizeLim,mean=muSize,sd=sigSize)
#propLegal <- 1-(1 / (1+exp(-1.7*((minSizeLim- muSize) / sigSize)))) # approximation of cummulative norm. probablity (may be faster for simulation???)

# Calculate the proportion of each stock that are sublegal (i.e., above min vulnerability to gear but below size limit)
# This is just 1-propNotLegal
propNotLegal<-pnorm(minSizeLim,mean=muSize,sd=sigSize)
#propNotLegal<- 1 / (1+exp(-1.7*((minSizeLim- muSize) / sigSize))) 


#Proportion too small to be caught in gear
propNotVuln<- A1<- pnorm(minVuln,mean=muSize,sd=sigSize)
#propNotVuln<- 1 / (1+exp(-1.7*((minVuln- muSize) / sigSize)))
propVuln <- 1-propNotVuln

#proportion of fish vulnerable, but also sublegal
propSubLegal <- A2 <-  propNotLegal - propNotVuln



# Do this 1000 times and take averages
CatchList <- list()
RelMortList <- list()
DropoffList <- list()
for( i in 1:1000) {
  ## Get fish -- guess age stock -- did it dropoff? -- get size -- thrown back -- die or not -- or catch
  Cohort<-matrix(c(2000,1000,500,200,6000,4000,2000,1000,950,500,250,150,
                   3000,1000,500,300,4000,3000,2000,1000),5,4,byrow=T)
  Landings <- 0
  Catch <- matrix(rep(0,20),nrow=5, ncol=4)
  RelMortMat <- matrix(rep(0,20),nrow=5, ncol=4)
  Dropoff <- matrix(rep(0,20),nrow=5, ncol=4)
  #index of fish encounters
  ind <- 1
  # Give each age, stock an id
  id <- matrix(1:20, nrow=5, ncol=4)
  while(Landings < TAC){
    A1<-   pnorm(minVuln,mean=muSize,sd=sigSize)
    A3 <-  1-pnorm(minSizeLim,mean=muSize,sd=sigSize)
    A2 <- 1-A1-A3
      
    # probabibility catching fish of each age stock
    p <- Cohort * ((A2+A3)/(A1+A2+A3)) / sum(Cohort * ((A2+A3)/(A1+A2+A3)))
    id_selected <- sample(x = id, 1, prob = p) 
    # Is it dropoff?
    Drop <- sample(x=c("Dead", "Alive"), 1, prob=c(dRate, 1-dRate))
    if(Drop=="Dead"){
      Cohort[id_selected] <- Cohort[id_selected]-1
      Dropoff[id_selected] <- Dropoff[id_selected]+1
      ind <- ind+1
    } else {
      # if fish brough on board alive, is it released or kept?
      # legal or sublegal?
      P_legal <- A3[id_selected]/(A2[id_selected]+A3[id_selected])
      Size <- sample(x=c("L", "SL"), 1, prob=c(P_legal, 1-P_legal ))
      if(Size=="L"){
        Catch[id_selected] <- Catch[id_selected]+1
        Landings <- Landings+1
      } else { # if sublegal
        # it is released, does it die?
        RelMort <- sample(x=c("Dead", "Alive"), 1, prob=c(relMortRate, (1-relMortRate)))
        if(RelMort == "Dead"){
          # if dead, remove from population
          Cohort[id_selected] <- Cohort[id_selected]-1
          RelMortMat[id_selected] <-  RelMortMat[id_selected] + 1
          ind <- ind+1
        } else { # else released alive, nothing happens
          ind <- ind+1
        } # end release else
      } # end sublegal else
    } # end not dropoff else
  
  } # end while loop
  CatchList[[i]] <- Catch
  RelMortList[[i]] <- RelMortMat
  DropoffList[[i]] <- Dropoff
} # end sim loop

AvgCatch <- apply(simplify2array(CatchList), 1:2, FUN=mean)
AvgRelMort <- apply(simplify2array(RelMortList), 1:2, FUN=mean)
AvgDropoff <- apply(simplify2array(DropoffList), 1:2, FUN=mean)
