#### REMOTE PASSIVE INTEGRATED TRANSPONDER (PIT) TAG READER = RePTaR
#### Brown treesnakes were PIT tagged and released into experimental trial arena to test the scanning success of this RePTaR reader prototype
#### Scanning success within a 2-in area around RePTaR antennas and distance from antenna gathered from videos of experimental trials
#### Analysis done via logistic regression in a Bayesian framework with subsequent model comparison via WAIC
#### Code written by Staci Amburgey
#### Involving data collected from August 1-August 30 2021

## The purpose of this script is to:
### 1) load and format data 
### 2) fit models to data
### 3) calculate WAIC for each model


rm(list=ls())

library(tidyverse); library(lubridate); library(dplyr); library(jagsUI); library(ggplot2)



#### PART ONE: LOAD AND FORMAT BEHAVIORAL TRIAL DATA.----

## Read in detailed trial info and format for analysis
alldat <- read_csv("Data/RePTaR_trials_AllData.csv",show_col_types = FALSE) %>%
  mutate(PITTAG = str_match(PITTAG, "\"(.*?)\"")[,2]) %>% # remove quotations used to maintain full 15-digit ID (instead of scientific notation issues)
  mutate(PITTAG = as.character(PITTAG)) %>%
  drop_na(TimeSide) %>%                              # drop instances where snake was too far away (> 2in)
  filter(TimeSide != 1) %>%                          # drop one instance where incorrect time entered (outside of range retained for analysis anyway)
  mutate(TimeSide2 = lubridate::hms(TimeSide)) %>%   # specify as time
  mutate(Date = dmy(Date)) %>%                       # specify as date
  filter((TimeSide2 >= c("19H 00M 00S") & TimeSide2 <= c("21H 00M 00S")) |       # limit to the subsetted trial times processed (7-9pm, 3-5am)
           (TimeSide2 >= c("03H 00M 00S") & TimeSide2 <= c("05H 00M 00S"))) %>%
  mutate(Antenna = `Side`) %>%                                                     # rename for ease
  mutate(Antenna2 = if_else(Antenna == 1 & ARENASIDE == 1, 1,                      # rename antennas to be 1 & 2 (for unit 1) and 3 & 4 (for unit 2)
          if_else(Antenna == 2 & ARENASIDE == 1, 2,
                  if_else(Antenna == 1 & ARENASIDE == 2, 3,
                          if_else(Antenna == 2 & ARENASIDE == 2, 4, -9999)))))


## Get expanded PITTAG IDs (full 15-digit ID) and reduced IDs (8-digit ID) for use in combining trait and scanning dataframes
tags <- as.data.frame(unique(alldat$PITTAG)); colnames(tags) <- "PITTAG"
tags$RFID <- substr(tags$PITTAG, 8, 15)

## Read in individual info and format for analysis
siz <- read_csv("Data/SnakeTraits.csv",show_col_types = FALSE) %>%
  select(DATETEST,TRIAL,PITTAG,SVL,TL,TAILBREAK,SEX,WEIGHT,BULGE,BATCHMARK,ARENASIDE) %>%  # subset to columns of interest
  mutate(PITTAG = str_match(PITTAG, "\"(.*?)\"")[,2]) %>% # remove quotations used to maintain full 15-digit ID (instead of scientific notation issues)
  mutate_at(vars(PITTAG), factor) %>%
  rename(Date = DATETEST) %>%                             # denote this day of trial as the date of interest
  mutate(Date = dmy(Date)) %>%                            # specify as date
  filter(PITTAG != 982091065198473) %>%                   # remove this tag as snake escaped when first trialed before scanning and there were no scans ever for this tag
  filter(!(PITTAG == 982091065198381 & TRIAL == 10)) %>%  # remove individual who lost tag before trial (but there other, earlier trials with this tag)
  mutate(SEX = ifelse(SEX == "F", 1, 0)) %>%              # specify sex as 1 (female), 0 (male)
  mutate(TagLoc = ifelse(TRIAL < 15, 1, 0))               # trials 1-14 had tags in neck (1) while 15 and 16 had tags in posterior (0)
siz$ID <- c(1:nrow(siz))                                  # denote individual ID

## Bind data to get scanning success by snake characteristics
alldat2 <- subset(alldat, !is.na(ApproxDist))  # remove instance where distance could not be approximated (e.g., snake blocked from view)
dscan <- inner_join(alldat2[,c("Date","PITTAG","Read","TimeSide","ApproxDist","Antenna2")],siz[,c("Date","TRIAL","PITTAG","SVL","SEX","WEIGHT","TagLoc","ID")], by=c("Date","PITTAG"))
colnames(dscan) <- c("Date","PITTAG","Scan","TimeSide","Dist","Antenna","TRIAL","SVL","SEX","WEIGHT","TagLoc","ID")

## Limit analysis to scans within 2 inches as failure to scan was not recorded beyond this distance
scan <- subset(dscan, Dist <= 2)

## Sort dataset by trial and snake ID
scan <- scan[order(scan$TRIAL,scan$ID),]

## Adjust ID of snakes to account for missing (unscanned during the 7-9pm or 3-5am times) snake 61 - allows model to still loop over snake ID
scan$ID <- ifelse(scan$ID < 61, scan$ID, scan$ID - 1)


#### CHECK FOR CORRELATION IN COVARIATES (r > 0.6 considered to be highly correlated).----
# Check the subset of variables used in models
cor(scan[,c(5:12)])
corrplot::corrplot(cor(scan[,c(5:12)]))
#weight and SVL, trial and ID, tag location and trial, and tag location and ID are all strongly correlated
#Weight and SVL are biologically related and should not be included in the same models
#However, tag location, trial, and ID are correlated due to the design of the study and can be included together


#### PREPARE MODEL INFO.----
N <- length(unique(scan$ID))                        # number of individuals
read <- as.vector(unlist(scan[,c("Scan")]))         # trial outcomes
ID <- as.vector(unlist(scan[,c("ID")]))             # ID of individuals
sex <- as.vector(unlist(scan[,c("SEX")])) + 1       # sex of individuals
trial <- as.vector(unlist(scan[,c("TRIAL")]))       # trial for each individual
size <- as.vector(unlist(scan[,c("SVL")]))          # SVL for each individual; used for deciding if retaining SVL or weight but not used further in this script file as weight better supported
dist <- as.vector(as.character(unlist(scan[,c("Dist")])))         # distance of individual from antenna
dist <- ifelse(dist == 0, 1,                                      # convert each distance to a whole number for modeling and plotting ease
               ifelse(dist == 0.5, 2,
                      ifelse(dist == 1, 3,
                             ifelse(dist == 1.5, 4, 5))))
weight <- as.vector(unlist(scan[,c("WEIGHT")]))     # weight of each individual
ant <- as.vector(unlist(scan[,c("Antenna")]))       # which antenna scanned
loc <- as.vector(unlist(scan[,c("TagLoc")])) + 1    # location of tag in snake




#### PART TWO: LOGISTIC REGRESSION IN JAGS.----

# MCMC settings 
nc = 3
ni = 1000
nb = 100
nthin = 1

inits <- function() {
  list()
}

########################################################

## MODEL ONE: Scanning success by global model (sex + dist + loc + weight + weight*loc) with random effect of ID

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  for(c in 1:C){
    b1[c] ~ dnorm(0,5)  # Female, male
    b3[c] ~ dnorm(0,5)  # Anterior, posterior
    b5[c] ~ dnorm(0,5)  # Anterior, posterior * Weight
  }
  for(d in 1:D){
    b2[d] ~ dnorm(0,5)  # Distance (0,0.5,1.0,1.5,2.0)
  }
  b4 ~ dnorm(0,5)       # Weight
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1[Sex[i]] + b2[Dist[i]] + b3[Loc[i]] + b4*Weight[i] + b5[Loc[i]]*Weight[i] + eta[ID[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCatCon_GlobalWeight.txt")

########################################################

## MODEL: Scan by Global + Weight.----
data <- list(read=read, ID=ID, Sex=sex, Dist=as.integer(dist), Loc=loc, Weight=weight, N=N, C=length(unique(sex)), D=length(unique(dist)))
modnam <- c("GlobalWeight")

## JAGS model details.----
parameters<-c('b0','b1','b2','b3','b4','b5','loglike','eta')

out <- jagsUI::jags(model.file ="Models/LogitCatCon_GlobalWeight.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL TWO: Scanning success by sex, loc, weight, loc*weight with random effect of ID

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  for(c in 1:C){
    b1[c] ~ dnorm(0,5)  # Female, male
    b2[c] ~ dnorm(0,5)  # Anterior, posterior
    b4[c] ~ dnorm(0,5)  # Anterior, posterior * Weight
  }
  b3 ~ dnorm(0,5)       # Weight
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1[Sex[i]] + b2[Loc[i]] + b3*Weight[i] + b4[Loc[i]]*Weight[i] + eta[ID[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCatCon_SexLocWeightInt.txt")

########################################################

## MODEL: Scan by Sex, Loc, Weight, Loc*Weight.----
data <- list(read=read, ID=ID, Sex=sex, Loc=loc, Weight=weight, N=N, C=length(unique(sex)))
modnam <- c("SexLocWeightInt")

## JAGS model details.----
parameters<-c('b0','b1','b2','b3','b4','loglike','eta')

out <- jagsUI::jags(model.file ="Models/LogitCatCon_SexLocWeightInt.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL THREE: Scanning success by sex, dist, loc, weight with random effect of ID

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  for(c in 1:C){
    b1[c] ~ dnorm(0,5)  # Female, male
    b3[c] ~ dnorm(0,5)  # Anterior, posterior
  }
  for(d in 1:D){
    b2[d] ~ dnorm(0,5)  # Distance (0,0.5,1.0,1.5,2.0)
  }
  b4 ~ dnorm(0,5)       # Weight
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1[Sex[i]] + b2[Dist[i]] + b3[Loc[i]] + b4*Weight[i] + eta[ID[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCatCon_SexDistLocWeight.txt")

########################################################

## MODEL: Scan by Sex, Dist, Loc, Weight.----
data <- list(read=read, ID=ID, Sex=sex, Dist=as.integer(dist), Loc=loc, Weight=weight, N=N, C=length(unique(sex)), D=length(unique(dist)))
modnam <- c("SexDistLocWeight")

## JAGS model details.----
parameters<-c('b0','b1','b2','b3','b4','loglike','eta')

out <- jagsUI::jags(model.file ="Models/LogitCatCon_SexDistLocWeight.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL FOUR: Scanning success by sex, dist, loc with random effect of ID

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  for(c in 1:C){
    b1[c] ~ dnorm(0,5)  # Female, male
    b3[c] ~ dnorm(0,5)  # Anterior, posterior
  }
  for(d in 1:D){
    b2[d] ~ dnorm(0,5)  # Distance (0,0.5,1.0,1.5,2.0)
  }
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1[Sex[i]] + b2[Dist[i]] + b3[Loc[i]] + eta[ID[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCatCon_SexDistLoc.txt")

########################################################

## MODEL: Scan by Sex, Dist, Loc.----
data <- list(read=read, ID=ID, Sex=sex, Dist=as.integer(dist), Loc=loc, N=N, C=length(unique(sex)), D=length(unique(dist)))
modnam <- c("SexDistLoc")

## JAGS model details.----
parameters<-c('b0','b1','b2','b3','loglike','eta')

out <- jagsUI::jags(model.file ="Models/LogitCatCon_SexDistLoc.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL FIVE: Scanning success by sex, dist, weight with random effect of ID

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  for(c in 1:C){
    b1[c] ~ dnorm(0,5)  # Female, male
  }
  for(d in 1:D){
    b2[d] ~ dnorm(0,5)  # Distance (0,0.5,1.0,1.5,2.0)
  }
  b3 ~ dnorm(0,5)       # Weight
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1[Sex[i]] + b2[Dist[i]] + b3*Weight[i] + eta[ID[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCatCon_SexDistWeight.txt")

########################################################

## MODEL: Scan by Sex, Dist, Weight.----
data <- list(read=read, ID=ID, Sex=sex, Dist=as.integer(dist), Weight=weight, N=N, C=length(unique(sex)), D=length(unique(dist)))
modnam <- c("SexDistWeight")

## JAGS model details.----
parameters<-c('b0','b1','b2','b3','loglike','eta')

out <- jagsUI::jags(model.file ="Models/LogitCatCon_SexDistWeight.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL SIX: Scanning success by dist, loc, weight with random effect of ID

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  for(c in 1:C){
    b1[c] ~ dnorm(0,5)  # Anterior, posterior
  }
  for(d in 1:D){
    b2[d] ~ dnorm(0,5)  # Distance (0,0.5,1.0,1.5,2.0)
  }
  b3 ~ dnorm(0,5)       # Weight
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1[Loc[i]] + b2[Dist[i]] + b3*Weight[i] + eta[ID[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCatCon_DistLocWeight.txt")

########################################################

## MODEL: Scan by Dist, Loc, Weight.----
data <- list(read=read, ID=ID, Loc=loc, Dist=as.integer(dist), Weight=weight, N=N, C=length(unique(sex)), D=length(unique(dist)))
modnam <- c("DistLocWeight")

## JAGS model details.----
parameters<-c('b0','b1','b2','b3','loglike','eta')

out <- jagsUI::jags(model.file ="Models/LogitCatCon_DistLocWeight.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL SEVEN: Scanning success by loc, weight, loc*weight with random effect of ID 

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  for(c in 1:C){
    b1[c] ~ dnorm(0,5)     # Anterior, posterior
    b3[c] ~ dnorm(0,5)     # Anterior, posterior * Weight
  }
  b2 ~ dnorm(0,5)          # Weight
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1[Loc[i]] + b2*Size[i] + b3[Loc[i]]*Weight[i] + eta[ID[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCatCon_LocWeightInt.txt")

########################################################

## MODEL: Scan by Loc, Weight, and Loc*Weight.----
data <- list(read=read, ID=ID, Loc=loc, Weight=weight, N=N, C=length(unique(loc)))
modnam <- c("LocWeightInt")

## JAGS model details.----
parameters<-c('b0','b1','b2',"b3","loglike",'eta')

out <- jagsUI::jags(model.file ="Models/LogitCatCon_LocWeightInt.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL EIGHT, NINE, TEN: Scanning success by two categorical covariate (e.g., sex and dist, sex and loc, dist and loc) with random effect of ID 

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  for(c in 1:C){
    b1[c] ~ dnorm(0,5)     # e.g., Female, male
  }
  for(d in 1:D){
    b2[d] ~ dnorm(0,5)  # e.g., Distance (0,0.5,1.0,1.5,2.0)
  }
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1[cov1[i]] + b2[cov2[i]] + eta[ID[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCatCat.txt")

########################################################

## MODEL: Scan by Sex and Dist.----
data <- list(read=read, ID=ID, cov1=sex, cov2=as.integer(dist), N=N, C=length(unique(loc)), D=length(unique(dist)))
modnam <- c("SexDist")

## JAGS model details.----
parameters<-c('b0','b1','b2',"loglike",'eta')

out <- jagsUI::jags(model.file ="Models/LogitCatCat.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))


## MODEL: Scan by Sex and Loc.----
data <- list(read=read, ID=ID, cov1=sex, cov2=loc, N=N, C=length(unique(sex)), D=length(unique(loc)))
modnam <- c("SexLoc")

## JAGS model details.----
parameters<-c('b0','b1','b2',"loglike",'eta')

out <- jagsUI::jags(model.file ="Models/LogitCatCat.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))


## MODEL: Scan by Dist and Loc.----
data <- list(read=read, ID=ID, cov1=as.integer(dist), cov2=loc, N=N, C=length(unique(dist)), D=length(unique(loc)))
modnam <- c("DistLoc")

## JAGS model details.----
parameters<-c('b0','b1','b2',"loglike",'eta')

out <- jagsUI::jags(model.file ="Models/LogitCatCat.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL ELEVEN, TWELVE, THIRTEEN: Scanning success by categorical covariate (e.g., sex or dist or loc) and continuous covariate (e.g., weight) with random effect of ID 

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  for(c in 1:C){
    b1[c] ~ dnorm(0,5)     # e.g., Female, male
  }
  b2 ~ dnorm(0,5)          # e.g., Weight
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1[cov1[i]] + b2*cov2[i] + eta[ID[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCatCon.txt")

########################################################

## MODEL: Scan by Sex and Weight.----
data <- list(read=read, ID=ID, cov1=sex, cov2=weight, N=N, C=length(unique(sex)))
modnam <- c("SexWeight")

## JAGS model details.----
parameters<-c('b0','b1','b2',"loglike",'eta')

out <- jagsUI::jags(model.file ="Models/LogitCatCon.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))


## MODEL: Scan by Dist and Weight.----
data <- list(read=read, ID=ID, cov1=as.integer(dist), cov2=weight, N=N, C=length(unique(dist)))
modnam <- c("DistWeight")

## JAGS model details.----
parameters<-c('b0','b1','b2',"loglike",'eta')

out <- jagsUI::jags(model.file ="Models/LogitCatCon.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))


## MODEL: Scan by Loc and Weight.----
data <- list(read=read, ID=ID, cov1=loc, cov2=weight, N=N, C=length(unique(loc)))
modnam <- c("LocWeight")

## JAGS model details.----
parameters<-c('b0','b1','b2',"loglike",'eta')

out <- jagsUI::jags(model.file ="Models/LogitCatCon.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL FOURTEEN, FIFTEEN, SIXTEEN: Scanning success by categorical covariate (e.g., sex) with random effect of ID 

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  for(c in 1:C){
    b1[c] ~ dnorm(0,5)      # e.g., Female, male
  }
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1[cov[i]] + eta[ID[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCat.txt")

########################################################

## MODEL: Scan by Sex.----
data <- list(read=read, ID=ID, cov=sex, N=N, C=length(unique(sex)))
modnam <- c("Sex")

## JAGS model details.----
parameters<-c('b0','b1','loglike','eta')

out <- jagsUI::jags(model.file ="Models/LogitCat.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))


## MODEL: Scan by Distance.----
data <- list(read=read, ID=ID, cov=as.integer(dist), N=N, C=length(unique(dist)))
modnam <- c("Dist")

## JAGS model details.----
parameters<-c('b0','b1','loglike','eta')

out <- jagsUI::jags(model.file ="Models/LogitCat.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))


## MODEL: Scan by Location.----
data <- list(read=read, ID=ID, cov=loc, N=N, C=length(unique(loc)))
modnam <- c("Loc")

## JAGS model details.----
parameters<-c('b0','b1','loglike','eta')

out <- jagsUI::jags(model.file ="Models/LogitCat.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL SEVENTEEN: Scanning success by continuous covariate (e.g., weight) with random effect of ID 

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  b1 ~ dnorm(0,5)     # Weight
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1*cov[i] + eta[ID[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCont.txt")

########################################################

## MODEL: Scan Weight.----
data <- list(read=read, ID=ID, cov=weight, N=N)
modnam <- c("Weight")

## JAGS model details.----
parameters<-c('b0','b1',"loglike","eta")

out <- jagsUI::jags(model.file ="Models/LogitCont.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL EIGHTEEN: Scanning success null model (just individual random effect)

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + eta[ID[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitNull.txt")

########################################################

## MODEL: Scan by Null.----
data <- list(read=read, ID=ID, N=N)
modnam <- c("Null")

## JAGS model details.----
parameters<-c('b0',"loglike","eta")

out <- jagsUI::jags(model.file ="Models/LogitNull.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))




## PART THREE: MODEL COMPARISON.----
# Watanabe-Akaike Information Criterion (WAIC)

## Modified from https://github.com/heathergaya/JAGS-NIMBLE-Tutorials/blob/master/Known_Fate/Known_Fate_Models.Rmd

#create a list of the files from your target directory
file_list <- list.files(path="Results")

calc.waic <- function(x){
  #find the output that relates to loglike
  like <- x$sims.list$loglike
  fbar <- colMeans(exp(like)) #mean likelihood 
  Pw <- sum(apply(like,2,var)) #mean variance in log-likelihood 
  WAIC<- -2*sum(log(fbar))+2*Pw
  return(WAIC)
}

#create matrix to hold WAIC
waicALL2 <- matrix(NA, nrow=length(file_list), ncol=1, dimnames = list(sub("\\..*","",file_list),c("waic")))

#calulcate WAIC
for(i in 1:length(file_list)){
  load(file=paste("Results/",file_list[i],sep=""))   # will throw an error when it runs through all files and hits any subfolders in here - just make sure all results are in this central folder
  waicALL2[i,1] <- calc.waic(out)
}
#format WAIC matrix 
waicALL3 <- as.data.frame(waicALL2[order(waicALL2),]); colnames(waicALL3) <- c("waic")
waicALL3$deltawaic <- waicALL3$waic - min(waicALL3$waic, na.rm = TRUE) # calculate deltaWAIC
waicALL3$deltawaic <- round(waicALL3$deltawaic, 2)

