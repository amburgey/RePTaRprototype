#### REMOTE PASSIVE INTEGRATED TRANSPONDER (PIT) TAG READER = RePTaR
#### Brown treesnakes were PIT tagged and released into experimental trial arena to test the scanning success of this RePTaR reader prototype
#### March 2022 - code written by Staci Amburgey
#### Involving data collected from August 1-August 30 2021

rm(list=ls())

library(tidyverse);library(dplyr);library(lubridate);library(jagsUI);library(ggplot2);library(lme4)


#### LOAD AND PREPARE DATA.----

### Read in detailed trial info
alldat <- read_csv("Data/RePTaR_trials_Aug2021_AllData_2.csv",show_col_types = FALSE) %>%
  mutate(PITTAG = as.character(PITTAG)) %>%
  drop_na(TimeSide) %>%
  filter(TimeSide != 1) %>%  #drop one instance where incorrect time entered (outside of range retained for analysis anyway)
  mutate(TimeSide2 = lubridate::hms(TimeSide)) %>%
  mutate(Date = dmy(Date)) %>%
  filter((TimeSide2 >= c("19H 00M 00S") & TimeSide2 <= c("21H 00M 00S")) |
           (TimeSide2 >= c("03H 00M 00S") & TimeSide2 <= c("05H 00M 00S"))) %>%
  mutate(Antenna = `Side (1/2)`) %>%
  mutate(Antenna2 = if_else(Antenna == 1 & Channel == 1, 1,
          if_else(Antenna == 2 & Channel == 1, 2,
                  if_else(Antenna == 1 & Channel == 2, 3,
                          if_else(Antenna == 2 & Channel == 2, 4, -9999)))))


### All animal PIT tags used
tags <- as.data.frame(unique(alldat$PITTAG)); colnames(tags) <- "PITTAG"
tags$RFID <- substr(tags$PITTAG, 8, 15)

### Read in individual info
siz <- read_csv("Data/SnakeSizes.csv",show_col_types = FALSE) %>%
  select(DATETEST,TRIAL,PITTAG,SVL,TL,TAILBREAK,SEX,WEIGHT,BULGE,BATCHMARK,ARENASIDE) %>%
  mutate_at(vars(PITTAG), factor) %>%
  rename(Date = DATETEST) %>%
  mutate(Date = dmy(Date)) %>%
  filter(PITTAG != 982091065198473) %>% # remove individual who lost tag before trial and there were no trials ever for this tag
  filter(!(PITTAG == 982091065198381 & TRIAL == 10)) %>%  # remove individual who lost tag before trial
  mutate(SEX = ifelse(SEX == "F", 1, 0)) %>%
  mutate(TagLoc = ifelse(TRIAL < 15, 1, 0)) # trials 1-14 had tags in neck (1) while 15 and 16 had tags in posterior (0)
siz$ID <- c(1:nrow(siz))

### Bind data to get scanning success by snake characteristics
alldat2 <- subset(alldat, !is.na(`ApproxDist(in)`))
dscan <- inner_join(alldat2[,c("Inits","Date","PITTAG","Read (1/0)","TimeSide","ApproxDist(in)","Antenna2")],siz[,c("Date","TRIAL","PITTAG","SVL","SEX","WEIGHT","TagLoc","ID")], by=c("Date","PITTAG"))
colnames(dscan) <- c("Inits","Date","PITTAG","Scan","TimeSide","Dist","Antenna","TRIAL","SVL","SEX","WEIGHT","TagLoc","ID")

### Limit analysis to scans within 2 inches as failure to scan was recorded beyond this
scan <- subset(dscan, Dist <= 2)

### Sort dataset by trial and snake ID
scan <- scan[order(scan$TRIAL,scan$ID),]

### Adjust ID of snakes to account for missing snake 61 - allows model to still loop over snake ID
scan$ID <- ifelse(scan$ID < 61, scan$ID, scan$ID - 1)


#### CHECK FOR CORRELATION IN COVARIATES (r > 0.6 considered to be highly correlated).----
cor(scan[,c(6:13)])
corrplot::corrplot(cor(scan[,c(6:13)]))
#weight and SVL, trial and ID, tag location and trial, and tag location and ID are all strongly correlated
#Weight and SVL are biologically related and should not be included in the same models
#However, tag location, trial, and ID are correlated due to the design of the study and can be included together


#### PREPARE MODEL INFO.----
N <- length(unique(scan$ID))                        # number of individuals
read <- as.vector(unlist(scan[,c("Scan")]))         # trial outcomes
ID <- as.vector(unlist(scan[,c("ID")]))             # ID of individuals
sex <- as.vector(unlist(scan[,c("SEX")])) + 1       # sex of individuals
trial <- as.vector(unlist(scan[,c("TRIAL")]))       # trial for each individual
size <- as.vector(unlist(scan[,c("SVL")]))          # svl for each individual
dist <- as.vector(unlist(scan[,c("Dist")]))         # distance of individual from antenna
dist <- as.vector(as.character(unlist(scan[,c("Dist")])))         # distance of individual from antenna
dist <- ifelse(dist == 0, 1,
               ifelse(dist == 0.5, 2,
                      ifelse(dist == 1, 3,
                             ifelse(dist == 1.5, 4, 5))))
weight <- as.vector(unlist(scan[,c("WEIGHT")]))     # weight of each individual
ant <- as.vector(unlist(scan[,c("Antenna")]))       # which antenna scanned
loc <- as.vector(unlist(scan[,c("TagLoc")])) + 1    # location of tag in snake



#### GLMER MODELS TO RUN.----

## MODEL ONE: Scan by Distance 
# dscan$Dist <- as.factor(dscan$Dist)
# m.dist <- glmer(Scan ~ Dist + (1|ID), data = dscan, family = binomial)
# summary(m.dist)
# anova(m.dist)

## MODEL TWO: Scan by Trial
# t.dist <- glm(Scan ~ TRIAL, data = dscan, family = binomial)
# summary(t.dist)

## MODEL THREE: Scan by SVL
# dscan$SVL <- scale(dscan$SVL)
# m.svl <- glmer(Scan ~ SVL + (1|ID), data = dscan, family = binomial)
# summary(m.svl)

## MODEL FOUR: Scan by Sex
# m.sex <- glmer(Scan ~ SEX + (1|ID), data = dscan, family = binomial)
# summary(m.sex)

## MODEL FIVE: Scan Weight
# m.wei <- glmer(Scan ~ WEIGHT + (1|ID), data = dscan, family = binomial)
# summary(m.wei)

## MODEL SIX: Scan by ArenaSide
# m.arena <- glmer(Scan ~ ARENASIDE + (1|ID), data = dscan, family = binomial)
# summary(m.arena)

## MODEL SEVEN: Scan by Location
# m.loc <- glmer(Scan ~ Location + (1|ID), data = dscan, family = binomial)
# summary(m.loc)

## MODEL EIGHT: Scan by Distance*SVL
# m.distsvl <- glmer(Scan ~ Dist + SVL + Dist * SVL + (1|ID), data = dscan, family = binomial)
# summary(m.distsvl)

## MODEL NINE: Scan by SVL and Location



#### LOGISTIC REGRESSION IN JAGS.----

########################################################

## MODEL ONE: Scanning success by global model (sex + dist + loc + size/weight + size/weight*loc) with random effect of ID

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  for(c in 1:C){
    b1[c] ~ dnorm(0,5)  # Female, male
    b3[c] ~ dnorm(0,5)  # Anterior, posterior
    b5[c] ~ dnorm(0,5)  # Anterior, posterior * Size
  }
  for(d in 1:D){
    b2[d] ~ dnorm(0,5)  # Distance (0,0.5,1.0,1.5,2.0)
  }
  b4 ~ dnorm(0,5)       # Size
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1[Sex[i]] + b2[Dist[i]] + b3[Loc[i]] + b4*Size[i] + b5[Loc[i]]*Size[i] + eta[ID[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCatCon_GlobalSize.txt")

########################################################

# MCMC settings 
nc = 3
ni = 1000
nb = 100
nthin = 1


## MODEL: Scan by Global + Size.----
# data <- list(read=read, ID=ID, Sex=sex, Dist=as.integer(dist), Loc=loc, Size=size, N=N, C=length(unique(sex)), D=length(unique(dist)))
# modnam <- c("GlobalSize")

## MODEL: Scan by Global + Weight.----
# data <- list(read=read, ID=ID, Sex=sex, Dist=as.integer(dist), Loc=loc, Size=weight, N=N, C=length(unique(sex)), D=length(unique(dist)))
# modnam <- c("GlobalWeight")

## All Models.----
parameters<-c('b0','b1','b2','b3','b4','b5','loglike','eta') #'p',

inits <- function() {
  list()
}

out <- jagsUI::jags(model.file ="Models/LogitCatCon_GlobalSize.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL TWO: Scanning success by sex, loc, size/weight, loc*size with random effect of ID

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  for(c in 1:C){
    b1[c] ~ dnorm(0,5)  # Female, male
    b2[c] ~ dnorm(0,5)  # Anterior, posterior
    b4[c] ~ dnorm(0,5)  # Anterior, posterior * size
  }
  b3 ~ dnorm(0,5)       # Size
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1[Sex[i]] + b2[Loc[i]] + b3*Size[i] + b4[Loc[i]]*Size[i] + eta[ID[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCatCon_SexLocSizeInt.txt")

########################################################

# MCMC settings 
nc = 3
ni = 1000
nb = 100
nthin = 1


## MODEL: Scan by Sex, Loc, Size, loc*size.----
# data <- list(read=read, ID=ID, Sex=sex, Loc=loc, Size=size, N=N, C=length(unique(sex)))
# modnam <- c("SexLocSizeInt")

## MODEL: Scan by Sex, Loc, Weight, Loc*Weight.----
# data <- list(read=read, ID=ID, Sex=sex, Loc=loc, Size=weight, N=N, C=length(unique(sex)))
# modnam <- c("SexLocWeightInt")

## All Models.----
parameters<-c('b0','b1','b2','b3','b4','loglike','eta') #'p',

inits <- function() {
  list()
}

out <- jagsUI::jags(model.file ="Models/LogitCatCon_SexLocSizeInt.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL THREE: Scanning success by sex, dist, loc, size/weight with random effect of ID

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
  b4 ~ dnorm(0,5)       # Size
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1[Sex[i]] + b2[Dist[i]] + b3[Loc[i]] + b4*Size[i] + eta[ID[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCatCon_SexDistLocSize.txt")

########################################################

# MCMC settings 
nc = 3
ni = 1000
nb = 100
nthin = 1


## MODEL: Scan by Sex, Dist, Loc, Size.----
# data <- list(read=read, ID=ID, Sex=sex, Dist=as.integer(dist), Loc=loc, Size=size, N=N, C=length(unique(sex)), D=length(unique(dist)))
# modnam <- c("SexDistLocSize")

## MODEL: Scan by Sex, Dist, Loc, Weight.----
# data <- list(read=read, ID=ID, Sex=sex, Dist=as.integer(dist), Loc=loc, Size=weight, N=N, C=length(unique(sex)), D=length(unique(dist)))
# modnam <- c("SexDistLocWeight")

## All Models.----
parameters<-c('b0','b1','b2','b3','b4','loglike','eta') #'p',

inits <- function() {
  list()
}

out <- jagsUI::jags(model.file ="Models/LogitCatCon_SexDistLocSize.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

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

# MCMC settings 
nc = 3
ni = 1000
nb = 100
nthin = 1


## MODEL: Scan by Sex, Dist, Loc.----
# data <- list(read=read, ID=ID, Sex=sex, Dist=as.integer(dist), Loc=loc, N=N, C=length(unique(sex)), D=length(unique(dist)))
# modnam <- c("SexDistLoc")


## All Models.----
parameters<-c('b0','b1','b2','b3','loglike','eta') #'p',

inits <- function() {
  list()
}

out <- jagsUI::jags(model.file ="Models/LogitCatCon_SexDistLoc.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL FIVE: Scanning success by sex, dist, size/weight with random effect of ID

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
  b3 ~ dnorm(0,5)       # Size
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1[Sex[i]] + b2[Dist[i]] + b3*Size[i] + eta[ID[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCatCon_SexDistSize.txt")

########################################################

# MCMC settings 
nc = 3
ni = 1000
nb = 100
nthin = 1


## MODEL: Scan by Sex, Dist, Size.----
# data <- list(read=read, ID=ID, Sex=sex, Dist=as.integer(dist), Size=size, N=N, C=length(unique(sex)), D=length(unique(dist)))
# modnam <- c("SexDistSize")

## MODEL: Scan by Sex, Dist, Weight.----
# data <- list(read=read, ID=ID, Sex=sex, Dist=as.integer(dist), Size=weight, N=N, C=length(unique(sex)), D=length(unique(dist)))
# modnam <- c("SexDistWeight")

## All Models.----
parameters<-c('b0','b1','b2','b3','loglike','eta') #'p',

inits <- function() {
  list()
}

out <- jagsUI::jags(model.file ="Models/LogitCatCon_SexDistSize.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL SIX: Scanning success by dist, loc, size/weight with random effect of ID

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
  b3 ~ dnorm(0,5)       # Size
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1[Loc[i]] + b2[Dist[i]] + b3*Size[i] + eta[ID[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCatCon_DistLocSize.txt")

########################################################

# MCMC settings 
nc = 3
ni = 1000
nb = 100
nthin = 1


## MODEL: Scan by Dist, Loc, Size.----
# data <- list(read=read, ID=ID, Loc=loc, Dist=as.integer(dist), Size=size, N=N, C=length(unique(sex)), D=length(unique(dist)))
# modnam <- c("DistLocSize")

## MODEL: Scan by Dist, Loc, Weight.----
# data <- list(read=read, ID=ID, Loc=loc, Dist=as.integer(dist), Size=weight, N=N, C=length(unique(sex)), D=length(unique(dist)))
# modnam <- c("DistLocWeight")

## All Models.----
parameters<-c('b0','b1','b2','b3','loglike','eta') #'p',

inits <- function() {
  list()
}

out <- jagsUI::jags(model.file ="Models/LogitCatCon_DistLocSize.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL SEVEN: Scanning success by loc, size/weight, loc*size/weight with random effect of ID 

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  for(c in 1:C){
    b1[c] ~ dnorm(0,5)     # Anterior, posterior
    b3[c] ~ dnorm(0,5)     # Anterior, posterior * Size
  }
  b2 ~ dnorm(0,5)          # Size
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1[Loc[i]] + b2*Size[i] + b3[Loc[i]]*Size[i] + eta[ID[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCatCon_LocSizeInt.txt")

########################################################

# MCMC settings 
nc = 3
ni = 1000
nb = 100
nthin = 1


## MODEL: Scan by Loc, Size, and Loc*Size.----
# data <- list(read=read, ID=ID, Loc=loc, Size=size, N=N, C=length(unique(loc)))
# modnam <- c("LocSizeInt")

## MODEL: Scan by Loc, Weight, and Loc*Weight.----
# data <- list(read=read, ID=ID, Loc=loc, Size=weight, N=N, C=length(unique(loc)))
# modnam <- c("LocWeightInt")


## All Models.----
parameters<-c('b0','b1','b2',"b3","loglike",'eta') #p

inits <- function() {
  list()
}

out <- jagsUI::jags(model.file ="Models/LogitCatCon_LocSizeInt.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

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

# MCMC settings 
nc = 3
ni = 1000
nb = 100
nthin = 1


## MODEL: Scan by Sex and Dist.----
# data <- list(read=read, ID=ID, cov1=sex, cov2=as.integer(dist), N=N, C=length(unique(loc)), D=length(unique(dist)))
# modnam <- c("SexDist")

## MODEL: Scan by Sex and Loc.----
# data <- list(read=read, ID=ID, cov1=sex, cov2=loc, N=N, C=length(unique(sex)), D=length(unique(loc)))
# modnam <- c("SexLoc")

## MODEL: Scan by Dist and Loc.----
# data <- list(read=read, ID=ID, cov1=as.integer(dist), cov2=loc, N=N, C=length(unique(dist)), D=length(unique(loc)))
# modnam <- c("DistLoc")


## All Models.----
parameters<-c('b0','b1','b2',"loglike",'eta') #p

inits <- function() {
  list()
}

out <- jagsUI::jags(model.file ="Models/LogitCatCat.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL ELEVEN, TWELVE, THIRTEEN: Scanning success by categorical covariate (e.g., sex or dist or loc) and continuous covariate (e.g., size/weight) with random effect of ID 

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  for(c in 1:C){
    b1[c] ~ dnorm(0,5)     # e.g., Female, male
  }
  b2 ~ dnorm(0,5)          # e.g., Size
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

# MCMC settings 
nc = 3
ni = 1000
nb = 100
nthin = 1


## MODEL: Scan by Sex and Size.----
# data <- list(read=read, ID=ID, cov1=sex, cov2=size, N=N, C=length(unique(sex)))
# modnam <- c("SexSize")

## MODEL: Scan by Dist and Size.----
# data <- list(read=read, ID=ID, cov1=as.integer(dist), cov2=size, N=N, C=length(unique(dist)))
# modnam <- c("DistSize")

## MODEL: Scan by Loc and Size.----
# data <- list(read=read, ID=ID, cov1=loc, cov2=size, N=N, C=length(unique(loc)))
# modnam <- c("LocSize")

## MODEL: Scan by Sex and Weight.----
# data <- list(read=read, ID=ID, cov1=sex, cov2=weight, N=N, C=length(unique(sex)))
# modnam <- c("SexWeight")

## MODEL: Scan by Dist and Weight.----
# data <- list(read=read, ID=ID, cov1=as.integer(dist), cov2=weight, N=N, C=length(unique(dist)))
# modnam <- c("DistWeight")

## MODEL: Scan by Loc and Weight.----
# data <- list(read=read, ID=ID, cov1=loc, cov2=weight, N=N, C=length(unique(loc)))
# modnam <- c("LocWeight")


## All Models.----
parameters<-c('b0','b1','b2',"loglike",'eta') #p

inits <- function() {
  list()
}

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

# MCMC settings 
nc = 3
ni = 1000
nb = 100
nthin = 1


## MODEL: Scan by Sex.----
# data <- list(read=read, ID=ID, cov=sex, N=N, C=length(unique(sex)))
# modnam <- c("Sex")

## MODEL: Scan by Distance.----
# data <- list(read=read, ID=ID, cov=as.integer(dist), N=N, C=length(unique(dist)))
# modnam <- c("Dist")

## MODEL: Scan by Location.----
# data <- list(read=read, ID=ID, cov=loc, N=N, C=length(unique(loc)))
# modnam <- c("Loc")

## All Models.----
parameters<-c('b0','b1','loglike','eta') #'p',

inits <- function() {
  list()
}

out <- jagsUI::jags(model.file ="Models/LogitCat.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL SEVENTEEN: Scanning success by continuous covariate (e.g., size/weight) with random effect of ID 

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  b1 ~ dnorm(0,5)     # Size
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

# MCMC settings 
nc = 3
ni = 1000
nb = 100
nthin = 1

## MODEL: Scan by SVL.----
# data <- list(read=read, ID=ID, cov=size, N=N)
# modnam <- c("Size")

## MODEL: Scan Weight.----
# data <- list(read=read, ID=ID, cov=weight, N=N)
# modnam <- c("Weight")

## All Models.----
parameters<-c('b0','b1',"loglike","eta") #p

inits <- function() {
  list()
}

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

# MCMC settings 
nc = 3
ni = 1000
nb = 100
nthin = 1

## MODEL: Scan by Null.----
# data <- list(read=read, ID=ID, N=N)
# modnam <- c("Null")

## All Models.----
parameters<-c('b0',"loglike","eta") #p

inits <- function() {
  list()
}

out <- jagsUI::jags(model.file ="Models/LogitNull.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))




## MODEL COMPARISON.----
# Watanabe-Akaike Information Criterion (WAIC)


### METHOD ONE? I DON'T THINK THIS IS CORRECT

# # Read in matrix of log-likelihood (LL) values for all models
# # Deviance = -2*LL
# waicvals <- matrix(NA, nrow = length(out$sims.list$deviance), ncol = length(file_list))
# 
# #create a list of the files from your target directory
# file_list <- list.files(path="Results")
# 
# library(loo)
# 
# waicALL <- matrix(NA, nrow=6, ncol=length(file_list), dimnames = list(c("elpd_waic","p_waic","waic","elpd_loo","p_loo","looic"), sub("\\..*","",file_list)))
# 
# for(i in 1:length(file_list)){
#   load(file=paste("Results/",file_list[i],sep=""))
#   waicvals <- cbind((out$sims.list$deviance[1:900])/(-2),(out$sims.list$deviance[901:1800])/(-2),(out$sims.list$deviance[1801:2700])/(-2))
#   int <- waic(waicvals, mc.cores = 3)
#   int2 <- loo(waicvals, mc.cores = 3)
#   waicALL[1,i] <- int$estimates[1,1]
#   waicALL[2,i] <- int$estimates[2,1]
#   waicALL[3,i] <- int$estimates[3,1]
#   waicALL[4,i] <- int2$estimates[1,1]
#   waicALL[5,i] <- int2$estimates[2,1]
#   waicALL[6,i] <- int2$estimates[3,1]
# }
# 
# test <- as.data.frame(waicALL[3,])
# test <- test[order(-test$`waicALL[3, ]`), , drop = FALSE]

### METHOD TWO?
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

waicALL2 <- matrix(NA, nrow=length(file_list), ncol=1, dimnames = list(sub("\\..*","",file_list),c("waic")))

for(i in 1:length(file_list)){
  load(file=paste("Results/",file_list[i],sep=""))
  waicALL2[i,1] <- calc.waic(out)
}
waicALL3 <- as.data.frame(waicALL2[order(waicALL2),]); colnames(waicALL3) <- c("waic")
waicALL3$deltawaic <- waicALL3$waic - min(waicALL3$waic, na.rm = TRUE)
waicALL3$deltawaic <- round(waicALL3$deltawaic, 2)


########################################################

## MODEL: Scanning success by categorical covariate (without random effect of ID) 

# cat("
# model {
# 
#   #Priors
#   #Intercept
#   b0 ~ dnorm(0,5)
#   #Covariate
#   for(c in 1:C){
#     b1[c] ~ dnorm(0,5)
#   }
# 
#   #Likelihood
#   for(i in 1:length(read)){
#     read[i] ~ dbern(p[i])
#     logit(p[i]) <- b0 + b1[cov[i]]
#     loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
#   }
# }
# ",file = "Models/LogitCat_noRE.txt")

########################################################

###### MODEL: Scan by Trial - not sure this is useful aside from being another random effect.----
# data <- list(read=read, ID=ID, cov=trial, C=length(unique(trial)))
# modnam <- c("Trial")


## MODEL: Scan by Antenna.----
# data <- list(read=read, ID=ID, cov=ant, N=N, C=length(unique(ant)))
# modnam <- c("Antenna")


## MODEL: Scan by Loc, Dist, and Loc*Dist.----
# data <- list(read=read, ID=ID, cov1=loc, cov2=dist, N=N, C=length(unique(loc)))
# modnam <- c("LocDist")

