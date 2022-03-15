library(tidyverse);library(dplyr);library(lubridate);library(jagsUI);library(ggplot2);library(lme4)

#### LOAD AND PREPARE DATA.----

### Read in detailed trial info
alldat <- read_csv("Data/RePTaR_trials_Aug2021_AllData.csv",show_col_types = FALSE) %>%
  mutate(PITTAG = as.character(PITTAG)) %>%
  drop_na(TimeSide) %>%
  mutate(TimeSide2 = lubridate::hms(TimeSide)) %>%   ## pointless warnings as lubridate is a POS
  mutate(Date = dmy(Date)) %>%
  # filter((TimeSide >= c("19:00:00") & TimeSide <= c("21:00:00")) |
  #          (TimeSide >= c("3:00:00") & TimeSide <= c("5:00:00")))
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
#weight and SVL, trial and ID, tag location and trial, and tag location and ID are all strongly correlated and should not be included in the same models


#### PREPARE MODEL INFO.----
N <- length(unique(scan$ID))                        # number of individuals
read <- as.vector(unlist(scan[,c("Scan")]))         # trial outcomes
ID <- as.vector(unlist(scan[,c("ID")]))             # ID of individuals
sex <- as.vector(unlist(scan[,c("SEX")])) + 1       # sex of individuals
trial <- as.vector(unlist(scan[,c("TRIAL")]))       # trial for each individual
size <- as.vector(unlist(scan[,c("SVL")]))          # svl for each individual
dist <- as.vector(unlist(scan[,c("Dist")]))         # distance of individual from antenna
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

## MODEL: Scanning success by categorical covariate with random effect of ID 

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  for(c in 1:C){
    b1[c] ~ dnorm(0,5)
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

## MODEL: Scan by Antenna.----
# data <- list(read=read, ID=ID, cov=ant, N=N, C=length(unique(ant)))
# modnam <- c("Antenna")

## All Models.----
parameters<-c('b0','b1','loglike','eta') #'p',

inits <- function() {
  list()
}

out <- jagsUI::jags(model.file ="Models/LogitCat.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL: Scanning success by categorical covariate (without random effect of ID) 

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  for(c in 1:C){
    b1[c] ~ dnorm(0,5)
  }

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1[cov[i]]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCat_noRE.txt")

########################################################

# MCMC settings 
nc = 3
ni = 1000
nb = 100
nthin = 1


###### MODEL: Scan by Trial - not sure this is useful aside from being another random effect.----
# data <- list(read=read, ID=ID, cov=trial, C=length(unique(trial)))
# modnam <- c("Trial")

## MODEL: Scan by Location.----
# data <- list(read=read, ID=ID, cov=loc, C=length(unique(loc)))
# modnam <- c("Loc")

## All Models.----
parameters<-c('b0','b1',"loglike") #'p',

inits <- function() {
  list()
}

out <- jagsUI::jags(model.file ="Models/LogitCat_noRE.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL: Scanning success by continuous covariate with random effect of ID 

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  b1 ~ dnorm(0,5)
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)
  
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }

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
# modnam <- c("SVL")

## MODEL: Scan by Distance.----
# data <- list(read=read, ID=ID, cov=dist, N=N)
# modnam <- c("Dist")

## MODEL: Scan Weight.----
# data <- list(read=read, ID=ID, cov=weight, N=N)
# modnam <- c("Weight")

## MODEL: Scan by Distance*SVL
## MODEL: Scan by SVL and Location

## All Models.----
parameters<-c('b0','b1',"loglike","eta") #p

inits <- function() {
  list()
}

out <- jagsUI::jags(model.file ="Models/LogitCont.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



########################################################

## MODEL: Scanning success by continuous and categorical covariate and interaction without random effect of ID 

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  for(c in 1:C){
    b1[c] ~ dnorm(0,5)
    b3[c] ~ dnorm(0,5)
  }
  b2 ~ dnorm(0,5)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1[cov1[i]] + b2*cov2[i] + b3[cov1[i]]*cov2[i]
    loglike[i] <- logdensity.bern(read[i],p[i])  # For WAIC
  }
}
",file = "Models/LogitCatCon_noRE.txt")

########################################################

# MCMC settings 
nc = 3
ni = 1000
nb = 100
nthin = 1


## MODEL: Scan by Loc, Size, and Loc*Size.----
# data <- list(read=read, ID=ID, cov1=loc, cov2=size, N=N, C=length(unique(loc)))
# modnam <- c("LocSize")

## MODEL: Scan by Loc, Dist, and Loc*Dist.----
# data <- list(read=read, ID=ID, cov1=loc, cov2=dist, N=N, C=length(unique(loc)))
# modnam <- c("LocDist")


## All Models.----
parameters<-c('b0','b1','b2',"b3","loglike") #p

inits <- function() {
  list()
}

out <- jagsUI::jags(model.file ="Models/LogitCatCon_noRE.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))




########################################################

## MODEL: Scanning success by continuous covariates and interaction with random effect of ID 

cat("
model {

  #Priors
  #Intercept
  b0 ~ dnorm(0,5)
  #Covariate
  b1 ~ dnorm(0,5)
  b2 ~ dnorm(0,5)
  b3 ~ dnorm(0,5)
  for(j in 1:N){
    #Random effect - ID
    eta[j] ~ dnorm(0, tau_p)
  }
  #Hyperprior random effect - ID
  tau_p ~ dunif(0,4)

  #Likelihood
  for(i in 1:length(read)){
    read[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1*cov1[i] + b2*cov2[i] + b3*cov1[i]*cov2[i] + eta[ID[i]]
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


## MODEL: Scan by Dist, Size, and Dist*Size.----
# data <- list(read=read, ID=ID, cov1=dist, cov2=size, N=N)
# modnam <- c("DistSize")


## All Models.----
parameters<-c('b0','b1','b2',"b3","loglike","eta") #p

inits <- function() {
  list()
}

out <- jagsUI::jags(model.file ="Models/LogitCatCat.txt", data, inits=inits, parameters.to.save = parameters, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin, parallel = TRUE, n.cores = 3)

save(out, file=paste("Results/ScanSuccess_",modnam,".Rdata",sep=""))



## MODEL COMPARISON.----
# Watanabe-Akaike Information Criterion (WAIC)


### METHOD ONE?

# Read in matrix of log-likelihood (LL) values for all models
# Deviance = -2*LL
waicvals <- matrix(NA, nrow = length(out$sims.list$deviance), ncol = 9)

#create a list of the files from your target directory
file_list <- list.files(path="Results")

library(loo)

waicALL <- matrix(NA, nrow=6, ncol=9, dimnames = list(c("elpd_waic","p_waic","waic","elpd_loo","p_loo","looic"), sub("\\..*","",file_list)))

for(i in 1:9){
  load(file=paste("Results/",file_list[i],sep=""))
  waicvals <- cbind((out$sims.list$deviance[1:900])/(-2),(out$sims.list$deviance[901:1800])/(-2),(out$sims.list$deviance[1801:2700])/(-2))
  int <- waic(waicvals, mc.cores = 3)
  int2 <- loo(waicvals, mc.cores = 3)
  waicALL[1,i] <- int$estimates[1,1]
  waicALL[2,i] <- int$estimates[2,1]
  waicALL[3,i] <- int$estimates[3,1]
  waicALL[4,i] <- int2$estimates[1,1]
  waicALL[5,i] <- int2$estimates[2,1]
  waicALL[6,i] <- int2$estimates[3,1]
}


### METHOD TWO?
## Modified from https://github.com/heathergaya/JAGS-NIMBLE-Tutorials/blob/master/Known_Fate/Known_Fate_Models.Rmd

calc.waic <- function(x){
  #find the output that relates to loglike
  like <- x$sims.list$loglike
  fbar <- colMeans(exp(like)) #mean likelihood 
  Pw <- sum(apply(like,2,var)) #mean variance in log-likelihood 
  WAIC<- -2*sum(log(fbar))+2*Pw
  return(WAIC)
}

waicALL2 <- matrix(NA, nrow=9, ncol=1, dimnames = list(sub("\\..*","",file_list),c("waic")))

for(i in 1:9){
  load(file=paste("Results/",file_list[i],sep=""))
  waicALL2[i,1] <- calc.waic(out)
}


