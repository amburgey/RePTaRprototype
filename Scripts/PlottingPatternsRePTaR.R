#### REMOTE PASSIVE INTEGRATED TRANSPONDER (PIT) TAG READER = RePTaR
#### Brown treesnakes were PIT tagged and released into experimental trial arena to test the scanning success of this RePTaR reader prototype
#### Code written by Staci Amburgey
#### Involving data collected from August 1-August 30 2021

rm(list = ls())

library(tidyr); library(dplyr); library(lubridate); library(tidyverse); library(ggplot2); library(LaCroixColoR)



#### DATA PREP FOR SCANNING DATA DOWNLOADED FROM REPTAR.----

# Read in scanning data files as items in a list and then create a combined data frame. 
# These are the files created by the RePTaR units that consist of records of when the snakes scan (so they do not contain failed scans).
file_names1 <- dir("Data/Scans/Station_1-2021/Station_1/")
data_frame1 <- as_tibble(do.call(rbind,lapply(paste("Data/Scans/Station_1-2021/Station_1/",file_names1, sep=""),read.csv)))
file_names2 <- dir("Data/Scans/Station_2-2021/Station_2/")
data_frame2 <- as_tibble(do.call(rbind,lapply(paste("Data/Scans/Station_2-2021/Station_2/",file_names2, sep=""),read.csv)))
data_frame <- rbind(data_frame1,data_frame2)

# Shape data into form for plotting
df <- data_frame %>%
  separate(Timestamp, sep = ", ", into = c("Date","Time")) %>%    # separate time stamp
  filter((Time >= c("18:00:00") & Time <= c("23:59:58")) |        # retain scans only when trials occurred (6pm to 6am)
           (Time <= c("06:00:00") & Time >= c("00:00:00"))) %>%
  unite(DateTime, Date:Time, sep = " ") %>%                       # recombine date and time for sorting
  mutate(DateTime = mdy_hms(DateTime)) %>%                        # specify column as date and time
  mutate(Day = mday(DateTime)) %>%                                # numeric day of the month (5-30)
  group_by(dt=cut(DateTime, breaks= seq(
    from=as.POSIXct("2021-08-05 18:00:00", tz="UTC"),             # create column of date-by-time categories and group scans by these
    to=as.POSIXct("2021-08-30 06:00:00", tz="UTC"),
    by="hour"
  )))

# Shape data to get frequency of scans by time category
dfsum <- df %>% 
  separate(dt, sep = " ", into = c("Date","Time")) %>%    # separate time stamp
  count(Time) %>%
  mutate(Time = as.factor(Time)) %>%
  mutate(Time = fct_relevel(Time, c("18:00:00","19:00:00","20:00:00","21:00:00","22:00:00","23:00:00",
                                    "00:00:00","01:00:00","02:00:00","03:00:00","04:00:00","05:00:00","06:00:00" )))


#### PLOT ONE: Figure of scanning frequency across 6pm to 6am trial.----

p <- lacroix_palette("Pamplemousse", n=20, type="continuous")
p[1:20]

labs <- c("6pm","7pm","8pm","9pm","10pm","11pm","12am","1am","2am","3am","4am","5am","6am")
# Specify which palette colors desired
colsLC <- c("#F7AA14","#F6B90B","#F5C903","#F5D523","#F6DE5F","#E4E292","#89CE9F","#2DBAAC","#14A7B3","#0C95BA","#0A7AAF","#10518C","#172869")
mid <- mean(dfsum$n) 

plotTimes <- ggplot(dfsum, aes(x = Time, y = n, color = n, fill = n)) +
  geom_bar(stat = "identity", alpha = 0.5) +
  scale_fill_gradient2(midpoint=mid, low=colsLC[9], mid=colsLC[4], high=colsLC[1]) +
  scale_color_gradient2(midpoint=mid, low=colsLC[9], mid=colsLC[4], high=colsLC[1]) +
  ylab("Scan Count") +
  scale_x_discrete(labels = labs) +
  theme(legend.position = "none", panel.background = element_rect(fill = "white", color = "darkgrey"), axis.title = element_text(size = 12), axis.text = element_text(size = 12))

# Save plot to file
# png(file="ScanCounts.png",width=8,height=5,units="in",res=600)
# plotTimes
# dev.off()



#### DATA PREP FOR DETAILED BEHAVIORAL TRIAL AND SNAKE INFO.----

## Read in detailed behavioral trial info and format for plotting
alldat <- read_csv("Data/RePTaR_trials_AllData.csv",show_col_types = FALSE) %>%
  mutate(PITTAG = as.character(PITTAG)) %>%
  mutate(PITTAG = str_match(PITTAG, "\"(.*?)\"")[,2]) %>% # remove quotations used to maintain full 15-digit ID (instead of scientific notation issues)
  drop_na(TimeSide) %>%                             # drop times when snake never scanned or didn't get close enough to the antenna to count as a scanning trial (<=2in)
  filter(TimeSide != 1) %>%                         # drop one instance where incorrect time entered (outside of time range retained for analysis anyway)
  mutate(TimeSide2 = lubridate::hms(TimeSide)) %>%  # specify as time
  mutate(Date = dmy(Date)) %>%                      # specify as date
  filter((TimeSide2 >= c("19H 00M 00S") & TimeSide2 <= c("21H 00M 00S")) |   # limit to the subsetted trial times processed (7-9pm, 3-5am)
           (TimeSide2 >= c("03H 00M 00S") & TimeSide2 <= c("05H 00M 00S")))

## Get expanded PITTAG IDs (full 15-digit ID) and reduced IDs (8-digit ID) for use in combining trait and scanning dataframes
tags <- as.data.frame(unique(alldat$PITTAG)); colnames(tags) <- "PITTAG"
tags$RFID <- substr(tags$PITTAG, 8, 15)

## Read in individual info and format for plotting
siz <- read_csv("Data/SnakeTraits.csv",show_col_types = FALSE) %>%
  select(DATETEST,TRIAL,PITTAG,SVL,TL,TAILBREAK,SEX,WEIGHT,BULGE,BATCHMARK,ARENASIDE) %>%  # subset to columns of interest
  mutate(PITTAG = str_match(PITTAG, "\"(.*?)\"")[,2]) %>% # remove quotations used to maintain full 15-digit ID (instead of scientific notation issues)
  mutate_at(vars(PITTAG), factor) %>% 
  rename(Date = DATETEST) %>%                             # denote this day of trial as the date of interest
  mutate(Date = dmy(Date)) %>%                            # specify as date
  filter(PITTAG != 982091065198473) %>%                   # remove this tag as snake escaped when first trialed and there were no scans ever for this tag
  filter(!(PITTAG == 982091065198381 & TRIAL == 10)) %>%  # remove individual who lost tag before trial (but there other, earlier trials with this tag)
  mutate(SEX = ifelse(SEX == "F", 1, 0))                  # specify sex as 1 (female), 0 (male)



#### FURTHER DATA PREP FOR SCANNING DATA DOWNLOADED FROM REPTAR AND JOINING WITH SNAKE INFO.----

# Read in scanning data files collected from RePTaR as items in a list
scan_df <- data_frame %>%
  mutate(RFID = as.character(RFID)) %>%
  inner_join(tags, by="RFID") %>%                                # expand tag ID to full 15-digits
  mutate(TimestampOG = Timestamp) %>%                            # create new column so not losing original timestamp
  separate(Timestamp, sep = ", ", into = c("Date","Time")) %>%   # split timestamp
  filter((Time >= c("18:00:00") & Time <= c("23:59:58")) |       # filter scans to 6pm to 6am trial period
           (Time <= c("06:00:00") & Time >= c("00:00:00"))) %>%
  mutate(Date = mdy(Date)) %>%                                   # specify as date
  arrange(TimestampOG) %>%                                       # sort times
  mutate(TimestampOG = str_replace(TimestampOG, ", ", " ")) %>%  # remove comma and replace with space
  mutate(TimestampOG = mdy_hms(TimestampOG)) %>%                 # specify date and time
  mutate(Timestamp2 = TimestampOG - hours(7)) %>%                # subtract 7 hours from times (to no longer straddle midnight) to be able to better identify trials
  mutate(Trial = as.numeric(factor(as.numeric(as.Date(Timestamp2))-18843)))   # subtract 18843 to get a numeric indicator of trial (1 through 16)

## Get number of scans per PITTAG ID and TRIAL (as PIT tags were re-used)
numscans <- scan_df %>%
  mutate_at(vars(PITTAG), factor) %>%
  rename(TRIAL = Trial) %>%
  group_by(TRIAL,PITTAG) %>%
  count()

## Join scanning info per snake per trial with snake information
nscans <- inner_join(numscans,siz, by=c("TRIAL","PITTAG")) %>%
  add_column(ID = 1:nrow(numscans))

## Use to determine colors available to use for plotting below (change n to get more or less colors)
p <- lacroix_palette("Pamplemousse", n=9, type="continuous")
p[1:9]
colsLC <- c("#EA7580","#F19097","#F6ABA2","#F7C79D","#89C1A5","#18B0B0","#0C95BA","#0D659E","#172869")



#### PLOT TWO: Two-panel figure of trial snake snout vent length (SVL) and weight.----

## Panel one. Snake sizes.

# Denote individual ID
siz$ID <- c(1:nrow(siz))
# Set up cut-off values for SVL categories
breaks <- c(634,727,820,913,1006,1099,1192,1285,1378,1434)
# Specify interval/bin labels
tags <- c("[634-727)", "[727-820)", "[820-913)", "[913-1006)", "[1006-1099)","[1099-1192)", "[1192-1285)","[1285-1378)", "[1378-1434)")
# Place values into bins
group_size <- cut(siz$SVL, 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)
# Inspect bins
summary(group_size)

mid <- mean(siz$SVL) 

plotSnake1 <- ggplot(data = as_tibble(group_size), aes(x=value)) + 
  geom_bar(fill="#18B0B0") + 
  labs(x='Snake SVL (mm)', y='Number of Snakes') +
  theme(panel.background = element_blank(),axis.text = element_text(size=12, angle = 50, vjust = 0.5, hjust = 0.8), axis.ticks = element_blank(), axis.title = element_text(size=14)) +
  guides(fill = "none")

## Panel two. Snake weights.

# Set up cut-off values for weight categories
breaks <- c(18,121,221,321,421,521,621,721,821,921,1093)
# Specify interval/bin labels
wts <- c("[18-121)", "[121-221)", "[221-321)", "[321-421)", "[421-521)","[521-621)", "[621-721)","[721-821)", "[821-921)","[921-1093)")
# Place values into bins
group_wt <- cut(siz$WEIGHT, 
                breaks=breaks, 
                include.lowest=TRUE, 
                right=FALSE, 
                labels=wts)
# Inspect bins
summary(group_wt)

plotSnake2 <- ggplot(data = as_tibble(group_wt), aes(x=value)) + 
  geom_bar(fill="#0C95BA") + 
  labs(x='Snake Weight (g)', y='Number of Snakes') +
  theme(panel.background = element_blank(),axis.text = element_text(size=12, angle = 50, vjust = 0.5, hjust = 0.8), axis.ticks = element_blank(), axis.title = element_text(size=14)) +
  guides(fill = "none")

# Save plot to file
# png(file="Figures/SnakeInfo.png",width=7,height=7,units="in",res=600)
# (plotSnake1 / plotSnake2)
# dev.off()



#### PLOT THREE: Figure of distance at which scanning occurred.----

# Subset data to just where scanning occurred
alldat2 <- subset(alldat, `Read (1/0)` == 1)
# Subset to instances where snakes were visible and a distance could be estimated
alldat2 <- subset(alldat2, !is.na(`ApproxDist(in)`))

# Join scanning info per snake per trial with snake information
dscans <- inner_join(alldat2,siz, by=c("Date","PITTAG"))
dscans$`ApproxDist(in)` <- as.factor(dscans$`ApproxDist(in)`)
dscans$NumDist <- as.numeric(as.character(droplevels(dscans$`ApproxDist(in)`)))
# Convert 1/0 sex to F/M for plotting
dscans$Sex2 <- ifelse(dscans$SEX == 0, "F", "M")

# Set up cut-off values for scanning distance bins
breaks <- c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)
# Specify interval/bin labels
tags <- c("0", "0.5", "1", "1.5", "2","2.5", "3","3.5", "4","4.5")
# Place values into bins
group_dist <- cut(as.numeric(as.character(alldat2$`ApproxDist(in)`)), 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)
# Inspect bins
summary(group_dist)


plotDist <- ggplot(data = as_tibble(group_dist), aes(x=value)) + 
  geom_bar(fill="#F7C69D") +
  labs(x='Distance of Scans (in)', y='Number of Scans') +
  theme(panel.background = element_blank(),axis.text = element_text(size=12),axis.title = element_text(size=15)) +
  guides(fill = "none")

# Save plot to file
# png(file="ScanDist.png",width=7,height=7,units="in",res=600)
# plotDist
# dev.off()



#### Qualitative description of scans and snake activity by trap and antenna.----
# Scanning by RePTaR device (unit 1 or 2)
scntrap <- table(alldat$ARENASIDE, alldat$`Read (1/0)`)
# Scanning by RePTaR device (unit 1 or 2) and antenna (side 1 or 2)
scnant <- table(alldat$ARENASIDE, alldat$`Side (1/2)`, alldat$`Read (1/0)`)
# Mouse watching by snakes
watch <- sum(alldat$WatchedMouse, na.rm = TRUE)/(nrow(subset(alldat, !is.na(WatchedMouse))))



#### OTHER POTENTIAL PLOTS (NOT IN MANUSCRIPT) ----

## Option. Number of scans by sex
scansex <- aggregate(nscans$n, by = list(Sex=nscans$SEX), FUN = sum)
# Add column of character M or F
scansex$Sex2 <- c("M","F")

plotScan1 <- ggplot(scansex, aes(y=x, x=Sex2, fill=as.factor(Sex2))) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values=c("#0B76AC","#172869")) +
  guides(fill=guide_legend(title="Sex")) +
  ylab(c("Number of Scans")) + xlab(c("Sex")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'darkgrey'), legend.position = NULL)


## Option. Count of scans by food bulge
scanbulge <- aggregate(nscans$n, by = list(Tail=nscans$BULGE), FUN = sum)
# Add column of Absent or Present
scanbulge$Meal <- c("Absent","Present")

plotScan2 <- ggplot(scanbulge, aes(y=x, x=Meal, fill=as.factor(Meal))) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values=c("#A6C4A3","#6CBEA8")) +
  guides(fill=guide_legend(title="Food bulge")) +
  ylab(c("Number of Scans")) + xlab(c("Food bulge")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'darkgrey'), legend.position = NULL)

