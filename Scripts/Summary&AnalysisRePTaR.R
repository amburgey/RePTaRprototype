library(tidyverse);library(dplyr);library(lubridate);library(lme4);library(performance);library(see);library(ggplot2);library(LaCroixColoR);library(patchwork);library(ordinal)

## Read in detailed behavioral trial info
alldat <- read_csv("Data/RePTaR_trials_Aug2021_AllData.csv",show_col_types = FALSE) %>%
  mutate(PITTAG = as.character(PITTAG)) %>%
  drop_na(TimeSide) %>%
  mutate(TimeSide2 = lubridate::hms(TimeSide)) %>%   ## warnings because some times are NA, ignore
  mutate(Date = dmy(Date)) %>%
  # filter((TimeSide >= c("19:00:00") & TimeSide <= c("21:00:00")) |
  #          (TimeSide >= c("3:00:00") & TimeSide <= c("5:00:00")))
  filter((TimeSide2 >= c("19H 00M 00S") & TimeSide2 <= c("21H 00M 00S")) |
           (TimeSide2 >= c("03H 00M 00S") & TimeSide2 <= c("05H 00M 00S")))

## Fully populate/expand PITTAG ID for use in combining later dataframes
tags <- as.data.frame(unique(alldat$PITTAG)); colnames(tags) <- "PITTAG"
tags$RFID <- substr(tags$PITTAG, 8, 15)

## Read in individual info, subset to columns of interest, convert PITTAG to factor, convert to recognized date, remove misc tags, and convert sex to 1/0 factor
siz <- read_csv("Data/SnakeSizes.csv",show_col_types = FALSE) %>%
  select(DATETEST,TRIAL,PITTAG,SVL,TL,TAILBREAK,SEX,WEIGHT,BULGE,BATCHMARK,ARENASIDE) %>%
  mutate_at(vars(PITTAG), factor) %>%
  rename(Date = DATETEST) %>%
  mutate(Date = dmy(Date)) %>%
  filter(PITTAG != 982091065198473) %>% # remove this tag as snake escaped when first used and there were no trials ever for this tag
  filter(!(PITTAG == 982091065198381 & TRIAL == 10)) %>%  # remove individual who lost tag before trial (but there other trials with this tag)
  mutate(SEX = ifelse(SEX == "F", 1, 0))
  

# Read in scanning data files collected from RePTaR as items in a list
file_names1 <- dir("Scans/Station_1-2021/Station_1/")
data_frame1 <- as_tibble(do.call(rbind,lapply(paste("Scans/Station_1-2021/Station_1/",file_names1, sep=""),read.csv)))  # hit error with read_csv
file_names2 <- dir("Scans/Station_2-2021/Station_2/")
data_frame2 <- as_tibble(do.call(rbind,lapply(paste("Scans/Station_2-2021/Station_2/",file_names2, sep=""),read.csv)))  # hit error with read_csv
data_frame <- rbind(data_frame1,data_frame2) %>%
  mutate(RFID = as.character(RFID)) %>%
  inner_join(tags, by="RFID") %>%
  mutate(TimestampOG = Timestamp) %>%
  separate(Timestamp, sep = ", ", into = c("Date","Time")) %>%
  filter((Time >= c("18:00:00") & Time <= c("23:59:58")) |
           (Time <= c("06:00:00") & Time >= c("00:00:00"))) %>%
  mutate(Date = mdy(Date)) %>%
  arrange(TimestampOG) %>%
  mutate(TimestampOG = str_replace(TimestampOG, ", ", " ")) %>%
  mutate(TimestampOG = mdy_hms(TimestampOG)) %>%
  mutate(Timestamp2 = TimestampOG - hours(7)) %>%
  mutate(Trial = as.numeric(factor(as.numeric(as.Date(Timestamp2))-18843)))
  
  
## Group individual info with scanning information
numscans <- data_frame %>%
  mutate_at(vars(PITTAG), factor) %>%
  rename(TRIAL = Trial) %>%
  group_by(TRIAL,PITTAG) %>%
  count()

## Create individual ID for random effect
nscans <- inner_join(numscans,siz, by=c("TRIAL","PITTAG")) %>%
  add_column(ID = 1:nrow(numscans))



#### CREATE PLOT OF SNAKE CHARACTERISTICS ----

## Panel one. Snake sizes.
siz$ID <- c(1:nrow(siz))

# set up cut-off values 
breaks <- c(634,727,820,913,1006,1099,1192,1285,1378,1434)
# specify interval/bin labels
tags <- c("[634-727)", "[727-820)", "[820-913)", "[913-1006)", "[1006-1099)","[1099-1192)", "[1192-1285)","[1285-1378)", "[1378-1434)")
# bucketing values into bins
group_size <- cut(siz$SVL, 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)
# inspect bins
summary(group_size)

p <- lacroix_palette("Pamplemousse", n=9, type="continuous")
colsLC <- c("#EA7580","#F19097","#F6ABA2","#F7C79D","#89C1A5","#18B0B0","#0C95BA","#0D659E","#172869")
mid <- mean(siz$SVL) 

plotSnake1 <- ggplot(data = as_tibble(group_size), aes(x=value)) + 
  geom_bar(fill="#18B0B0") + 
  labs(x='Snake SVL (mm)', y='Number of Snakes') +
  theme(panel.background = element_blank(),axis.text = element_text(size=12, angle = 50, vjust = 0.5, hjust = 0.8), axis.ticks = element_blank(), axis.title = element_text(size=14)) +
  guides(fill = "none")

## Panel two. Snake weights.

# set up cut-off values 
breaks <- c(18,121,221,321,421,521,621,721,821,921,1093)
# specify interval/bin labels
wts <- c("[18-121)", "[121-221)", "[221-321)", "[321-421)", "[421-521)","[521-621)", "[621-721)","[721-821)", "[821-921)","[921-1093)")
# bucketing values into bins
group_wt <- cut(siz$WEIGHT, 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=wts)
# inspect bins
summary(group_wt)

plotSnake2 <- ggplot(data = as_tibble(group_wt), aes(x=value)) + 
  geom_bar(fill="#0C95BA") + 
  labs(x='Snake Weight (g)', y='Number of Snakes') +
  theme(panel.background = element_blank(),axis.text = element_text(size=12, angle = 50, vjust = 0.5, hjust = 0.8), axis.ticks = element_blank(), axis.title = element_text(size=14)) +
  guides(fill = "none")


png(file="Figures/SnakeInfo.png",width=7,height=7,units="in",res=600)
(plotSnake1 / plotSnake2)
dev.off()


#### CREATE PLOTS OF NUMSCANS ----

## Option One. Count of scans by sex
scansex <- aggregate(nscans$n, by = list(Sex=nscans$SEX), FUN = sum)

# Not doing model selection because we're doing hypothesis testing: NumScans ~ Sex + (1|ID)
# model.SCAN <- lm(log(n) ~ SEX, data=nscans)
# summary(model.SCAN)
# plot(model.SCAN)

# Test assumptions
# hist(log(nscans$n))
# check_model(model.SCAN)
# shapiro.test(log(nscans$n))

# Plot scans by snake sex
p <- lacroix_palette("Pamplemousse", n=20, type="continuous")
# p[1:20]
# colsLC <- c("#0B76AC","#172869")
scansex$Sex2 <- c("M","F")

plotScan1 <- ggplot(scansex, aes(y=x, x=Sex2, fill=as.factor(Sex2))) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values=c("#0B76AC","#172869")) +
  guides(fill=guide_legend(title="Sex")) +
  ylab(c("Number of Scans")) + xlab(c("Sex")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'darkgrey'), legend.position = NULL)



## Option Two. Count of scans by tail break
scantail <- aggregate(nscans$n, by = list(Tail=nscans$TAILBREAK), FUN = sum)

# Plot scans by snake sex
scantail$Break <- c("Not broken","Broken")

plotScan2 <- ggplot(scantail, aes(y=x, x=Break, fill=as.factor(Break))) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values=c("#EA7580","#F7BA9F")) +
  guides(fill=guide_legend(title="Tail")) +
  ylab(c("Number of Scans")) + xlab(c("Tail")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'darkgrey'), legend.position = NULL)
plotScan2



# Get count of scans by food bulge
scanbulge <- aggregate(nscans$n, by = list(Tail=nscans$BULGE), FUN = sum)

# Plot scans by snake sex
scanbulge$Meal <- c("Absent","Present")

plotScan2 <- ggplot(scanbulge, aes(y=x, x=Meal, fill=as.factor(Meal))) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values=c("#A6C4A3","#6CBEA8")) +
  guides(fill=guide_legend(title="Food bulge")) +
  ylab(c("Number of Scans")) + xlab(c("Food bulge")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'darkgrey'), legend.position = NULL)
plotScan2



# Get count of scans by snake size
scansize <- aggregate(nscans$n, by = list(Size=nscans$SVL), FUN = sum)

# set up cut-off values 
breaks <- c(634,727,820,913,1006,1099,1192,1285,1378,1434)
# specify interval/bin labels
tags <- c("[634-727)", "[727-820)", "[820-913)", "[913-1006)", "[1006-1099)","[1099-1192)", "[1192-1285)","[1285-1378)", "[1378-1434)")
# bucketing values into bins
group_size <- cut(scansize$Size, 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)
# inspect bins
summary(group_size)

# Plot scans by snake sex
scanbulge$Meal <- c("Absent","Present")

plotScan2 <- ggplot(scanbulge, aes(y=x, x=Meal, fill=as.factor(Meal))) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values=c("#A6C4A3","#6CBEA8")) +
  guides(fill=guide_legend(title="Food bulge")) +
  ylab(c("Number of Scans")) + xlab(c("Food bulge")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'darkgrey'), legend.position = NULL)
plotScan2



#### DISTANCE SCANNED ----

alldat2 <- subset(alldat, `Read (1/0)` == 1)
alldat2 <- subset(alldat2, !is.na(`ApproxDist(in)`))

dscans <- inner_join(alldat2[,c(5:7,14:15,18:19)],siz[,c(1:3,7,11:12)], by=c("Date","PITTAG"))
dscans$`ApproxDist(in)` <- as.factor(dscans$`ApproxDist(in)`)
dscans$NumDist <- as.numeric(as.character(droplevels(dscans$`ApproxDist(in)`)))
dscans$Sex2 <- ifelse(dscans$SEX == 0, "F", "M")

# Not doing model selection because we're doing hypothesis testing: Distance ~ Sex + (1|ID)
# Do ordered logistic regression because the distance value (e.g., 0, 0.5, 1, etc) mean something by the order they are in
model.DIST <- clmm(`ApproxDist(in)` ~ SEX + (1|ID), data=dscans, Hess = TRUE)
summary(model.DIST)


# set up cut-off values 
breaks <- c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)
# specify interval/bin labels
tags <- c("0", "0.5", "1", "1.5", "2","2.5", "3","3.5", "4","4.5")
# bucketing values into bins
group_dist <- cut(as.numeric(as.character(alldat$`ApproxDist(in)`)), 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)
# inspect bins
summary(group_dist)


plotDist <- ggplot(data = as_tibble(group_dist), aes(x=value)) + 
  geom_bar(fill="#F7C69D") + 
  # scale_fill_manual(values = p) +
  # scale_color_manual(values = p) +
  labs(x='Snake Size (mm)', y='Distance of Scan') +
  theme(panel.background = element_blank(),axis.text = element_text(size=50),axis.title = element_text(size=55)) +
  guides(fill = "none")


png(file="ScanDist.png",width=7,height=7,units="in",res=600)
plotDist
dev.off()


plotDistSex <- ggplot(dscans, aes(x = as.factor(Sex2), y = NumDist, fill = Sex2)) +
  geom_boxplot(size = .75) +
  scale_fill_manual(values = c("#0B76AC","#172869")) +
  xlab("Sex") + ylab("Distance Scanned (in)") +
  theme(panel.background = element_blank(),axis.text = element_text(size=50), axis.title = element_text(size=55), legend.position = "none")

png(file="ScanDistSex.png",width=7,height=7,units="in",res=600)
plotDistSex
dev.off()
