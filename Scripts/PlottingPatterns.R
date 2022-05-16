#### REMOTE PASSIVE INTEGRATED TRANSPONDER (PIT) TAG READER = RePTaR
#### Brown treesnakes were PIT tagged and released into experimental trial arena to test the scanning success of this RePTaR reader prototype
#### Nov 11, 2021 - code written by Staci Amburgey
#### Involving data collected from August 1-August 30 2021

rm(list = ls())

library(ggplot2);library(tidyr);library(tidyverse);library(lubridate);library(LaCroixColoR)

# Read in scanning data files as items in a list
file_names1 <- dir("Scans/Station_1-2021/Station_1/")
your_data_frame1 <- as_tibble(do.call(rbind,lapply(paste("Scans/Station_1-2021/Station_1/",file_names1, sep=""),read.csv)))  # hit error with read_csv
file_names2 <- dir("Scans/Station_2-2021/Station_2/")
your_data_frame2 <- as_tibble(do.call(rbind,lapply(paste("Scans/Station_2-2021/Station_2/",file_names2, sep=""),read.csv)))  # hit error with read_csv
your_data_frame <- rbind(your_data_frame1,your_data_frame2)

# Shape data into useful form
df <- your_data_frame %>%
  separate(Timestamp, sep = ", ", into = c("Date","Time")) %>%
  filter((Time >= c("18:00:00") & Time <= c("23:59:58")) |
           (Time <= c("06:00:00") & Time >= c("00:00:00"))) %>%
  unite(DateTime, Date:Time, sep = " ") %>%
  mutate(DateTime = mdy_hms(DateTime)) %>%
  mutate(Day = mday(DateTime)) %>%
  group_by(dt=cut(DateTime, breaks= seq(
    from=as.POSIXct("2021-08-05 18:00:00", tz="UTC"),
    to=as.POSIXct("2021-08-30 06:00:00", tz="UTC"),
    by="hour"
  )))

dfsum <- df %>% 
  separate(dt, sep = " ", into = c("Date","Time")) %>%
  count(Time) %>%
  mutate(Time = as.factor(Time)) %>%
  mutate(Time = fct_relevel(Time, c("18:00:00","19:00:00","20:00:00","21:00:00","22:00:00","23:00:00",
                                    "00:00:00","01:00:00","02:00:00","03:00:00","04:00:00","05:00:00","06:00:00" )))



# p <- lacroix_palette("Lemon", n=13, type="continuous")
# p[1:13]
# 
# labs <- c("6pm","7pm","8pm","9pm","10pm","11pm","12am","1am","2am","3am","4am","5am","6am")
# colsLC <- c("#F7AA14","#F6B90B","#F5C903","#F5D523","#F6DE5F","#E4E292","#89CE9F","#2DBAAC","#14A7B3","#0C95BA","#0A7AAF","#10518C","#172869")
# mid <- mean(dfsum$n)           
# 
# # Number of scans at different times of the day
# p1 <- ggplot(dfsum, aes(x = Time, y = n, color = n, fill = n)) +
#   geom_bar(stat = "identity", alpha = 0.4) +
#   scale_fill_gradient2(midpoint=mid, low=colsLC[11], mid=colsLC[9], high=colsLC[4]) +
#   scale_color_gradient2(midpoint=mid, low=colsLC[11], mid=colsLC[9], high=colsLC[4]) +
#   ylab("Scan Count") +
#   scale_x_discrete(labels = labs) +
#   theme(legend.position = "none", panel.background = element_rect(fill = "white", color = "darkgrey"), axis.title = element_text(size = 16), axis.text = element_text(size = 12))
# 
# p1


p <- lacroix_palette("Pamplemousse", n=9, type="continuous")

labs <- c("6pm","7pm","8pm","9pm","10pm","11pm","12am","1am","2am","3am","4am","5am","6am")
colsLC <- c("#F7AA14","#F6B90B","#F5C903","#F5D523","#F6DE5F","#E4E292","#89CE9F","#2DBAAC","#14A7B3","#0C95BA","#0A7AAF","#10518C","#172869")
mid <- mean(dfsum$n) 

p1 <- ggplot(dfsum, aes(x = Time, y = n, color = n, fill = n)) +
  geom_bar(stat = "identity", alpha = 0.5) +
  scale_fill_gradient2(midpoint=mid, low=colsLC[9], mid=colsLC[4], high=colsLC[1]) +
  scale_color_gradient2(midpoint=mid, low=colsLC[9], mid=colsLC[4], high=colsLC[1]) +
  ylab("Scan Count") +
  scale_x_discrete(labels = labs) +
  theme(legend.position = "none", panel.background = element_rect(fill = "white", color = "darkgrey"), axis.title = element_text(size = 12), axis.text = element_text(size = 12))

png(file="ScanCounts.png",width=8,height=5,units="in",res=600)
p1
dev.off()



### Plotting beta estimates from JAGS models.----

modResults <- read_csv()

pSex <- ggplot()

