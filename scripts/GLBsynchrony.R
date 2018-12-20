### THIS CODE IS FOR CALCULATING PAIRWISE 'SYNCHRONY' (CROSS-CORRELATION) 
### AMONG BIRD POPULATION GROWTH RATES IN US STATES (BREEDING BIRD SURVEY DATA)
### FOR NOW I AM ONLY USING IT FOR GRASSLAND BIRDS AT THE STATE LEVEL
### SPECIES CAN BE CHANGED BY ALTERING THE AOU SPECIES CODE ON LINE 37
### SOME DAY I WILL ALSO MAKE A SCRIPT FOR EXAMINING SYNCHRONY AT THE ROUTE (OR AGGREGATED ROUTE) LEVEL

library(ggplot2)
library(utils)
library(sp)
library(rgeos)
library(downloader)
library(maptools)
library(mapmisc)
library(sgeostat)
library(psych)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(rgdal)
library(sp)

# READ IN DATA
dstate=read.csv("data/dstate.csv")
bbs=read.csv("data/BBS_Annual_Indices_Estimates_2015_7-29-2016.csv") # read in BBS data
colnames(bbs)<-c("sp","state","year","ind","cred")
levels(bbs$state)
levels(bbs$state)=c("AL","CAN","AZ","AR","CAN","CA","X","X","CO","CT", "DE", "X", "FL", "GA", "IA",
                   "ID", "IL", "IN", "KS", "KY", "LA", "CAN", "MA", "MD", "ME", "MI", "MN", "MS", 
                   "MO", "MT", "CAN", "NC", "ND", "NE", "NV", "NH", "NJ", "NM", "CAN", "NY", "OH",
                   "OK", "CAN", "OR", "PA", "CAN", "CAN", "RI", "X", "X", "X", "X", "X", "X", "X",
                   "X", "X", "X", "X", "X", "X", "X", "X", "X", "X", "X", "X", "X", "X", "X", "X",
                   "X", "X", "X", "X", "X", "X", "X", "X", "X", "X", "CAN", "SC", "SD", "X", "TN", 
                   "TX", "X", "UT", "VA", "VT", "WA", "X", "WI", "WV", "WY")

# census shapefile for US states and territories
states <- readOGR("data/cb_2014_us_state_5m.shp")
# downloaded from: https://www2.census.gov/geo/tiger/GENZ2014/shp/cb_2014_us_state_5m.zip
states.latlong = spTransform(states,CRS("+proj=longlat +datum=WGS84"))

### SUBSET BBS DATA TO THE DESIRED SPECIES 
# (for now, look up species codes on BBS website and change the sp=="" to desired species)
bird=subset(bbs,sp=="s05460" & state!= "X" &state != "CAN"); head(bird)

# creating population growth rate (% change) time series from index values
# (note this only works because years are in order within states)
bird$r = 100*c(NA,diff(bird$ind))/lag(bird$ind)
bird$r = as.numeric(ifelse(bird$year==1966,"NA",bird$r))
bird=droplevels(bird)
head(bird,10)


# checking all time series data graphically (both index and growth rate)

ggplot(bird) + 
  geom_line(aes(x=year,y=ind),size = 1) +
  facet_wrap(~state, scales = "free_y")
 
ggplot(bird) + 
  geom_line(aes(x=year,y=r),size = 1) +
  facet_wrap(~state, scales = "free_y")


# check average index values to identify states with very low populations 
# working criteria: exclude states with mean index < 0.1 
# also check number of years where index < 0.1

mean.index=data.frame(mean.index = tapply(bird$ind,bird$state,mean, na.rm=T))
  mean.index$state = rownames(mean.index)
  mean.index$exclude = ifelse(mean.index$mean.index < 0.1, "Y", "N")
  
tapply(bird$ind<0.1,bird$state,sum, na.rm=T)

# reshaping data to wide format for correlations
bird.r.matrix=reshape(bird[c("state","year", "r")],
                      idvar="year",timevar = "state", direction = "wide")

colnames(bird.r.matrix) <- substr(colnames(bird.r.matrix),3,4)
bird.r.matrix = bird.r.matrix[2:ncol(bird.r.matrix)]
head(bird.r.matrix); tail(bird.r.matrix); ncol(bird.r.matrix)

bird.corr.test=corr.test(bird.r.matrix, adjust = "none")
bird.cor = data.frame(reshape2::melt(bird.corr.test$r), reshape2::melt(bird.corr.test$p)[3]); head(bird.cor); nrow(bird.cor)
colnames(bird.cor) = c("state1","state2","cor","pval"); head(bird.cor)
bird.cor$pairid = paste(as.character(bird.cor$state1),"-",as.character(bird.cor$state2),sep=""); head(bird.cor); max(bird.cor$cor);nrow(bird.cor)
bird.cor=bird.cor[bird.cor$cor!=1,]; head(bird.cor); nrow(bird.cor)
bird.cor[230:235,]

### merge with the low-population exclusion criteria table
bird.cor2 = merge(bird.cor,mean.index,by.x=c("state1"),by.y="state"); head(bird.cor2); nrow(bird.cor2)
bird.cor3 = merge(bird.cor2,mean.index,by.x=c("state2"),by.y="state"); head(bird.cor3); nrow(bird.cor3)
bird.cor3$exclude = ifelse(bird.cor3$exclude.x=="Y"|bird.cor3$exclude.y=="Y","Y","N")
head(bird.cor3)

#### combining GRSP correlation and interstate distances
head(bird.cor3); nrow(bird.cor3)
head(dstate); nrow(dstate)

cordist = merge(bird.cor3,dstate,by='pairid',all.x=T); head(cordist); nrow(cordist)
colnames(cordist)[c(2,3)] <- c("state1","state2"); head(cordist)
cordist = cordist[,c("pairid","state1","state2","cor","pval","dist","dcat","exclude","reg")]; head(cordist)
cordist$sig = ifelse(cordist$pval<0.05,"Y","N"); head(cordist); nrow(cordist)

### delete duplicate state pairs
cordist = distinct(cordist, cor, .keep_all = TRUE); head(cordist); nrow(cordist)

## exclude states with very low populations and 
  ## ## make subsets of the two main geographic areas (east and west of the Mississippi River)
cordist = filter(cordist, exclude=="N"); nrow(cordist)
cordist.east = filter(cordist,reg=="east"); head(cordist.east); nrow(cordist.east)
  #cordist.east = droplevels(cordist.east)
cordist.west = filter(cordist,reg=="west"); head(cordist.west); nrow(cordist.west)
  #cordist.west = droplevels(cordist.west) 

# to change geographic region plotted below:
# change following code to cordist.east (for east), 
  # cordist.west (for west), 
  # or cordist (for all)
region = cordist.east

### SCATTERPLOT OF SYNCHRONY VS. DISTANCE
#X11(13,9)
ggplot(region) + 
  geom_point(aes(x=dist,y=cor,color=sig),size=4) + 
  geom_smooth(aes(x=dist,y=cor),method="lm",color="black") 
#ggsave("figures/GRSP_cor_all_1966-2015.png")
#ggsave("figures/GRSP_cor_east_1966-2015.png")
#ggsave("figures/GRSP_cor_west_1966-2015.png")

# with significant pairs labeled
ggplot(region) + 
  geom_point(aes(x=dist,y=cor,color=sig),size=4) + 
  geom_smooth(aes(x=dist,y=cor),method="lm",color="black") +
  geom_label(data=subset(region,sig=="Y"),aes(x=dist,y=cor,label=pairid))


###### plotting % signficant by distance bins 

region$sig1=ifelse(region$sig=="Y",1,0)
head(region)
binlabels = c("<100","100-300","301-500","501-700","701-900",
              "901-1100","1101-1300","1301-1500",">1500")
corbin = cbind.data.frame(pctsig=100*tapply(region$sig1,region$dcat,mean,na.rm=T),
      n=tapply(region$sig1,region$dcat,length))
corbin$bin=factor(rownames(corbin),levels=binlabels); corbin

#x11(13,8.5)
ggplot(corbin) + 
  geom_col(aes(y=pctsig,x=bin)) + 
  labs(x="distance (km)",y="% of correlations that are signifcant") 
#ggsave("figures/GRSP_bincor_all_1966-2015.png")
#ggsave("figures/GRSP_bincor_east_1966-2015.png")
#ggsave("figures/GRSP_bincor_west_1966-2015.png")


## NEXT: need to create n = 1000 null (random) distributions for each bin like Martin et al. 
## to test spatial extent of synchrony



##### MAPPING AVERAGE LOCAL SYNCHRONY

region$corres = resid(lm(cor~dist,data=region)) 

state.abr = unique(dstate$state1)

# sometime I could make all this into a for loop?
# for(i in 1:length(state.abr))  {

region$AL = ifelse(region$state1 %in% "AL"|region$state2 %in% "AL","AL",""); head(region)
region$DE = ifelse(region$state1 %in% "DE"|region$state2 %in% "DE","DE",""); head(region)
region$GA = ifelse(region$state1 %in% "GA"|region$state2 %in% "GA","GA",""); head(region)
region$IL = ifelse(region$state1 %in% "IL"|region$state2 %in% "IL","IL",""); head(region)
region$IN = ifelse(region$state1 %in% "IN"|region$state2 %in% "IN","IN",""); head(region)
region$KY = ifelse(region$state1 %in% "KY"|region$state2 %in% "KY","KY",""); head(region)
region$MD = ifelse(region$state1 %in% "MD"|region$state2 %in% "MD","MD",""); head(region)
region$MI = ifelse(region$state1 %in% "MI"|region$state2 %in% "MI","MI",""); head(region)
region$NC = ifelse(region$state1 %in% "NC"|region$state2 %in% "NC","NC",""); head(region)
region$NJ = ifelse(region$state1 %in% "NJ"|region$state2 %in% "NJ","NJ",""); head(region)
region$NY = ifelse(region$state1 %in% "NY"|region$state2 %in% "NY","NY",""); head(region)
region$OH = ifelse(region$state1 %in% "OH"|region$state2 %in% "OH","OH",""); head(region)
region$PA = ifelse(region$state1 %in% "PA"|region$state2 %in% "PA","PA",""); head(region)
region$SC = ifelse(region$state1 %in% "SC"|region$state2 %in% "SC","SC",""); head(region)
region$TN = ifelse(region$state1 %in% "TN"|region$state2 %in% "TN","TN",""); head(region)
region$VA = ifelse(region$state1 %in% "VA"|region$state2 %in% "VA","VA",""); head(region)
region$WI = ifelse(region$state1 %in% "WI"|region$state2 %in% "WI","WI",""); head(region)
region$WV = ifelse(region$state1 %in% "WV"|region$state2 %in% "WV","WV",""); head(region)

meansynch = 
rbind.data.frame(
region %>%
  filter(dist<500,AL=="AL") %>%
  summarise(avg = mean(corres)) %>%
  arrange(avg) %>%
  cbind(state="AL")
,
region %>%
  filter(dist<500,DE=="DE") %>%
  summarise(avg = mean(corres)) %>%
  arrange(avg) %>%
  cbind(state="DE")
,
region %>%
  filter(dist<500,GA=="GA") %>%
  summarise(avg = mean(corres)) %>%
  arrange(avg) %>%
  cbind(state="GA")
,
region %>%
  filter(dist<500,IL=="IL") %>%
  summarise(avg = mean(corres)) %>%
  arrange(avg) %>%
  cbind(state="IL")
,
region %>%
  filter(dist<500,IN=="IN") %>%
  summarise(avg = mean(corres)) %>%
  arrange(avg) %>%
  cbind(state="IN")
,
region %>%
  filter(dist<500,KY=="KY") %>%
  summarise(avg = mean(corres)) %>%
  arrange(avg) %>%
  cbind(state="KY")
,
region %>%
  filter(dist<500,MD=="MD") %>%
  summarise(avg = mean(corres)) %>%
  arrange(avg) %>%
  cbind(state="MD")
,
region %>%
  filter(dist<500,MI=="MI") %>%
  summarise(avg = mean(corres)) %>%
  arrange(avg) %>%
  cbind(state="MI")
,
region %>%
  filter(dist<500,NC=="NC") %>%
  summarise(avg = mean(corres)) %>%
  arrange(avg) %>%
  cbind(state="NC")
,
region %>%
  filter(dist<500,NJ=="NJ") %>%
  summarise(avg = mean(corres)) %>%
  arrange(avg) %>%
  cbind(state="NJ")
,
region %>%
  filter(dist<500,NY=="NY") %>%
  summarise(avg = mean(corres)) %>%
  arrange(avg) %>%
  cbind(state="NY")
,
region %>%
  filter(dist<500,OH=="OH") %>%
  summarise(avg = mean(corres)) %>%
  arrange(avg) %>%
  cbind(state="OH")
,
region %>%
  filter(dist<500,PA=="PA") %>%
  summarise(avg = mean(corres)) %>%
  arrange(avg) %>%
  cbind(state="PA")
,
region %>%
  filter(dist<500,SC=="SC") %>%
  summarise(avg = mean(corres)) %>%
  arrange(avg) %>%
  cbind(state="SC")
,
region %>%
  filter(dist<500,TN=="TN") %>%
  summarise(avg = mean(corres)) %>%
  arrange(avg) %>%
  cbind(state="TN")
,
region %>%
  filter(dist<500,VA=="VA") %>%
  summarise(avg = mean(corres)) %>%
  arrange(avg) %>%
  cbind(state="VA")
,
region %>%
  filter(dist<500,WI=="WI") %>%
  summarise(avg = mean(corres)) %>%
  arrange(avg) %>%
  cbind(state="WI")
,
region %>%
  filter(dist<500,WV=="WV") %>%
  summarise(avg = mean(corres)) %>%
  arrange(avg) %>%
  cbind(state="WV")
)
head(meansynch); nrow(meansynch)

mean.synch.map = merge(states.latlong,meansynch, by.y="state",by.x="STUSPS")


my.palette <- brewer.pal(n = 7, name = "YlOrRd")

#x11(11,11)
mean.synch.map %>%
  subset(STUSPS %in% c("AL","DE","GA","IL","IN","KY","MD","MI","NC","NJ","NY",
                       "OH","PA","SC","TN","VA","WI","WV")) %>%
  spplot(zcol="avg",col.regions = my.palette, cuts = 6, 
         main = "GRSP Synchrony w/in 500 km, 1966-2015", col="transparent")
         


       