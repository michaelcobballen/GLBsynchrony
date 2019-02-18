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
library(tidyr)
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
# make a look-up table used to join tables based on state codes
states.lookup = 
  data.frame(states@data) %>%
  rename(state = STUSPS) %>%
  mutate(State.ANSI = as.numeric(as.character(STATEFP)))

#####
### A FUNCTION TO GENERATE PAIRWISE CORRELATIONS IN DETRENDED YIELD AMONG STATES
### FOR A GIVEN RANGE OF YEARS
#####

birdsynch = function(yearstart, yearend, spcode) {

# yearstart should be the year before the first pop change (r) value, e.g., use 1966 when r starts at 1967
# spcode is the AOU species code which you can look up here: https://www.pwrc.usgs.gov/bbl/manual/speclist.cfm
  #  spcode = s05460 for GRSP, s05010 for EAME, s04940 for BOBO, s05420 for SAVS
  # (for now: look up other species codes on BBS website and change the sp=="" to desired species)
# filtering and joining BBS data to state code data; creating "population growth rate" (r) variable
#yearstart=1990; yearend=2015; spcode="s05460" # for testing the function

bird =
  bbs %>%
  filter(sp==spcode & state!= "X" & state != "CAN" & year > yearstart-1 & year < yearend+1) %>%
  arrange(state,year) %>%
  mutate(r = ifelse(year==yearstart,NA,100*c(NA,diff(ind))/lag(ind))) %>%
  left_join(states.lookup, by = "state") %>%
  dplyr::select(sp, state, year, ind, cred, r, State.ANSI = STATEFP, name = NAME) 
bird=droplevels(bird)

# check average BBS index values to exclude states with very low populations (working criteria: mean index < 0.1)

mean.index=data.frame(mean.index = tapply(bird$ind,bird$state,mean, na.rm=T))
mean.index$state = rownames(mean.index)
mean.index$exclude = ifelse(mean.index$mean.index < 0.1, "Y", "N")
  
# put data into a 'wide format' matrix, with a column for each state and row for each year
bird.mat = 
  bird %>%
  dplyr::select(state,year,r) %>%
  spread(state,r) %>%
  filter(year > yearstart)
bird.mat = bird.mat[,2:ncol(bird.mat)]

# create matrix of Pearson correlations among states
bird.test = 
  bird.mat %>%
  corr.test(adjust="none")

bird.cor = 
  data.frame(reshape2::melt(bird.test$r), 
             reshape2::melt(bird.test$p)[3])
colnames(bird.cor) = c("state1","state2","cor","pval")
bird.cor$pairid = paste(as.character(bird.cor$state1),"-",as.character(bird.cor$state2),sep="")
bird.cor=bird.cor[bird.cor$cor!=1,]
bird.cor[230:235,]

### merge with the low-population exclusion criteria table
bird.cor2 = merge(bird.cor,mean.index,by.x=c("state1"),by.y="state")
bird.cor3 = merge(bird.cor2,mean.index,by.x=c("state2"),by.y="state")
bird.cor3$exclude = ifelse(bird.cor3$exclude.x=="Y"|bird.cor3$exclude.y=="Y","Y","N")

#### combining GRSP correlation and interstate distances
### delete duplicate state pairs & exclude low population states
  # note: 'distinct' this assumes that no two pairs of states have the same 9-digit correlation value, should check this manually to be sure; will write into code eventually
bird.cordist = bird.cor3 %>%
  left_join(dstate,by='pairid') %>%
  select(pairid, state1=state1.x, state2=state2.x, cor,pval,dist,dcat,exclude,reg) %>%
  mutate(sig = ifelse(pval < 0.05,"Y","N"), sig1 = as.numeric(ifelse(pval < 0.05,"1","0"))) %>%
  distinct(cor, .keep_all = TRUE) %>%
  filter(exclude=="N")

return(bird.cordist)
}

### NOTE: NEED TO CHECK COR VALUES TO MAKE SURE THE FUNCTION IS WORKING...
### ALSO NEED TO EXCLUDE LIKE I DID IN BIRD-YIELD-CROSS-COR (CRITERIA SEPARATE BY TIME PERIOD)
### PLOTS LOOK DIFFERENT - INCLUDING THE BINNED ONE

### SCATTERPLOT OF SYNCHRONY VS. DISTANCE
#  spcode = s05460 for GRSP, s05010 for EAME, s04940 for BOBO, s05420 for SAVS

bird.67.90 = cbind.data.frame(birdsynch(1966,1990,"s04940"), per = "1967-1990")
bird.91.15 = cbind.data.frame(birdsynch(1990,2015,"s04940"), per = "1991-2015")
bird.67.15 = cbind.data.frame(birdsynch(1966,2015,"s04940"), per = "1967-2015")

X11(13,9)
ggplot(bird.67.90) + 
  geom_point(aes(x=dist,y=cor,color=sig),size=4,alpha=0.5) + 
  geom_smooth(aes(x=dist,y=cor),method="loess",color="black") 
#ggsave("figures/sync_scatter/SAVS_cor_all_1967-1990.png")

X11(13,9)
ggplot(bird.91.15) + 
  geom_point(aes(x=dist,y=cor,color=sig),size=4,alpha=0.5) + 
  geom_smooth(aes(x=dist,y=cor),method="loess",color="black") 
#ggsave("figures/sync_scatter/SAVS_cor_all_1991-2015.png")

X11(13,9)
ggplot(bird.67.15) + 
  geom_point(aes(x=dist,y=cor,color=sig),size=4,alpha=0.5) + 
  geom_smooth(aes(x=dist,y=cor),method="loess",color="black") #+
  #xlim(c(0,3000))
#ggsave("figures/sync_scatter/zSAVS_cor_all_1967-2015.png")

### a better way to do it: rbinding together to make it easier to plot in ggplot
bird2time = rbind(bird.67.90,bird.91.15)

X11(13,9)
ggplot(bird2time) +
  #  geom_point(aes(x = dist, y = cor, color=per, shape = sig), size=2, alpha = 0.5) +
  scale_shape_manual(values = c(16,1)) +
  scale_color_manual(values = c("#de2d26","#2b8cbe")) +
  geom_smooth(aes(x = dist, y = cor, color=per), method = "loess", size = 2) +
  labs(title = "Spatial Synchrony in SAVS Pop. Growth Rate")
  #xlim(c(0,2500))
#ggsave("figures/temp_change_sync/SAVS_synchrony.1967.vs.2015.nopts.png")


######
###### plotting % signficant correlations by distance bins 
######

binlabels = c("<100","100-300","301-500","501-700","701-900",
              "901-1100","1101-1300","1301-1500",">1500")
corbin = cbind.data.frame(pctsig=100*tapply(bird.67.90$sig1,bird.67.90$dcat,mean,na.rm=T),
      n=tapply(bird.67.90$sig1,bird.67.90$dcat,length))
corbin$bin=factor(rownames(corbin),levels=binlabels)

#x11(13,8.5)
ggplot(corbin) + 
  geom_col(aes(y=pctsig,x=bin)) + 
  labs(x="distance (km)",y="% of correlations with P < 0.05") +
  ylim(c(0,35))
#ggsave("figures/bin_cor/zBOBO_bincor_all_1967-2015.png")
#ggsave("figures/bin_cor/BOBO_bincor_all_1967-1990.png")
#ggsave("figures/bin_cor/BOBO_bincor_all_1991-2015.png")


## NEXT: need to create n = 1000 null (random) distributions for each bin like Martin et al. 
## to test spatial extent of synchrony




###########################################################################
#####... CODE BELOW HERE NO LONGER ACTIVE
###########################################################################


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
         

# NEXT: REPEAT FOR 3 DIFF TIME STEPS (1966-1981, 1981-1996, 1997-2015), AND PUT
# AND PUT LOESS LINES ON THE SAME GRAPH



# checking all time series data graphically (both index and growth rate)

ggplot(bird) + 
  geom_line(aes(x=year,y=ind),size = 1) +
  facet_wrap(~state, scales = "free_y")

ggplot(bird) + 
  geom_line(aes(x=year,y=r),size = 1) +
  facet_wrap(~state, scales = "free_y")
       