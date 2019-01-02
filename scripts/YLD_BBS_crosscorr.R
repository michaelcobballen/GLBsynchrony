### THIS CODE IS FOR CALCULATING RELATIONSHIP BETWEEN DETRENDED HAY YIELD (TONS/ACRE/YEAR; USDA NASS DATA) 
### AND GRASSLAND BIRD POPULATION GROWTH RATE (BREEDING BIRD SURVEY) IN US STATES 
### STATE-LEVEL DATA FOR NOW; COULD BE DONE AT ~COUNTY AND ROUTE-LEVEL

library(ggplot2)
library(dplyr)
library(tidyr)
library(psych)
library(cowplot)
library(RColorBrewer)
library(rgdal)
library(sp)
library(broom)
library(geofacet)

### READ IN AND PROCESS DATA

# census shapefile for US states and territories
states <- readOGR("data/cb_2014_us_state_5m.shp")
# downloaded from: https://www2.census.gov/geo/tiger/GENZ2014/shp/cb_2014_us_state_5m.zip
states.latlong = spTransform(states,CRS("+proj=longlat +datum=WGS84"))
# make a look-up table used to join tables based on state codes
states.lookup = 
  data.frame(states@data) %>%
  rename(state = STUSPS) %>%
  mutate(State.ANSI = as.numeric(as.character(STATEFP)))

# Breeding Bird Survey data (annual state-level population index)
bbs = read.csv("data/BBS_Annual_Indices_Estimates_2015_7-29-2016.csv") # read in BBS data
colnames(bbs)<-c("sp","state","year","ind","cred") # rename columns
#rename region codes to match state abbreviations
levels(bbs$state)=c("AL","CAN","AZ","AR","CAN","CA","X","X","CO","CT", "DE", "X", "FL", "GA", "IA",
                    "ID", "IL", "IN", "KS", "KY", "LA", "CAN", "MA", "MD", "ME", "MI", "MN", "MS", 
                    "MO", "MT", "CAN", "NC", "ND", "NE", "NV", "NH", "NJ", "NM", "CAN", "NY", "OH",
                    "OK", "CAN", "OR", "PA", "CAN", "CAN", "RI", "X", "X", "X", "X", "X", "X", "X",
                    "X", "X", "X", "X", "X", "X", "X", "X", "X", "X", "X", "X", "X", "X", "X", "X",
                    "X", "X", "X", "X", "X", "X", "X", "X", "X", "X", "CAN", "SC", "SD", "X", "TN", 
                    "TX", "X", "UT", "VA", "VT", "WA", "X", "WI", "WV", "WY")

# filtering and joining BBS data to state code data; creating "population growth rate" (r) variable
### TO CHANGE SPECIES AND YEARS, ALTER THE CODE IN THE "filter" COMMAND HERE
# (for now: look up species codes on BBS website and change the sp=="" to desired species)
# create "statelagyear" variable to join with yield data in the proper time lag structure (see below)
bird =
  bbs %>%
  filter(sp=="s05460" & state!= "X" & state != "CAN" & year>1965 & year<2016) %>%
  arrange(state,year) %>%
  mutate(r = ifelse(year==1966,NA,100*c(NA,diff(ind))/lag(ind))) %>%
  left_join(states.lookup, by = "state") %>%
  select(sp, state, year, ind,cred, r, State.ANSI = STATEFP, name = NAME) %>%
  mutate(statelagyear = paste(state,year,sep=""))
bird=droplevels(bird)

# Hay yield data for US states from USDA NASS ('all hay') in tons / acre
yld = read.csv("data/AllHay_Yield_survey_allstates.csv")
# filtering data and joining with state code data
# creating detrended yield variable (ryld) by regressing against year (within states)
# subsetted 1966-2014 to relate "lag 1" yield with "lag 0" population growth rate  

### SOMETHING IS TERRIBLY WRONG WITH DETR.YLD RIGHT NOW!!!!
# a joining problem
# fix states.lookup
detr.yld = 
  yld %>% 
  filter(Year>1965 & Year<2015 & Period=="YEAR") %>%
  group_by(State) %>%
  mutate(ryld = resid(lm(Value~Year))) %>%
  arrange(State,Year) %>%
  select(State,State.ANSI,Year,Value,ryld)

detr.yld = 
  detr.yld %>%
  left_join(states.lookup, by = "State.ANSI") %>%
  select(State, state, year = Year, yield = Value, ryld, State.ANSI = STATEFP, name = NAME) %>%
  mutate(lagyear = year + 1, statelagyear = paste(state,lagyear,sep=""), State.ANSI = as.numeric(State.ANSI))

# joining the bird data with the yield data (lagged 1 year)

r.ryld = 
  left_join(bird,detr.yld,by = "statelagyear") %>%
  select(sp,State, state = state.x, ansi = State.ANSI.x, name = name.x, 
         yearbird = year.x, yearyld = year.y,ind,r,yield.lag1 = yield,ryld.lag1 = ryld)


# checking all time series data graphically: BBS index, growth rate, yield, and detrended yield

#X11(26,18)
ggplot(r.ryld) + 
  geom_line(aes(x=yearbird,y=ind)) +
  facet_geo(~state, scales = "free_y") +
  theme_bw() +
  xlab("Year") + ylab("BBS Index") +
  theme(axis.text=element_text(size=3.25))
#ggsave("figures/GRSP_all_states_BBSindex_US.png")

#X11(26,18)
ggplot(r.ryld) + 
  geom_line(aes(x=yearbird,y=r,group=state)) +
  facet_wrap(~state, scales = "free_y") +
  theme_bw() +
  xlab("Year") + ylab("BBS Index Growth Rate") +
  theme(axis.text=element_text(size=5))
#ggsave("figures/GRSP_all_states_r.png")

#X11(26,18)
ggplot(r.ryld) + 
  geom_line(aes(x=yearyld,y=yield.lag1,group=state)) + 
  geom_line(aes(x=yearyld,y=ryld.lag1)) +
  facet_wrap(~state, scales = "free_y") +
  theme_bw() + 
  theme(axis.text=element_text(size=5))
#ggsave("figures/YLD_all_states_detrended.png")


# check average BBS index values to identify states with very low populations
# working criteria: exclude states with mean index < 0.1 

mean.index=data.frame(mean.index = tapply(bird$ind,bird$state,mean, na.rm=T))
mean.index$state = rownames(mean.index)
mean.index$exclude = ifelse(mean.index$mean.index < 0.1, "Y", "N")

# calculate lm for r ~ ryld in each state
# from: https://stackoverflow.com/questions/22713325/fitting-several-regression-models-with-dplyr
head(all)
allsimple = all %>% filter(is.na(r)!=T,is.na(ryld.lag1)!=T) %>% select(state,yearbird,r,ryld.lag1)
allsimple
allstate = allsimple %>% group_by(state) %>% do(fitstate = lm(r ~ ryld.lag1, na.action=na.exclude, data = .))

(r.yld.lm = tidy(allstate, fitstate))



# put data into a 'wide format' matrix, with a column for each state and row for each year
yld.mat = 
  detr.yld %>%
  select(State,Year,ryld) %>%
  spread(State,ryld)

# name state columns by USPS postal abbreviations
colnames(yld.mat) = c("year","AL","AK","AZ","AR","CA","CO","CT","DE","FL","GA","ID","IL",
                     "IN","IA","KS","KY","LA","ME","MD","MA","MI","MN","MS","MO","MT",
                     "NE","NV","NH","NJ","NM","NY","NC","ND","OH","OK","OR","PA","RI",
                     "SC","SD","TN","TX","UT","VT","VA","WA","WV","WI","WY")

yld.names = reshape2::melt(yld.mat[,2:ncol(yld.mat)])
yld.names$year = rep(2016:1966,49); head(yld.names)

ggplot(yld.names) + 
  geom_line(aes(x=year,y=value)) +
  facet_wrap(~variable, scales = "free_y") 
# check this against the facet plot above to make sure it worked (years were in order)



# progress ...




# create matrix of Pearson correlations among states
yld.test = 
  yld.mat %>%
  corr.test(adjust="none")
  
yld.cor = 
  data.frame(reshape2::melt(yld.test$r), 
             reshape2::melt(yld.test$p)[3])
colnames(yld.cor) = c("state1","state2","cor","pval"); head(yld.cor)
yld.cor$pairid = paste(as.character(yld.cor$state1),"-",as.character(yld.cor$state2),sep=""); head(yld.cor); nrow(yld.cor)
yld.cor=yld.cor[yld.cor$cor!=1,]; head(yld.cor); nrow(yld.cor); max(yld.cor$cor)
yld.cor[230:235,] 

#### combining yield correlations and interstate distances
head(yld.cor); head(dstate)

yldcordist = merge(yld.cor,dstate,by='pairid',all.x=T); head(yldcordist); nrow(yldcordist)
colnames(yldcordist)[c(2,3)] <- c("state1","state2"); head(yldcordist)
yldcordist = yldcordist[,c("pairid","state1","state2","cor","pval","dist","dcat","reg")]; head(yldcordist)
yldcordist$sig = ifelse(yldcordist$pval<0.05,"Y","N"); head(yldcordist); nrow(yldcordist)

### delete duplicate state pairs
nrow(yldcordist) == length(unique(yldcordist$cor))*2 # this tests that all values of "cor" are unique
  # if that returns "TRUE", then proceed (if not you need to find another way of removing duplicates)
yldcordist = distinct(yldcordist, cor, .keep_all = TRUE); head(yldcordist); nrow(yldcordist)

## exclude Alaska and... 
## make subsets of the two main geographic areas (east and west of the Mississippi River)
yldcordist = filter(yldcordist, state1!="AK" & state2!="AK"); nrow(yldcordist)
yldcordist.east = filter(yldcordist,reg=="east") # 325 rows
  # cordist.east = droplevels(cordist.east)
yldcordist.west = filter(yldcordist,reg=="west")
  # cordist.west = droplevels(cordist.west) # 231 rows

# to change geographic region plotted below:
# change following code to cordist.east (for east), 
# cordist.west (for west), 
# or cordist (for all)
yldreg = yldcordist.east


### SCATTERPLOT OF SYNCHRONY VS. DISTANCE
#X11(13,9)
ggplot(yldreg) + 
  geom_point(aes(x=dist,y=cor,color=sig),size=4) + 
  geom_smooth(aes(x=dist,y=cor),method="loess",color="black") 
#ggsave("figures/YLD_cor_all_1966-2015_loess.png")
#ggsave("figures/YLD_cor_east_1966-2015.png")
#ggsave("figures/YLD_cor_west_1966-2015.png")

# with significant pairs labeled
ggplot(yldreg) + 
  geom_point(aes(x=dist,y=cor,color=sig),size=4) + 
  geom_smooth(aes(x=dist,y=cor),method="lm",color="black") +
  geom_label(data=subset(yldreg,sig=="Y"),aes(x=dist,y=cor,label=pairid))


###### plotting % signficant by distance bins 

yldreg$sig1=ifelse(yldreg$sig=="Y",1,0)
binlabels = c("<100","100-300","301-500","501-700","701-900",
              "901-1100","1101-1300","1301-1500",">1500")
corbin = cbind.data.frame(pctsig=100*tapply(yldreg$sig1,yldreg$dcat,mean,na.rm=T),
                          n=tapply(yldreg$sig1,yldreg$dcat,length))
corbin$bin=factor(rownames(corbin),levels=binlabels); corbin

#x11(13,8.5)
ggplot(corbin) + 
  geom_col(aes(y=pctsig,x=bin)) + 
  labs(x="distance (km)",y="% of correlations that are signifcant") 
#ggsave("figures/YLD_bincor_all_1966-2015.png")
#ggsave("figures/YLD_bincor_east_1966-2015.png")
#ggsave("figures/YLD_bincor_west_1966-2015.png")


## NEXT: need to create n = 1000 null (random) distributions for each bin like Martin et al. 
## to test spatial extent of synchrony



##### MAPPING AVERAGE LOCAL SYNCHRONY (EAST ONLY FOR NOW)

yldreg$corres = resid(lm(cor~dist,data=yldreg)) 

state.abr = unique(dstate$state1)

# could make all this into a for loop?
# something like: for(i in 1:length(state.abr))  {

yldreg$AL = ifelse(yldreg$state1 %in% "AL"|yldreg$state2 %in% "AL","AL",""); head(yldreg)
yldreg$DE = ifelse(yldreg$state1 %in% "DE"|yldreg$state2 %in% "DE","DE",""); head(yldreg)
yldreg$GA = ifelse(yldreg$state1 %in% "GA"|yldreg$state2 %in% "GA","GA",""); head(yldreg)
yldreg$IL = ifelse(yldreg$state1 %in% "IL"|yldreg$state2 %in% "IL","IL",""); head(yldreg)
yldreg$IN = ifelse(yldreg$state1 %in% "IN"|yldreg$state2 %in% "IN","IN",""); head(yldreg)
yldreg$KY = ifelse(yldreg$state1 %in% "KY"|yldreg$state2 %in% "KY","KY",""); head(yldreg)
yldreg$MD = ifelse(yldreg$state1 %in% "MD"|yldreg$state2 %in% "MD","MD",""); head(yldreg)
yldreg$MI = ifelse(yldreg$state1 %in% "MI"|yldreg$state2 %in% "MI","MI",""); head(yldreg)
yldreg$NC = ifelse(yldreg$state1 %in% "NC"|yldreg$state2 %in% "NC","NC",""); head(yldreg)
yldreg$NJ = ifelse(yldreg$state1 %in% "NJ"|yldreg$state2 %in% "NJ","NJ",""); head(yldreg)
yldreg$NY = ifelse(yldreg$state1 %in% "NY"|yldreg$state2 %in% "NY","NY",""); head(yldreg)
yldreg$OH = ifelse(yldreg$state1 %in% "OH"|yldreg$state2 %in% "OH","OH",""); head(yldreg)
yldreg$PA = ifelse(yldreg$state1 %in% "PA"|yldreg$state2 %in% "PA","PA",""); head(yldreg)
yldreg$SC = ifelse(yldreg$state1 %in% "SC"|yldreg$state2 %in% "SC","SC",""); head(yldreg)
yldreg$TN = ifelse(yldreg$state1 %in% "TN"|yldreg$state2 %in% "TN","TN",""); head(yldreg)
yldreg$VA = ifelse(yldreg$state1 %in% "VA"|yldreg$state2 %in% "VA","VA",""); head(yldreg)
yldreg$WI = ifelse(yldreg$state1 %in% "WI"|yldreg$state2 %in% "WI","WI",""); head(yldreg)
yldreg$WV = ifelse(yldreg$state1 %in% "WV"|yldreg$state2 %in% "WV","WV",""); head(yldreg)

yldmeansynch = 
  rbind.data.frame(
    yldreg %>%
      filter(dist<500,AL=="AL") %>%
      summarise(avg = mean(corres)) %>%
      arrange(avg) %>%
      cbind(state="AL")
    ,
    yldreg %>%
      filter(dist<500,DE=="DE") %>%
      summarise(avg = mean(corres)) %>%
      arrange(avg) %>%
      cbind(state="DE")
    ,
    yldreg %>%
      filter(dist<500,GA=="GA") %>%
      summarise(avg = mean(corres)) %>%
      arrange(avg) %>%
      cbind(state="GA")
    ,
    yldreg %>%
      filter(dist<500,IL=="IL") %>%
      summarise(avg = mean(corres)) %>%
      arrange(avg) %>%
      cbind(state="IL")
    ,
    yldreg %>%
      filter(dist<500,IN=="IN") %>%
      summarise(avg = mean(corres)) %>%
      arrange(avg) %>%
      cbind(state="IN")
    ,
    yldreg %>%
      filter(dist<500,KY=="KY") %>%
      summarise(avg = mean(corres)) %>%
      arrange(avg) %>%
      cbind(state="KY")
    ,
    yldreg %>%
      filter(dist<500,MD=="MD") %>%
      summarise(avg = mean(corres)) %>%
      arrange(avg) %>%
      cbind(state="MD")
    ,
    yldreg %>%
      filter(dist<500,MI=="MI") %>%
      summarise(avg = mean(corres)) %>%
      arrange(avg) %>%
      cbind(state="MI")
    ,
    yldreg %>%
      filter(dist<500,NC=="NC") %>%
      summarise(avg = mean(corres)) %>%
      arrange(avg) %>%
      cbind(state="NC")
    ,
    yldreg %>%
      filter(dist<500,NJ=="NJ") %>%
      summarise(avg = mean(corres)) %>%
      arrange(avg) %>%
      cbind(state="NJ")
    ,
    yldreg %>%
      filter(dist<500,NY=="NY") %>%
      summarise(avg = mean(corres)) %>%
      arrange(avg) %>%
      cbind(state="NY")
    ,
    yldreg %>%
      filter(dist<500,OH=="OH") %>%
      summarise(avg = mean(corres)) %>%
      arrange(avg) %>%
      cbind(state="OH")
    ,
    yldreg %>%
      filter(dist<500,PA=="PA") %>%
      summarise(avg = mean(corres)) %>%
      arrange(avg) %>%
      cbind(state="PA")
    ,
    yldreg %>%
      filter(dist<500,SC=="SC") %>%
      summarise(avg = mean(corres)) %>%
      arrange(avg) %>%
      cbind(state="SC")
    ,
    yldreg %>%
      filter(dist<500,TN=="TN") %>%
      summarise(avg = mean(corres)) %>%
      arrange(avg) %>%
      cbind(state="TN")
    ,
    yldreg %>%
      filter(dist<500,VA=="VA") %>%
      summarise(avg = mean(corres)) %>%
      arrange(avg) %>%
      cbind(state="VA")
    ,
    yldreg %>%
      filter(dist<500,WI=="WI") %>%
      summarise(avg = mean(corres)) %>%
      arrange(avg) %>%
      cbind(state="WI")
    ,
    yldreg %>%
      filter(dist<500,WV=="WV") %>%
      summarise(avg = mean(corres)) %>%
      arrange(avg) %>%
      cbind(state="WV")
  )
head(yldmeansynch); nrow(yldmeansynch)

yld.mean.synch.map = merge(states.latlong,yldmeansynch, by.y="state",by.x="STUSPS")


my.palette <- brewer.pal(n = 7, name = "YlOrRd")

#x11(11,11)
yld.mean.synch.map %>%
  subset(STUSPS %in% c("AL","DE","GA","IL","IN","KY","MD","MI","NC","NJ","NY",
                       "OH","PA","SC","TN","VA","WI","WV")) %>%
  spplot(zcol="avg",col.regions = my.palette, cuts = 6, 
         main = "Hay Yield Synchrony w/in 500 km, 1966-2015", col="transparent")


# NEXT: REPEAT FOR 3 DIFF TIME STEPS (1966-1981, 1981-1996, 1997-2015), AND PUT
  # AND PUT LOESS LINES ON THE SAME GRAPH