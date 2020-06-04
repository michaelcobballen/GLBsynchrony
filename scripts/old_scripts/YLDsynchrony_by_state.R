### THIS CODE IS FOR CALCULATING PAIRWISE 'SYNCHRONY' (CROSS-CORRELATION) 
### AMONG DETRENDED HAY YIELD (TONS/ACRE/YEAR ) IN US STATES (USDA NASS DATA)
### FOR NOW I AM ONLY USING YIELD DATA AT THE STATE LEVEL
### CODE COULD BE MODIFIED FOR OTHER SPATIAL, 'TREND-STATIONARY' TIME SERIES
### SOME DAY I WILL ALSO MAKE A SCRIPT FOR EXAMINING SYNCHRONY AT THE COUNTY LEVEL

library(ggplot2)
library(dplyr)
library(tidyr)
library(rgdal)
library(sp)
library(psych)
#library(cowplot)
library(RColorBrewer)

### SET YEAR RANGE
yearstart = 1952
yearend = 1972


### READ IN DATA
h = subset(read.csv("data/AllHay_Yield_survey_allstates.csv"),Period=="YEAR"); head(h) # 5241
dstate = read.csv("data/dstate.csv")
# census shapefile for US states and territories
states <- readOGR("data/cb_2014_us_state_5m.shp")
# downloaded from: https://www2.census.gov/geo/tiger/GENZ2014/shp/cb_2014_us_state_5m.zip
states.latlong = spTransform(states,CRS("+proj=longlat +datum=WGS84"))

### A FUNCTION TO GENERATE PAIRWISE CORRELATIONS IN DETRENDED YIELD AMONG STATES
### FOR A GIVEN RANGE OF YEARS
yldsynch = function(yearstart, yearend) {

### detrend yield data by regressing against year (within state)
detr.yld = 
  h %>% 
  filter(Year > yearstart - 1 & Year < yearend+1 & State !="ALASKA") %>%
  # exclude Alaska because no data (not a state) before 1959
  group_by(State) %>%
  mutate(ryld = resid(lm(Value~Year)))

# put data into a 'wide format' matrix, with a column for each state and row for each year
yld.mat = 
  detr.yld %>%
  select(State,Year,ryld) %>%
  spread(State,ryld)

# name state columns by USPS postal abbreviations
colnames(yld.mat) = c("year","AL","AZ","AR","CA","CO","CT","DE","FL","GA","ID","IL",
                     "IN","IA","KS","KY","LA","ME","MD","MA","MI","MN","MS","MO","MT",
                     "NE","NV","NH","NJ","NM","NY","NC","ND","OH","OK","OR","PA","RI",
                     "SC","SD","TN","TX","UT","VT","VA","WA","WV","WI","WY")
yld.mat=yld.mat[,2:ncol(yld.mat)]

# create matrix of Pearson correlations among states
yld.test = 
  yld.mat %>%
  corr.test(adjust="none")
  
yld.cor = 
  data.frame(reshape2::melt(yld.test$r), 
             reshape2::melt(yld.test$p)[3])
colnames(yld.cor) = c("state1","state2","cor","pval")
yld.cor$pairid = paste(as.character(yld.cor$state1),"-",as.character(yld.cor$state2),sep="")
yld.cor=yld.cor[yld.cor$cor!=1,]

#### combining yield correlations and interstate distances
yldcordist = merge(yld.cor,dstate,by='pairid',all.x=T)
colnames(yldcordist)[c(2,3)] <- c("state1","state2")
yldcordist = yldcordist[,c("pairid","state1","state2","cor","pval","dist","dcat","reg")]
yldcordist$sig = ifelse(yldcordist$pval<0.05,"Y","N")

### delete duplicate state pairs
nrow(yldcordist) == length(unique(yldcordist$cor))*2 # this tests that all values of "cor" are unique
  # if that returns "TRUE", then proceed (if not you need to find another way of removing duplicates)
yldcordist = distinct(yldcordist, cor, .keep_all = TRUE); head(yldcordist); nrow(yldcordist)

return(yldcordist)
}

yldreg = yldsynch(1966,2015) # this is for general plotting in period of BBS
yldreg.09.30 = yldsynch(1909,1930)
yldreg.31.51 = yldsynch(1931,1951)
yldreg.52.72 = yldsynch(1952,1972)
yldreg.73.94 = yldsynch(1973,1994)
yldreg.95.16 = yldsynch(1995,2016)


### FUCNTION TO PLOT ALL RAW AND DETRENDED YIELD TIME SERIES FOR ALL STATES TO VISUALLY CHECK STATIONARITY (CONSTANT MEAN AND VARIANCE OVER TIME)
yldplot = function(yearstart, yearend){
  detr.yld2 = 
    h %>% 
    filter(Year > yearstart - 1 & Year < yearend+1) %>%
    group_by(State) %>%
    mutate(ryld = resid(lm(Value~Year)))
  ggplot(detr.yld2) + 
    geom_line(aes(x=Year,y=Value)) + 
    geom_line(aes(x=Year,y=ryld)) +
    facet_wrap(~State, scales = "free_y") + 
    theme_bw() +
    labs(title = paste("Years: ",yearstart,"-",yearend,sep=""))
}

# plotting the time series using the yldplot function
X11(26,18);(yldplot.09.16 = yldplot(1909,2016)) # entire 1909-2016 time series for fun - don't expect it to be linear
X11(26,18);(yldplot.66.15 = yldplot(1966,2015)) # same period of the Breeding Bird Survey data
X11(26,18);(yldplot.09.30 = yldplot(1909,1930))
X11(26,18);(yldplot.31.51 = yldplot(1931,1951))
X11(26,18);(yldplot.52.72 = yldplot(1952,1972))
X11(26,18);(yldplot.73.94 = yldplot(1973,1994))
X11(26,18);(yldplot.95.16 = yldplot(1995,2016))
#ggsave("figures/YLD_all_states_detrended.png")

# creating dataframe that has pairwise synchrony for each time period (in separate columns)
# need to change to make it "long" format
yldtime = cbind.data.frame(all.09.30 = yldsynch(1909,1930), all.31.51 = yldsynch(1931,1951), 
                     all.52.72 = yldsynch(1952,1972), all.73.94 = yldsynch(1973,1994),
                      all.95.16 = yldsynch(1995,2016)) %>%
    select(c(pairid=all.09.30.pairid,state1=all.09.30.state1,state2=all.09.30.state2,
            dist=all.09.30.dist, dcat=all.09.30.dcat,reg=all.09.30.reg,
                      all.09.30.cor,all.09.30.pval,all.09.30.sig,
                      all.31.51.cor,all.31.51.pval,all.31.51.sig,
                      all.52.72.cor,all.52.72.pval,all.52.72.sig,
                      all.73.94.cor,all.73.94.pval,all.73.94.sig,
                      all.95.16.cor,all.95.16.pval,all.95.16.sig))

### SCATTERPLOT OF SYNCHRONY VS. DISTANCE (FILTER BY REGION IF DESIRED)
X11(13,9)
ggplot(filter(yldtime,reg %in% c("east","west","both"))) + 
  geom_point(aes(x=dist,y=all.31.51.cor,color=all.31.51.sig),size=4, alpha = 0.5) + 
  geom_smooth(aes(x=dist,y=all.31.51.cor),method="loess",color="black")
#ggsave("figures/YLD_cor_all_1931-1951.loess.png")
#ggsave("figures/YLD_cor_east_1966-2015.png")
#ggsave("figures/YLD_cor_west_1966-2015.png")


### a better way to do it: rbinding together to make it easier to plot in ggplot
yld.66.90 = cbind.data.frame(yldsynch(1966,1990), per = "1966-1990")
yld.91.15 = cbind.data.frame(yldsynch(1991,2015), per = "1991-2015")
yld2time = rbind(yld.66.90,yld.91.15)

ggplot(yld2time) +
#  geom_point(aes(x = dist, y = cor, color=per, shape = sig), size=2, alpha = 0.5) +
  scale_shape_manual(values = c(16,1)) +
  scale_color_manual(values = c("#de2d26","#2b8cbe")) +
  geom_smooth(aes(x = dist, y = cor, color=per), method = "loess", size = 2) +
  labs(title = "Spatial Synchrony in Hay Yield") +
  xlim(c(0,3000))
ggsave(figures/"yield_synchrony.1966.vs.2015.nopts.png")


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

