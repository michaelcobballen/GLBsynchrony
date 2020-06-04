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

# Hay yield data for US states from USDA NASS ('all hay') in tons / acre
yld = read.csv("data/AllHay_Yield_survey_allstates.csv")


### Make a function to conduct bird r vs. yield regressions by state

birdyldlm = function(spcode, birdyearstart,birdyearend) { 
  # birdyearstart should be the year before the first pop change (r) value, e.g., 1966 when r starts at 1967
  # use birdcode = s05460 for GRSP
  # (for now: look up other species codes on BBS website and change the sp=="" to desired species)
# filtering and joining BBS data to state code data; creating "population growth rate" (r) variable
# create "statelagyear" variable to join with yield data in the proper time lag structure (see below)

bird =
  bbs %>%
  filter(sp==spcode & state!= "X" & state != "CAN" & year > birdyearstart-1 & year < birdyearend+1) %>%
  arrange(state,year) %>%
  mutate(r = ifelse(year==birdyearstart,NA,100*c(NA,diff(ind))/lag(ind))) %>%
  left_join(states.lookup, by = "state") %>%
  select(sp, state, year, ind,cred, r, State.ANSI = STATEFP, name = NAME) %>%
  mutate(statelagyear = paste(state,year,sep=""))
bird=droplevels(bird)

# filtering data and joining with state code data
# creating detrended yield variable (ryld) by regressing against year (within states)
# subsetted to relate "lag 1" yield with "lag 0" population growth rate  

detr.yld = 
  yld %>% 
  filter(Year > birdyearstart-1 & Year < birdyearend & Period=="YEAR") %>%
  group_by(State) %>%
  mutate(ryld = resid(lm(Value~Year))) %>%
  arrange(State,Year) %>%
  select(State,State.ANSI,Year,Value,ryld)

detr.yld2 = 
  detr.yld %>%
  left_join(states.lookup, by = "State.ANSI") %>%
  select(State, state, year = Year, yield = Value, ryld, State.ANSI = STATEFP, name = NAME) %>%
  mutate(lagyear = year + 1, statelagyear = paste(state,lagyear,sep=""), State.ANSI = as.numeric(State.ANSI))

# joining the bird data with the yield data (lagged 1 year)

r.ryld1 = 
  full_join(bird,detr.yld2,by = "statelagyear") %>%
  select(stateyearbird = statelagyear, sp,State, state = state.x, ansi = State.ANSI.x, name = name.x, 
         yearbird = year.x, yearyld = year.y,ind,r,yield.lag1 = yield,ryld.lag1 = ryld)

# check average BBS index values to exclude states with very low populations (working criteria: mean index < 0.1)

mean.index=data.frame(mean.index = tapply(bird$ind,bird$state,mean, na.rm=T))
mean.index$state = rownames(mean.index)
mean.index$exclude = ifelse(mean.index$mean.index < 0.1, "Y", "N")

# calculate lm for r ~ ryld in each state
# from: https://stackoverflow.com/questions/22713325/fitting-several-regression-models-with-dplyr

r.ryld = r.ryld1 %>% 
  filter(is.na(r)!=T & is.na(ryld.lag1)!=T) %>% 
  select(state,yearbird,r,ryld.lag1) 

r.ryld.lm = r.ryld %>%  
  group_by(state) %>% 
  do(fit.r.ryld = lm(r ~ ryld.lag1, na.action=na.exclude, data = .)) 

r.ryld.lm.tidy = tidy(r.ryld.lm, fit.r.ryld)

r.ryld.lm.beta = r.ryld.lm.tidy %>%
    left_join(mean.index,by="state") %>%
    filter(term == "ryld.lag1") # add this to exclude these before spatializing: & exclude =="N"
return(r.ryld.lm.beta)
}

r.ryld.lm.beta.66.90 = birdyldlm("s05460",1966,1990)
r.ryld.lm.beta.90.15 = birdyldlm("s05460",1990,2015)

r.ryld.lm.beta.states1.66.90 = merge(states.latlong,r.ryld.lm.beta.66.90, by.y="state",by.x="STUSPS")
r.ryld.lm.beta.states.66.90 = subset(r.ryld.lm.beta.states1.66.90,as.numeric(as.character(GEOID))<57 & STUSPS != "HI" 
              & STUSPS != "AK")
# write.csv(r.ryld.lm.beta.states.66.90@data,"r.ryld.lm.beta.66.90-function.csv")


r.ryld.lm.beta.states1.90.15 = merge(states.latlong,r.ryld.lm.beta.90.15, by.y="state",by.x="STUSPS")
r.ryld.lm.beta.states.90.15 = subset(r.ryld.lm.beta.states1.90.15,as.numeric(as.character(GEOID))<57 & STUSPS != "HI" 
                                     & STUSPS != "AK")
# write.csv(r.ryld.lm.beta.states.90.15@data,"r.ryld.lm.beta.90.15-function.csv")


####
#### Map of GRSP v. Yield Relationship, 1990-2015
####

r.ryld.lm.beta.states.fort.90.15 =
  r.ryld.lm.beta.states.90.15 %>%
  fortify(region = "STUSPS")

r.ryld.lm.beta.ggmap.90.15 = merge(r.ryld.lm.beta.states.fort.90.15,
                                   r.ryld.lm.beta.states.90.15@data, by.x ="id",by.y="STUSPS") %>%
  mutate(cat = cut(estimate,breaks = c(-5000,-40,-20,-10,0,10,20,40,5000),
                   labels = c("< -40","-40 to -20","-20 to -10","-10 to 0","0 to 10","10 to 20","20 to 40","> 40"))) 

r.ryld.lm.sig.90.15 = cbind.data.frame(r.ryld.lm.beta.states.90.15@data,
                                       lon = coordinates(r.ryld.lm.beta.states.90.15)[,1],
                                       lat = coordinates(r.ryld.lm.beta.states.90.15)[,2])
r.ryld.lm.sig.90.15 = filter(r.ryld.lm.sig.90.15, is.na(p.value) != T)


x11(17,10)
ggplot(r.ryld.lm.beta.ggmap.90.15) +
  borders("usa",fill=gray(1)) +
  geom_polygon(aes(long, lat, group=group, fill = cat)) +
  coord_equal()  +
  scale_fill_manual(values=c("#8c2d04","#cc4c02","#ec7014","#fe9929",gray(.8),gray(.7),gray(.6)),
                    labels = c("< -40","-40 to -20","-20 to -10","-10 to 0","0 to 10","10-20","20-40")) +
  labs(fill="slope", title = "Grasshopper Sparrow Pop. Growth Rate vs. Previous Year's Hay Yield", subtitle = "1990-2015") +
  geom_point(data=r.ryld.lm.sig.90.15, aes(x = lon,y= lat,color = ifelse(p.value<0.1,"Y","N")),
             color = ifelse(r.ryld.lm.sig.90.15$p.value<0.1,"black","white"), size=2.5)

ggsave("figures/r.ryld.lm.beta.states.1990-2015.png")


####
#### Map of GRSP v. Yield Relationship, 1966-1990
####

r.ryld.lm.beta.states.fort.66.90 =
  r.ryld.lm.beta.states.66.90 %>%
  fortify(region = "STUSPS")

r.ryld.lm.beta.ggmap.66.90 = merge(r.ryld.lm.beta.states.fort.66.90,
                                   r.ryld.lm.beta.states.66.90@data, by.x ="id",by.y="STUSPS") %>%
  mutate(cat = cut(estimate,breaks = c(-5000,-40,-20,-10,0,10,20,40,5000),
                   labels = c("< -40","-40 to -20","-20 to -10","-10 to 0","0 to 10","10 to 20","20 to 40","> 40"))) 

r.ryld.lm.sig.66.90 = cbind.data.frame(r.ryld.lm.beta.states.66.90@data,
                                       lon = coordinates(r.ryld.lm.beta.states.66.90)[,1],
                                       lat = coordinates(r.ryld.lm.beta.states.66.90)[,2])
r.ryld.lm.sig.66.90 = filter(r.ryld.lm.sig.66.90, is.na(p.value) != T)


x11(17,10)
ggplot(r.ryld.lm.beta.ggmap.66.90) +
  borders("usa",fill=gray(1)) +
  geom_polygon(aes(long, lat, group=group, fill = cat)) +
  coord_equal()  +
  scale_fill_manual(values=c("#8c2d04","#cc4c02","#ec7014","#fe9929",gray(.8),gray(.7),gray(.6),gray(.5)),
                    labels = c("< -40","-40 to -20","-20 to -10","-10 to 0","0 to 10","10-20","20-40","> 40")) +
  labs(fill="slope", title = "Grasshopper Sparrow Pop. Growth Rate vs. Previous Year's Hay Yield", subtitle = "1966-1990") +
  geom_point(data=r.ryld.lm.sig.66.90, aes(x = lon,y= lat,color = ifelse(p.value<0.05,"Y","N")),
             color = ifelse(r.ryld.lm.sig.66.90$p.value<0.05,"black","white"), size=2.5)

ggsave("figures/r.ryld.lm.beta.states.1966-1990.png")






#### a test to make sure lm's worked (top code borrowed from inside the function above - see comments there)

birdyearstart = 1966
birdyearend = 1990
spcode = "s05460"

bird =
  bbs %>%
  filter(sp==spcode & state!= "X" & state != "CAN" & year > birdyearstart-1 & year < birdyearend+1) %>%
  arrange(state,year) %>%
  mutate(r = ifelse(year==birdyearstart,NA,100*c(NA,diff(ind))/lag(ind))) %>%
  left_join(states.lookup, by = "state") %>%
  select(sp, state, year, ind,cred, r, State.ANSI = STATEFP, name = NAME) %>%
  mutate(statelagyear = paste(state,year,sep=""))
bird=droplevels(bird)

detr.yld = 
  yld %>% 
  filter(Year > birdyearstart-1 & Year < birdyearend & Period=="YEAR") %>%
  group_by(State) %>%
  mutate(ryld = resid(lm(Value~Year))) %>%
  arrange(State,Year) %>%
  select(State,State.ANSI,Year,Value,ryld)

detr.yld2 = 
  detr.yld %>%
  left_join(states.lookup, by = "State.ANSI") %>%
  select(State, state, year = Year, yield = Value, ryld, State.ANSI = STATEFP, name = NAME) %>%
  mutate(lagyear = year + 1, statelagyear = paste(state,lagyear,sep=""), State.ANSI = as.numeric(State.ANSI))

r.ryld1 = 
  full_join(bird,detr.yld2,by = "statelagyear") %>%
  select(stateyearbird = statelagyear, sp,State, state = state.x, ansi = State.ANSI.x, name = name.x, 
         yearbird = year.x, yearyld = year.y,ind,r,yield.lag1 = yield,ryld.lag1 = ryld)

r.ryld = r.ryld1 %>% 
  filter(is.na(r)!=T & is.na(ryld.lag1)!=T) %>% 
  select(state,yearbird,r,ryld.lag1)
head(r.ryld)

### this code checks that lm's worked - they did
test = r.ryld %>%
  filter(state=="NV") 
test1 = test %>%
  lm(r~ryld.lag1, na.action=na.exclude, data = .)
tidy(test1)  
