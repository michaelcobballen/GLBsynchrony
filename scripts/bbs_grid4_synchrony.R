###############################################################################################
###
### Summarize route-level BBS data by 2x2 degree grid-cell for spatial synchrony analysis
### Eventually will use methods to correct for observer effects, e.g., the conditional autoregressive (CAR) spatial model described with code in Smith et al. PlosOne (2015) and used by Michel et al. Ecography (2016)
###
###############################################################################################

# see Links at bottom for data sources, papers, etc.
# to do:
# 1) include raw "inclusion" variables (e.g., num zeros) into the final data 'outputs' 
      # that way grid-cells can be excluded as needed on the fly in the final synchrony step
      # without re-running entire code
# 2) incorporate a re-sampling / montecarlo test for the binned bar graphs and/or maps

#######################################
### load libraries
#######################################
library(dplyr)
library(stringr)
library(sp)
library(rgdal)
library(maptools)
library(rgeos)
library(geosphere)
library(ggplot2)
library(psych)
library(tidyr)
library(RColorBrewer)
library(ggthemes)
library(ggridges)
library(viridis)
library(forcats)
library(lme4)
library(Hmisc)
load("~/Research/Chapter 3 - Synchrony/GLBsynchrony/bbsync.RData")

#######################################
### Read in BBS data (route, weather, and separate bird data files by state)
#######################################
path <- "data//states//"
files <- list.files(path=path, pattern="*.csv")
for(file in files)
{
  perpos <- which(strsplit(file, "")[[1]]==".")
  assign(
    gsub(" ","",substr(file, 1, perpos-1)), 
    read.csv(paste(path,file,sep=""),header=T,skip=0))
}

# combine data files for birds (by state), and read in route and weather data
# BBS fields: ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/NABBS_DataFieldDefinitions_1996-2017.xml
b = rbind.data.frame(Alabama,Alaska,Alberta,Arizona,Arkansa,BritCol,Califor,Colorad,Connect,Delawar,
                   Florida,Georgia,Idaho,Illinoi,Indiana,Iowa,Kansas,Kentuck, Louisia,Maine,
                   Manitob,Marylan,Massach,Michiga,Minneso,Mississ,Missour,Montana,NBrunsw,
                   NCaroli,NDakota,Nebrask,Nevada,Newfoun,NHampsh,NJersey,NMexico,NovaSco,Nunavut,
                   NWTerri,NYork,Ohio,Oklahom,Ontario,Oregon,PEI,Pennsyl,Quebec,RhodeIs,
                   Saskatc,SCaroli,SDakota,Tenness,Texas,Utah,Vermont,Virgini,W_Virgi,Washing,
                   Wiscons,Wyoming,Yukon) %>% 
        mutate(rid = paste(StateNum,".",Route,sep=""))
r = read.csv("data/routes.csv") %>% 
        mutate(statenum=as.character(statenum),rid = paste(statenum,".",Route,sep=""))
w = read.csv("data/weather.csv") %>% 
        mutate(rid = paste(StateNum,".",Route,sep=""))

#######################################
### BBSync function
#######################################

bbsync = function(spcode, yearstart, yearend, method, output){

# if running these as part of the function, then put "#" before them; if not, then remove the "#"
#spcode = "5460" # AOU species number; #  spcode = s05460 for GRSP, s05010 for EAME, s04940 for BOBO, s05420 for SAVS
#yearstart = 1992
#yearend = 2017
#method = "AR" # "AR" = residuals from 1st order autoregressive model, "change" = first difference of logged counts
minyearsdata = 10 # each grid cell contains at least one route with at least this many years of data (10 is based on criteria of Michel et al. Ecography 2016)
minsproutes = 2 # each grid cell contains at least this many routes at which the species was detected (2 is based on criteria of Michel et al. Ecography 2016)
maxzeros = 5 # maximum number of years allowable where the mean count for the species is zero (I chose five)
minyearscor = 20 # minimum number of years two grid cells must have in common to calculate a correlation coefficient between them (Hanski used 7; I chose 20, more conservative)
#output = "grid.map" #or... CAR.data, sync.dist, grid.year, route.year

# functions for identifying even and odd numbers
is.even <- function(x) x %% 2 == 0
is.odd <- function(x) x %% 2 != 0

# combining the bird, route, and weather data
# assigning routes to named grid cells 
# creating a variable for abundance of target species
brw = b %>%
  left_join(r,by = 'rid') %>%
  left_join(w,by = 'RouteDataID') %>%
  filter(RunType==1,RPID.x ==101) %>%   
  dplyr::select(RouteDataID, CountryNum = CountryNum.x, StateNum = StateNum.x, Route = Route.x,
         Year = Year.x, AOU, StopTotal, rid = rid.x, RouteName, lat = Latitude, lon = Longitude, 
         Stratum, BCR, ObsN) %>%
  mutate(sp.number = as.numeric(ifelse(str_detect(AOU,spcode),StopTotal,0))) %>%
  mutate(grid4_lat = ifelse(is.odd(trunc(lat)),trunc(lat)+1,trunc(lat))) %>%
  mutate(grid4_lon = ifelse(is.odd(trunc(lon)),trunc(lon)-1,trunc(lon))) %>%
  mutate(grid4 = paste(grid4_lat,grid4_lon,sep=""))

# summarizing data for species of interest by route-year
brw.route.year = brw %>%
  filter(Year>yearstart-1, Year<yearend+1) %>%
  group_by(RouteDataID) %>%
  summarise(sp.number = sum(sp.number)) %>%
  left_join(distinct(dplyr::select(brw,-AOU,-sp.number,-StopTotal)),by="RouteDataID") %>%
  mutate(ObsN.rid = paste(ObsN,".",rid,sep="")) %>%
  group_by(ObsN.rid) %>%
  arrange(ObsN.rid,Year) %>%
  mutate(ObsN.rid.year = row_number(), yr = Year-1965) %>%  
  mutate(ObsN.rid.first = as.numeric(ifelse(ObsN.rid.year==1,1,0))) %>%
  ungroup() %>%
  filter(is.na(lat)==F) %>% # this excludes ~ 20 routes without attribute data
  group_by(grid4) %>%
  mutate(ObsN.rid.count = as.numeric(as.factor(ObsN.rid))) %>%
  ungroup()

# summarizing route-year data by route 
# including # years species present & # years data on each route
brw.route = brw.route.year %>%
  group_by(rid) %>%
  summarise(sp.yrs = sum(sp.number>0), yrs.data = length(sp.number), sp.mean = mean(sp.number),
            nobservers = max(ObsN.rid.count)) %>%
  left_join(distinct(dplyr::select(brw.route.year,-sp.number,-Year,-RouteDataID,-ObsN,-ObsN.rid,-ObsN.rid.year,-yr,-ObsN.rid.first,-ObsN.rid.count)),by="rid") %>%
  mutate(sp.present = as.numeric(ifelse(sp.yrs>0,1,0))) %>%
  arrange(lat) %>%
  ungroup() 

# summarize data by 2x2 grid cell
# including whether cells meet the "minyearsdata" and "minsproutes" criteria of Michel et al. (Ecography)
grid4.summary = brw.route %>%
  group_by(grid4) %>%
  summarise(grid4_lat=mean(grid4_lat), grid4_lon=mean(grid4_lon), num.routes = length(sp.present), 
            sp.mean = mean(sp.mean), num_sp_rtes = sum(sp.present), 
            min_yr_rtes = sum(yrs.data>=minyearsdata), nobservers = max(nobservers)) %>%
  mutate(sp.min.rtes.yes = as.numeric(ifelse(num_sp_rtes >= minsproutes,1,0))) %>%
  mutate(min.years.yes = as.numeric(ifelse(min_yr_rtes > 0,1,0))) %>%
  mutate(include_strat = as.numeric(ifelse(sp.min.rtes.yes + min.years.yes == 2,1,0))) %>%
  mutate(nonzeroweight = num_sp_rtes/num.routes) %>%
  dplyr::select(-sp.min.rtes.yes,-min.years.yes)

# add grid-level information to the route-year data
# remove routes in cells excluded based on criteria in Michel et al. (Ecography)
route.year.grid4 = brw.route.year %>%
  left_join(dplyr::select(grid4.summary, grid4, include_strat), by = "grid4") %>%
  filter(include_strat==1) %>%
  dplyr::select(count = sp.number, grid4, obser = ObsN.rid.count, firstyr = ObsN.rid.first, 
         year = yr, -include_strat, grid4_lat, grid4_lon) %>%
  mutate(grid4 = gsub("-","_",grid4))

# summarize route-year data by grid-year
grid4.year = route.year.grid4 %>%
  group_by(grid4, year, grid4_lat, grid4_lon) %>%
  summarise(mean.count = mean(count)) %>%
  ungroup() %>%
  mutate(year = year+1965, grid4.year = paste(grid4,year,sep="."))

### create criteria for including only grid cells with < XX years of zeros
grid4.include2 = grid4.year %>%
  group_by(grid4) %>%
  summarise(num.zeros = sum(mean.count==0)) %>%
  mutate(include2 = ifelse(num.zeros >maxzeros-1,0,1)) %>%
  ungroup()

# subsetting route.year.grid4 based on the second inclusion criteria (number of years with zeros)
    # cells already removed based on the Michel et al. criteria
    # this is the data file to be used in the OpenBUGS conditional autoregressive (CAR) model via Smith et al. PlosONE 2015 code
route.year.grid4b = route.year.grid4 %>%
  left_join(grid4.include2, by = "grid4") %>%
  filter(include2 == 1)

# subsetting grid4.summary based on both criteria (Michel et al. and number of years with zeros) 
grid4.summary2 = grid4.summary %>%
  mutate(grid4 = gsub("-","_",grid4)) %>%
  left_join(grid4.include2, by = "grid4") %>%
  filter(include_strat==1 & include2 == 1)

### subsetting the grid-year data based on 2nd criteria (number of years with zeros)
  # this had already been subsetted to exclude based on Michel et al. criteria 
grid4.year = grid4.year %>%
  left_join(grid4.include2,by = "grid4") %>%
  filter(include2 == 1)

### calculating inter-grid distances for 2 X 2 degree lat/long grid
grid4.spatial = SpatialPointsDataFrame(distinct(grid4.year[,c("grid4_lon","grid4_lat")]),
                                             distinct(grid4.year[,c("grid4_lon","grid4_lat","grid4")]))
proj4string(grid4.spatial) = CRS("+init=epsg:4326")
#plot(grid4.spatial)

df.grid4 = data.frame(name = grid4.spatial$grid4, 
                                  lat = coordinates(grid4.spatial)[,2],
                                  lon = coordinates(grid4.spatial)[,1])

#calculate Haversine distances between grid cells
dmat.grid4 = round(distm(df.grid4[,c("lon","lat")])/1000)
colnames(dmat.grid4) <- df.grid4$name; head(dmat.grid4)
dmat.grid4B = data.frame(grid4B=df.grid4$name, dmat.grid4)

dgrid4 = data.frame(reshape2::melt(dmat.grid4B)) %>%
  filter(value!=0) %>%
  dplyr::select(grid4B,grid4A = variable, dist = value) %>%
  mutate(grid4A = gsub("X","",grid4A)) %>%
  mutate(pairid = paste(grid4B,".",grid4A,sep=""))

### synchrony data formatting
g4 = dgrid4

# make a data frame of all possible grid/year combinations
grid4.list = unique(grid4.year$grid4)
year.list = (yearstart:yearend)
grid4.year.list = expand.grid(grid4 = grid4.list,year = year.list) %>%
  arrange(grid4,year) %>%
  mutate(grid4.year = paste(grid4,year,sep="."))

# attach the bird data to the complete list and format for correlation analysis
grid4.year.full = grid4.year.list %>%
  left_join(grid4.year, by = "grid4.year") %>%
  dplyr::select(grid4=grid4.x, year=year.x, grid4.year, mean.count) %>%
  mutate(log.count = log(mean.count+1), logcountlag1 = lag(log.count)) %>% # log+1 transform count ala Hanski and create lag for AR and differencing analyses
  mutate(logcountlag1 = as.numeric(ifelse(year==yearstart,"NA",logcountlag1))) %>% # assign "NA" to first year which should have no lagged value
  filter(is.na(log.count) == F & is.na(logcountlag1) == F) %>% # temporarily remove NAs to get residuals and first diffs of logged counts
  mutate(change = log.count - logcountlag1) %>% # first diff of logged counts, ~ % change
  group_by(grid4) %>%
  mutate(AR1.R = resid(lm(log.count~logcountlag1))) %>% # calculate residuals from first order autoregression (Hanski-style)
  ungroup() %>%
  right_join(grid4.year.list, by = "grid4.year") %>% # merge back with the full grid/year list to allow matching in correlation analyses
  dplyr::select(grid4 = grid4.y, year = year.y, grid4.year, mean.count, log.count,logcountlag1, change, AR = AR1.R)

### put bird population data into 'wide' format for correlation analysis (rows = years, columns = grid cells)
grid4.bird.mat = 
  grid4.year.full %>%
  dplyr::select(grid4,year,method) %>%
  spread(grid4,method) %>%
  filter(year > yearstart & year <= yearend) %>%
  # removes first year which have NAs for all calculations involving lagged counts
  dplyr::select(-year)

# create matrix of Pearson correlations among grid cells
grid4.bird.test = 
  grid4.bird.mat %>%
  corr.test(adjust="none", ci = F)

# ntest checks if the n column of correlation test has length 1
# allows us to deal with it differently in the next section
    # n has length 1 when all pairwise correlations have a sample size = yearend-yearstart (i.e., a quirk of the cortest function)
    # ...e.g., in Eastern Meadowlark, 1992-2017
ntest = nrow(reshape2::melt(grid4.bird.test$n)) 

# reshaping the giant correlation matrix into 'long' format for easier graphing and analysis
grid4.bird.cor = 
  data.frame(reshape2::melt(grid4.bird.test$r), 
             reshape2::melt(grid4.bird.test$p)[3],
             n = ifelse(ntest > 1, reshape2::melt(grid4.bird.test$n)[3],
                        as.list(yearend-yearstart))[[1]]) %>%
  dplyr::select(grid4A=Var1, grid4B=Var2, cor=value, pval=value.1, n) %>%
  mutate(pairid = paste(grid4A,".",grid4B,sep="")) %>%
  filter(n>minyearscor-1) %>% # excludes grid pairs with less than XX years in common (Hanski's criteria was 7)
  filter(grid4A!=grid4B) %>% # excludes self-correlations
  filter(is.na(pval)==F) %>% # excludes NA correlations resulting from too many zeros if these weren't removed earlier - they were! So nothing changes.
  distinct(pval, .keep_all = TRUE) # should cut number in half, and yes it does! (eventually add check and warning message if it doesn't)

# adding in the distance between grid cells again
grid4.bird.cordist = grid4.bird.cor %>%
  left_join(g4,by='pairid') %>%
  dplyr::select(pairid, grid4A=grid4A.x, grid4B=grid4B.x, cor,pval,n,dist) %>%
  mutate(sig = ifelse(pval < 0.05,"Y","N"), sig1 = as.numeric(ifelse(pval < 0.05,"1","0"))) %>%
  distinct(cor, .keep_all = TRUE) 

### making spatial grid cell polygons to extract spatial data for use in the OpenBUGS CAR model via Smith et al. PlosONE 2015 code
if(output == "grid.map" | output == "CAR.data"){
  
# make a grid of spatial polygons to match up with data
library(raster)
rr <- raster::raster(ext = extent(-177, -52, 24, 69), res=c(2,2))
values(rr) <- 1:ncell(rr)
p <- rasterToPolygons(rr)
attributes = tibble(grid4_lon = coordinates(p)[,1], grid4_lat = coordinates(p)[,2]) %>%
  mutate(join = 1, grid4 = paste(trunc(grid4_lat),trunc(grid4_lon),sep="")) %>%
  mutate(grid4 = gsub("-","_",grid4)) %>%
  left_join(grid4.summary2, by="grid4") 

p$grid4 = attributes$grid4
p$include_strat = attributes$include_strat
p$include2 = attributes$include2
p$nonzeroweight = attributes$nonzeroweight
p$nobservers = attributes$nobservers
p$strata = log(attributes$sp.mean)

grid4.bird = p[p$include_strat==1 & is.na(p$include_strat)==F,]
grid4.bird@data = grid4.bird@data %>%
  dplyr::select(layer, grid4, include_strat, include2, nonzeroweight, nobservers, strata) %>%
  mutate(strat = 1:nrow(grid4.bird@data))
head(grid4.bird@data); nrow(grid4.bird@data)

# project into Lambert conformal conic to allow plotting in OpenBUGS (with meters scale) if needed
proj = "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=10000000 +y_0=10000000 +datum=NAD83 +units=m +no_defs"
grid4.bird.lamb = spTransform(grid4.bird, proj)
plot(grid4.bird.lamb, col="gray")
head(grid4.bird.lamb); nrow(grid4.bird.lamb)
  
# changing names of the ID variable of the polygons to grid labels so names are visible in OpenBUGS map
# https://eco-data-science.github.io/spatial_analysis2_R/
sapply(slot(grid4.bird.lamb, "polygons"), function(x) slot(x, "ID"))
for (i in 1:length(slot(grid4.bird.lamb, "polygons"))){
  slot(slot(grid4.bird.lamb, "polygons")[[i]], "ID") = grid4.bird.lamb$grid4[i]
}
sapply(slot(grid4.bird.lamb, "polygons"), function(x) slot(x, "ID"))
}

# conditional statements telling the function what to provide as output (with code for making some)
if(output == "grid.map"){
return(grid4.bird.lamb)
}

if(output == "CAR.model"){
# converting spatial data to geobugs format for CAR model
sp2WB(grid4.bird.lamb, "grid4.bird.lamb.txt")

# creating BUGS data for CAR model 
count = route.year.grid4b$count
strat = route.year.grid4b$strat
obser = route.year.grid4b$obser
firstyr = route.year.grid4b$firstyr
year = route.year.grid4b$year

ncounts = c(nrow(route.year.grid4b))
nstrata = nrow(grid4.bird.lamb@data)
ymin = 1
ymax = 2017-1965

nonzeroweight = grid4.bird.lamb@data$nonzeroweight
strata = grid4.bird.lamb@data$strata
nobservers = grid4.bird.lamb@data$nobservers

require(R2WinBUGS)
bugs.data(list("count"=count, "strat"=strat, "obser"=obser, "firstyr"=firstyr,
               "year"=year,"ncounts"=ncounts,"nstrata"=nstrata,"ymin"=ymin,
               "ymax"=ymax,"nonzeroweight"=nonzeroweight,"nobservers"=nobservers, 
               "strata"=strata), 
          getwd(),12,data.file="bugs.data.txt")

# note: additional "adj", "num", & "sumNumNeigh" are calculated in OpenBUGs from map file
}

if(output == "sync.dist"){ 
  return(grid4.bird.cordist)
}

if(output == "grid.year"){ 
  return(grid4.year.full %>% filter(year > yearstart))
}

if(output == "route.year"){ 
  return(brw.route.year)
}

}


#######################################
###### make maps and shapefiles of bird grids
#######################################
#grsp.map = bbsync("5460",1992,2017,"AR","grid.map"); x11(13,9); plot(grsp.map,col="gray")
# writeOGR(grsp.map, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "grsp_grid4_poly_lamb.92.17", driver="ESRI Shapefile")
#bobo.map = bbsync("4940",1992,2017,"AR","grid.map"); x11(13,9); plot(bobo.map,col="gray")
# writeOGR(bobo.map, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "bobo_grid4_poly_lamb.92.17", driver="ESRI Shapefile")
#eame.map = bbsync("5010",1992,2017,"AR","grid.map"); x11(13,9); plot(eame.map,col="gray")
# writeOGR(eame.map, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "eame_grid4_poly_lamb.92.17", driver="ESRI Shapefile")
#savs.map = bbsync("5420",1992,2017,"AR","grid.map"); x11(13,9); plot(savs.map,col="gray")
# writeOGR(savs.map, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "savs_grid4_poly_lamb.92.17", driver="ESRI Shapefile")

#######################################
###### summarizing function results in dataframe to map / graph more easily
#######################################
# AOU species number; #  spcode = s05460 for GRSP, s05010 for EAME, s04940 for BOBO, s05420 for SAVS
                      # s07550 for WOTH, s06740 for OVEN, 
                      # https://www.pwrc.usgs.gov/bbl/manual/speclist.cfm
allsync1 = rbind(
grsp.66.91 = cbind.data.frame(bbsync("5460",1966,1991,"AR","sync.dist"),sp="GRSP",per="1966-1991"),
grsp.92.17 = cbind.data.frame(bbsync("5460",1992,2017,"AR","sync.dist"),sp="GRSP",per="1992-2017"),
bobo.66.91 = cbind.data.frame(bbsync("4940",1966,1991,"AR","sync.dist"),sp="BOBO",per="1966-1991"),
bobo.92.17 = cbind.data.frame(bbsync("4940",1992,2017,"AR","sync.dist"),sp="BOBO",per="1992-2017"),
eame.66.91 = cbind.data.frame(bbsync("5010",1966,1991,"AR","sync.dist"),sp="EAME",per="1966-1991"),
eame.92.17 = cbind.data.frame(bbsync("5010",1992,2017,"AR","sync.dist"),sp="EAME",per="1992-2017")
)

# in batches because memory couldn't handle it (ctrl + shift + F10 resets R and allows next to run)
allsync2 = rbind(
savs.66.91 = cbind.data.frame(bbsync("5420",1966,1991,"AR","sync.dist"),sp="SAVS",per="1966-1991"),
savs.92.17 = cbind.data.frame(bbsync("5420",1992,2017,"AR","sync.dist"),sp="SAVS",per="1992-2017"),
noha.66.91 = cbind.data.frame(bbsync("3310",1966,1991,"AR","sync.dist"),sp="NOHA",per="1966-1991"),
noha.92.17 = cbind.data.frame(bbsync("3310",1992,2017,"AR","sync.dist"),sp="NOHA",per="1992-2017"),
upsa.66.91 = cbind.data.frame(bbsync("2610",1966,1991,"AR","sync.dist"),sp="UPSA",per="1966-1991"),
upsa.92.17 = cbind.data.frame(bbsync("2610",1992,2017,"AR","sync.dist"),sp="UPSA",per="1992-2017")
)

# in batches because memory couldn't handle it
allsync3 = rbind(
lbcu.66.91 = cbind.data.frame(bbsync("2640",1966,1991,"AR","sync.dist"),sp="LBCU",per="1966-1991"),
lbcu.92.17 = cbind.data.frame(bbsync("2640",1992,2017,"AR","sync.dist"),sp="LBCU",per="1992-2017"),
hola.66.91 = cbind.data.frame(bbsync("4740",1966,1991,"AR","sync.dist"),sp="HOLA",per="1966-1991"),
hola.92.17 = cbind.data.frame(bbsync("4740",1992,2017,"AR","sync.dist"),sp="HOLA",per="1992-2017"),
hesp.66.91 = cbind.data.frame(bbsync("5470",1966,1991,"AR","sync.dist"),sp="HESP",per="1966-1991"),
hesp.92.17 = cbind.data.frame(bbsync("5470",1992,2017,"AR","sync.dist"),sp="HESP",per="1992-2017")
)

# in batches because memory couldn't handle it
allsync4 = rbind(
sewr.66.91 = cbind.data.frame(bbsync("7240",1966,1991,"AR","sync.dist"),sp="SEWR",per="1966-1991"),
sewr.92.17 = cbind.data.frame(bbsync("7240",1992,2017,"AR","sync.dist"),sp="SEWR",per="1992-2017"),
sppi.66.91 = cbind.data.frame(bbsync("7000",1966,1991,"AR","sync.dist"),sp="SPPI",per="1966-1991"),
sppi.92.17 = cbind.data.frame(bbsync("7000",1992,2017,"AR","sync.dist"),sp="SPPI",per="1992-2017"),
dick.66.91 = cbind.data.frame(bbsync("6040",1966,1991,"AR","sync.dist"),sp="DICK",per="1966-1991"),
dick.92.17 = cbind.data.frame(bbsync("6040",1992,2017,"AR","sync.dist"),sp="DICK",per="1992-2017")
)

# in batches because memory couldn't handle it
allsync5 = rbind(
vesp.66.91 = cbind.data.frame(bbsync("5400",1966,1991,"AR","sync.dist"),sp="VESP",per="1966-1991"),
vesp.92.17 = cbind.data.frame(bbsync("5400",1992,2017,"AR","sync.dist"),sp="VESP",per="1992-2017"),
larb.66.91 = cbind.data.frame(bbsync("6050",1966,1991,"AR","sync.dist"),sp="LARB",per="1966-1991"),
larb.92.17 = cbind.data.frame(bbsync("6050",1992,2017,"AR","sync.dist"),sp="LARB",per="1992-2017"),
bais.66.91 = cbind.data.frame(bbsync("5450",1966,1991,"AR","sync.dist"),sp="BAIS",per="1966-1991"),
bais.92.17 = cbind.data.frame(bbsync("5450",1992,2017,"AR","sync.dist"),sp="BAIS",per="1992-2017")
)

# in batches because memory couldn't handle it
allsync6 = rbind(
casp.66.91 = cbind.data.frame(bbsync("5780",1966,1991,"AR","sync.dist"),sp="CASP",per="1966-1991"),
casp.92.17 = cbind.data.frame(bbsync("5780",1992,2017,"AR","sync.dist"),sp="CASP",per="1992-2017"),
lesp.66.91 = cbind.data.frame(bbsync("5480",1966,1991,"AR","sync.dist"),sp="LESP",per="1966-1991"),
lesp.92.17 = cbind.data.frame(bbsync("5480",1992,2017,"AR","sync.dist"),sp="LESP",per="1992-2017"),
weme.66.91 = cbind.data.frame(bbsync("5011",1966,1991,"AR","sync.dist"),sp="WEME",per="1966-1991"),
weme.92.17 = cbind.data.frame(bbsync("5011",1992,2017,"AR","sync.dist"),sp="WEME",per="1992-2017")
)

# in batches because memory couldn't handle it
allsync7 = rbind(
rnep.66.91 = cbind.data.frame(bbsync("3091",1966,1991,"AR","sync.dist"),sp="RNEP",per="1966-1991"),
rnep.92.17 = cbind.data.frame(bbsync("3091",1992,2017,"AR","sync.dist"),sp="RNEP",per="1992-2017"),
mclo.66.91 = cbind.data.frame(bbsync("5390",1966,1991,"AR","sync.dist"),sp="MCLO",per="1966-1991"),
mclo.92.17 = cbind.data.frame(bbsync("5390",1992,2017,"AR","sync.dist"),sp="MCLO",per="1992-2017"),
cclo.66.91 = cbind.data.frame(bbsync("5380",1966,1991,"AR","sync.dist"),sp="CCLO",per="1966-1991"),
cclo.92.17 = cbind.data.frame(bbsync("5380",1992,2017,"AR","sync.dist"),sp="CCLO",per="1992-2017")
)

# in batches because memory couldn't handle it
# FEHA didn't work - no qualifying cells?
#allsync8 = rbind(
#feha.66.91 = cbind.data.frame(bbsync("3480",1966,1991,"AR","sync.dist"),sp="FEHA",per="1966-1991"),
#feha.92.17 = cbind.data.frame(bbsync("3480",1992,2017,"AR","sync.dist"),sp="FEHA",per="1992-2017")
#)

# GRPC didn't work - no qualifying cells?
#allsync9 = rbind(
#grpc.66.91 = cbind.data.frame(bbsync("3050",1966,1991,"AR","sync.dist"),sp="GRPC",per="1966-1991"),
#grpc.92.17 = cbind.data.frame(bbsync("3050",1992,2017,"AR","sync.dist"),sp="GRPC",per="1992-2017")
#)

allsync10 = rbind(
stgr.66.91 = cbind.data.frame(bbsync("3080",1966,1991,"AR","sync.dist"),sp="STGR",per="1966-1991"),
stgr.92.17 = cbind.data.frame(bbsync("3080",1992,2017,"AR","sync.dist"),sp="STGR",per="1992-2017")
)

# in batches because memory couldn't handle it
allsync11 = rbind(
mopl.66.91 = cbind.data.frame(bbsync("2810",1966,1991,"AR","sync.dist"),sp="MOPL",per="1966-1991"),
mopl.92.17 = cbind.data.frame(bbsync("2810",1992,2017,"AR","sync.dist"),sp="MOPL",per="1992-2017")
)

# BNOW didn't work - no qualifying cells?
#allsync12 = rbind(
#bnow.66.91 = cbind.data.frame(bbsync("3650",1966,1991,"AR","sync.dist"),sp="BNOW",per="1966-1991"),
#bnow.92.17 = cbind.data.frame(bbsync("3650",1992,2017,"AR","sync.dist"),sp="BNOW",per="1992-2017")
#)

# SEOW didn't work - no qualifying cells?
#allsync13 = rbind(
#seow.66.91 = cbind.data.frame(bbsync("3670",1966,1991,"AR","sync.dist"),sp="SEOW",per="1966-1991"),
#seow.92.17 = cbind.data.frame(bbsync("3670",1992,2017,"AR","sync.dist"),sp="SEOW",per="1992-2017")
#)

allsync = rbind(allsync1,allsync2,allsync3,allsync4,allsync5,allsync6,allsync7,allsync10,allsync11) 

allsync = allsync %>%
  group_by(sp, per) %>%
  mutate(pairs = length(sp)) %>%
  ungroup()

table(allsync$sp, allsync$per)
# number of grid-cell-pairs per species per period
#     1966-1991 1992-2017
#GRSP      3566      6441
#BOBO      3383      4750
#EAME      7138      7875
#SAVS      4992     13659
#NOHA       969      3001
#UPSA       703      1953
#LBCU        65       561
#HOLA      9481     16457
#HESP         6         6
#SEWR       136       375
#SPPI        15        55
#DICK      1890      2701
#VESP      5015     10132
#LARB       276       378
#BAIS        28        78
#CASP       125       300
#LESP        10       245
#WEME      6515     11898
#RNEP      4269      6441
#MCLO         1        28
#CCLO        91       136
#STGR         1        36
#MOPL         1         6

#######################################
###### graphing synchrony - continous / smoothed
#######################################

X11(13,9)
allsync %>%
  filter(sp %in% c("WEME","EAME","VESP")) %>%
  ggplot() + 
  geom_smooth(aes(x=dist, y=cor, color=sp, linetype=per), method="auto", se=F) +
  xlim(c(0,2000)) +
  theme_bw() +
  theme(text = element_text(size=15)) +
  labs(x = "Distance (km)", y = "Pearson's correlation", linetype = "Period", color = "Species")
#ggsave("figures/grid4_sync_scatter/EAME.GRSP.BOBO.SAVS_cor_grid4_66.17_no.se_AR_5zeros_n20.jpg")
    # suggestion to use method="auto" as "loess" didn't have enough memory: https://groups.google.com/forum/#!topic/ggplot2/enavD18MmkY

#######################################
###### graphing synchrony - % signficant correlations by distance bins 
#######################################

binlabels = c("<100","101-200","201-300","301-400","401-500","501-600",
              "601-700","701-800","801-900","901-1000","1001-1100","1101-1200",
              "1201-1300","1301-1400","1401-1500","1501-1600","1601-1700",
              "1701-1800","1801-1900","1901-2000",">2000")
binlabels2 = c("<200","201-400","401-600","601-800","801-1000","1001-1200",
               "1201-1400","1401-1600","1601-1800","1801-2000",">2000")

allsync1 = allsync %>%
  mutate(dcat=cut(allsync$dist,c(0,100,200,300,400,500,600,700,800,900,
                                                  1000,1100,1200,1300,1400,1500,1600,1700,
                                                  1800,1900,2000,100000),
                  labels=binlabels)) %>%
  mutate(dcat2 = cut(allsync$dist,c(0,200,400,600,800,1000,
                                                     1200,1400,1600,1800,2000,100000),
                     labels=binlabels2))

corbin1 = filter(allsync1,per=="1992-2017") %>%
  group_by(dcat, sp) %>%
  summarise(pctsig = mean(sig1)*100, n = length(sig1)) %>%
  ungroup()

corbin2 = filter(allsync1,per=="1992-2017") %>%
  group_by(dcat2,sp) %>%
  summarise(pctsig = mean(sig1)*100, n = length(sig1)) %>%
  ungroup()

x11(18,8.5)
ggplot(corbin1) + 
  geom_col(aes(y=pctsig,x=dcat,fill=sp), position = "dodge") + 
  labs(x="Distance (km)",y="% of correlations with P < 0.05", fill="Species") +
  theme_bw() +
  theme(text = element_text(size=15))
#ggsave("figures/grid4_bin_cor/GRSPgrid4_100bin_92.17_AR_5zeros_n20.png")

x11(13,8.5)
ggplot(corbin2) + 
  geom_col(aes(y=pctsig,x=dcat2,fill=sp), position = "dodge") + 
  labs(x="Distance (km)",y="% of correlations with P < 0.05", fill="Species") +
  theme_bw() +
  theme(text = element_text(size=15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#ggsave("figures/grid4_bin_cor/allsp.grid4_200bin_92.17_AR_5zeros_n20.png")

## NEXT: need to create n = 1000 null (random) distributions for each bin like Martin et al. 
## to test spatial extent of synchrony

#######################################
###### mapping mean synchrony
#######################################

# calculating mean synchrony within 400 m
mapsync = 
  rbind(data.frame(allsync, grid = allsync$grid4A), data.frame(allsync, grid = allsync$grid4B)) %>%
  filter(dist<400) %>%
  group_by(grid, per, sp) %>%
  summarise(meancor = mean(cor), n = length(cor)) %>%
  ungroup() %>%
  mutate(grid_lat = as.numeric(substr(grid,1,2)), grid_lon = as.numeric(paste("-",substr(grid,4,10),sep=""))) %>%
  mutate(corcat = cut(meancor, breaks = c(-1,0,0.1,.2,.3,.4,.5,.6,.7))) %>%
  group_by(sp, per) %>%
  mutate(cells = length(sp)) %>%
  ungroup()

# making North America map
can = readOGR("data/province.shp")
canada.fort = ggplot2::fortify(can) %>%
  filter(lat <60.0001)

na <- ggplot() +
  borders("world", xlim = c(-150, -60), ylim = c(35, 40), colour = "gray85", fill = "gray80")  +
  borders("state", colour = "gray85", fill = "gray80") +
  geom_polygon(data=canada.fort, aes(x = long, y = lat, group=group), colour = "gray85", fill = "gray80") +
  theme_map() + 
  coord_map('albers', lat0=30, lat1=60) #+
  #xlim(c(-150,-60))

# mapping mean synchrony 1992-2017 
x11(30,20)
na +
  geom_point(aes(x = grid_lon, y = grid_lat, color = corcat), size = 1, alpha = .7, data = filter(mapsync,per=="1992-2017",grid_lat<=59)) + 
  scale_colour_manual(values = c("#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026", "#800026")) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  facet_wrap(~sp) +
  labs(colour="Mean\nsynchrony", x="", y="", title="1992-2017") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59))
#ggsave("figures/grid4_maps/ALLsp.grid4_92.17_AR_5zeros_n20_60lat.png")

# mapping mean synchrony 1966-1991 
x11(30,20)
na +
  geom_point(aes(x = grid_lon, y = grid_lat, color = corcat), size = 1, alpha = .7, data = filter(mapsync,per=="1966-1991",grid_lat<=59)) + 
  scale_colour_manual(values = c("#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026", "#800026")) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  facet_wrap(~sp) +
  labs(colour="Mean\nsynchrony", x="", y="", title="1966-1991") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59))
#ggsave("figures/grid4_maps/ALLsp.grid4_66.91_AR_5zeros_n20_60lat.png")

#######################################
###### mapping change in synchrony
#######################################

mapsync.early = mapsync %>% filter(per=="1966-1991") %>% mutate(grid.sp = paste(grid,sp,sep="."))
mapsync.late = mapsync %>% filter(per=="1992-2017") %>% mutate(grid.sp = paste(grid,sp,sep="."))
syncdiff = full_join(mapsync.early, mapsync.late, by = "grid.sp") %>%
  filter(is.na(meancor.x)==F, is.na(meancor.y)==F) %>%
  mutate(syncdiff = meancor.y - meancor.x) %>%
  mutate(diffcat = cut(syncdiff, breaks = c(-.6,-.5,-.4,-.3,-.2,-.1,0,.1,.2,.3,.4,.5,.6))) %>%
  mutate(incdec = ifelse(syncdiff>0,"Increasing","Decreasing")) %>%
  dplyr::select(-grid.y,-grid_lat.y, -grid_lat.y)

# mapping change in synchrony - magnitude
x11(30,20)
na +
  geom_point(aes(x = grid_lon.x, y = grid_lat.x, color = diffcat), size = 1, alpha = .7, data = syncdiff) + 
  scale_colour_manual(values = c("#053061","#5e4fa2","#3288bd","#66c2a5","#abdda4","#e6f598",
                                 "#fee08b","#fdae61","#f46d43","#d53e4f","#9e0142", "#67001f")) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  facet_wrap(~sp.x) +
  labs(colour="Change in\nmean\nsynchrony", x="", y="", title="1966-1991 vs. 1992-2017") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59))
#ggsave("figures/grid4_maps/ALLsp.grid4_change_AR_5zeros_n20_60lat.png")

# mapping change in synchrony - increase/decrease
x11(30,20)
na +
  geom_point(aes(x = grid_lon.x, y = grid_lat.x, color = incdec), size = 1, alpha = .7, data = syncdiff) + 
  scale_colour_manual(values = c("blue","red")) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  facet_wrap(~sp.x) +
  labs(colour="Change in\nmean\nsynchrony", x="", y="", title="1966-1991 vs. 1992-2017") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59))
#ggsave("figures/grid4_maps/ALLsp.grid4_incdec_AR_5zeros_n20_60lat.png")

# mapping change in synchrony - magnitude, continous with viridis color scheme
x11(30,20)
na +
  geom_point(aes(x = grid_lon.x, y = grid_lat.x, color = syncdiff), size = 1, alpha = .7, data = syncdiff) + 
  scale_color_viridis(name = "change", option = "viridis") +
  theme_classic() +
  theme(text = element_text(size=15)) +
  facet_wrap(~sp.x) +
  labs(colour="Change in\nmean\nsynchrony", x="", y="", title="1966-1991 vs. 1992-2017") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59))
#ggsave("figures/grid4_maps/ALLsp.grid4_change.viridis_AR_5zeros_n20_60lat.png")

#######################################
###### mapping MEAN change in synchrony across all species
#######################################

meansyncdiff = syncdiff %>%
  group_by(grid.x, grid_lat.x, grid_lon.x) %>%
  summarise(meansyncdiff = mean(syncdiff), n = length(sp.x)) %>%
  ungroup() %>%
  filter(n >= 3) %>%
  mutate(meandiffcat = cut(meansyncdiff, breaks = c(-.3,-.25,-.2,-.15,-.1,-.05,0,.05,.1,.15,.2,.25,.3))) %>%
  mutate(meanincdec = ifelse(meansyncdiff>0,"Increasing","Decreasing")) %>%
  mutate(species = "All species")
  
# mapping change in synchrony - magnitude
x11(13,9)
na +
  geom_point(aes(x = grid_lon.x, y = grid_lat.x, color = meandiffcat), size = 4, alpha = 1, data = meansyncdiff) + 
  scale_colour_manual(values = c("#053061","#5e4fa2","#3288bd","#66c2a5","#abdda4","#e6f598",
                                 "#fee08b","#fdae61","#f46d43","#d53e4f","#9e0142", "#67001f")) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  labs(colour="Change in\nmean spatial\nsynchrony", x="", y="", title="1966-1991 vs. 1992-2017") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59))
#ggsave("figures/grid4_maps/Mean.ALLsp.grid4_change_AR_5zeros_n20_60lat.png")

# mapping change in synchrony - increase/decrease
x11(13,9)
na +
  geom_point(aes(x = grid_lon.x, y = grid_lat.x, color = meanincdec), size = 4, alpha = .7, data = meansyncdiff) + 
  scale_colour_manual(values = c("blue","red")) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  labs(colour="Change in\nmean spatial\nsynchrony", x="", y="", title="1966-1991 vs. 1992-2017") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59))
#ggsave("figures/grid4_maps/Mean.ALLsp.grid4_incdec_AR_5zeros_n20_60lat.png")

# mapping change in synchrony - magnitude - continuous with viridis color ramp
x11(13,9)
na +
  geom_point(aes(x = grid_lon.x, y = grid_lat.x, color = meansyncdiff), size = 4, alpha = 1, data = meansyncdiff) + 
  scale_color_viridis(name = "Mean\nchange", option = "viridis") +
  theme_classic() +
  theme(text = element_text(size=15)) +
  labs(colour="Change in\nmean spatial\nsynchrony", x="", y="", title="1966-1991 vs. 1992-2017") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59))
#ggsave("figures/grid4_maps/Mean.ALLsp.grid4_change.viridis_AR_5zeros_n20_60lat.png")

#######################################
###### ridgeline plot change in synchrony across all species
#######################################

ridgeplotdata = syncdiff %>%
  select(species = sp.x, syncdiff, incdec, grid_lon.x) %>%
  bind_rows(select(meansyncdiff, species, syncdiff = meansyncdiff, incdec = meanincdec, grid_lon.x)) %>%
  filter(species != "STGR" & species != "MCLO") %>%
  mutate(region = ifelse(grid_lon.x <= -90, "West", "East")) %>%
  group_by(species) %>%
  mutate(medmeandiff = median(syncdiff), n = length(species), 
         meanmeandiff = mean(syncdiff), sdmeandiff = sd(syncdiff), 
         wtmeandiff = 1/(sdmeandiff^2)) %>%
  ungroup() %>%
  mutate(medmeandiff = as.numeric(ifelse(species == "All species","1",medmeandiff))) %>%
  mutate(species = fct_reorder(species, medmeandiff, .desc = T))

# rename codes to be full species names
levels(ridgeplotdata$species) =  
         rev(c("E. Meadowlark","Sprague's Pipit",
            "N. Harrier","Chestnut-col. Longspur",
            "Grasshopper Sparrow","Upland Sandpiper",
            "Ring-necked Pheasant", "Vesper Sparrow",
            "Bobolink", "LeConte's Sparrow", "Horned Lark",
            "Savannah Sparrow", "W. Meadowlark",
            "Lark Bunting", "Dickcissel",
            "Baird's Sparrow","Sedge Wren", 
            "Cassin's Sparrow", "Long-billed Curlew","All species"))

# summary by species for plotting sample sizes and calculating weighted mean/CI
ridgeplot_sum = ridgeplotdata %>%
  arrange(desc(medmeandiff)) %>%
  select(species,medmeandiff,n,meanmeandiff,sdmeandiff,wtmeandiff) %>%
  distinct() %>%
  mutate(n = as.numeric(ifelse(species == "All species", "", n)))

# linear mixed-effects model to find grand mean and CI
library(lme4)

lmedata = filter(ridgeplotdata, species != "All species")
lmemod = lmer(syncdiff ~ 1 + (1|species), data = lmedata)
allspmean = data.frame(mean = fixef(lmemod), 
           se = sqrt(diag(vcov(lmemod, useScale = FALSE)))) %>%
  mutate(ucl = mean + 2*se, lcl = mean - 2*se)

x11(10,9)
ridgeplotdata %>%
  mutate(syncdiff = as.numeric(ifelse(species=="All species", NA, syncdiff))) %>%
ggplot(aes(x = syncdiff, y = species, fill = ..x..), fill="gray") +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, quantile_lines = TRUE, quantiles=2) +
  scale_fill_viridis(name = "Mean\nchange", option = "viridis") +
  labs(x = 'Change in spatial synchrony, 1966-1991 vs. 1992-2017', y = "") +
  annotate("text", x = 0.754, y = (2:21)+0.5, 
           label = c(paste0(ridgeplot_sum$n[2:20]),"grid cells (n) ="), 
           hjust = 1, color = "gray25", size=3) +
  theme_bw() +
  xlim(c(-0.754,0.754)) +
  geom_vline(aes(xintercept = 0), color="red", linetype = 2) +
  theme(axis.text=element_text(size=12, color = "black"), axis.title=element_text(size=12)) +
  geom_segment(aes(x = allspmean$lcl, y = 1.4, xend = allspmean$ucl, yend = 1.4), color="black", size = 1) +
  geom_point(data=allspmean, aes(x = mean, y = 1.4), color="red", shape = 18, size = 3.5)
#ggsave("figures/grid4_maps/ALLsp.grid4_change.ridgeline_AR_5zeros_n20_60latB.png")

#######################################
###### forest plot change in synchrony across all species
#######################################

# means for each species and bootstrapped confidence intervals
# https://stackoverflow.com/questions/38554383/bootstrapped-confidence-intervals-with-dplyr
# note: should ultimately use "spatial block bootstrapping": https://sites.google.com/site/halkingwang/programming/r/spatiallycorrelatederrorsspatialblockbootstrappingbagoflittlebootstrapblb
dotplotdata = ridgeplotdata %>%
#  filter(species != "All species") %>%
  group_by(species) %>%
  do(data.frame(rbind(Hmisc::smean.cl.boot(.$syncdiff))))  %>%
  ungroup() %>%
  mutate(mean.order = as.numeric(ifelse(species == "All species","-1",Mean))) %>%
  mutate(species = fct_reorder(species, mean.order, .desc = F)) %>%
  # input overall mean and CI from the linear mixed-effect model intercept
  mutate(Mean = as.numeric(ifelse(species == "All species",allspmean$mean,Mean))) %>%
  mutate(Lower = as.numeric(ifelse(species == "All species",allspmean$lcl,Lower))) %>%
  mutate(Upper = as.numeric(ifelse(species == "All species",allspmean$ucl,Upper)))

dotplot_n = dotplotdata %>%
  left_join(ridgeplot_sum, by = "species") %>%
  select(species, Mean, Lower, Upper, n) %>%
  mutate(order = ifelse(species == "All species", -1, Mean)) %>%
  arrange(desc(order)) %>%
  mutate(n = ifelse(species == "All species", "", n))
  
x11(10,9)
ggplot(data = dotplotdata) +
  geom_vline(aes(xintercept = 0), color="red", linetype = 2, size = 1.5) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper, y = species, height = .2)) +
  geom_point(aes(x = Mean, y = species, color = Mean), size = 4.5) +
  scale_color_viridis(name = "Mean\nchange", option = "viridis",limits=c(-.1,.21)) +
  labs(x = 'Change in spatial synchrony, 1966-1991 vs. 1992-2017', y = "") +
  annotate("text", x = 0.48, y = c(20.4,20:1), 
           label = c("",paste0(dotplot_n$n[1:20])), 
           hjust = 1, color = "gray25", size=4) +
  annotate("text", x = 0.43, y = 20, label = "n = ", color = "gray25") +
  xlim(c(-0.3,0.48)) +
  theme_bw() +
  theme(axis.text=element_text(size=12, color = "black"), axis.title=element_text(size=12,color="black"))
#ggsave("figures/grid4_maps/ALLsp.grid4_change.forestplot_AR_5zeros_n20.png")

#######################################
###### traits vs. syncrhony
#######################################
# migratory status based on BBS guild classifications: 
# https://www.mbr-pwrc.usgs.gov/bbs/guild/guildlst.html (also used by Zuckerberg et al. 2009)
    # note: coding Ring-necked Pheasant as "Short" (it is actually Res) to combine those categories
# trend info from # https://www.mbr-pwrc.usgs.gov/cgi-bin/tf15.pl
traits = dotplot_n %>%
  filter(species != "All species") %>%
  mutate(mig = c("Neo","Short","Short","Neo","Neo","Neo","Short","Short","Short","Short",
                 "Neo","Neo","Short","Short","Neo","Short","Short","Short","Short")) %>%
  mutate(ord = c("non",rep("pas",10),"non","non",rep("pas",3),"non","pas","pas")) %>%
  mutate(firstegg = c("Apr 1-15","Jul 15-31","Jun 1-15","Jun 1-15","May 15-31","May 15-31","May 1-15",
                       "May 15-31","Jun 1-15","Mar 15-31","May 15-31","May 1-15","Apr 15-30",
                       "May 1-15","Jun 1-15","May 1-15","Apr 15-30","May 1-15","May 15-31")) %>%
  mutate(phen = c("early",rep("late",5),"early",rep("late",2),"early","late",rep("early",3),"late",
                  rep("early",3),"late")) %>%
  mutate(massL = c(640.1, 18.3, 7.73, 17.8, 35.9, 25.2, 89.4, 13, 19.7, 104.4, 28.7, 151, 
                    917, 20.3, 17.34, 24.9, 336, 100.1, 23.6)) %>% # avg. mass of the smaller sex (largest n on breeding grounds, BNA); or pooled avg. if sample size larger
  mutate(massU = c(758.6, 18.3, 8.25, 19.1, 38.5, 28.5, 106, 13.2, 19.7, 111.3, 33, 164, 
                    1263, 20.3, 18.75, 26.5, 513, 123.2, 23.9)) %>% # avg. mass of the larger sex (largest n on breeding grounds, BNA); or pooled avg. if sample size larger
  mutate(mass = (massL + massU)/2) %>%
  dplyr::select(-massL, -massU) %>%
  mutate(trend.66.15 = c(0.17,-0.75,0.40,-2.17,-2.90,-0.36,-1.29,-2.59,-1.36,-2.46,-2.06,0.40,
                         -0.64,-4.19,-2.52,-0.85,-1.21,-3.28,-3.10)) %>%
  mutate(trend.66.91 = c(-0.40,-0.98,0.07,-1.92,-4.74,-1.04,-1.51,-1.93,-1.63,-2.47,-3.15,0.08,
                         -1.31,-4.56,-2.93,-1.06,-1.83,-3.48,-3.20)) %>%
  mutate(trend.92.15 = c(1.17,-1.21,0.62,-2.39,-1.46,0.18,-1.10,-2.39,-1.27,-2.35,-0.80,0.70,0.02,
                         -4.18,-2.01,-0.67,-1.10,-3.24,-2.18)) %>%
  mutate(trend.change = trend.92.15 - trend.66.91) %>%
  mutate(n = as.numeric(n))

# model set for species-level traits vs. change in synchrony
mass.mod = lm(Mean ~ log(mass), data = traits, weights=n); summary(mass.mod)
trend.mod = lm(Mean ~ trend.66.15, data = traits, weights=n); summary(trend.mod)
mig.mod = lm(Mean ~ mig, data = traits, weights=n); summary(mig.mod)
phen.mod = lm(Mean ~ phen, data = traits, weights=n); summary(phen.mod)
glob.mod = lm(Mean ~ log(mass) + trend.66.15 + mig + phen, data = traits, weights = n); summary(glob.mod)

MuMIn::AICc(mass.mod, trend.mod,mig.mod,phen.mod,glob.mod)

x11(9,9)
ggplot(traits) +
  geom_point(aes(x=trend.66.15, y = Mean, size = as.numeric(n))) +
  geom_smooth(aes(x = trend.66.15, y = Mean), method = "lm", se=T, color="black") +
  cowplot::theme_cowplot() +
  labs(size = "n", x = "BBS Trend (1966-2015)", y = "Change in spatial synchrony")
#ggsave("species.traits.vs.change.synchrony.jpg")

# conserve-o-gram of all grassland species
# note: soon add other species: McCownan's LS: -7.60 to -5.20
x11(12,12)
traits %>%
  mutate(species = fct_reorder(species, trend.92.15)) %>%
ggplot() +
  geom_segment(aes(x = trend.66.91, y = species, xend = trend.92.15, yend = species), size = 1, 
               color = "gray50", arrow = arrow(length = unit(.5,"cm"))) +
  geom_segment(aes(x = 0, y = 1, xend = 0, yend = 19), size = 2, color = "black", linetype = 2) +
  geom_point(aes(x=trend.66.91, y = species), size = 4, color = "gray50") +
  geom_point(aes(x=trend.92.15, y = species), size = 4, color = "firebrick") +
  labs(x = "BBS Pop. Trends (% change per year)", y = "") +
  labs(title = "Change in trend, 1966-1991 vs. 1992-2015") +
  theme(axis.text = element_text(size = 12), title = element_text(size = 12)) +
  annotate("text", x = c(-.25, .25), y = 1, 
         label = c("decreasing","increasing"),
         hjust = c(1,0), color = "gray25", size=4)
# ggsave("conservatogram_grassland_birds.jpg")




#######################################
###### making polygon file with synchrony change data, for cluster analysis
#######################################

# make a grid of spatial polygons to match up with data
library(raster)
rr <- raster::raster(ext = extent(-177, -52, 24, 69), res=c(2,2))
values(rr) <- 1:ncell(rr)
p <- rasterToPolygons(rr)
attributes = tibble(grid4_lon = coordinates(p)[,1], grid4_lat = coordinates(p)[,2]) %>%
  mutate(join = 1, grid.x = paste(trunc(grid4_lat),trunc(grid4_lon),sep="")) %>%
  mutate(grid.x = gsub("-","_",grid.x)) %>%
  left_join(meansyncdiff, by="grid.x") 

p$grid = attributes$grid.x
p$meansyncdiff = attributes$meansyncdiff
p$n = attributes$n
p$meandiffcat = attributes$meandiffcat

meansyncdiff.poly = p[is.na(p$n)!=T,]

plot(meansyncdiff.poly, col="gray")
head(meansyncdiff.poly)
# writeOGR(meansyncdiff.poly, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "meansyncdiff_poly", driver="ESRI Shapefile")

# plot meansyncdiff polygon grid in ggplot
meansyncdiff.data = fortify(meansyncdiff.poly, region = "grid") %>%
  left_join(rename(meansyncdiff.poly@data, id=grid),by = "id")

# read in cluster shapefile (local Moran's I analysis performed in GeoDa)
clusters = readOGR("data/meansyncdiff_cluster.shp")
cluster_fort = fortify(clusters, region = "grid") %>% left_join(rename(clusters@data,id=grid), by = "id")

x11(13,9)
na + geom_polygon(data=meansyncdiff.data, aes(x=long, y = lat, group=group, fill = meansyncdiff)) + 
  scale_fill_viridis(name = "Mean\nchange", option = "viridis") +
  theme_classic() +
  theme(text = element_text(size=15)) +
  labs(colour="Change in\nmean spatial\nsynchrony", x="", y="", title="1966-1991 vs. 1992-2017") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59)) +
  geom_polygon(data = filter(cluster_fort, LISA_CL==1), aes(x = long, y = lat, group=group), color = "red",fill="transparent") +
  geom_polygon(data = filter(cluster_fort, LISA_CL==2), aes(x = long, y = lat, group=group), color = "white",fill="transparent")
#ggsave("figures/grid4_maps/Mean.ALLsp.grid4_change.viridis.poly_AR_5zeros_n20_60lat.png")

# who contributes to main increasing synchrony cluster?
# make a list of the 13 grid cells in that cluster
grid_inc = data.frame(grid=unique(filter(cluster_fort, LISA_CL==1)$id)) %>% 
  filter(grid != "50_96")
# count the number of grid cells (out of 13) occupied by each species
who_inc = syncdiff %>%
  filter(grid.x %in% grid_inc$grid) %>%
  group_by(sp.x) %>%
  summarise(contribution = length(sp.x)) %>%
  arrange(desc(contribution))

# who contributes to main decreasing synchrony cluster?
# make a list of the 13 grid cells in that cluster
grid_dec = data.frame(grid=unique(filter(cluster_fort, LISA_CL==2)$id)) %>% 
  filter(grid != "44_66")
# count the number of grid cells (out of 12) occupied by each species
who_dec = syncdiff %>%
  filter(grid.x %in% grid_dec$grid) %>%
  group_by(sp.x) %>%
  summarise(contribution = length(sp.x)) %>%
  arrange(desc(contribution))

#######################################
###### Links to data sources, code, papers, etc.
#######################################
# BBS Route-level data fields: ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/Summary.txt
# BBS RouteID info: ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/RegionCodes.txt
# BBS phys strata map: https://www.pwrc.usgs.gov/bbs/StrataNames/strata_map_small.htm
# BBS Route info txt with BCR and strata codes: ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/RouteInf.txt
# BBS BCR shapefile: http://www.pwrc.usgs.gov/bba/index.cfm?fa=bba.getdata
# Smith et al. sup data for CAR model and indices of aerial insectivores: 
#       https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0130768#sec016
#       https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0130768.s005&type=supplementary
# Michel et al. README file: https://datadryad.org/bitstream/handle/10255/dryad.96773/README.txt?sequence=2
# Michel et al. paper with link to sup code files: https://onlinelibrary.wiley.com/doi/abs/10.1111/ecog.01798

