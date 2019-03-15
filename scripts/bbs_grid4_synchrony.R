###############################################################################################
###
### An attempt to generate the BBS indices myself
### for 2x2 degree grid-cell strata via methods in Michel et al. Ecography (2016)
### to use in the conditional autoregressive (CAR) spatial model described with code in Smith et al. PlosOne (2015)
###
###############################################################################################

# see Links at bottom for data sources, papers, etc.
# to do:
# 1) include raw inclusion variables into the final data products 
      # so they can be excluded as needed on the fly in the final synchrony step
      # without re-running data
# 2) move grid-cell polygon creation code to separate script file

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

#######################################
### BBSync function
#######################################

bbsync = function(spcode, yearstart, yearend, method){

# if not a function, then input these parameters:
#spcode = "5460" # AOU species number; #  spcode = s05460 for GRSP, s05010 for EAME, s04940 for BOBO, s05420 for SAVS
#yearstart = 1992
#yearend = 2017
#method = "change" # "AR" = autoregressive, "change" = first difference
minyearsdata = 10 # each grid cell contains at least one route with at least this many years of data (10 is based on criteria of Michel et al. Ecography 2016)
minsproutes = 2 # each grid cell contains at least this many routes at which the species was detected (2 is based on criteria of Michel et al. Ecography 2016)
maxzeros = 5 # maximum number of years allowable where the mean count for the species is zero (I chose five)
minyearscor = 20 # minimum number of years two grid cells must have in common to calculate a correlation coefficient between them (Hanski used 7; I chose 20, more conservative)

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
  select(RouteDataID, CountryNum = CountryNum.x, StateNum = StateNum.x, Route = Route.x,
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
  left_join(distinct(select(brw.route.year,-sp.number,-Year,-RouteDataID,-ObsN,-ObsN.rid,-ObsN.rid.year,-yr,-ObsN.rid.first,-ObsN.rid.count)),by="rid") %>%
  mutate(sp.present = as.numeric(ifelse(sp.yrs>0,1,0))) %>%
  arrange(lat) %>%
  ungroup() 

# use this to make a shapefile of all BBS routes if desired
#brw.route.spatial = sp::SpatialPointsDataFrame(brw.route[,c("lon","lat")],brw.route)
#proj4string(brw.route.spatial) = CRS("+init=epsg:4326")
#plot(brw.route.spatial)
# writeOGR(brw.route.spatial, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "bbs_routes1", driver="ESRI Shapefile")

# summarize data by 2x2 grid cell
# including whether cells meet the "minyearsdata" and "minsproutes" criteria of Michel et al. (Ecography)
grid4.summary = brw.route %>%
  group_by(grid4) %>%
  summarise(grid4_lat=mean(grid4_lat), grid4_lon=mean(grid4_lon), num.routes = length(sp.present), 
            sp.mean = mean(sp.mean), num_sp_rtes = sum(sp.present), 
            min_yr_rtes = sum(yrs.data>minyearsdata-1), nobservers = max(nobservers)) %>%
  mutate(sp.min.rtes.yes = as.numeric(ifelse(num_sp_rtes > minsproutes-1,1,0))) %>%
  mutate(min.years.yes = as.numeric(ifelse(min_yr_rtes > 0,1,0))) %>%
  mutate(include_strat = as.numeric(ifelse(sp.min.rtes.yes + min.years.yes == 2,1,0))) %>%
  mutate(nonzeroweight = num_sp_rtes/num.routes) %>%
  select(-sp.min.rtes.yes,-min.years.yes)

# use this to make a point shapefile of centroids for all grid cells included based on criteria in Michel et al. (Ecography)
    # this is to compare with the map with exclusions based on # of zeros below
#grid4.summary.spatial = sp::SpatialPointsDataFrame(grid4.summary[grid4.summary$include_strat==1,c("grid4_lon","grid4_lat")],grid4.summary[grid4.summary$include_strat==1,])
#proj4string(grid4.summary.spatial) = CRS("+init=epsg:4326")
#plot(grid4.summary.spatial); nrow(grid4.summary.spatial@data)
# writeOGR(grid4.summary.spatial, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "bbs_grid4", driver="ESRI Shapefile")

# add grid-level information to the route-year data
# remove routes in cells excluded based on criteria in Michel et al. (Ecography)
route.year.grid4 = brw.route.year %>%
  left_join(select(grid4.summary, grid4, include_strat), by = "grid4") %>%
  filter(include_strat==1) %>%
  select(count = sp.number, grid4, obser = ObsN.rid, firstyr = ObsN.rid.first, 
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

# plot or make shapefile of grid cell centroids included based on all criteria (if needed)
#grid4.grsp.spatial = SpatialPointsDataFrame(grid4.summary2[,c("grid4_lon","grid4_lat")],grid4.summary2)
#proj4string(grid4.grsp.spatial) = CRS("+init=epsg:4326")
#plot(grid4.grsp.spatial); nrow(grid4.grsp.spatial@data)
# writeOGR(grid.grsp.spatial, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "grsp_grid", driver="ESRI Shapefile")

### subsetting the grid-year data based on 2nd criteria (number of years with zeros)
  # this had already been subsetted to exclude based on Michel et al. criteria 
grid4.year = grid4.year %>%
  left_join(grid4.include2,by = "grid4") %>%
  filter(include2 == 1)


################# 
### calculating inter-grid distances for 2 X 2 degree lat/long grid
################# 

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
  select(grid4B,grid4A = variable, dist = value) %>%
  mutate(grid4A = gsub("X","",grid4A)) %>%
  mutate(pairid = paste(grid4B,".",grid4A,sep=""))


#################
# synchrony data formatting
#################

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
  select(grid4=grid4.x, year=year.x, grid4.year, mean.count) %>%
  mutate(log.count = log(mean.count+1), logcountlag1 = lag(log.count)) %>% # log+1 transform count ala Hanski and create lag for AR and differencing analyses
  mutate(logcountlag1 = as.numeric(ifelse(year==yearstart,"NA",logcountlag1))) %>% # assign "NA" to first year which should have no lagged value
  filter(is.na(log.count) == F & is.na(logcountlag1) == F) %>% # temporarily remove NAs to get residuals and first diffs of logged counts
  mutate(change = log.count - logcountlag1) %>% # first diff of logged counts, ~ % change
  group_by(grid4) %>%
  mutate(AR1.R = resid(lm(log.count~logcountlag1))) %>% # calculate residuals from first order autoregression (Hanski-style)
  ungroup() %>%
  right_join(grid4.year.list, by = "grid4.year") %>% # merge back with the full grid/year list to allow matching in correlation analyses
  select(grid4 = grid4.y, year = year.y, grid4.year, mean.count, log.count,logcountlag1, change, AR = AR1.R)

### put bird population data into 'wide' format for correlation analysis (rows = years, columns = grid cells)
grid4.bird.mat = 
  grid4.year.full %>%
  dplyr::select(grid4,year,method) %>%
  spread(grid4,method) %>%
  filter(year > yearstart & year < yearend+1) %>%
  select(-year)

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
  select(grid4A=Var1, grid4B=Var2, cor=value, pval=value.1, n) %>%
  mutate(pairid = paste(grid4A,".",grid4B,sep="")) %>%
  filter(n>minyearscor-1) %>% # excludes grid pairs with less than XX years in common (Hanski's criteria was 7)
  filter(grid4A!=grid4B) %>% # excludes self-correlations
  filter(is.na(pval)==F) %>% # excludes NA correlations resulting from too many zeros if these weren't removed earlier - they were! So nothing changes.
  distinct(pval, .keep_all = TRUE) # should cut number in half, and yes it does! (eventually add check and warning message if it doesn't)

grid4.bird.cordist = grid4.bird.cor %>%
  left_join(g4,by='pairid') %>%
  select(pairid, grid4A=grid4A.x, grid4B=grid4B.x, cor,pval,n,dist) %>%
  mutate(sig = ifelse(pval < 0.05,"Y","N"), sig1 = as.numeric(ifelse(pval < 0.05,"1","0"))) %>%
  distinct(cor, .keep_all = TRUE) 

return(grid4.bird.cordist)
}






#######################################
###### summarizing function results in dataframe to graph more easily
#######################################
# AOU species number; #  spcode = s05460 for GRSP, s05010 for EAME, s04940 for BOBO, s05420 for SAVS
                      # s07550 for WOTH, s06740 for OVEN, 
                      # https://www.pwrc.usgs.gov/bbl/manual/speclist.cfm
allsync = rbind(
grsp.66.91 = cbind.data.frame(bbsync("5460",1966,1991,"AR"),sp="GRSP",per="1966-1991"),
bobo.66.91 = cbind.data.frame(bbsync("4940",1966,1991,"AR"),sp="BOBO",per="1966-1991"),
eame.66.91 = cbind.data.frame(bbsync("5010",1966,1991,"AR"),sp="EAME",per="1966-1991"),
grsp.92.17 = cbind.data.frame(bbsync("5460",1992,2017,"AR"),sp="GRSP",per="1992-2017"),
bobo.92.17 = cbind.data.frame(bbsync("4940",1992,2017,"AR"),sp="BOBO",per="1992-2017"),
eame.92.17 = cbind.data.frame(bbsync("5010",1992,2017,"AR"),sp="EAME",per="1992-2017")
)

allsync1 = rbind(allsync,
savs.66.91 = cbind.data.frame(bbsync("5420",1966,1991,"AR"),sp="SAVS",per="1966-1991"),
savs.92.17 = cbind.data.frame(bbsync("5420",1992,2017,"AR"),sp="SAVS",per="1992-2017")
)

#######################################
###### graphing synchrony - continous / smoothed
#######################################

X11(13,9)
ggplot() + 
  geom_smooth(data= allsync1, aes(x=dist,y=cor, color=sp, linetype=per),method="auto", se=F) +
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
###### mapping synchrony
#######################################

mapsync = 
  rbind(data.frame(allsync1, grid = allsync1$grid4A), data.frame(allsync1, grid = allsync1$grid4B)) %>%
  filter(dist<400) %>%
  group_by(grid, per, sp) %>%
  summarise(meancor = mean(cor), n = length(cor)) %>%
  ungroup() %>%
  mutate(grid_lat = as.numeric(substr(grid,1,2)), grid_lon = as.numeric(paste("-",substr(grid,4,10),sep=""))) %>%
  mutate(corcat = cut(meancor, breaks = c(-1,0,0.1,.2,.3,.4,.5,.6)))

can = readOGR("data/province.shp")

canada.fort = ggplot2::fortify(can) %>%
  filter(lat <60.0001)

na <- ggplot() +
  borders("world", xlim = c(-150, -60), ylim = c(30, 40), colour = "gray85", fill = "gray80")  +
  #borders("state", colour = "gray85", fill = "gray80") +
  #geom_polygon(data=canada.fort, aes(x = long, y = lat, group=group), colour = "gray85", fill = "gray80") +
  theme_map() + 
  coord_map('albers', lat0=30, lat1=60)

x11(21,7)
ggplot(filter(mapsync,per=="1992-2017")) +
  geom_point(aes(x = grid_lon, y = grid_lat, color = corcat), size = 3, alpha = .7) + 
  scale_colour_manual(values = c("#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")) +
  theme_bw() +
  theme(text = element_text(size=15)) +
  facet_grid(~sp)

x11(21,7)
ggplot(filter(mapsync,per=="1966-1991")) +
  geom_point(aes(x = grid_lon, y = grid_lat, color = corcat), size = 3, alpha = .7) + 
  scale_colour_manual(values = c("#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")) +
  theme_bw() +
  theme(text = element_text(size=15)) +
  facet_grid(~sp)



#################
# below is for making a spatial grid file to extract spatial data for use in the OpenBUGS conditional autoregressive (CAR) model via Smith et al. PlosONE 2015 code
#################

# make a grid of spatial polygons to match up with data
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

grid4.grsp = p[p$include_strat==1 & is.na(p$include_strat)==F,]
grid4.grsp@data = grid4.grsp@data %>%
  select(layer, grid4, include_strat, include2, nonzeroweight, nobservers, strata) %>%
  mutate(strat = 1:nrow(grid4.grsp@data))
head(grid4.grsp@data); nrow(grid4.grsp@data)

# plot or save shapefile of included grid cells if you want
plot(grid4.grsp, col="gray")
#writeOGR(grid4.grsp, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "grsp_grid4_poly", driver="ESRI Shapefile")

# project into Lambert conformal conic to allow plotting in OpenBUGS (with meters scale) if needed
#proj = "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=10000000 +y_0=10000000 +datum=NAD83 +units=m +no_defs"
grid4.grsp.lamb = spTransform(grid4.grsp, proj)
plot(grid4.grsp.lamb, col="gray")
head(grid4.grsp.lamb); nrow(grid4.grsp.lamb)
# writeOGR(grid4.grsp.lamb, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "grsp_grid4_poly_lamb", driver="ESRI Shapefile")

# changing names of the ID variable of the polygons to grid labels so names are visible in OpenBUGS map
# https://eco-data-science.github.io/spatial_analysis2_R/
sapply(slot(grid4.grsp.lamb0, "polygons"), function(x) slot(x, "ID"))
for (i in 1:length(slot(grid4.grsp.lamb0, "polygons"))){
  slot(slot(grid4.grsp.lamb0, "polygons")[[i]], "ID") = grid4.grsp.lamb0$grid4[i]
}
sapply(slot(grid4.grsp.lamb0, "polygons"), function(x) slot(x, "ID"))


### converting spatial data to geobugs format if needed
sp2WB(grid.grsp.lamb0, "grid.grsp.lamb.txt")




### Links to data sources, code, papers, etc.
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
