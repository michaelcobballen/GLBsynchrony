### THIS CODE IS FOR CALCULATING PAIRWISE HAVERSINE DISTANCES BETWEEN US STATE CENTROIDS
    ### THIS WAS SAVED AS DSTATE.CSV AND USED IN SEPARATE ANALYSES OF SPATIAL SYNCHRONY 
    ### IN BIRD POPULATIONS AND AGRICULTURE
    ### CAN BE ADAPTED FOR ANY PURPOSE WHERE YOU WANT TO CALCULATE INTER-POLYGON DISTANCES

library(rgdal)
library(rgeos)
library(geosphere)


# Get the census shapefile for  US states and territories
# download("https://www2.census.gov/geo/tiger/GENZ2014/shp/cb_2014_us_state_5m.zip", dest="states.zip", mode="wb") 
# unzip("states.zip", exdir = getwd())
states <- readOGR("data/cb_2014_us_state_5m.shp")
states.proj <- CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
#proj4string(states) <- wgs1984.proj
states.equal <- spTransform(states, states.proj)

plot(states.equal,col="gray")

# get centroids of all states
# NOTE: the 'coordinates' function produces identical centroids to packages such as rgeos
plot(states.equal)
points(coordinates(states.equal),pch=1)

states.latlong = spTransform(states.equal,CRS("+proj=longlat +datum=WGS84"))

plot(states.latlong)
points(coordinates(states.latlong))

(df.states=data.frame(name = states.latlong$STUSPS, 
                      lat = coordinates(states.latlong)[,2],
                      lon = coordinates(states.latlong)[,1]))

#calculate Haversine distances between states
dmat = round(distm(df.states[,c("lon","lat")])/1000)
colnames(dmat) <- df.states$name; dmat
dmat2 = data.frame(state2=df.states$name, dmat); dmat2

dstate= data.frame(reshape2::melt(dmat2)); dstate
dstate = dstate[dstate$value!=0,]
colnames(dstate)[2:3] = c("state1","dist"); head(dstate); nrow(dstate)
dstate$pairid = paste(dstate$state2,"-",dstate$state1,sep=""); head(dstate); nrow(dstate)

binlabels=c("<100","100-300","301-500","501-700","701-900","901-1100","1101-1300",
            "1301-1500",">1500")

dstate$dcat=cut(dstate$dist,c(0,100,300,500,700,900,1100,1300,1500,100000),
                  labels=binlabels)
head(dstate)

### creating region column: east vs. west of the Mississippi River
dstate$reg1 = ifelse(dstate$state1 %in% c("AL","CT","DE","FL","GA","IL","IN",
  "KY","MA","MD","ME","MI","MS","NC","NH","NJ","NY","OH","PA","RI","SC","TN", 
  "VA","VT","WI","WV"), "east", "west")

dstate$reg2 = ifelse(dstate$state2 %in% c("AL","CT","DE","FL","GA","IL","IN",
  "KY","MA","MD","ME","MI","MS","NC","NH","NJ","NY","OH","PA","RI","SC","TN", 
  "VA","VT","WI","WV"), "east", "west")

dstate$reg = ifelse(dstate$reg1=="east"& dstate$reg2=="east","east","west")
dstate$reg = ifelse(dstate$reg1!=dstate$reg2,"both",dstate$reg)
head(dstate); nrow(dstate)

dstate = select(dstate, pairid,state1,state2,dist,dcat,reg)
head(dstate)

#write.csv(dstate,"data/dstate.csv")



