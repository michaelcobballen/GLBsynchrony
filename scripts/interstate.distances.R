### GETTING STATE CENTROIDS AND DISTANCES BETWEEN STATES

#download("https://www2.census.gov/geo/tiger/GENZ2014/shp/cb_2014_us_state_5m.zip", dest="states.zip", mode="wb") 
#unzip ("states.zip", exdir = getwd())
states <- rgdal::readOGR("cb_2014_us_state_5m.shp")
states.proj <- CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
#proj4string(states) <- wgs1984.proj
states.equal <- spTransform(states, states.proj)
reg = subset(states.equal,STUSPS %in% c("NH","NY","PA","NJ"
                                        ,"DE","MD","WV","VA","NC","KY","TN", "OH"))
plot(states.equal,col="gray",main="The USA")
plot(reg,col="gray",main="Study region")

head(states.equal)

# get centroids of all states
cents = rgeos::gCentroid(reg,byid=TRUE)
head(cents)
plot(states.equal)
points(coordinates(reg),pch=1)
points(cents,pch=2)

reg.latlong = spTransform(reg,CRS("+proj=longlat +datum=WGS84"))
states.latlong = spTransform(states.equal,CRS("+proj=longlat +datum=WGS84"))

plot(reg.latlong)
plot(states.latlong)
points(coordinates(reg.latlong))
points(coordinates(states.latlong))
head(reg.latlong)
head(states.latlong)

(df.states=data.frame(name = states.latlong$STUSPS, 
                      lat = coordinates(states.latlong)[,2],
                      lon = coordinates(states.latlong)[,1]))

dmat = round(distm(df.states[,c("lon","lat")])/1000)
colnames(dmat) <- df.states$name; dmat
dmat2 = data.frame(state2=df.states$name, dmat); dmat2

dstate= data.frame(reshape2::melt(dmat2)); dstate
dstate = dstate[dstate$value!=0,]
colnames(dstate)[2:3] = c("state1","dist"); head(dstate); nrow(dstate)
dstate$pairid = paste(dstate$state2,"-",dstate$state1,sep=""); head(dstate); nrow(dstate)

#write.csv(dstate,"dstate.csv")
