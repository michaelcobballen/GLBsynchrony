### GETTING STATE CENTROIDS AND DISTANCES BETWEEN STATES

#download("https://www2.census.gov/geo/tiger/GENZ2014/shp/cb_2014_us_state_5m.zip", dest="states.zip", mode="wb") 
#unzip ("states.zip", exdir = getwd())
states <- rgdal::readOGR("data/cb_2014_us_state_5m.shp")
states.proj <- CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
#proj4string(states) <- wgs1984.proj
states.equal <- spTransform(states, states.proj)
reg = subset(states.equal,STUSPS %in% c("NH","NY","PA","NJ","ME","CT","RI","MA","VT",
                                        "DE","MD","WV","VA","NC","KY","TN", "OH","AL",
                                        "MS","SC","VA","NC", "MI","GA","IN","IL","WI"))
plot(states.equal,col="gray")
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

binlabels=c("100-300","301-500","501-700","701-900","901-1100","1101-1300",
            "1301-1500",">1500")

dstate$dcat=cut(dstate$dist,c(0,300,500,700,900,1100,1300,1500,10000),
                  labels=binlabels)
    head(dstate)

gbin = cbind.data.frame(pctsig=100*tapply(gcordist$sig1,gcordist$dcat,mean),
                        n=tapply(gcordist$sig1,gcordist$dcat,length))
gbin$bin=factor(rownames(gbin),levels=binlabels); gbin


#write.csv(dstate,"dstate.csv")
