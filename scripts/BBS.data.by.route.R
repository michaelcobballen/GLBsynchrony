# practice attempts using code from blog post in "Spatial Ecology and R"
# http://spatialecology.weebly.com/r-code--data/31
library(sp)
library(downloader)
library(maptools)
library(mapmisc)

#Get all NJ sightings from the BBS (10 route summaries and Stop Totals)
temp <- tempfile()
download.file("ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/States/NJersey.zip",destfile=temp,method="libcurl")
nj <- read.csv(unzip(temp))
unlink(temp)

#Get data on BBS route locations
temp <- tempfile()
download.file("ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/Archivefiles/Version2016v0/Routes.zip",destfile=temp,method="libcurl")
routes <- read.csv(unzip(temp))
unlink(temp)

#Get weather and quality data on BBS route locations
temp <- tempfile()
download.file("ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/Weather.zip",destfile=temp,method="libcurl")
weather <- read.csv(unzip(temp))
unlink(temp)


#Merge the two data sets to combine sightings with route locations and plot the sightings
data = merge(nj, routes); head(data)
#data = data[data$CountryNum == 840,]  #Occurrences in United States only
head(data)
data = data[data$AOU == 5460, c("Year","StopTotal","Longitude", "Latitude","AOU")]; head(data)  #Sp X occurrences only
coordinates(data) = c("Longitude","Latitude")
wgs1984.proj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4string(data) <- wgs1984.proj
download("http://www2.census.gov/geo/tiger/GENZ2014/shp/cb_2014_us_nation_20m.zip", dest="us.zip", mode="wb") 
unzip ("us.zip", exdir = getwd())
us <- rgdal::readOGR("cb_2014_us_nation_20m.shp")
us.atlas.proj <- CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
proj4string(us) <- wgs1984.proj
us.equal <- spTransform(us, us.atlas.proj)
data.equal <- spTransform(data, us.atlas.proj)
plot(us.equal,col="gray",main=paste("GRSP Sightings in the United States during ",range(data$Year)[1],"-",range(data$Year)[2]," (n = ",nrow(data),")",sep=""))
points(data.equal,pch=16,cex=0.6,col="red")
par(xpd=TRUE)
legend("bottomleft",inset=c(-0.0425,0.55),cex=1,legend=c("Occurrence"),pch=c(16),col=c("red"),bty="n")
scaleBar(data.equal,pos="bottomleft",bg="white",inset=c(-0.0425,0.45),bty="n")
legend("bottomleft",inset=-0.1,legend="Note: data are displayed in a Lambert Azimuthal Equal Area projection.",bty="n")
