# GIS IN R - Introduction to GIS (Section 2)
========================================================

# Here are some of the stuff we learned in Introduction to GIS course.
# GIS is probably better for data collection, data manipulation/editting,
# and creating a visually pleasing map layout. But, R might be a good 
# option for sharing & publishing your analysis, including both processes 
# and results.

# what does a package do?
### library(help = "sp") 
### help(package = "sp") 
### packageDescription("sp") 

# get ready for doing something spatial - by installing packages/libraries
# install.packages("rgdal")
# install.packages("spdep")
# install.packages("maptools")
# install.packages("raster")
library(rgdal)
library(maptools)
library(spdep)
library(raster)

# http://stackoverflow.com/questions/15155814/check-if-r-package-is-installed-then-load-library
install_load <- function (package1, ...)  {   
  # convert arguments to vector
  packages <- c(package1, ...)
  # start loop to determine if each package is installed
  for(package in packages){
    # if package is installed locally, load
    if(package %in% rownames(installed.packages()))
      do.call('library', list(package))
    # if package is not installed locally, download, then load
    else {
      install.packages(package)
      do.call("library", list(package))
    }
  } 
}

install_load(c("sp","rgdal","rgeos","maptools","raster","spdep"))

# change the fowlling working directory accordingly
setwd("F:/Geog_Intro2GIS/2013/Misc/R")

# download and unzip the exercise data
download.file("http://home.uchicago.edu/~cmaene/data.zip",destfile="data.zip")
unzip("data.zip")
setwd("data")

##### MAPPING POINTS, LINES, POLYGONS & RASTER ALL AT ONCE #####

# read (upload) data
roads<-readShapeSpatial("roads.shp")

# what exactly did we do?
help(readShapeSpatial)

# take a look at the spatial object
class(roads)
summary(roads) # data summary
plot(roads) # look at the spatial object graphically
dim(roads) # check the number of obsevations & variables
names(roads) # variable names
head(roads) # view the first 6 observations
tail(roads) # view the last 6 observations
str(roads@data) # another way of looking at the data (data only - as opposed to "lines" object)

# add more data
studyarea<-readShapeSpatial("studyarea.shp")
elevation<-raster("elevation.asc")
lilac<-read.csv("LilacTreesGPS.csv",sep=",")

# lay data over the elevation raster
plot(elevation)
lines(roads)
points(y=lilac$y, x=lilac$x, col="purple",cex=0.8) # add the GPS points
plot(studyarea, border="red", pch=10, add=T)

##### DATUM & PROJECTIONS #####

# read data
wgs84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
graticule<-readShapeSpatial("latlong.shp", proj4string=wgs84)
countries<-readShapeSpatial("country.shp", proj4string=wgs84)
cities<-readShapeSpatial("cities.shp", proj4string=wgs84)
route<-readShapeSpatial("route.shp", proj4string=wgs84)

# plot data - map
par(mar = c(1, 1, 1, 1))
plot(graticule, col="grey", cex=.1)
plot(countries, add=T)
citiesxy<-coordinates(cities)
plot(route, col="red", add=T)
text(citiesxy, labels=cities@data$NAME, col="red", cex=.5)

# use "maps" library - comes with data and functions
# install.packages("maps")
library(maps)
map("world")
wgs84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
cities<-readShapeSpatial("cities.shp", proj4string=wgs84)
route<-readShapeSpatial("route.shp", proj4string=wgs84)
citiesxy<-coordinates(cities)
plot(route, col="red", add=T)
text(citiesxy, labels=cities@data$NAME, col="red", cex=.5)

# change the map view
proj.type <- "azequidistant"
proj.orient <- c(90,0,0)
map("world", proj=proj.type, orient=proj.orient, resolution=0, wrap=T)
map.grid(col="black", labels=F, lty=1)
proj.coords <- mapproject(citiesxy[,1], citiesxy[,2], proj=proj.type, orient=proj.orient)
text(proj.coords, labels=cities@data$NAME, col="red", cex=0.5)

# another method - reproject all data
azequaldistant<-CRS("+proj=aeqd +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 
                    +datum=WGS84 +units=m +no_defs")
graticule2<-spTransform(graticule,CRS=azequaldistant)
countries2<-spTransform(countries,CRS=azequaldistant)
route2<-spTransform(route,CRS=azequaldistant)
plot(graticule2, col="grey")
plot(countries2, add=T)
plot(route2, col="red", add=T)
text(proj.coords, labels=cities@data$NAME, col="red", cex=0.5)

##### THEMATIC MAPPING #####

# install.packages("RColorBrewer")
# install.packages("classInt") 
library(RColorBrewer)
library(classInt) 
par(mar = c(1, 3, 1, 1))
display.brewer.all() # see the color palette

# get the spatial data ready
stateplane<-CRS("+proj=tmerc +lat_0=36.66666666666666 +lon_0=-88.33333333333333 
                +k=0.999975 +x_0=300000 +y_0=0 +ellps=GRS80 +datum=NAD83 
                +to_meter=0.3048006096012192")
tracts<-readShapeSpatial("tracts.shp", ID="CTIDFP00", proj4string=stateplane)

# equal-frequency class intervals
age_under5 <- tracts@data$AGEU5
nclr <- 8
colors <- brewer.pal(nclr,"YlOrBr")
class <- classIntervals(age_under5, nclr, style="quantile")
colorcode <- findColours(class, colors)

par(mar=c(1,1,2,2))
plot(tracts, col=colorcode)
title(main="Age Under 5: Quantile (Equal-Frequency)")
legend("topright", legend=names(attr(colorcode, "table")), fill=attr(colorcode, "palette"), cex=0.6, bty="n")

# dot density map
hispanic <- tracts@data$HISPANIC
dotper<-hispanic/500
dothispanic<-dotsInPolys(tracts, as.integer(dotper), f="random")
plot(tracts)
plot(dothispanic, col="brown", pch=19, cex=0.2, add=T)
title(main="Dot Density Map: dot=500 persons")

# add dots for black & a legend
black <- tracts@data$BLK
dotper2<-black/500
dotblack<-dotsInPolys(tracts, as.integer(dotper2), f="random")
plot(dotblack, col="blue", pch=19, cex=0.2, add=T)
legend("bottomleft", legend=c("Hispanic", "Black"), fill = c("brown","blue"), bty = "n")

##### OVERLAY ANALYSIS #####

# read data
wgs84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
stateplane<-CRS("+proj=tmerc +lat_0=36.66666666666666 +lon_0=-88.33333333333333 
                +k=0.999975 +x_0=300000 +y_0=0 +ellps=GRS80 +datum=NAD83 
                +to_meter=0.3048006096012192")
apt<-readShapeSpatial("Overlay_apt.shp", ID="APTID", proj4string=wgs84)
apartments<-spTransform(apt,CRS=stateplane)
buffer<-readShapeSpatial("Buffer.shp", ID="OBJECTID", proj4string=stateplane)
ChicagoTracts<-readShapeSpatial("tracts.shp", ID="CTIDFP00", proj4string=stateplane)

# run simple overlay analysis
overlay<-overlay(buffer,apartments) # appears to be old but still works
apartments@data<-cbind(apartments,overlay)

# subset the apartments where OBJECTID is not NULL (i.e. intersected with buffer)
within<-!is.na(apartments@data[,c("OBJECTID")])

# see the result
plot(ChicagoTracts)
plot(buffer, col="gray", add=T)
plot(apartments, col="blue", add=T)
plot(apartments[within,], col="red", pch=19, add=T)
title(main="Apartments (red) near the CTA Stations")

##### CALCULATE DISTANCE #####

# install.packages("SDMTools")
library(SDMTools)

# read origin-destination data
distancedata<-read.csv("od25.csv",sep=",",header=T)

# distance tool returns a data frame
head(distancedata) # check the first 6 rows
result<-data.frame(lat1=numeric(),lon1=numeric(),lat2=numeric(),lon2=numeric(),distance=numeric())

# run the function repeatedly
for(i in 1:nrow(distancedata)){
  y1<-distancedata[i,2]
  x1<-distancedata[i,1]
  y2<-distancedata[i,4]
  x2<-distancedata[i,3]
  rdf<-distance(lat1=y1, lon1=x1, lat2=y2, lon2=x2) # add variable names (header)
  result<-rbind(result,rdf)
}
# distance is in meters
head(result) # check the first 6 rows

##### MEAN CENTER #####

# read data
snowmap<-raster("snowmap.asc")
death<-readShapeSpatial("JohnSnow_CholeraDeaths.shp")
pumps<-readShapeSpatial("JohnSnow_Pumps.shp")

# get mean center
meanx<-mean(death$X)
meany<-mean(death$Y)
plot(snowmap)
plot(death, pch=19, add=T)
plot(pumps, col="red", add=T)
points(meanx, meany, pch=23, col="red", bg="red")
title(main="John Snow's Map Analysis (mean center)")

##### INTERPOLATION #####

#install.packages("gstat")
library(gstat)

# read data and plot
par(mar = c(3, 3, 3, 3))
states48 <- readShapeSpatial("States48.shp")
weatherst <- readShapeSpatial("WeatherSt.shp")
plot(weatherst, pch=20, cex=0.5, axes=T)
plot(states48, add=T)

# create a regularly spaced cell grid
n <- 1 # grid cell increment by 1 degree
stcoord<-coordinates(weatherst)
xmin <- min(stcoord[,1], na.rm=TRUE)
xmax <- max(stcoord[,1], na.rm=TRUE)
x <- seq(xmin, xmax, by=n)
ymin <- min(stcoord[,2], na.rm=TRUE)
ymax <- max(stcoord[,2], na.rm=TRUE)
y <- seq(ymin, ymax, by=n)
idwgrid<-expand.grid(x,y)
names(idwgrid)<-c("x","y")

# Interpolate "annual total precipitation" value
idwdata<-data.frame(stcoord[,1], stcoord[,2], weatherst$PYR)
names(idwdata)<-c("x","y","precipitation")
coordinates(idwdata) <- c("x", "y")
coordinates(idwgrid) <- c("x", "y")

# interpolate (IDW, power=2) the total precipitation value
idw<-idw(precipitation ~ 1, idwdata, idwgrid, idp=2)
idwresult <- data.frame(x=idwgrid$x, y=idwgrid$y, z=idw$var1.pred)
idwraster <- rasterFromXYZ(as.data.frame(idwresult)[, c("x", "y", "z")])

# see the result - 25 rainbow color, half-transparent
plot(idwraster, col=rainbow(25, alpha=0.5), axes=T)
plot(states48, add=T)
plot(weatherst, pch=20, cex=0.5, add=T)

##### SPATIAL REGRESSION MODELS #####

# read data
stateplane<-CRS("+proj=tmerc +lat_0=36.66666666666666 +lon_0=-88.33333333333333 
                +k=0.999975 +x_0=300000 +y_0=0 +ellps=GRS80 +datum=NAD83 
                +to_meter=0.3048006096012192")
tracts<-readShapeSpatial("Spatial_tracts.shp", ID="OBJECTID", proj4string=stateplane)

# plot our dependent variable - average price. 
# Note - many tracts have average price=0 (i.e. no data problem)
tracts2<-subset(tracts,tracts@data$AVPRCSQFT>0)
averageprice <- tracts2@data$AVPRCSQFT
nclr <- 8
colors <- brewer.pal(nclr, "RdYlGn")
class <- classIntervals(averageprice, nclr, style = "jenks")
colorcode <- findColours(class, colors)
par(mar = c(1, 1, 2, 2.5))
plot(tracts)
plot(tracts2, col=colorcode, add=T)
title(main="Average Property Price (per sq. foot)")
legend("topright", legend = names(attr(colorcode, "table")), fill = attr(colorcode, 
    "palette"), cex = 0.6, bty = "n")

# regular data checking.
# Dependent var <- AVPRCSQFT "average house price per square feet"
par(mfrow = c(1, 1))
histg <- hist(tracts$AVPRCSQFT)
xfit <- seq(min(tracts$AVPRCSQFT), max(tracts$AVPRCSQFT), length = 25)
yfit <- dnorm(xfit, mean = mean(tracts$AVPRCSQFT), sd = sd(tracts$AVPRCSQFT))
yfit <- yfit * diff(histg$mids[1:2]) * length(tracts$AVPRCSQFT)
lines(xfit, yfit, col = "blue", lwd = 2) # not a normal distribution but we'll ignore it for now

# Lets find neighbors - create a spatial weights matrix using 1st neighbor queen contiquity
tractsxy <- coordinates(tracts)
tracts.queen <- poly2nb(tracts, queen = T)
tracts.queenw <- nb2listw(tracts.queen, style = "W")

# see how it looks
par(mfrow = c(1, 2))
plot(tractsxy)
title(main = "centroid of tracts")
plot(tracts, border = "grey")
plot(tracts.queen, tractsxy, add = T)
title(main = "1st order queen contiquity")

# Check autocorrelation - and plot the moran's I using the spatial weights. 
par(mfrow = c(1, 1))
moran.test(tracts$AVPRCSQFT, listw = tracts.queenw)
# Outliers with heavier influence will be labeled with tracts name:
moran.plot(tracts$AVPRCSQFT, listw = tracts.queenw, labels = as.character(tracts$TRACT))

# though data needs transformation, we'll ignore it for now.
# Try linear regression
reg <-lm(AVPRCSQFT ~ AVG_HOUSES + PUBLICTRAN + MEDIAN_HHI
         + PCT_COLLGR + RATIO_TOTC, data = tracts@data)
summary(reg)

# good to check autocorrelation in residuals - show some sp. autocorrelation
moran.test(reg$residuals, tracts.queenw, alternative = "two.sided")

# try spatial lag model - remember data needs some work
spatiallag <- lagsarlm(AVPRCSQFT ~ AVG_HOUSES + PUBLICTRAN + MEDIAN_HHI
                       + PCT_COLLGR + RATIO_TOTC, data = tracts@data, 
                       tracts.queenw, method = "eigen", quiet = FALSE, tol.solve=1e-14)
summary(spatiallag)
# AIC actually decreases in comparison to linear regression (above)
# but, spatial autocorrelation in residuals deminishes
moran.test(spatiallag$residuals, tracts.queenw, alternative = "two.sided")

# try spatial error model - remember data needs some work
spatialerror <- errorsarlm(AVPRCSQFT ~ AVG_HOUSES + PUBLICTRAN + MEDIAN_HHI
                       + PCT_COLLGR + RATIO_TOTC, data = tracts@data, 
                       tracts.queenw, method = "eigen", quiet = FALSE, tol.solve=1e-14)
summary(spatialerror)
# AIC actually decreases in comparison to linear regression (above)
# but, spatial autocorrelation in residuals deminishes
moran.test(spatialerror$residuals, tracts.queenw, alternative = "two.sided")
