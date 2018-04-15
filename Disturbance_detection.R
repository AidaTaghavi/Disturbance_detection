# load packages
require(RStoolbox)
require(sp)				
require(rgdal)			
require(raster)			
require(bfast)
require(bfastSpatial)	
require(maptools)	
require(devtools)	
require(rasterVis)
require(igraph)
require(rts)
#download MODIS product
x=1 #x="MOD13Q1"
#Downloads selected tiles, and mosaic, reproject them in UTM_WGS84, zone 30 projection and convert all bands into Geotif format
ModisDownload(x=x,h=c(21,22),v=5,dates=c('2000.01.01','2014.01.01'),MRTpath='.', mosaic=T,proj=T,proj_type="UTM",utm_zone=39,datum="WGS84",pixel_size=250)
# define here your working directory
path <- "F:/R_project/"		# working directory
setwd(path)
# read data from GeoTiff
list <- list.files(path=".",pattern='.tif$',full.names=TRUE)
# stack files together and read                 
r <- stack(list)
#---------------------------------
# read some shapefiles for mapping
#---------------------------------
setwd(path.shp)
# read vector data 
# readOGR is a generic function from package rgdal to read vector data
region9 <- readOGR(".", "khilikhoshk_21")
region6 <- readOGR(".", "khilimartobt_2")

#---------------------------------
# subset data for area of interest
#---------------------------------
verydry_ndvi <- crop(r,region9)
vrdry_ndvi <- mask(verydry_ndvi, region9)

# transform data to lat/lon projection
#-------------------------------------

# check the projection 
projection(vrdry_ndvi)		# data is in Lambert Conformal Conic projection 
# transform data from Lambert projection to Latitude/Longitude standard projection
change_pro_ras<-projectRaster(vrdry_ndvi,crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
change_pro<-spTransform(region9,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# write raster
writeRaster(change_pro_ras,file.path(".",names(vrdry_ndvi)), bylayer=TRUE,format='GTiff')
#---------------------------------
# Display the studyarea
#---------------------------------
plot(r[[1]],main="Iran",legend=FALSE)
plot(change_pro6,pch=18, col="red", add=TRUE)
plot(change_pro9,pch=18, col="yellow", add=TRUE)
legend("bottomright", title="Regions", c("very Arid","very Humid"), fill=c("red","yellow"), horiz=FALSE, cex=0.7)

#-------------
# Temporal analyze NDVI data
#-------------
setwd(path)
#####Read NDVI data set 2008 and explore the structure
vrdry_ndvi.2008 <- subset(r, grep(paste("MOD13Q1.A2008", sep =""),names(r)))
vrdry_ndvi.2008.median <- calc(vrdry_ndvi.2008, fun = function(x) median(x, na.rm =TRUE))
vrdry_ndvi.2008.var <- calc(vrdry_ndvi.2008, fun = function(x) var(x, na.rm =TRUE))
vrdry_ndvi.2008.mean <- calc(vrdry_ndvi.2008, fun = function(x) mean(x, na.rm =TRUE))
vrdry_ndvi.2008.sd <- calc(vrdry_ndvi.2008, fun = function(x) sd(x, na.rm =TRUE))
#--------------------
# plot graphics
#--------------------
cols <- brewer.pal(5, "RdYlBu")
plot(vrdry_ndvi.2008.median/10000,main = "Median NDVI values 2008 in veryArid regions-Iran",col=cols)
plot(vrdry_ndvi.2008.mean/10000,main = "Mean NDVI values 2008 in veryArid regions-Iran",col=cols)
plot(vrdry_ndvi.2008.var/10000,main = "variance NDVI values 2008 in veryArid regions-Iran",col=cols)
plot(vrdry_ndvi.2008.sd/10000,main = "standard deviation NDVI values 2008 in veryArid regions-Iran",col=cols)
#######Read NDVI data set 2003 and explore the structure
vrdry_ndvi.2003 <- subset(r, grep(paste("MOD13Q1.A2003", sep =""),names(r)))
vrdry_ndvi.2003
vrdry_ndvi.2003.median <- calc(vrdry_ndvi.2003, fun = function(x) median(x, na.rm =TRUE))
vrdry_ndvi.2003.var <- calc(vrdry_ndvi.2003, fun = function(x) var(x, na.rm =TRUE))
vrdry_ndvi.2003.mean <- calc(vrdry_ndvi.2003, fun = function(x) mean(x, na.rm =TRUE))
vrdry_ndvi.2003.sd <- calc(vrdry_ndvi.2003, fun = function(x) sd(x, na.rm =TRUE))
#--------------------
# plot graphics
#--------------------
cols <- brewer.pal(5, "YlOrRd")
plot(vrdry_ndvi.2003.median/10000,main = "Median NDVI values 2003 in veryArid regions-Iran",col=cols)
plot(vrdry_ndvi.2003.mean/10000,main = "Mean NDVI values 2003 in veryArid regions-Iran",col=cols)
plot(vrdry_ndvi.2003.var/10000,main = "variance NDVI values 2003 in veryArid regions-Iran",col=cols)
plot(vrdry_ndvi.2003.sd/10000,main = "standard deviation NDVI values 2003 in veryArid regions-Iran",col=cols)
#####analyze monthly Modis NDVI
doy2date <- function(year, doy) {
  as.Date(doy - 1, origin = paste0(year, "-01-01"))
}
data.dates <- as.character(doy2date(2003, doy = as.numeric(substr(names(vrdry_ndvi.2008), 14, 16))))
data.dates
modis.stack <- stack()
for (i in 1:12) {
bands_i <- which(sapply(strsplit(data.dates, "-"),
function(x) as.numeric(x[2])) == i)
modis.stack <- stack(modis.stack, calc(subset(vrdry_ndvi.2003,
bands_i),
fun = mean))
}
names(modis.stack) <- month.abb
levelplot(modis.stack, par.settings = PuOrTheme)
bwplot(modis.stack)
#----------------------------
# BFASTSpatial for very Arid region
#----------------------------
#####create raster stacks of MODIS layers
list <- list.files(path=".",pattern='.tif$',full.names=TRUE)
verydry_ndvi<-timeStackMODIS(list)
######2003
#####Spatial BFASTMonitor_Breaks For Additive Season and Trend
bfm2003<- bfmSpatial(verydry_ndvi, start=c(2003, 1),returnLayers=c("breakpoint", "magnitude", "error", "history", "r.squared", "adj.r.squared", "coefficients"))
# extract change raster
change_2003 <- raster(bfm2003, 1)
months <- changeMonth(change_2003)
# set up labels and colourmap for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, pch=16,cex=0.7, fill=cols, ncol=3 )
# extract magnitude raster
magn_2003 <- raster(bfm2003, 2)/10000
# make a version showing only breakpoing pixels
magn_bkp2003 <- magn_2003
magn_bkp2003[is.na(change_2003)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp2003, main="Magnitude: breakpoints")
plot(magn_2003, main="Magnitude: all pixels")
par(op)
magn_2003 <- raster(bfm2003, 2) / 10000
# extract and rescale magnitude and apply a -500 threshold
magn03thresh <- magn_2003
magn03thresh[magn_2003 > -0.05] <- NA
magn03_areasieve <- areaSieve(magn03thresh)
changeSize <- clumpSize(magn03_areasieve, f=62500/10000)
plot(changeSize, col=bpy.colors(50), main="Size of change(hectares)")
changeSize <- clumpSize(magn03_areasieve, f=62500/10000, stats=TRUE)
print(changeSize$stats)
#######2008
#####Spatial BFASTMonitor_Breaks For Additive Season and Trend
bfm2008_10<- bfmSpatial(verydry_ndvi, start=c(2008, 1), monend=c(2010, 1),order=1,returnLayers=c("breakpoint", "magnitude", "error", "history", "r.squared", "adj.r.squared", "coefficients"))
change_2008 <- raster(bfm2008_10, 1)
cols <- brewer.pal(5, "Spectral")
plot(change_2008, col=cols, breaks=c(1:12), legend=FALSE)
months <- changeMonth(change_2008)
# set up labels and colourmap for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, pch=16,cex=0.7, fill=cols, ncol=3 )
# extract magnitude raster
magn_2008 <- raster(bfm2008_10, 2)/10000
# make a version showing only breakpoing pixels
magn_bkp2008 <- magn_2008
magn_bkp2008[is.na(change_2008)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp2008, main="Magnitude: breakpoints",col=cols)
plot(magn_2008, main="Magnitude: all pixels",col=cols)
par(op)
magn_2008 <- raster(bfm2008_10, 2) / 10000
# extract and rescale magnitude and apply a -500 threshold
magn08thresh <- magn_2008
magn08thresh[magn_2008 > -0.05] <- NA
magn08_areasieve <- areaSieve(magn08thresh)
changeSize_2008 <- clumpSize(magn08_areasieve, f=62500/10000)
plot(changeSize_2008, col=bpy.colors(50), main="Size of change(hectares)")
changeSize_2008 <- clumpSize(magn08_areasieve, f=62500/10000, stats=TRUE)
print(changeSize_2008$stats) 
#################Working with pixels: bfmPixel
#####In 2008
#####create raster stacks of MODIS layers
list <- list.files(path=".",pattern='.tif$',full.names=TRUE)
verydry_ndvi<-timeStackMODIS(list)
# plot the 22nd layer
plot(verydry_ndvi, 22)
# run bfmPixel() in interactive mode with a monitoring period
# choose the pixel whose time series you want to see by clicking on the map you plotted a moment ago
# starting from the 1st day in 2008
bfm08 <- bfmPixel(verydry_ndvi, start=c(2008, 1), interactive=TRUE)
targcell <- 76315
bfm08 <- bfmPixel(verydry_ndvi, cell=76315, start=c(2008, 1))
#to get the model and plot it
bfm08$bfm
plot(bfm08$bfm)
####
# starting from the 1st day in 2003
bfm03 <- bfmPixel(verydry_ndvi, start=c(2003, 1), interactive=TRUE)
targcell <- 76315
bfm03 <- bfmPixel(verydry_ndvi, cell=targcell, start=c(2003, 1))
#to get the model and plot it
bfm03$bfm
plot(bfm03$bfm)
##################
#########Region shapefile:very humid
region6 <- readOGR("E:/NDVI_EVI_MOD - Copy/shp", "khilimartobt_2")
vryhumid_ndvi <- crop(r, extent(region6))
vrhumid_ndvi <- mask(vryhumid_ndvi, region6)/10000
change_pro<-spTransform(region6,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
change_pro_ras<-projectRaster(vrhumid_ndvi,crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
writeRaster(change_pro_ras, file.path('D:/NDVI_EVI_MOD - Copy/ndvi-test.MOD13Q1', names(r)), bylayer=TRUE,format='GTiff')
setwd("E:/new_humid")
list <- list.files(path=".",pattern='.tif$',full.names=TRUE)
r_humid <- stack(list)
vrhu_ndvi.2008 <- subset(r_humid, grep(paste("MOD13Q1.A2008", sep =""),names(r_humid)))
vrhu_ndvi.2008
vrhu_ndvi.2008.median <- calc(vrhu_ndvi.2008, fun = function(x) median(x, na.rm =TRUE))
vrhu_ndvi.2008.var <- calc(vrhu_ndvi.2008, fun = function(x) var(x, na.rm =TRUE))
vrhu_ndvi.2008.mean <- calc(vrhu_ndvi.2008, fun = function(x) mean(x, na.rm =TRUE))
cols <- brewer.pal(9, "Greens")
cols <- brewer.pal(9, "YlGn")
plot(vrhu_ndvi.2008.median,main = "Median NDVI values 2008 in Humid regions-Iran",col=cols)
plot(vrhu_ndvi.2008.mean,main = "Mean NDVI values 2008 in Humid regions-Iran",col=cols)
doy2date <- function(year, doy) {
  as.Date(doy - 1, origin = paste0(year, "-01-01"))
}
data.dates <- as.character(doy2date(2008, doy = as.numeric(substr(names(vrhu_ndvi.2008), 14, 16))))
data.dates
modis.stack <- stack()
for (i in 1:12) {
  bands_i <- which(sapply(strsplit(data.dates, "-"),
                          function(x) as.numeric(x[2])) == i)
  modis.stack <- stack(modis.stack, calc(subset(vrhu_ndvi.2008,
  bands_i),
 fun = mean))
}
names(modis.stack) <- month.abb
levelplot(modis.stack, par.settings = RdBuTheme)
bwplot(modis.stack)
veryhumid_region<-timeStackMODIS(list)
bfm8<- bfmSpatial(veryhumid_region, start=c(2008, 1), monend=c(2009, 1),order=1)
bfm8 <- bfmSpatial(veryhumid_region, start=c(2008, 1))
bfm8_layers <- bfmSpatial(veryhumid_region, start=c(2008, 1),returnLayers=c("breakpoint", "magnitude", "error", "history", "r.squared", "adj.r.squared", "coefficients"))
bfm3_layers <- bfmSpatial(veryhumid_region, start=c(2003, 1),returnLayers=c("breakpoint", "magnitude", "error", "history", "r.squared", "adj.r.squared", "coefficients"))
change <- raster(bfm8, 1)
months <- changeMonth(change)
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun",
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
legend("bottomright", legend=monthlabs, pch=16,cex=0.7, fill=cols, ncol=1 )
####legend(0.03, 0.015,legend=monthlabs,cex=0.8, fill=cols, ncol=1,inset =0.01)
magn <- raster(bfm8, 2)
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")
magn08thresh <- magn
magn08thresh[magn > -0.05] <- NA
magn08_areasieve <- areaSieve(magn08thresh)
changeSize <- clumpSize(magn08_areasieve, f=62500/10000)
changeSize <- clumpSize(magn08thresh, f=62500/10000)
plot(changeSize, col=bpy.colors(50), main="size of changes in Humid regions during 2008 to 2009(hectares)")
changeSize <- clumpSize(magn08thresh, f=62500/10000,stats=TRUE)
print(changeSize$stats)
#################Working with pixels: bfmPixel
#####In 2008
#####create raster stacks of MODIS layers
list <- list.files(path=".",pattern='.tif$',full.names=TRUE)
humid_ndvi<-timeStackMODIS(list)
# plot the 22nd layer
plot(humid_ndvi, 22)
# run bfmPixel() in interactive mode with a monitoring period
# choose the pixel whose time series you want to see by clicking on the map you plotted a moment ago
# starting from the 1st day in 2008
bfm08 <- bfmPixel(humid_ndvi, start=c(2008, 1), interactive=TRUE)
targcell <- 56315
bfm08 <- bfmPixel(humid_ndvi, cell=56315, start=c(2008, 1))
#to get the model and plot it
bfm08$bfm
plot(bfm08$bfm)
####
# starting from the 1st day in 2003
bfm03 <- bfmPixel(humid_ndvi, start=c(2003, 1), interactive=TRUE)
targcell <- 34315
bfm03 <- bfmPixel(humid_ndvi, cell=targcell, start=c(2003, 1))
#to get the model and plot it
bfm03$bfm
plot(bfm03$bfm)
