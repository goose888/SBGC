# ===================================================
# Script for reproject permafrost extent data
# to WGS84 half by half degree geographical coord
# So called climate modeling grid (CMG)
# Author: Shijie Shu
# ===================================================
library(RNetCDF)
library(sp)
library(raster)
library(rgdal)

## Functions
# EPSG:3410, NSIDC EASE-Grid Global
convertGrid <- function(gridfile, header, name, inCRS="+init=epsg:3410"){

#  inCRS="+init=epsg:3408"

  d <- as.matrix(read.table(gridfile))
  pf <- matrix(0, nrow=1441, ncol=1441)
  for (i in 1:1441) {
    for (j in 1:1441) {
      pf[j, i] <- d[((i-1)*1441+j)]
    }
  }

 # finfo <- readLines(header)
  # Get the x-north, y 
#  x_westing <- as.numeric(substr(finfo[13], 47, 60))
#  y_northing <- as.numeric(substr(finfo[13], 63, 74))
#  x_res <- as.numeric(substr(finfo[13], 77, 93))
#  y_res <- as.numeric(substr(finfo[13], 96, 112))

  x_westing <- -9024309
  y_northing <- 9024309
  x_res <- 12533.7625
  y_res <- 12533.7625

#  sm$x <- seq((x_westing+x_res/2), (-x_westing-x_res/2), x_res)
#  sm$y <- seq((y_northing-y_res/2), (-y_northing+y_res/2), -y_res)

#  sm$x <- seq(1, length(sm$x), 1)
#  sm$y <- seq(1, length(sm$y), 1)

#  xp <- data.frame(x=sm$x,y=0)
#  coordinates(xp)=~x+y
#  proj4string(xp)=CRS(inCRS)
#  xp=spTransform(xp,CRS("+init=epsg:3410"))

#  yp <- data.frame(x=0,y=sm$y)
#  coordinates(yp)=~x+y
#  proj4string(yp)=CRS(inCRS)
#  yp=spTransform(yp,CRS("+init=epsg:3410"))

 # sm$z[is.na(sm$z)] <- -9999.

#  sm$xp <- coordinates(xp)[,1]
#  sm$yp <- coordinates(yp)[,2]
#  sm$xp <- round(sm$xp, digits=4)
#  sm$yp <- round(sm$yp, digits=4)

  b<-extent(c(x_westing, -x_westing, -y_northing, y_northing))
  smr <- raster(b, ncol=dim(pf)[1], nrow=dim(pf)[2], crs="+init=epsg:3408")
  smr[] <- t(pf)
  return(smr)
}

transformTo <- function(r1){
   ## 0.5*0.5 degree resolution and extent -180, 180, -90, 90
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90,
    nrows=180*2,ncols=360*2,crs="+init=epsg:4326")
  r <- projectRaster(r1,r,method='ngb')
  return(r)
}

## Main Program
# User defined variable
header <- 'FW_ML_DOMAIN.HDR'   # Share the same header
varname <- 'pfcode'
output_nc <- TRUE
output_gtiff <- FALSE

# Do not edit below unless you have confidence

# Read in latitude and longitude

fname <- 'pf.ascii'
smr <- convertGrid(fname,header,varname)
# plot(smr)

## anything below zero is NA (-1 is missing data, soil moisture is +ve)
smr[smr < -0.1] <- NA
smrp = transformTo(smr) # takes a short while

if(output_nc) {
   longvector <- seq(-179.75, 179.75, 0.5)
   latvector <- seq(89.75, -89.75, -0.5)
   avgdata <- t(as.matrix(smrp))

   nc <- create.nc(paste(fname, '_proj.nc',sep=''))

   # Define the dimensions
   dim.def.nc(nc, "Lon", length(longvector))
   dim.def.nc(nc, "Lat", length(latvector))

   # Define the coordinate variables
   var.def.nc(nc, "Lon", "NC_FLOAT", "Lon")
   var.def.nc(nc, "Lat",  "NC_FLOAT", "Lat")

   var.def.nc(nc, "PFCODE", "NC_FLOAT", c("Lon", "Lat"))
   # Put attributes
   att.put.nc(nc, "PFCODE", "missing_value", "NC_FLOAT", -9999.)
   att.put.nc(nc, "PFCODE", "long_name", "NC_CHAR", "Permafrost extent and ground ice conditions.")

   # Put the data
   var.put.nc(nc, "Lon", longvector, na.mode=0)
   var.put.nc(nc, "Lat", latvector, na.mode=0)
   var.put.nc(nc, "PFCODE", avgdata, na.mode=0)

   close.nc(nc)
}

if(output_gtiff) {
   rf <- writeRaster(smrp, filename=paste(fname,'_proj.tif',sep=''), format="GTiff", overwrite=TRUE)
}
