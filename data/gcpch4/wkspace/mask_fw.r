# ------------------------------------------------------------
# R script for calculating US total CH4 oxidation and emission
# Author: Shijie Shu
# Date: 06/25/2017
# ------------------------------------------------------------

# Load Libraries
library(RNetCDF)

# User options
args = commandArgs(trailingOnly = TRUE)
fnc <- as.character(args[1])
fmasknc <- 'fw_frac_mask.nc'
maskvar <- 'FW'
fwvar <- 'fw'
nlon <- 720
nlat <- 360
res <- 0.5

## DO NOT EDIT BELOW
# Open and read netcdf file
nc <- open.nc(fmasknc, write=FALSE)
# get the size of lat and lon
londim <- dim.inq.nc(nc, "Lon")
latdim <- dim.inq.nc(nc, "Lat")
# get mask
maskval <- var.get.nc(nc, maskvar)
maskval[maskval>=0.] <- 1
maskval[maskval<0.] <- 0
close.nc(nc)

# # Reverse the latitude
# maskval = maskval[,360:1]

# get FW
nc <- open.nc(fnc, write=TRUE)
fw <- var.get.nc(nc, fwvar)
lat_new <- var.get.nc(nc, 'lat')
lon_new <- var.get.nc(nc, 'lon')
print(dim(fw))

# Adjust lon
tmp <- fw[1:360,]
fw[1:360,] <- fw[361:720,]
fw[361:720,] <- tmp
print("TTT")

# Reverse the latitude and change the longitude
fw <- fw[,360:1]
lat_new <- lat_new[360:1]
lon_new <- lon_new + 180.

# Mask out values outside the US
maskval[fw>0.] <- 2
fw[maskval<=0] <- -999.
fw[maskval==1] <- 0.0
print(dim(fw))

#new_var2d <- aperm(fw, c(2,1))
#print(dim(new_var2d))

var.put.nc(nc, 'lon', lon_new, na.mode=0)
var.put.nc(nc, 'lat', lat_new, na.mode=0)
var.put.nc(nc, fwvar, fw, na.mode=0)
close.nc(nc)

