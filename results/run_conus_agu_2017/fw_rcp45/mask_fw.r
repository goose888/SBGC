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
# fnc <- 'FW_2013_3.nc'
fmasknc <- 'fw_frac_max.nc'
maskvar <- 'FW'
fwvar <- 'FW'
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
# maskval[maskval>=0] <- 1
close.nc(nc)

# get FW
nc <- open.nc(fnc, write=TRUE)
fw <- var.get.nc(nc, fwvar)

# Mask out values outside the US
fw[fw>0.8] <- 0.8
fw[is.na(fw)] <- 0.
fw[is.na(maskval)] <- -9999.

var.put.nc(nc, fwvar, fw, na.mode=0)
close.nc(nc)

