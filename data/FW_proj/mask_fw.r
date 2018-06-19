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
londim <- dim.inq.nc(nc, "lon")
latdim <- dim.inq.nc(nc, "lat")
# get mask
maskval <- var.get.nc(nc, maskvar)
maskval[maskval>=0] <- 1
maskval[maskval>0.8] <- 2
close.nc(nc)

# get FW
nc <- open.nc(fnc, write=TRUE)
fw <- var.get.nc(nc, fwvar)

# Mask out values outside the US
fw[maskval==2] <- 0.8
fw[maskval==1 && fw<0.0] <- 0.0

var.put.nc(nc, fwvar, fw, na.mode=0)
close.nc(nc)

