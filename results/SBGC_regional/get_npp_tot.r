# ------------------------------------------------------------
# R script for calculating US total CH4 oxidation and emission
# Author: Shijie Shu
# Date: 06/25/2017
# ------------------------------------------------------------

# Load Libraries
library(RNetCDF)

# User options
args = commandArgs(trailingOnly = TRUE)
fnppnc <- as.character(args[1])
fmasknc <- 'surfdata_05x05.nc'
maskvar <- 'PFMASK'
nppvar <- 'anpp_avg_yr'
respvar <- 'resp'
nlon <- 720
nlat <- 360
res <- 0.5

## DO NOT EDIT BELOW
# Open and read netcdf file
nc <- open.nc(fmasknc, write=FALSE)
# get the size of lat and lon
londim <- dim.inq.nc(nc, "lon")
latdim <- dim.inq.nc(nc, "lat")
# get permafrost mask
maskval <- var.get.nc(nc, maskvar)
close.nc(nc)

nc <- open.nc(fnppnc, write=FALSE)
# get SOC
npp <- var.get.nc(nc, nppvar)
close.nc(nc)

# Mask out values outside the permafrost region
npp[maskval < 1]  <- NA

# Get the grid area
grid_area <- matrix(-9999., nrow=nlon, ncol=nlat)
EARTH_AREA =  5.096e14;
lat <- seq(-89.75, 89.75, 0.5)
res = 0.5;

for (i in 1:nlat) {
  for (j in 1:nlon) {
     grid_area[j,i] <- (EARTH_AREA/2)*abs(sin((lat[i] - res/2)*pi/180) -
              sin((lat[i] + res/2)*pi/180))/(360/res)
  }
}

##### CH4 for the separated categories only #####
npp_a <- npp*grid_area

# Summarize the total number, TgCH4 m-2 yr-1
npp_tot <- sum(colSums(npp_a, na.rm = TRUE), na.rm = TRUE)/1e15

print(npp_tot)
