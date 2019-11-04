# ------------------------------------------------------------
# R script for calculating US total CH4 oxidation and emission
# Author: Shijie Shu
# Date: 06/25/2017
# ------------------------------------------------------------

# Load Libraries
library(RNetCDF)

# User options
args = commandArgs(trailingOnly = TRUE)
fch4nc <- as.character(args[1])
ffwnc <- as.character(args[2])
fmasknc <- 'surfdata_05x05_13reg.nc'
fwetnc <- 'surfdata_05x05.nc'
maskvar <- 'REGION_MASK_CRU_NCEP'
ch4var <- 'ch4_aere'
ch4var_dry <- 'ch4_aere_dry'
ch4var_wet <- 'ch4_aere_wet'
ch4var_inund <- 'ch4_aere_inund'
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
close.nc(nc)

# Open and get the wetland fraction
nc <- open.nc(fwetnc, write=FALSE)
frac_wet <- var.get.nc(nc, 'FW')

nc <- open.nc(fch4nc, write=FALSE)
# get flux
ch4_avg <- var.get.nc(nc, ch4var)
ch4_wet <- var.get.nc(nc, ch4var_wet)
ch4_inund <- var.get.nc(nc, ch4var_inund)
ch4_dry <- var.get.nc(nc, ch4var_dry)
close.nc(nc)
#oxid <- var.get.nc(nc, oxidvar)
#resp <- var.get.nc(nc, respvar)

# get fractional water
nc <- open.nc(ffwnc, write=FALSE)
fw <- var.get.nc(nc, fwvar)
close.nc(nc)
#print('tag3')

frac_inund <- fw - frac_wet

# Now the output already considered wetland, inundated land and 
# dryland separately

## CH4 for all non-wetland #####
#ch4 <- ch4_dry
#ch4[ch4>0] = NA     # The condition to only get the results of dry soil
## CH4 for wetland #####
ch4 <- ch4_wet + ch4_inund
## CH4 for all #####
#ch4 <- ch4_avg

# Mask out values outside the US
ch4[maskval < 12]  <- NA
ch4[maskval > 12]  <- NA
fw[maskval < 12]  <- NA
fw[maskval > 12]  <- NA

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
ch4_m <- ch4*grid_area

# Summarize the total number, TgCH4 m-2 yr-1
ch4_tot <- sum(colSums(ch4_m, na.rm = TRUE), na.rm = TRUE)/1e12

print(ch4_tot)
