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
fmasknc <- 'surfdata_05x05_13reg.nc'
maskvar <- 'REGION_MASK_CRU_NCEP'
oxidvar <- 'wetland_CH4_emissions'
ch4var <- 'wetland_CH4_emissions'
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

nc <- open.nc(fch4nc, write=FALSE)
# get oxid
oxid_ori <- var.get.nc(nc, oxidvar)
ch4_ori <- var.get.nc(nc, ch4var)
close.nc(nc)

#oxid <- oxid_ori[,,1]
#ch4 <- ch4_ori[,,1]
oxid <- oxid_ori
ch4 <- ch4_ori

temp = maskval[1:360,]
maskval[1:360,] = maskval[361:720,]
maskval[361:720,] = temp[1:360,]

# Mask out values outside the US
oxid[maskval < 12] <- NA
oxid[maskval > 12] <- NA
#ch4[maskval < 13]  <- NA
#ch4[maskval > 13]  <- NA

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
oxid_m <- oxid*grid_area
ch4_m <- ch4*grid_area

# Summarize the total number 
oxid_tot <- sum(colSums(oxid_m, na.rm = TRUE), na.rm = TRUE)/1e12
ch4_tot <- sum(colSums(ch4_m, na.rm = TRUE), na.rm = TRUE)/1e12

#print(oxid_tot)
print(ch4_tot)
