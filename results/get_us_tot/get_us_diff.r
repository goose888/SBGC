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
maskvar <- 'REGION_MASK_CRU_NCEP'
oxidvar <- 'ch4_oxid'
ch4var <- 'ch4_flux'
ch4var_2 <- 'ch4_ebul'
ch4var_3 <- 'ch4_aere'
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
#print('tag1')

# get sink
nc <- open.nc(fch4nc, write=FALSE)
# get oxid
oxid <- var.get.nc(nc, oxidvar)
ch4 <- var.get.nc(nc, ch4var)
close.nc(nc)
#print('tag2')

# get fractional water
nc <- open.nc(ffwnc, write=FALSE)
fw <- var.get.nc(nc, fwvar)
close.nc(nc)
#print('tag3')

#ch4[ch4>0]=0
fw[fw<=0.01] <- NA
fw[fw>0.01] <- 1.0

# Mask out values outside the US
oxid[maskval < 12] <- NA
oxid[maskval > 12] <- NA
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
oxid_m <- oxid*grid_area
ch4_m <- ch4*grid_area*fw

# Summarize the total number 
oxid_tot <- sum(colSums(oxid_m, na.rm = TRUE), na.rm = TRUE)/1e12/2.2
ch4_tot <- sum(colSums(ch4_m, na.rm = TRUE), na.rm = TRUE)/1e12/2.2

#print(oxid_tot)
print(ch4_tot)
