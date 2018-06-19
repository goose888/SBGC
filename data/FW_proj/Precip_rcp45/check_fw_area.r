# ------------------------------------------------------------
# R script for convert NC into ASCII
# Author: Shijie Shu
# Date: 06/25/2017
# ------------------------------------------------------------

# Load Libraries
library(RNetCDF)
#library(caTools)

# User defined variable
fheader <- 'FW_ML_DOMAIN.HDR'

prefix <- 'FW_frac_'
suffix <- '.nc'
year <- 2006:2015
prd <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', 
        '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', 
        '30', '31', '32', '33', '34', '35', '36' )
varname <- 'FW'

fmask <- 'surfdata_05x05_13reg.nc'
varmask <- 'REGION_MASK_CRU_NCEP'

lenyear <- length(year)
lenprd <- length(prd)

nlat <- 360
nlon <- 720

# Open and read multi-year minimum as mask
d <- open.nc(fmask)

#lon <- var.get.nc(d,"Lon")
#lat <- var.get.nc(d,"Lat")
mask <- var.get.nc(d,varmask)

close.nc(d)

mask[mask<=11] <- NA
mask[mask>=13] <- NA
#mask = mask[,ncol(mask):1]

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


# Open each file and calculate the total area
for (i in 1:lenyear) {
   datapt <- 0
   for (k in 1:lenprd) {
      # Counter
      datapt <- datapt + 1

      # Open and read data
      fname <- paste(prefix, year[i], '_', prd[k], suffix, sep='')

      d <- open.nc(fname)

      #lon <- var.get.nc(d,"Lon")
      lat <- var.get.nc(d,"Lat")
      fw <- var.get.nc(d,varname)

      fw[fw<0] <- NA

      close.nc(d)

      fw[is.na(mask)] <- 0.

      fwa <- fw*grid_area
      fwa[fwa<0] <- NA

      # Calculate the total area
      area_fw <- sum(sum(fwa, na.rm=TRUE), na.rm=TRUE)

      print(area_fw)
   }
}
