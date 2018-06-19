# ===================================================
# Script for reproject permafrost extent data
# to WGS84 half by half degree geographical coord
# So called climate modeling grid (CMG)
# Author: Shijie Shu
# ===================================================
library(RNetCDF)

## Main Program
# User defined variable
fsurfnc <- 'surfdata_05x05.nc'   # Share the same header
fnc <- 'permafrost_mask.nc'
varname <- 'PFCODE'

# Open and read in NC file
nc <- open.nc(fnc, write=FALSE)
# get the size of lat and lon
londim <- dim.inq.nc(nc, "Lon")
latdim <- dim.inq.nc(nc, "Lat")
var_raw <- var.get.nc(nc, varname)
close.nc(nc)

# Swap the longitude
temp = var_raw[1:360,]
var_raw[1:360,] = var_raw[361:720,]
var_raw[361:720,] = temp
var_raw <- var_raw[,ncol(var_raw):1]

# Open and read in NC file
nc <- open.nc(fsurfnc, write=TRUE)
# get the size of lat and lon
londim <- dim.inq.nc(nc, "lon")
latdim <- dim.inq.nc(nc, "lat")
mask <- var.get.nc(nc, 'MASK')

# Check and extract only continuous and discontinuous permafrost
for (i in 1:720) {
   for (j in 1:360) {
      if(mask[i,j]==1) {
         if(var_raw[i,j] %in% c(1, 2, 5, 6, 9, 10, 13, 14, 17, 18)) {
            var_raw[i,j] = 1
         } else {
            var_raw[i,j] = 0
         }
      } else {
         var_raw[i,j] = 0
      }
   }
}

# Add new variable into the NC file
var.def.nc(nc, 'PFMASK', 'NC_INT', c("lon", "lat"))
att.put.nc(nc, 'PFMASK', "missing_value", "NC_FLOAT", -9999.)
att.put.nc(nc, 'PFMASK', "long_name", "NC_CHAR", 'Mask for the permafrost')

# Put value in
var.put.nc(nc, 'PFMASK', var_raw)

# Close NC file
close.nc(nc)

