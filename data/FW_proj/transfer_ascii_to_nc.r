# ------------------------------------------------------------
# R script for transfering ascii file into netcdf file for 
# visualization and other calculation purpose 
# Author: Shijie Shu
# Date: 06/15/2017
# ------------------------------------------------------------

# Load Libraries
library(RNetCDF)

# User options
fdat <- './IP_DATA/HYDE32/urbanarea_hyde32_isam_cesm.dat'
fldmaskdat <- 'mask_isam_cesm.dat'
iyr <- 246
yr <- 1765:2010
nlat <- 360
nlon <- 720
tpoint <- 90672

# ----- Options ------
d3 <- FALSE
use_mask <- FALSE
save_netcdf <- TRUE
to_fraction <- TRUE
# --------------------

varname <- 'urban'
var_description <- 'Crop area (in fraction) in 1850 based on SSP2 model.'
fname_out <- 'urban_1765_hyde32'

## DO NOT EDIT BELOW
# Open and read in the dat file
var <- as.matrix(read.table(fdat, header=FALSE, sep=""))
if(d3) {

   # 3D file
   var3d <- array(-9999., dim=c(iyr, nlat, nlon))

   for (k in 1:iyr) {
      if(use_mask) {
         igpmap <- as.matrix(read.table(fldmaskdat, header=FALSE, sep=""))
         for (i in 1:tpoint) {
             var3d[k,igpmap[i,3],igpmap[i,2]] <- var[(k-1)*tpoint+i]
         }
      } else {
         for (i in 1:nlat) {
            for (j in 1:nlon) {
                var3d[k,i,j] <- var[(k-1)*nlon*nlat+(i-1)*nlon+j]
            }
         }
      }
      
      if(to_fraction) {
          # Calculate grid area
          grid_area <- matrix(-9999., nrow=nlat, ncol=nlon)
          EARTH_AREA =  5.096e14;
          lat <- seq(-89.75, 89.75, 0.5)
          res = 0.5;
      
          for (i in 1:nlat) {
            for (j in 1:nlon) {
               grid_area[i,j] <- (EARTH_AREA/2)*abs(sin((lat[i] - res/2)*pi/180) -
                        sin((lat[i] + res/2)*pi/180))/(360/res)
            }
          }
          var3d[k,,] <- var3d[k,,]/grid_area
          var3d[var3d>1.] <- 1.
      }
   }
   var3d[var3d<0] <- 0.
   
   if(save_netcdf) {
      fname <- paste(fname_out, '.nc', sep='')
   
      longvector <- seq(-179.75, 179.75, 0.5)
      latvector <- seq(-89.75, 89.75, 0.5)
   
      nc <- create.nc(fname)
   
      # Define the dimensions
      dim.def.nc(nc, "lon", length(longvector))
      dim.def.nc(nc, "lat", length(latvector))
      dim.def.nc(nc, "year", iyr)
   
      # Define the coordinate variables
      var.def.nc(nc, "lon", "NC_FLOAT", "lon")
      var.def.nc(nc, "lat",  "NC_FLOAT", "lat")
      var.def.nc(nc, "year",  "NC_INT", "year")
   
      var.def.nc(nc, varname, "NC_FLOAT", c("lon", "lat", "year"))
      # Put attributes
      att.put.nc(nc, varname, "missing_value", "NC_FLOAT", -9999.)
      att.put.nc(nc, varname, "long_name", "NC_CHAR", var_description)
   
      # Put the data
      var.put.nc(nc, "lon", longvector, na.mode=0)
      var.put.nc(nc, "lat", latvector, na.mode=0)
      var.put.nc(nc, "year", yr, na.mode=0)

      # re-order the dimension first
      new_var3d <- aperm(var3d, c(3,2,1))
      var.put.nc(nc, varname, new_var3d, na.mode=0)
   
      close.nc(nc)
   
   }

} else {

   # 2D file
   var2d <- matrix(-9999., nrow=nlat, ncol=nlon)
   if(use_mask) {
       igpmap <- as.matrix(read.table(fldmaskdat, header=FALSE, sep=""))
       for (i in 1:tpoint) {
           var2d[igpmap[i,3],igpmap[i,2]] <- var[i]
       }
   } else {
      for (i in 1:nlat) {
         for (j in 1:nlon) {
             var2d[i,j] <- var[(i-1)*nlon+j]
         }
      }
   }
   
   if(to_fraction) {
       # Calculate grid area
       grid_area <- matrix(-9999., nrow=nlat, ncol=nlon)
       EARTH_AREA =  5.096e14;
       lat <- seq(-89.75, 89.75, 0.5)
       res = 0.5;
   
       for (i in 1:nlat) {
         for (j in 1:nlon) {
            grid_area[i,j] <- (EARTH_AREA/2)*abs(sin((lat[i] - res/2)*pi/180) -
                     sin((lat[i] + res/2)*pi/180))/(360/res)
         }
       }
       var2d <- var2d/grid_area
       var2d[var2d>1.] <- 1.
   }
   
   if(save_netcdf) {
      fname <- paste(fname_out, '.nc', sep='')
   
      longvector <- seq(-179.75, 179.75, 0.5)
      latvector <- seq(-89.75, 89.75, 0.5)
   
      nc <- create.nc(fname)
   
      # Define the dimensions
      dim.def.nc(nc, "lon", length(longvector))
      dim.def.nc(nc, "lat", length(latvector))
   
      # Define the coordinate variables
      var.def.nc(nc, "lon", "NC_FLOAT", "lon")
      var.def.nc(nc, "lat",  "NC_FLOAT", "lat")
   
      var.def.nc(nc, varname, "NC_FLOAT", c("lon", "lat"))
      # Put attributes
      att.put.nc(nc, varname, "missing_value", "NC_FLOAT", -9999.)
      att.put.nc(nc, varname, "long_name", "NC_CHAR", var_description)
   
      # Put the data
      var.put.nc(nc, "lon", longvector, na.mode=0)
      var.put.nc(nc, "lat", latvector, na.mode=0)
      var.put.nc(nc, varname, t(var2d), na.mode=0)
   
      close.nc(nc)
   
   }

}
