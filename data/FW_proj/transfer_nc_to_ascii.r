# ------------------------------------------------------------
# R script for adjusting the FW projection.
# Author: Shijie Shu
# Date: 06/25/2017
# ------------------------------------------------------------

# Load Libraries
library(RNetCDF)

# User options
fcropnc <- './IP_DATA/SSP2/GLANDCOVER_rcp37_crops.nc'
fpasturenc <- './IP_DATA/SSP2/GLANDCOVER_rcp37_pasture.nc'
furbannc <- './IP_DATA/SSP2/GLANDCOVER_rcp37_urban.nc'
fforestnc <- './IP_DATA/SSP2/GLANDCOVER_rcp37_forest.nc'
fothernc <- './IP_DATA/SSP2/GLANDCOVER_rcp37_other.nc'
joblist <- list(fcropnc, fpasturenc, furbannc, fforestnc, fothernc)
fldmaskdat <- 'mask_isam_cesm.dat'
iyr <- 251
yr <- 1850:2100
nlat <- 360
nlon <- 720
tpoint <- 90672

# <time, lat, lon>
varname_crops <- 'crops'
varname_pasture <- 'pasture'
varname_urban <- 'urban'
varname_forests <- 'forests'
varname_other <- 'other_natural'
varlist <- list(varname_crops, varname_pasture, varname_urban, varname_forests, varname_other)

#fcrops_o <-   './IP_DATA/SSP2/crop_ssp2_igp_isam_cesm.dat'
#fpasture_o <- './IP_DATA/SSP2/pasture_ssp2_igp_isam_cesm.dat'
#furban_o <-   './IP_DATA/SSP2/urban_ssp2_igp_isam_cesm.dat'
#fforests_o <- './IP_DATA/SSP2/forest_ssp2_igp_isam_cesm.dat'
#fother_o <-   './IP_DATA/SSP2/other_ssp2_igp_isam_cesm.dat'
#
fcrops_o <-   './IP_DATA/SSP2/croparea_1850_ssp2_isam_cesm.dat'
fpasture_o <- './IP_DATA/SSP2/pastarea_1850_ssp2_isam_cesm.dat'
furban_o <-   './IP_DATA/SSP2/urbanarea_ssp2_isam_cesm.dat'
fforests_o <- './IP_DATA/SSP2/forestarea_1850_ssp2_isam_cesm.dat'
fother_o <-   './IP_DATA/SSP2/otherarea_1850_ssp2_isam_cesm.dat'
outlist <- list(fcrops_o, fpasture_o, furban_o, fforests_o, fother_o)

# ----- Options ------
d3 <- FALSE
use_mask <- FALSE
save_ascii <- TRUE
as_fraction <- FALSE
# --------------------

## DO NOT EDIT BELOW

for (f in 1:5) {

   # Open and read netcdf file
   nc <- open.nc(as.character(joblist[f]), write=FALSE)
   # get the size of lat and lon
   londim <- dim.inq.nc(nc, "longitude")
   latdim <- dim.inq.nc(nc, "latitude")
   timedim <- dim.inq.nc(nc, "time")
 
   # start from 1982 to 2005
   var_raw <- var.get.nc(nc, as.character(varlist[f]))
 
   close.nc(nc)

   if(d3) {
      # 3-D
      var <- var_raw
      # 1. Consider writing area or directly the fraction?
      if(as_fraction) {
         # Do nothing
         return
      } else {
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
         for (i in 1:iyr) {
            var[,,i] <- var[,,i]*grid_area
         }
      }
      # 2. Transfer matrix into vector
      if(use_mask) {
         # Quick initialization
         var_dat <- 1:(tpoint*iyr)
         # Read in the mask
         igpmap <- as.matrix(read.table(fldmaskdat, header=FALSE, sep=""))
         for (j in 1:iyr) {
            for (i in 1:tpoint) {
               # For SSP2
                var_dat[((j-1)*tpoint+i)] <- var[igpmap[i,2], (nlat - igpmap[i,3] + 1), j]
               # For Hyde
               # var_dat[((j-1)*tpoint+i)] <- var[igpmap[i,2], igpmap[i,3], j]
            }
         }
      } else {
         # Quick initialization
         var_dat <- 1:(nlat*nlon*iyr)
         # See if this loop is too slow.
         for (k in 1:iyr) {
            for (i in 1:nlat) {
               for (j in 1:nlon) {
                 # For SSP2
                  var_dat[((k-1)*nlat*nlon+(i-1)*nlon+j)] <- var[j, (nlat - i + 1), k]
                 # For Hyde
                 # var_dat[((k-1)*nlat*nlon+(i-1)*nlon+j)] <- var[j, i, k]
               }
            }
         }
      }
      # Clear all NAs
      var_dat[is.na(var_dat)] <- 0.
      # 3. Write the vector into ASCII file
      if(save_ascii) {
         write.table(var_dat, file = as.character(outlist[f]), eol = "\n", dec = ".", row.names = FALSE,
                    col.names = FALSE)
      }

   } else {
      # 2-D
      var <- var_raw[,,1]
      # 1. Consider writing area or directly the fraction?
      if(as_fraction) {
         # Do nothing
         return
      } else {
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
         var <- var*grid_area
      }
      # 2. Transfer matrix into vector
      if(use_mask) {
         # Quick initialization
         var_dat <- 1:tpoint
         # Read in the mask
         igpmap <- as.matrix(read.table(fldmaskdat, header=FALSE, sep=""))
         for (i in 1:tpoint) {
             # For SSP2
             var_dat[i] <- var[igpmap[i,2],(nlat - igpmap[i,3] + 1)]
             # For HYDE
            # var_dat[i] <- var[igpmap[i,2], igpmap[i,3]]
         }
      } else {
         # Quick initialization
         var_dat <- 1:(nlat*nlon)
         # See if this loop is too slow.
         for (i in 1:nlat) {
            for (j in 1:nlon) {
              # For SSP2
               var_dat[((i-1)*nlon+j)] <- var[j,(nlat - i + 1)]
              # For HYDE
              # var_dat[((i-1)*nlon+j)] <- var[j, i]
            }
         }
      }
      # Clear all NAs
      var_dat[is.na(var_dat)] <- 0.
      # 3. Write the vector into ASCII file
      if(save_ascii) {
         write.table(var_dat, file = as.character(outlist[f]), eol = "\n", dec = ".", row.names = FALSE,
                    col.names = FALSE)
      }
   }  # END OF 2/3-D
}
