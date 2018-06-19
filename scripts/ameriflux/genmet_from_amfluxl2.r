library(ncdf4)

# Read in Ameriflux netcdf file
readnc <- function(fname) {
   ncdata <- nc_open(fname)
 
   ta <- ncvar_get(ncdata, 'TA')
   rh <- ncvar_get(ncdata, 'RH')
   ws <- ncvar_get(ncdata, 'WS')
   prec <- ncvar_get(ncdata, 'PREC')
   press <- ncvar_get(ncdata, 'PRESS')
   rg <- ncvar_get(ncdata, 'Rg')
   rgl <- ncvar_get(ncdata, 'Rgl')
   co2 <- ncvar_get(ncdata, 'CO2')
   flux_data <- rbind(ta, rh, ws,prec,press,rg,rgl,co2)

   return(flux_data)
 
}

# Begin of the code
# Read in data ever year
prr2011 <- readnc('AMF_USPrr_2011_L2_GF_V001.nc')
prr2012 <- readnc('AMF_USPrr_2012_L2_GF_V001.nc')
#prr2013 <- readnc('AMF_USPrr_2013_L2_GF_V001.nc')
#prr2014 <- readnc('AMF_USPrr_2014_L2_GF_V001.nc')

prrdata <- prr2011
#prrdata <- cbind(prr2011, prr2012, prr2013, prr2014)
siz = dim(prrdata)

#print(prrdata[c('prec'), 9000:14000])

# Tranform unit when necessary
# TA: degreeC -> K
# PREC: kg/m2/30min -> kg/m2/s
# PRESS: kPa -> Pa
# Other variables are in the consistent unit
prrdata[c('ta'),]    <- prrdata[c('ta'),] + 273.16
prrdata[c('prec'),]  <- prrdata[c('prec'),] / 1800.
prrdata[c('press'),] <- prrdata[c('press'),] * 1000
for(i in 1:siz[2]) {
   if(prrdata[c('ta'), i] < 0 ){
     prrdata[c('ta'), i] <- prrdata[c('ta'), i+48]
   }
   if(prrdata[c('rh'), i] < 0. ){
     if(i > 1 ){
        prrdata[c('rh'), i] <- prrdata[c('rh'), i-1]
     } else {
        prrdata[c('rh'), i] <- prrdata[c('rh'), i+1]
     }
   }
   if(prrdata[c('ws'), i] < 0 ){
     prrdata[c('ws'), i] <- prrdata[c('ws'), i-48]
   }
   if(prrdata[c('press'), i] < 0 ){
     prrdata[c('press'), i] <- prrdata[c('press'), i+48]
   }
   if(prrdata[c('prec'), i] < 0 ){
     prrdata[c('prec'), i] <- 0.
   }
   if(prrdata[c('rg'), i] < 0 ){
     prrdata[c('rg'), i] <- 0.
   }
   if(prrdata[c('rgl'), i] < 0 ){
     if(i > 1){
        prrdata[c('rgl'), i] <- prrdata[c('rgl'), i-1]
     } else {
        prrdata[c('rgl'), i] <- prrdata[c('rgl'), i+1]
     }
   }
   if(prrdata[c('co2'), i] < 0 ){
     prrdata[c('co2'), i] <- prrdata[c('co2'), i+48]
   }
}

# Write NetCDF file 

t <- seq(1,siz[2],1)
x <- 1
y <- 1
z <- 1

dim4 <- ncdim_def('x', 'longitude', as.double(x))
dim3 <- ncdim_def('y', 'latitude', as.double(y))
dim2 <- ncdim_def('z', 'vertical metres', as.double(z))
dim1 <- ncdim_def('t', '30min', as.double(t))
 
ta <- ncvar_def('Tair', 'Kelvin', list(dim1, dim2, dim3, dim4), -9999., longname = 'Near surface air temperature')
rh <- ncvar_def('RH', '100%', list(dim1, dim2, dim3, dim4), -9999., longname = 'Near surface relative humidity')
ws <- ncvar_def('Wind', 'm/s', list(dim1, dim2, dim3, dim4), -9999., longname = 'Near surface module of the wind')
prec <- ncvar_def('Rainf', 'kg/m2/s', list(dim1, dim2, dim3, dim4), -9999., longname = 'Rainfall rate')
press <- ncvar_def('Psurf', 'Pa', list(dim1, dim2, dim3, dim4), -9999., longname = 'Surface pressure')
rg <- ncvar_def('SWdown', 'W/m2', list(dim1, dim2, dim3, dim4), -9999., longname = 'Surface incident shortwave radiation')
rgl <- ncvar_def('LWdown', 'W/m2', list(dim1, dim2, dim3, dim4), -9999., longname = 'Surface incident longwave radiation')
co2 <- ncvar_def('CO2air', 'ppm', list(dim1, dim2, dim3, dim4), -9999., longname = 'Near surface CO2 concentration')
 
outnc <- nc_create('US-Prrforcing.nc', list(ta, rh, ws, prec, press, rg, rgl, co2))
 
ncvar_put(outnc, ta, prrdata[c('ta'),])
ncvar_put(outnc, rh, prrdata[c('rh'),])
ncvar_put(outnc, ws, prrdata[c('ws'),])
ncvar_put(outnc, prec, prrdata[c('prec'),])
ncvar_put(outnc, press, prrdata[c('press'),])
ncvar_put(outnc, rg, prrdata[c('rg'),])
ncvar_put(outnc, rgl, prrdata[c('rgl'),])
ncvar_put(outnc, co2, prrdata[c('co2'),])
 
#ncatt_put(outnc, 't', 'units', '30 mins steps since 2011-01-01')
#ncatt_put(outnc, 'SN', 'axis', 'To read as Y axis data, from South to North')
#ncatt_put(outnc, 0, 'Author', 'The DataJoy Team')
 
nc_close(outnc)
