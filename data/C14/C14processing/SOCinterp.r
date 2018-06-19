## Load packages
library(aqp)         # Algorithms for Quantitative Pedology
library(plyr)        # Tools for Splitting, Applying and Combining Data 
library(sp)          # Classes and Methods for Spatial Data
library(rgdal)
library(GSIF)        # Global Soil Information Facilities

### Main Progame Start Here
# Fetch command line arguments
# myArgs <- commandArgs(trailingOnly = TRUE)

# Convert to numerics
#fname <- myArgs[1]

## User predefined variables
plot_ensemble <- FALSE

## Read in all circumpolar data from comma delimited CSV
# Table fields: ProfileID, Veg Class, Lat, Lon, top, bottom, Layer depth, Soil Order, Suborder, Internal,
#               AL_depth (cm), Horizon type, pH, bulk density (kg/m3), C Density (kg/m3), N Density (kg/m3),
#               C/N
#mishradata <- read.table(fname, header=TRUE, sep=",")
# mishradata <- read.table("sel_sites_bd.csv", header=TRUE, sep=",", quote="\"")
mishradata <- read.table("sel_sites_soc.csv", header=TRUE, sep=",", quote="\"")

## Also calculate soil node depth in mishra data for plotting purpose
mishradata$idx <- seq(1,nrow(mishradata),1)
mishradata$top <- mishradata$Basal_depth - mishradata$Layer_depth
mishradata$bottom <- mishradata$Basal_depth
mishradata$nodedepth <- mishradata$top + (mishradata$bottom - mishradata$top)/2
# C.Den [g cm-3], depth [cm], Ccontent [g cm-2] -> 10 [kg m-2]
mishradata$Ccontent <- mishradata$C_Density * mishradata$Layer_depth * 10.   # From g/cm2 to kg/m2

## Extract Site and Soil Order information. Should do it before transfering to SPC class
info <- mishradata[, c("ProfileID","Veg_Class", "Lat","Lon","top","bottom","nodedepth","Layer_depth", "Soil_Order",
                        "Suborder","Horizon_type","bulk_density","C_Density","Ccontent")]

## Collapse the data frame by deleting redundant records
info <- unique(info)

## Quality control before starting interpolation

## Promote Mishra dataset to SPC class
depths(mishradata) <- ProfileID ~ top + bottom   # Promote normal data.frame to SoilProfileCollection data frame
mishradata$idcol <- mishradata$ProfileID

## Fit MP spline by profile
try(m<-mpspline(mishradata, 'C_Density', lam = 0.1))
#try(m<-mpspline(mishradata, 'bulk_density', lam = 0.1))
#try(m<-mpspline(mishradata, 'C_Density', lam = 0.5))
#try(m<-mpspline(mishradata, 'Ccontent', lam = 0.1))

## Sould have these profiles fitting failed: 
# 45, 46, 47, 169, 186, 221, 249

## Sorting interpolated profiles by ID
m$idcol <- as.numeric(m$idcol)
m$interpolated <- t(m$var.1cm[,order(m$idcol)])
m$fitted <- m$var.fitted[order(m$idcol),]
m$idcol_sorted <- m$idcol[order(m$idcol)]

# Write out interpolated SOC profile
write.csv(m$fitted, file='SOCfitted.csv', na="", row.names = TRUE)
write.csv(m$interpolated, file='SOCprofile.csv', na="", row.names = TRUE)
# write.csv(m$fitted, file='BDfitted.csv', na="", row.names = TRUE)
# write.csv(m$interpolated, file='BDprofile.csv', na="", row.names = TRUE)


#DONE
