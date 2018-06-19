## Load packages
library(aqp)         # Algorithms for Quantitative Pedology
library(plyr)        # Tools for Splitting, Applying and Combining Data 
library(sp)          # Classes and Methods for Spatial Data
library(GSIF)        # Global Soil Information Facilities

## Define functions
# 1. Extract function
extract <- function(m, info, key, type) {
  # Check if arguments are missing
  if(missing(m) || missing(info) || missing(key) || missing(type)) {
    cat("Error, missing necessary input!\n")
    cat("Usage: \n")
    cat("result <- extract(m, info, key, type) \n")
    cat("m: List contains the value needed to be extracted. Must have at least one value column \n")
    cat("'interpolated' and another id column 'idcol_sorted'. \n")
    cat("info: Data frame contains the auxilliary information needed to be extracted. Must have at least one factor column \n")
    cat("as specified by type. \n")
    cat("key: Specific value in the factor that will be extracted. \n")
    cat("type: Corresponding factor field in the 'info' data frame that contains the value of 'key'. \n")
    return(NA)
  }
  # Validate necessary columns needed for the function
  if(!(type %in% colnames(info)))
  {
    cat("Error, the type you specified is not in info list!\n")
    cat("Usage: \n")
    cat("result <- extract(m, info, key, type) \n")
    cat("m: List contains the value needed to be extracted. Must have at least one value column \n")
    cat("'interpolated' and another id column 'idcol_sorted'. \n")
    cat("info: Data frame contains the auxilliary information needed to be extracted. Must have at least one factor column \n")
    cat("as specified by type. \n")
    cat("key: Specific value in the factor that will be extracted. \n")
    cat("type: Corresponding factor field in the 'info' data frame that contains the value of 'key'. \n")
    return(NA)
  }
  if(!('interpolated' %in% names(m)))
  {
    cat("Error, could not find 'interpolated' in the list provided!\n")
    cat("Usage: \n")
    cat("result <- extract(m, info, key, type) \n")
    cat("m: List contains the value needed to be extracted. Must have at least one value column \n")
    cat("'interpolated' and another id column 'idcol_sorted'. \n")
    cat("info: Data frame contains the auxilliary information needed to be extracted. Must have at least one factor column \n")
    cat("as specified by type. \n")
    cat("key: Specific value in the factor that will be extracted. \n")
    cat("type: Corresponding factor field in the 'info' data frame that contains the value of 'key'. \n")
    return(NA)
  }
  if(!('idcol_sorted' %in% names(m)))
  {
    cat("Error, could not find 'idcol_sorted' in the list provided!\n")
    cat("Usage: \n")
    cat("result <- extract(m, info, key, type) \n")
    cat("m: List contains the value needed to be extracted. Must have at least one value column \n")
    cat("'interpolated' and another id column 'idcol_sorted'. \n")
    cat("info: Data frame contains the auxilliary information needed to be extracted. Must have at least one factor column \n")
    cat("as specified by type. \n")
    cat("key: Specific value in the factor that will be extracted. \n")
    cat("type: Corresponding factor field in the 'info' data frame that contains the value of 'key'. \n")
    return(NA)
  }
  # Extract corresponding interpolated field according to specified field and value
  # Build the expression
  expr <- paste('m$interpolated[,m$idcol_sorted[info$', type, ' == key]]', sep='')
  # Execute
  result = eval(parse(text=expr))
  return(result)

} # END OF THE FUNCTION

# 2. Get ID Function
getid <- function(m, info, key, type) {
  # Check if arguments are missing
  if(missing(m) || missing(info) || missing(key) || missing(type)) {
    cat("Error, missing necessary input!\n")
    cat("Usage: \n")
    cat("result <- getid(m, info, key, type) \n")
    cat("m: List contains the corresponding ID needed to be extracted. Must have at least one value column \n")
    cat("'interpolated' and another id column 'idcol_sorted'. \n")
    cat("info: Data frame contains the auxilliary information needed to be extracted. Must have at least one factor column \n")
    cat("as specified by type. \n")
    cat("key: Specific value in the factor that will be extracted. \n")
    cat("type: Corresponding factor field in the 'info' data frame that contains the value of 'key'. \n")
    return(NA)
  }
  # Validate necessary columns needed for the function
  if(!(type %in% colnames(info)))
  {
    cat("Error, the type you specified is not in info list!\n")
    cat("Usage: \n")
    cat("result <- getid(m, info, key, type) \n")
    cat("m: List contains the corresponding ID needed to be extracted. Must have at least one value column \n")
    cat("'interpolated' and another id column 'idcol_sorted'. \n")
    cat("info: Data frame contains the auxilliary information needed to be extracted. Must have at least one factor column \n")
    cat("as specified by type. \n")
    cat("key: Specific value in the factor that will be extracted. \n")
    cat("type: Corresponding factor field in the 'info' data frame that contains the value of 'key'. \n")
    return(NA)
  }
  if(!('interpolated' %in% names(m)))
  {
    cat("Error, could not find 'interpolated' in the list provided!\n")
    cat("Usage: \n")
    cat("result <- getid(m, info, key, type) \n")
    cat("m: List contains the corresponding ID needed to be extracted. Must have at least one value column \n")
    cat("'interpolated' and another id column 'idcol_sorted'. \n")
    cat("info: Data frame contains the auxilliary information needed to be extracted. Must have at least one factor column \n")
    cat("as specified by type. \n")
    cat("key: Specific value in the factor that will be extracted. \n")
    cat("type: Corresponding factor field in the 'info' data frame that contains the value of 'key'. \n")
    return(NA)
  }
  if(!('idcol_sorted' %in% names(m)))
  {
    cat("Error, could not find 'idcol_sorted' in the list provided!\n")
    cat("Usage: \n")
    cat("result <- getid(m, info, key, type) \n")
    cat("m: List contains the corresponding ID needed to be extracted. Must have at least one value column \n")
    cat("'interpolated' and another id column 'idcol_sorted'. \n")
    cat("info: Data frame contains the auxilliary information needed to be extracted. Must have at least one factor column \n")
    cat("as specified by type. \n")
    cat("key: Specific value in the factor that will be extracted. \n")
    cat("type: Corresponding factor field in the 'info' data frame that contains the value of 'key'. \n")
    return(NA)
  }
  # Extract corresponding interpolated field according to specified field and value
  # Build the expression
  expr <- paste('m$idcol_sorted[info$', type, ' == key]', sep='')
  # Execute
  result = eval(parse(text=expr))
  return(result)

} # END OF THE FUNCTION


### Main Progame Start Here

## User predefined variables
plot_ensemble <- FALSE
classification <- 'SoilEcoregions'  # SoilOrder, Ecoregions, Landcover, SoilEcoregions, SoilLandcover

## Read in all circumpolar data from comma delimited CSV
# Table fields: ID, Profile ID, Lat, Lon, top, bottom, Layer depth, Soil Order, Suborder, Internal,
#               AL_depth (cm), Horizon type, pH, bulk density (kg/m3), C Density (kg/m3), N Density (kg/m3),
#               C/N
mishradata <- read.table("mishradata.csv", header=TRUE, sep=",")
#mishradata <- read.table("circumpolar_data.csv", header=TRUE, sep=",")

## Also calculate soil node depth in mishra data for plotting purpose
mishradata$idx <- seq(1,nrow(mishradata),1)
mishradata$nodedepth <- mishradata$top + (mishradata$bottom - mishradata$top)/2
mishradata$Ccontent <- mishradata$C.Density * mishradata$Layer.depth / 100   # From kg/m3 to kg/m2

## Extract Site and Soil Order information. Should do it before transfering to SPC class
info <- mishradata[, c("ID","Lat","Lon","Soil.Order","Suborder","Ecoregions","Landcover")]

## Shrink the data frame by deleting redundant records
info <- unique(info)

## Quality control before starting interpolation

## Promote Mishra dataset to SPC class
depths(mishradata) <- ID ~ top + bottom   # Promote normal data.frame to SoilProfileCollection data frame

## Fit MP spline by profile
#try(m<-mpspline(mishradata, 'C.Density', lam = 0.1))
try(m<-mpspline(mishradata, 'Ccontent', lam = 0.1))

## Sould have these profiles fitting failed: 
# 45, 46, 47, 169, 186, 221, 249

## Sorting interpolated profiles by ID
m$idcol <- as.numeric(m$idcol)
m$interpolated <- m$var.1cm[,order(m$idcol)]
m$idcol_sorted <- m$idcol[order(m$idcol)]

## Calculate mean soil C profile based on defined classification: Soil order, Veg type, Soil order + veg type

if(classification == 'SoilOrder') {
   source("meanprof_soil.r")
}
if(classification == 'Landcover') {
   source("meanprof_landcover.r")
}
if(classification == 'Ecoregions') {
   source("meanprof_ecoregion.r")
}
if(classification == 'SoilEcoregions') {
   source("meanprof_soil_ecoregion.r")
}
if(classification == 'SoilLandcover') {
   source("meanprof_soil_lc.r")
}

## Write out interpolated SOC profile
# write.csv(m$var.fitted, file='SOCfitted.csv', na="", row.names = TRUE)
# write.csv(m$var.1cm, file='SOCprofile.csv', na="", row.names = TRUE)

################# DONE

## subplot two
###par(mar=c(0.6,0.6,8,8)+0.1)
###plot(bdecid$Measured, bdecid$depth, xlab = "", ylab = "", 
###      xlim=c(0.0,1.0), ylim=rev(c(0,100)), pch=17, cex = 4, xaxt='n', yaxt='n')
###axis(3,at=seq(0.0,1.0,by=0.2), cex.axis=3)
###axis(3,at=seq(0.0,1.0,by=0.2),labels=FALSE,tcl=-0.25)
###grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
###lines(bdecid$Dynamic, bdecid$depth.1, lwd = 4, xaxt='n', yaxt='n')
###lines(bdecid$Static, bdecid$depth.2, lwd = 4, xaxt='n', yaxt='n', col="gray50")
###text(0.03, 10, "(b)", cex = 3)
##### subplot three
###par(mar=c(8,8,0.6,0.6)+0.1)
###plot(grass$Measured, grass$depth, xlab = "", ylab = "", 
###      xlim=c(0.0,1.0), ylim=rev(c(0,100)), pch=17, cex = 4, xaxt='n', yaxt='n')
###axis(2,at=seq(100,0,by=-20), cex.axis=3, las=1)
###axis(2,at=seq(100,0,by=-20),labels=FALSE,tcl=-0.25)
###grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
###lines(grass$Dynamic, grass$depth.1, lwd = 4, xaxt='n', yaxt='n')
###lines(grass$Static, grass$depth.2, lwd = 4, xaxt='n', yaxt='n', col="gray50")
###text(0.03, 10, "(c)", cex = 3)
#### subplot four
###par(mar=c(8,0.6,0.6,8)+0.1)
###plot(bever$Measured, bever$depth, xlab = "", ylab = "", 
###      xlim=c(0.0,1.0), ylim=rev(c(0,100)), pch=17, cex = 4, cex.lab=1, cex.axis = 1, xaxt='n', yaxt='n')
###grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
###lines(bever$Dynamic, bever$depth.1, lwd = 4, xaxt='n', yaxt='n')
###lines(bever$Static, bever$depth.2, lwd = 4, xaxt='n', yaxt='n', col="gray50")
#### Add subplot labels
###text(0.03, 10, "(d)", cex = 3)

