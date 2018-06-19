## Load packages
library(aqp)         # Algorithms for Quantitative Pedology
library(plyr)        # Tools for Splitting, Applying and Combining Data 
library(sp)          # Classes and Methods for Spatial Data
library(GSIF)        # Global Soil Information Facilities
library(ncdf4)       # NetCDF 4 library
library(Hmisc)

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

## User options
plot_ensemble <- FALSE
extract_suborder <- FALSE
inerpolate_model <- FALSE
write_csv <- FALSE
group_data <- FALSE
group_figure <- FALSE
statistical_analysis <- FALSE
readin_model <- FALSE
classification <- 'Ecoregions'  # SoilOrder, Ecoregions, Landcover, SoilEcoregions, SoilLandcover


if(readin_model) {
  model_fname <- 'socprofeq.dat'
  model_idfname <- 'caselist'
  model_ndepth <- c(0.0071006, 0.0279250, 0.0622586, 0.1188651, 0.2121934, 0.3660658, 0.6197585, 1.0380271, 1.7276353, 2.8646071)
  model_ldepth <- c(0.017513, 0.027579, 0.045470, 0.074967, 0.123600, 0.203783, 0.335981, 0.553938, 0.913290, 1.136972)
  group_model <- TRUE
  group_veg <- FALSE
} else{
  group_model <- FALSE
  group_veg <- FALSE
}
if(group_veg){
   veg_fname <- 'vegcover.dat'
}
 

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
# For Mishra data only
info <- mishradata[, c("ID","Profile.ID","Lat","Lon","Soil.Order","Suborder","Ecoregions","Landcover")]
# For Circumpolar data
# info <- mishradata[, c("ID","Lat","Lon","Soil.Order","Suborder")]

## Shrink the data frame by deleting redundant records
info <- unique(info)

## Quality control before starting interpolation

## Promote Mishra dataset to SPC class
depths(mishradata) <- ID ~ top + bottom   # Promote normal data.frame to SoilProfileCollection data frame

## Fit MP spline by profile
try(m<-mpspline(mishradata, 'C.Density', lam = 0.1))    # kgC/m3
#try(m<-mpspline(mishradata, 'Ccontent', lam = 0.1))     # kgC/m2

## Sould have these profiles fitting failed: 
# 45, 46, 47, 169, 186, 221, 249

## Sorting interpolated profiles by ID
m$idcol <- as.numeric(m$idcol)
m$interpolated <- m$var.1cm[,order(m$idcol)]
m$idcol_sorted <- m$idcol[order(m$idcol)]

## Extract Gelisols based on suborders by looking at data
if(extract_suborder) {
   mishradf <- as(mishradata, 'data.frame')
   c_turbel <- mishradf[mishradf$Suborder == 'Turbel',]
   c_orthel <- mishradf[mishradf$Suborder == 'Orthel',]
   c_histel <- mishradf[mishradf$Suborder == 'Histel',]
   write.csv(c_turbel, file='Turbel.csv', na="", row.names = FALSE)
   write.csv(c_orthel, file='Orthel.csv', na="", row.names = FALSE)
   write.csv(c_histel, file='Histel.csv', na="", row.names = FALSE)
   write.csv(info[info$Soil.Order == 'Gelisol',], file='gelisol_info.csv', na="", row.names = FALSE)
}

# Read in model results
if(readin_model) {
  soc_modeled <- t(read.table(model_fname, sep="", na.strings='-9999.', fill=TRUE))
  soc_case <- read.table(model_idfname, sep="/", na.strings='-9999.', fill=TRUE)
  if(group_veg) {
     vegcase <- read.table(veg_fname, sep="", na.strings='-9999.', fill=TRUE)
     veg_case <- vegcase[,"V1"]
  }
  soc_pname <- soc_case[,dim(soc_case)[2]]
  soc_den_modeled <- soc_modeled / model_ldepth
  pid_modeled_idx <- info$ID[info$Profile.ID %in% soc_pname]
  eid_modeled <- info$Ecoregions[pid_modeled_idx]
  rm(soc_case)
  # Mean SOC density profile grouped by ecoregion
  if(group_model){
     if(group_veg){
        veglevel <- levels(as.factor(veg_case))
        model_mean <- matrix(0., dim(soc_den_modeled)[1], length(veglevel))
        model_sd <- model_mean
        model_ci <- model_mean
        cnt <- 1
        for (i in as.numeric(veglevel)) {
           model_subset <- soc_den_modeled[,veg_case == i]
           if(is.null(dim(model_subset))) {
             # model_subset[is.na(model_subset)] <- 0.
              model_mean[,cnt] <- model_subset
           } else {
              model_mean[,cnt] <- rowMeans(model_subset, na.rm=TRUE)
              model_sd[,cnt] <- apply(model_subset, 1, sd, na.rm=TRUE)
              model_ci[,cnt] <- 1.96 * model_sd[,cnt] / sqrt(rowSums(!is.na(model_subset)))
           }
           rm(model_subset)
           cnt <- cnt + 1
        }
     } else {
        ecolevel <- levels(as.factor(eid_modeled))
        model_mean <- matrix(0., dim(soc_den_modeled)[1], length(ecolevel))
        model_sd <- model_mean
        model_ci <- model_mean
        cnt <- 1
        for (i in as.numeric(ecolevel)) {
           model_subset <- soc_den_modeled[,eid_modeled == i]
           if(is.null(dim(model_subset))) {
             # model_subset[is.na(model_subset)] <- 0.
              model_mean[,cnt] <- model_subset
           } else {
              model_mean[,cnt] <- rowMeans(model_subset, na.rm=TRUE)
              model_sd[,cnt] <- apply(model_subset, 1, sd, na.rm=TRUE)
              model_ci[,cnt] <- 1.96 * model_sd[,cnt] / sqrt(rowSums(!is.na(model_subset)))
           }
           rm(model_subset)
           cnt <- cnt + 1
        }
     }
  }
}

# Group Gelisol sample profile by ecoregion
if(group_data) {
   # obtain interpolated SOC density data for Gelisols
   pid_turbel_idx <- getid(m, info, c('Turbel'), 'Suborder')
   pid_orthel_idx <- getid(m, info, c('Orthel'), 'Suborder')
   pid_histel_idx <- getid(m, info, c('Histel'), 'Suborder')
   pid_gelisol_idx <- info$ID[info$Profile.ID %in% soc_pname] #c(pid_turbel_idx, pid_orthel_idx)#, pid_histel_idx)
#   pid_gelisol_idx <- pid_turbel_idx

   pid_turbel <- info$Profile.ID[pid_turbel_idx]
   pid_orthel <- info$Profile.ID[pid_orthel_idx]
   pid_histel <- info$Profile.ID[pid_histel_idx]
   pid_gelisol <- info$Profile.ID[pid_gelisol_idx]
#   pid_gelisol <- info$Profile.ID[pid_gelisol_idx]

   eid_turbel <- info$Ecoregions[pid_turbel_idx]
   eid_orthel <- info$Ecoregions[pid_orthel_idx]
   eid_histel <- info$Ecoregions[pid_histel_idx]
   eid_gelisol <- info$Ecoregions[pid_gelisol_idx]
#   eid_gelisol <- info$Ecoregions[pid_gelisol_idx]

   selected_interpolated_turbel <- extract(m ,info, c('Turbel'), 'Suborder')
   selected_interpolated_orthel <- extract(m ,info, c('Orthel'), 'Suborder')
   selected_interpolated_histel <- extract(m ,info, c('Histel'), 'Suborder')
   selected_interpolated_gelisol <- m$interpolated[,pid_gelisol_idx] #cbind(selected_interpolated_turbel, selected_interpolated_orthel)#, selected_interpolated_histel)
# Totally 40 samples
#   selected_interpolated_gelisol <- selected_interpolated_turbel

   # Group them by ecoregion and calculate statistics
   if(group_veg){
      veglevel <- levels(as.factor(veg_case))
      veg_case_model <- veg_case[soc_pname %in% info$Profile.ID[pid_gelisol_idx]]
      gelisol_mean <- matrix(0., dim(selected_interpolated_gelisol)[1], length(veglevel))
      gelisol_sd <- gelisol_mean
      gelisol_ci <- gelisol_mean
      cnt <- 1
      for (i in as.numeric(veglevel)) {
         gelisol_subset <- selected_interpolated_gelisol[,veg_case_model == i]
         if(is.null(dim(gelisol_subset))) {
           # gelisol_subset[is.na(gelisol_subset)] <- 0.
            gelisol_mean[,cnt] <- gelisol_subset
         } else {
            gelisol_mean[,cnt] <- rowMeans(gelisol_subset, na.rm=TRUE)
            gelisol_sd[,cnt] <- apply(gelisol_subset, 1, sd, na.rm=TRUE)
            gelisol_ci[,cnt] <- 1.96 * gelisol_sd[,cnt] / sqrt(rowSums(!is.na(gelisol_subset)))
         }
         rm(gelisol_subset)
         cnt <- cnt + 1
      }
   } else{
      ecolevel <- levels(as.factor(eid_gelisol))
      gelisol_mean <- matrix(0., dim(selected_interpolated_gelisol)[1], length(ecolevel))
      gelisol_sd <- gelisol_mean
      gelisol_ci <- gelisol_mean
      cnt <- 1
      for (i in as.numeric(ecolevel)) {
         gelisol_subset <- selected_interpolated_gelisol[,eid_gelisol == i]
         if(is.null(dim(gelisol_subset))) {
           # gelisol_subset[is.na(gelisol_subset)] <- 0.
            gelisol_mean[,cnt] <- gelisol_subset
         } else {
            gelisol_mean[,cnt] <- rowMeans(gelisol_subset, na.rm=TRUE)
            gelisol_sd[,cnt] <- apply(gelisol_subset, 1, sd, na.rm=TRUE)
            gelisol_ci[,cnt] <- 1.96 * gelisol_sd[,cnt] / sqrt(rowSums(!is.na(gelisol_subset)))
         }
         rm(gelisol_subset)
         cnt <- cnt + 1
      }
   }

   # Calculate accumulated SOC content
   
}

## Produce figures of Mean SOC density profiles and accumulated SOC profiles for Gelisol grouped by ecoregion or vegcover
if(group_data & group_model & group_figure){
   # Setting up
   depth_obs <- seq(1,200,1)
   depth_model <- model_ndepth*100

   if(group_veg) {
     # pdf(file="Grouped_ecoregion.pdf", width=17, height=10)
      pdf(file="Grouped_veg.pdf", width=17, height=10)
      ### Failed to capture SOC density
      # Set margin
      par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page
      # 
      plot(gelisol_mean[,1], depth_obs, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
                  xlim=c(0,200), ylim=rev(c(0,280)), pch=17, cex= 0, xaxt='n', yaxt='n')
      lines(gelisol_mean[,1], depth_obs, lty = 1, lwd = 4, col="red", xaxt='n', yaxt='n')    # Boreal
      lines(gelisol_mean[,1]+gelisol_sd[,1], depth_obs, lwd = 4, xaxt='n', yaxt='n', col="indianred")
      lines(gelisol_mean[,1]-gelisol_sd[,1], depth_obs, lwd = 4, xaxt='n', yaxt='n', col="indianred")
      lines(gelisol_mean[,2], depth_obs, lty = 1, lwd = 4, col="blue", xaxt='n', yaxt='n')    # Tundra
      lines(gelisol_mean[,2]+gelisol_sd[,2], depth_obs, lwd = 4, xaxt='n', yaxt='n', col="deepskyblue")
      lines(gelisol_mean[,2]-gelisol_sd[,2], depth_obs, lwd = 4, xaxt='n', yaxt='n', col="deepskyblue")
 
      lines(model_mean[,1], depth_model, lty = 6, lwd = 4, col="red", xaxt='n', yaxt='n')    # Boreal
      #lines(c_1_mean+c_1_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
      #lines(c_1_mean-c_1_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
      lines(model_mean[,2], depth_model, lty = 6, lwd = 4, col="blue", xaxt='n', yaxt='n')    # Tundra
 
      # Move axis to top and left
      axis(3,at=seq(0.0,200.0,by=20),cex.axis=3)
      axis(3,at=seq(0.0,200.0,by=20),labels=FALSE,tcl=-0.25)
      axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
      axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
      grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
    #  legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
      legend(60, 80, legend = c("Obs - Boreal", "Obs - Tundra"), cex=3, lty=1, col=c("red","blue"), lwd=3, bty = "n")#, lty = 1:7,
      legend(60, 160, legend = c("Model - Boreal", "Model - Tundra"), cex=3, lty=6, col=c("red","blue"), lwd=3, bty = "n")#, lty = 1:7,
           #  pch = "*", ncol = 4, cex = 0.8)
      # Put Label for X and Y axis
      mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
      mtext("Organic Carbon Density (kgC/m3)", side=3, outer=TRUE,line=-3.7,cex=3)
     # text(120, 120, "Mean Profiles", cex = 3)
      dev.off()
   } else {
     # pdf(file="Grouped_ecoregion.pdf", width=17, height=10)
      pdf(file="Grouped_ecoregion.pdf", width=17, height=10)
      ### Failed to capture SOC density
      # Set margin
      par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page
      # 
      plot(gelisol_mean[,2], depth_obs, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
                  xlim=c(0,200), ylim=rev(c(0,160)), pch=17, cex= 0, xaxt='n', yaxt='n')
      collect <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120)
      Dotplot(c("","","","","","","","","","") ~ Cbind(gelisol_mean[collect,2],
                                                       gelisol_mean[collect,2]-gelisol_sd[collect,2],
                                                       gelisol_mean[collect,2]+gelisol_sd[collect,2]), col="blue", pch=20, xlab="n",ylab="n")
      lines(gelisol_mean[,2], depth_obs, lty = 1, lwd = 4, col="red", xaxt='n', yaxt='n')    # Ecoregion 5= Arctic Tundra
      lines(gelisol_mean[,2]+gelisol_sd[,2], depth_obs, lwd = 4, xaxt='n', yaxt='n', col="indianred")
      lines(gelisol_mean[,2]-gelisol_sd[,2], depth_obs, lwd = 4, xaxt='n', yaxt='n', col="indianred")
    #  lines(gelisol_mean[,5], depth_obs, lty = 1, lwd = 4, col="blue", xaxt='n', yaxt='n')    # Ecoregion 6 = Bering Tundra
      lines(gelisol_mean[,3], depth_obs, lty = 1, lwd = 4, col="green", xaxt='n', yaxt='n')    # Ecoregion 8 =Intermontane Boreal
      lines(gelisol_mean[,3]+gelisol_sd[,3], depth_obs, lwd = 4, xaxt='n', yaxt='n', col="deepskyblue")
      lines(gelisol_mean[,3]-gelisol_sd[,3], depth_obs, lwd = 4, xaxt='n', yaxt='n', col="deepskyblue")
   
      lines(model_mean[,2], depth_model, lty = 6, lwd = 4, col="red", xaxt='n', yaxt='n')    # Ecoregion 5= Arctic Tundra
      #lines(c_1_mean+c_1_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
      #lines(c_1_mean-c_1_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
    #  lines(model_mean[,4], depth_model, lty = 6, lwd = 4, col="blue", xaxt='n', yaxt='n')    # Ecoregion 6 = Bering Tundra
      lines(model_mean[,3], depth_model, lty = 6, lwd = 4, col="green", xaxt='n', yaxt='n')    # Ecoregion 8 =Intermontane Boreal
   
      # Move axis to top and left
      axis(3,at=seq(0.0,160.0,by=20),cex.axis=3)
      axis(3,at=seq(0.0,160.0,by=20),labels=FALSE,tcl=-0.25)
      axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
      axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
      grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
    #  legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
   #   legend(60, 80, legend = c("Obs - Arctic Tundra", "Obs - Bering Tundra", "Obs - Intermontane Boreal"), cex=3, lty=1, col=c("red","blue","green"), lwd=3, bty = "n")#, lty = 1:7,
   #   legend(60, 160, legend = c("Model - Arctic Tundra", "Model - Bering Tundra", "Model - Intermontane Boreal"), cex=3, lty=6, col=c("red","blue","green"), lwd=3, bty = "n")#, lty = 1:7,
      legend(60, 80, legend = c("Obs - Arctic Tundra", "Obs - Intermontane Boreal"), cex=3, lty=1, col=c("red","green"), lwd=3, bty = "n")#, lty = 1:7,
      legend(60, 120, legend = c("Model - Arctic Tundra", "Model - Intermontane Boreal"), cex=3, lty=6, col=c("red","green"), lwd=3, bty = "n")#, lty = 1:7,
           #  pch = "*", ncol = 4, cex = 0.8)
      # Put Label for X and Y axis
      mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
      mtext("Organic Carbon Density (kgC/m3)", side=3, outer=TRUE,line=-3.7,cex=3)
     # text(120, 120, "Mean Profiles", cex = 3)
      dev.off()
   }
}

# Statistical analysis
if(statistical_analysis) {
}

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
if(write_csv) {
   write.table(m$interpolated, file='SOCfitted.csv', na="", sep=",", col.names=as.character(info$Profile.ID))
}

# Averaged to model depth
if(inerpolate_model) {
   avg_len <- c(2, 3, 4, 8, 12, 21, 33, 55, 92, 113)
   dz <- c(0.017513, 0.027579, 0.045470, 0.074967, 0.123600, 0.203783, 0.335981, 0.553938, 0.913290, 1.136972)
   prof <- matrix(0., 10, dim(m$interpolated)[2])
   
   for (entry in 1:dim(m$interpolated)[2]){
      for (i in 1:10) {
         if(i == 1){
            prof[i,entry] <- mean(m$interpolated[1:avg_len[1],entry], na.rm = TRUE)
         }
         if((i >= 2) & (i <= 8)){
            prof[i,entry] <- mean(m$interpolated[(sum(avg_len[1:i-1])+1):sum(avg_len[1:i]), entry], na.rm = TRUE)
         }
         if(i == 9){
            prof[i,entry] <- mean(m$interpolated[(sum(avg_len[1:8])+1):200, entry], na.rm = TRUE)
         }
         if(i == 10){
            prof[i,entry] <- prof[i-1,entry]
         }
         # Let's assume the layer which all obs missing having the same C density as the layer above 
         if(prof[i,entry] < 0.01 | is.na(prof[i,entry])) {
            prof[i,entry] = prof[i-1,entry]
         }
         # From kgC/m3 to kgC/m2
         prof[i,entry] = prof[i,entry] * dz[i]
      }
   }
   
   # Get Profile ID
   
   ## Write out NC file contained
   profnum <- seq(1, dim(m$interpolated)[2], 1)
   nlevsoc <- seq(1, 10, 1)
   
   dim1 <- ncdim_def('profnum', 'Profile ID', as.double(profnum))
   dim2 <- ncdim_def('nlevsoc', 'Vertical Levels for SOC', as.double(nlevsoc))
   
   profID <- ncvar_def('ProfID', 'ID of the profile', list(dim1))
   soc <- ncvar_def('SOC', 'kgC/m2', list(dim2, dim1), -9999., longname = 'Soil organic carbon conetnt')
   
   outnc <- nc_create('SOC_profile.nc', list(profID, soc))
   
   ncvar_put(outnc, profID, profnum)
   ncvar_put(outnc, soc, prof)
   
   nc_close(outnc)
}


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

