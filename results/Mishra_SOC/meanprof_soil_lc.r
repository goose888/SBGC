# Mean profiles classified first by soil order then by ecoregion
# must be called by MishraSOCspline.r
############## Landcover ######################
# 9 different landcover types are involved
# 11 = open water
# 12 = perrenial ice
# 21 = Developed
# 31 = Barren land
# 41 = Forrest
# 51 = Scrub
# 71 = Herbaceous vegetation
# 81 = Pasture
# 90 = wetland
############## Soil Orders ####################
# 9 different landcover types are involved
# 1 = Alfisol
# 2 = Andisol
# 3 = Entisol
# 4 = Histosol
# 5 = Inceptisol
# 6 = Spodosol
# 7 = Histel
# 8 = Orthel
# 9 = Turbel
title_lc = c('Open water', 'Perrenial ice', 'Developed', 'Barren land',
                'Forest', 'Scrub', 'Herbaceous vegetation', 'Pasture', 'Wetland')
code_lc = c(11, 12, 21, 31, 41, 51, 71, 81, 90)

# Alfisol
#c_alfisol = extract(m ,info, 'Alfisol', 'Soil.Order')
#c_alfisol_mean = rowMeans(c_alfisol, na.rm=TRUE)
#c_alfisol_sd = apply(c_alfisol, 1, sd, na.rm=TRUE)
#c_alfisol_se = c_alfisol_sd / rowSums(!is.na(c_alfisol))
# Andisol
# Overall mean profile
c_andisol = extract(m ,info, 'Andisol', 'Soil.Order')
c_andisol_mean = rowMeans(c_andisol, na.rm=TRUE)
c_andisol_sd = apply(c_andisol, 1, sd, na.rm=TRUE)
c_andisol_se = c_andisol_sd / rowSums(!is.na(c_andisol))
# calculate mean and std for each veg type
c_andisol_id = getid(m, info, 'Andisol', 'Soil.Order')
c_andisol_ecoregion = info$Landcover[c_andisol_id]
s_ecoregion <- summary(as.factor(c_andisol_ecoregion))   # Get the occurence frequency by different ecoregion
s_andisol_ecoregion <- s_ecoregion[s_ecoregion>=2]
lenclass = length(s_andisol_ecoregion)
andisol_mean = matrix(0, nrow=200, ncol=lenclass)
andisol_sd   = matrix(0, nrow=200, ncol=lenclass)
andisol_se   = matrix(0, nrow=200, ncol=lenclass)
for(k in 1:lenclass){
   ecoid = as.numeric(names(s_andisol_ecoregion)[k])
   if(s_andisol_ecoregion[k] >= 2){
      andisol_mean[,k] = rowMeans(c_andisol[,c_andisol_ecoregion == ecoid], na.rm=TRUE)
      andisol_sd[,k] = apply(c_andisol[,c_andisol_ecoregion == ecoid], 1, sd, na.rm=TRUE)
      andisol_se[,k] = andisol_sd[,k] / rowSums(!is.na(c_andisol[,c_andisol_ecoregion == ecoid]))
   }
}
# Entisol
# Overall mean profile
c_entisol = extract(m ,info, 'Entisol', 'Soil.Order')
c_entisol_mean = rowMeans(c_entisol, na.rm=TRUE)
c_entisol_sd = apply(c_entisol, 1, sd, na.rm=TRUE)
c_entisol_se = c_entisol_sd / rowSums(!is.na(c_entisol))
# calculate mean and std for each veg type
c_entisol_id = getid(m, info, 'Entisol', 'Soil.Order')
c_entisol_ecoregion = info$Landcover[c_entisol_id]
s_ecoregion <- summary(as.factor(c_entisol_ecoregion))   # Get the occurence frequency by different ecoregion
s_entisol_ecoregion <- s_ecoregion[s_ecoregion>=2]
lenclass = length(s_entisol_ecoregion)
entisol_mean = matrix(0, nrow=200, ncol=lenclass)
entisol_sd   = matrix(0, nrow=200, ncol=lenclass)
entisol_se   = matrix(0, nrow=200, ncol=lenclass)
for(k in 1:lenclass){
   ecoid = as.numeric(names(s_entisol_ecoregion)[k])
   if(s_entisol_ecoregion[k] >= 2){
      entisol_mean[,k] = rowMeans(c_entisol[,c_entisol_ecoregion == ecoid], na.rm=TRUE)
      entisol_sd[,k] = apply(c_entisol[,c_entisol_ecoregion == ecoid], 1, sd, na.rm=TRUE)
      entisol_se[,k] = entisol_sd[,k] / rowSums(!is.na(c_entisol[,c_entisol_ecoregion == ecoid]))
   }
}
# Histosol
#c_histosol = extract(m ,info, 'Histosol', 'Soil.Order')
#c_histosol_mean = rowMeans(c_histosol, na.rm=TRUE)
#c_histosol_sd = apply(c_histosol, 1, sd, na.rm=TRUE)
#c_histosol_se = c_histosol_sd / rowSums(!is.na(c_histosol))
# Inceptisol
# Overall mean profile
c_inceptisol = extract(m ,info, 'Inceptisol', 'Soil.Order')
c_inceptisol_mean = rowMeans(c_inceptisol, na.rm=TRUE)
c_inceptisol_sd = apply(c_inceptisol, 1, sd, na.rm=TRUE)
c_inceptisol_se = c_inceptisol_sd / rowSums(!is.na(c_inceptisol))
# calculate mean and std for each veg type
c_inceptisol_id = getid(m, info, 'Inceptisol', 'Soil.Order')
c_inceptisol_ecoregion = info$Landcover[c_inceptisol_id]
s_ecoregion <- summary(as.factor(c_inceptisol_ecoregion))   # Get the occurence frequency by different ecoregion
s_inceptisol_ecoregion <- s_ecoregion[s_ecoregion>=2]
lenclass = length(s_inceptisol_ecoregion)
inceptisol_mean = matrix(0, nrow=200, ncol=lenclass)
inceptisol_sd   = matrix(0, nrow=200, ncol=lenclass)
inceptisol_se   = matrix(0, nrow=200, ncol=lenclass)
for(k in 1:lenclass){
   ecoid = as.numeric(names(s_inceptisol_ecoregion)[k])
   if(s_inceptisol_ecoregion[k] >= 2){
      inceptisol_mean[,k] = rowMeans(c_inceptisol[,c_inceptisol_ecoregion == ecoid], na.rm=TRUE)
      inceptisol_sd[,k] = apply(c_inceptisol[,c_inceptisol_ecoregion == ecoid], 1, sd, na.rm=TRUE)
      inceptisol_se[,k] = inceptisol_sd[,k] / rowSums(!is.na(c_inceptisol[,c_inceptisol_ecoregion == ecoid]))
   }
}
# Spodosol
# Overall mean profile
c_spodosol = extract(m ,info, 'Spodosol', 'Soil.Order')
c_spodosol_mean = rowMeans(c_spodosol, na.rm=TRUE)
c_spodosol_sd = apply(c_spodosol, 1, sd, na.rm=TRUE)
c_spodosol_se = c_spodosol_sd / rowSums(!is.na(c_spodosol))
# calculate mean and std for each veg type
c_spodosol_id = getid(m, info, 'Spodosol', 'Soil.Order')
c_spodosol_ecoregion = info$Landcover[c_spodosol_id]
s_ecoregion <- summary(as.factor(c_spodosol_ecoregion))   # Get the occurence frequency by different ecoregion
s_spodosol_ecoregion <- s_ecoregion[s_ecoregion>=2]
lenclass = length(s_spodosol_ecoregion)
spodosol_mean = matrix(0, nrow=200, ncol=lenclass)
spodosol_sd   = matrix(0, nrow=200, ncol=lenclass)
spodosol_se   = matrix(0, nrow=200, ncol=lenclass)
for(k in 1:lenclass){
   ecoid = as.numeric(names(s_spodosol_ecoregion)[k])
   if(s_spodosol_ecoregion[k] >= 2){
      spodosol_mean[,k] = rowMeans(c_spodosol[,c_spodosol_ecoregion == ecoid], na.rm=TRUE)
      spodosol_sd[,k] = apply(c_spodosol[,c_spodosol_ecoregion == ecoid], 1, sd, na.rm=TRUE)
      spodosol_se[,k] = spodosol_sd[,k] / rowSums(!is.na(c_spodosol[,c_spodosol_ecoregion == ecoid]))
   }
}
# Gelisol is seperated into three suborder
# Histel
# Overall mean profile
c_gelisol_histel = extract(m ,info, 'Histel', 'Suborder')
c_gelisol_histel_mean = rowMeans(c_gelisol_histel, na.rm=TRUE)
c_gelisol_histel_sd = apply(c_gelisol_histel, 1, sd, na.rm=TRUE)
c_gelisol_histel_se = c_gelisol_histel_sd / rowSums(!is.na(c_gelisol_histel))
# calculate mean and std for each veg type
c_gelisol_histel_id = getid(m, info, 'Histel', 'Suborder')
c_gelisol_histel_ecoregion = info$Landcover[c_gelisol_histel_id]
s_ecoregion <- summary(as.factor(c_gelisol_histel_ecoregion))   # Get the occurence frequency by different ecoregion
s_gelisol_histel_ecoregion <- s_ecoregion[s_ecoregion>=2]
lenclass = length(s_gelisol_histel_ecoregion)
gelisol_histel_mean = matrix(0, nrow=200, ncol=lenclass)
gelisol_histel_sd   = matrix(0, nrow=200, ncol=lenclass)
gelisol_histel_se   = matrix(0, nrow=200, ncol=lenclass)
for(k in 1:lenclass){
   ecoid = as.numeric(names(s_gelisol_histel_ecoregion)[k])
   if(s_gelisol_histel_ecoregion[k] >= 2){
      gelisol_histel_mean[,k] = rowMeans(c_gelisol_histel[,c_gelisol_histel_ecoregion == ecoid], na.rm=TRUE)
      gelisol_histel_sd[,k] = apply(c_gelisol_histel[,c_gelisol_histel_ecoregion == ecoid], 1, sd, na.rm=TRUE)
      gelisol_histel_se[,k] = gelisol_histel_sd[,k] / rowSums(!is.na(c_gelisol_histel[,c_gelisol_histel_ecoregion == ecoid]))
   }
}
# Orthel
# Overall mean profile
c_gelisol_orthel = extract(m ,info, 'Orthel', 'Suborder')
c_gelisol_orthel_mean = rowMeans(c_gelisol_orthel, na.rm=TRUE)
c_gelisol_orthel_sd = apply(c_gelisol_orthel, 1, sd, na.rm=TRUE)
c_gelisol_orthel_se = c_gelisol_orthel_sd / rowSums(!is.na(c_gelisol_orthel))
# calculate mean and std for each veg type
c_gelisol_orthel_id = getid(m, info, 'Orthel', 'Suborder')
c_gelisol_orthel_ecoregion = info$Landcover[c_gelisol_orthel_id]
s_ecoregion <- summary(as.factor(c_gelisol_orthel_ecoregion))   # Get the occurence frequency by different ecoregion
s_gelisol_orthel_ecoregion <- s_ecoregion[s_ecoregion>=2]
lenclass = length(s_gelisol_orthel_ecoregion)
gelisol_orthel_mean = matrix(0, nrow=200, ncol=lenclass)
gelisol_orthel_sd   = matrix(0, nrow=200, ncol=lenclass)
gelisol_orthel_se   = matrix(0, nrow=200, ncol=lenclass)
for(k in 1:lenclass){
   ecoid = as.numeric(names(s_gelisol_orthel_ecoregion)[k])
   if(s_gelisol_orthel_ecoregion[k] >= 2){
      gelisol_orthel_mean[,k] = rowMeans(c_gelisol_orthel[,c_gelisol_orthel_ecoregion == ecoid], na.rm=TRUE)
      gelisol_orthel_sd[,k] = apply(c_gelisol_orthel[,c_gelisol_orthel_ecoregion == ecoid], 1, sd, na.rm=TRUE)
      gelisol_orthel_se[,k] = gelisol_orthel_sd[,k] / rowSums(!is.na(c_gelisol_orthel[,c_gelisol_orthel_ecoregion == ecoid]))
   }
}
# Turbel
# Overall mean profile
c_gelisol_turbel = extract(m ,info, 'Turbel', 'Suborder')
c_gelisol_turbel_mean = rowMeans(c_gelisol_turbel, na.rm=TRUE)
c_gelisol_turbel_sd = apply(c_gelisol_turbel, 1, sd, na.rm=TRUE)
c_gelisol_turbel_se = c_gelisol_turbel_sd / rowSums(!is.na(c_gelisol_turbel))
# calculate mean and std for each veg type
c_gelisol_turbel_id = getid(m, info, 'Turbel', 'Suborder')
c_gelisol_turbel_ecoregion = info$Landcover[c_gelisol_turbel_id]
s_ecoregion <- summary(as.factor(c_gelisol_turbel_ecoregion))   # Get the occurence frequency by different ecoregion
s_gelisol_turbel_ecoregion <- s_ecoregion[s_ecoregion>=2]
lenclass = length(s_gelisol_turbel_ecoregion)
gelisol_turbel_mean = matrix(0, nrow=200, ncol=lenclass)
gelisol_turbel_sd   = matrix(0, nrow=200, ncol=lenclass)
gelisol_turbel_se   = matrix(0, nrow=200, ncol=lenclass)
for(k in 1:lenclass){
   ecoid = as.numeric(names(s_gelisol_turbel_ecoregion)[k])
   if(s_gelisol_turbel_ecoregion[k] >= 2){
      gelisol_turbel_mean[,k] = rowMeans(c_gelisol_turbel[,c_gelisol_turbel_ecoregion == ecoid], na.rm=TRUE)
      gelisol_turbel_sd[,k] = apply(c_gelisol_turbel[,c_gelisol_turbel_ecoregion == ecoid], 1, sd, na.rm=TRUE)
      gelisol_turbel_se[,k] = gelisol_turbel_sd[,k] / rowSums(!is.na(c_gelisol_turbel[,c_gelisol_turbel_ecoregion == ecoid]))
   }
}

## Raw profile data
# 1. Bering Taiga
depth_1 <- mishradata$nodedepth[mishradata$idx[mishradata$Ecoregions == "1"]]
#c_data_1 <- mishradata$C.Density[mishradata$idx[mishradata$Ecoregions == "1"]]
c_data_1 <- mishradata$Ccontent[mishradata$idx[mishradata$Ecoregions == "1"]]
# 2. Aleutian Meadows
#depth_2 <- mishradata$nodedepth[mishradata$idx[mishradata$Ecoregions == "2"]]
##c_data_2 <- mishradata$C.Density[mishradata$idx[mishradata$Ecoregions == "2"]]
#c_data_2 <- mishradata$Ccontent[mishradata$idx[mishradata$Ecoregions == "2"]]
# 3. Alaska Range Transition
depth_3 <- mishradata$nodedepth[mishradata$idx[mishradata$Ecoregions == "3"]]
#c_data_3 <- mishradata$C.Density[mishradata$idx[mishradata$Ecoregions == "3"]]
c_data_3 <- mishradata$Ccontent[mishradata$idx[mishradata$Ecoregions == "3"]]
# 4. Coastal Rainforest
depth_4 <- mishradata$nodedepth[mishradata$idx[mishradata$Ecoregions == "4"]]
#c_data_4 <- mishradata$C.Density[mishradata$idx[mishradata$Ecoregions == "4"]]
c_data_4 <- mishradata$Ccontent[mishradata$idx[mishradata$Ecoregions == "4"]]
# 5. Arctic Tundra
depth_5 <- mishradata$nodedepth[mishradata$idx[mishradata$Ecoregions == "5"]]
#c_data_5 <- mishradata$C.Density[mishradata$idx[mishradata$Ecoregions == "5"]]
c_data_5 <- mishradata$Ccontent[mishradata$idx[mishradata$Ecoregions == "5"]]
# 6. Bering Tundra
depth_6 <- mishradata$nodedepth[mishradata$idx[mishradata$Ecoregions == "6"]]
#c_data_6 <- mishradata$C.Density[mishradata$idx[mishradata$Ecoregions == "6"]]
c_data_6 <- mishradata$Ccontent[mishradata$idx[mishradata$Ecoregions == "6"]]
# 7. Pacific Mountain Transition
depth_7 <- mishradata$nodedepth[mishradata$idx[mishradata$Ecoregions == "7"]]
#c_data_7 <- mishradata$C.Density[mishradata$idx[mishradata$Ecoregions == "7"]]
c_data_7 <- mishradata$Ccontent[mishradata$idx[mishradata$Ecoregions == "7"]]
# 8. Intermontane Boreal
depth_8 <- mishradata$nodedepth[mishradata$idx[mishradata$Ecoregions == "8"]]
#c_data_8 <- mishradata$C.Density[mishradata$idx[mishradata$Ecoregions == "8"]]
c_data_8 <- mishradata$Ccontent[mishradata$idx[mishradata$Ecoregions == "8"]]
# 9. Coast Mountain Transition
#depth_9 <- mishradata$nodedepth[mishradata$idx[mishradata$Ecoregions == "9"]]
##c_data_9 <- mishradata$C.Density[mishradata$idx[mishradata$Ecoregions == "9"]]
#c_data_9 <- mishradata$Ccontent[mishradata$idx[mishradata$Ecoregions == "9"]]

## Plot profiles
# Assign depths
depth <- seq(1,200,1)

# 1. Alfisol
## Open PDF file
#if(plot_ensemble){
#   pdf(file="Alfisol_ensemble.pdf", width=17, height=10)
#} else{
#   pdf(file="Alfisol_ecoregion.pdf", width=17, height=10)
#}
## Set margin
#par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
#if(plot_ensemble){
#   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_1)[2]), ncol=dim(c_1)[2], byrow=TRUE)
#   matplot(c_1, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
#} else{
#   depth_ens = matrix(rep(seq(1,200,1),each=dim(alfisol_mean)[2]), ncol=dim(alfisol_mean)[2], byrow=TRUE)
#   matplot(alfisol_mean, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
##   lines(c_1_mean, depth, lwd = 4, xaxt='n', yaxt='n')
##   lines(c_1_mean+c_1_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
##   lines(c_1_mean-c_1_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
#}
## Move axis to top and left
#axis(3,at=seq(0.0,200.0,by=20),cex.axis=3)
#axis(3,at=seq(0.0,200.0,by=20),labels=FALSE,tcl=-0.25)
#axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
#axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
#grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
#legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
#legend(100, 150, legend = c("Fitted profile", "Fitted profile +- 1 STD"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
#     #  pch = "*", ncol = 4, cex = 0.8)
## Put Label for X and Y axis
#mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
#mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
#text(140, 120, "Alfisol", cex = 3)
#dev.off()

# 2. Andisol
# Open PDF file
if(plot_ensemble){
   pdf(file="Andisol_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Andisol_lc.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_2)[2]), ncol=dim(c_2)[2], byrow=TRUE)
   matplot(c_2, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   depth_ens = matrix(rep(seq(1,200,1),each=dim(andisol_mean)[2]), ncol=dim(andisol_mean)[2], byrow=TRUE)
   matplot(andisol_mean, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
  # lines(c_2_mean, depth, lwd = 4, xaxt='n', yaxt='n')
  # lines(c_2_mean+c_2_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
  # lines(c_2_mean-c_2_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,40.0,by=5),cex.axis=3)
axis(3,at=seq(0.0,40.0,by=5),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
legend("bottomright", legend = title_lc[code_lc %in% as.numeric(names(s_andisol_ecoregion))], lty = 1, cex = 3, col = 1:length(s_andisol_ecoregion)) #cex=3, pch = 17, bty = "n")
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
text(98, 120, "Andisol", cex = 3)
dev.off()

# 3. Entisol
# Open PDF file
if(plot_ensemble){
   pdf(file="Entisol_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Entisol_lc.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_3)[2]), ncol=dim(c_3)[2], byrow=TRUE)
   matplot(c_3, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   depth_ens = matrix(rep(seq(1,200,1),each=dim(entisol_mean)[2]), ncol=dim(entisol_mean)[2], byrow=TRUE)
   matplot(entisol_mean, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
  # lines(c_3_mean, depth, lwd = 4, xaxt='n', yaxt='n')
  # lines(c_3_mean+c_3_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
  # lines(c_3_mean-c_3_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,40.0,by=5),cex.axis=3)
axis(3,at=seq(0.0,40.0,by=5),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
legend("bottomright", legend = title_lc[code_lc %in% as.numeric(names(s_entisol_ecoregion))], lty = 1, cex = 3, col = 1:length(s_entisol_ecoregion)) #cex=3, pch = 17, bty = "n")
#legend(100, 150, legend = c("Fitted profile", "Fitted profile +- 1 STD"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
text(140, 120, "Entisol", cex = 3)
dev.off()

# 4. Histosol
# Open PDF file
#if(plot_ensemble){
#   pdf(file="Histosol_ensemble.pdf", width=17, height=10)
#} else{
#   pdf(file="Histosol_lc.pdf", width=17, height=10)
#}
## Set margin
#par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
#if(plot_ensemble){
#   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_4)[2]), ncol=dim(c_4)[2], byrow=TRUE)
#   matplot(c_4, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
#} else{
#   depth_ens = matrix(rep(seq(1,200,1),each=dim(histosol_mean)[2]), ncol=dim(histosol_mean)[2], byrow=TRUE)
#   matplot(histosol_mean, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
#  # lines(c_4_mean, depth, lwd = 4, xaxt='n', yaxt='n')
#  # lines(c_4_mean+c_4_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
#  # lines(c_4_mean-c_4_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
#}
## Move axis to top and left
#axis(3,at=seq(0.0,250.0,by=20),cex.axis=3)
#axis(3,at=seq(0.0,250.0,by=20),labels=FALSE,tcl=-0.25)
#axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
#axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
#grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
#legend("bottomright", legend = title_lc[code_lc %in% as.numeric(names(s_histosol_ecoregion))], lty = 1, cex = 3, col = 1:length(s_histosol_ecoregion)) #cex=3, pch = 17, bty = "n")
#     #  pch = "*", ncol = 4, cex = 0.8)
## Put Label for X and Y axis
#mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
#mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
#text(155, 120, "Histosol", cex = 3)
#dev.off()

# 5. Inceptisol
# Open PDF file
if(plot_ensemble){
   pdf(file="Inceptisol_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Inceptisol_lc.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_5)[2]), ncol=dim(c_5)[2], byrow=TRUE)
   matplot(c_5, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   depth_ens = matrix(rep(seq(1,200,1),each=dim(inceptisol_mean)[2]), ncol=dim(inceptisol_mean)[2], byrow=TRUE)
   matplot(inceptisol_mean, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
  # lines(c_5_mean, depth, lwd = 4, xaxt='n', yaxt='n')
  # lines(c_5_mean+c_5_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
  # lines(c_5_mean-c_5_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,40.0,by=5),cex.axis=3)
axis(3,at=seq(0.0,40.0,by=5),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
legend("bottomright", legend = title_lc[code_lc %in% as.numeric(names(s_inceptisol_ecoregion))], lty = 1, cex = 3, col = 1:length(s_inceptisol_ecoregion)) #cex=3, pch = 17, bty = "n")
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
text(140, 120, "Inceptisol", cex = 3)
dev.off()

# 6. Spodosol
# Open PDF file
if(plot_ensemble){
   pdf(file="Spodosol_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Spodosol_lc.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_6)[2]), ncol=dim(c_6)[2], byrow=TRUE)
   matplot(c_6, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   depth_ens = matrix(rep(seq(1,200,1),each=dim(spodosol_mean)[2]), ncol=dim(spodosol_mean)[2], byrow=TRUE)
   matplot(spodosol_mean, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
  # lines(c_6_mean, depth, lwd = 4, xaxt='n', yaxt='n')
  # lines(c_6_mean+c_6_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
  # lines(c_6_mean-c_6_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,50.0,by=5),cex.axis=3)
axis(3,at=seq(0.0,50.0,by=5),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
legend("bottomright", legend = title_lc[code_lc %in% as.numeric(names(s_spodosol_ecoregion))], lty = 1, cex = 3, col = 1:length(s_spodosol_ecoregion)) #cex=3, pch = 17, bty = "n")
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
text(140, 120, "Spodosol", cex = 3)
dev.off()

# 7. Histel
# Open PDF file
if(plot_ensemble){
   pdf(file="Histel_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Histel_lc.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_7)[2]), ncol=dim(c_7)[2], byrow=TRUE)
   matplot(c_7, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   depth_ens = matrix(rep(seq(1,200,1),each=dim(gelisol_histel_mean)[2]), ncol=dim(gelisol_histel_mean)[2], byrow=TRUE)
   matplot(gelisol_histel_mean, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
  # lines(c_7_mean, depth, lwd = 4, xaxt='n', yaxt='n')
  # lines(c_7_mean+c_7_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
  # lines(c_7_mean-c_7_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,70.0,by=10),cex.axis=3)
axis(3,at=seq(0.0,70.0,by=10),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
legend("bottomright", legend = title_lc[code_lc %in% as.numeric(names(s_gelisol_histel_ecoregion))], lty = 1, cex = 3, col = 1:length(s_gelisol_histel_ecoregion)) #cex=3, pch = 17, bty = "n")
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
text(140, 120, "Gelisol Histel", cex = 3)
dev.off()

# 8. Orthel
# Open PDF file
if(plot_ensemble){
   pdf(file="Orthel_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Orthel_lc.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_8)[2]), ncol=dim(c_8)[2], byrow=TRUE)
   matplot(c_8, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   depth_ens = matrix(rep(seq(1,200,1),each=dim(gelisol_orthel_mean)[2]), ncol=dim(gelisol_orthel_mean)[2], byrow=TRUE)
   matplot(gelisol_orthel_mean, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
  # lines(c_8_mean, depth, lwd = 4, xaxt='n', yaxt='n')
  # lines(c_8_mean+c_8_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
  # lines(c_8_mean-c_8_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,100.0,by=10),cex.axis=3)
axis(3,at=seq(0.0,100.0,by=10),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
legend("bottomright", legend = title_lc[code_lc %in% as.numeric(names(s_gelisol_orthel_ecoregion))], lty = 1, cex = 3, col = 1:length(s_gelisol_orthel_ecoregion)) #cex=3, pch = 17, bty = "n")
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
text(140, 120, "Gelisol Orthel", cex = 3)
dev.off()

# 9. Turbel
# Open PDF file
if(plot_ensemble){
   pdf(file="Turbel_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Turbel_lc.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_9)[2]), ncol=dim(c_9)[2], byrow=TRUE)
   matplot(c_9, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   depth_ens = matrix(rep(seq(1,200,1),each=dim(gelisol_turbel_mean)[2]), ncol=dim(gelisol_turbel_mean)[2], byrow=TRUE)
   matplot(gelisol_turbel_mean, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
  # lines(c_9_mean, depth, lwd = 4, xaxt='n', yaxt='n')
  # lines(c_9_mean+c_9_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
  # lines(c_9_mean-c_9_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,60.0,by=10),cex.axis=3)
axis(3,at=seq(0.0,60.0,by=10),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
legend("bottomright", legend = title_lc[code_lc %in% as.numeric(names(s_gelisol_turbel_ecoregion))], lty = 1, cex = 3, col = 1:length(s_gelisol_turbel_ecoregion)) #cex=3, pch = 17, bty = "n")
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
text(140, 120, "Gelisol Turbel", cex = 3)
dev.off()

