# Mean profiles classified by Landcover type
# must be called by MishraSOCspline.r 
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

# 11. Open Water
c_11 = extract(m ,info, '11', 'Landcover')
c_11_mean = rowMeans(c_11, na.rm=TRUE)
c_11_sd = apply(c_11, 1, sd, na.rm=TRUE)
c_11_se = c_11_sd / rowSums(!is.na(c_11))
# 12. Perrenial Ice
c_12 = extract(m ,info, '12', 'Landcover')
c_12_mean = rowMeans(c_12, na.rm=TRUE)
c_12_sd = apply(c_12, 1, sd, na.rm=TRUE)
c_12_se = c_12_sd / rowSums(!is.na(c_12))
# 21. Developed
c_21 = extract(m ,info, '21', 'Landcover')
c_21_mean = rowMeans(c_21, na.rm=TRUE)
c_21_sd = apply(c_21, 1, sd, na.rm=TRUE)
c_21_se = c_21_sd / rowSums(!is.na(c_21))
# 31. Barren Land
c_31 = extract(m ,info, '31', 'Landcover')
c_31_mean = rowMeans(c_31, na.rm=TRUE)
c_31_sd = apply(c_31, 1, sd, na.rm=TRUE)
c_31_se = c_31_sd / rowSums(!is.na(c_31))
# 41. Forrest
c_41 = extract(m ,info, '41', 'Landcover')
c_41_mean = rowMeans(c_41, na.rm=TRUE)
c_41_sd = apply(c_41, 1, sd, na.rm=TRUE)
c_41_se = c_41_sd / rowSums(!is.na(c_41))
# 51. Scrub
c_51 = extract(m ,info, '51', 'Landcover')
c_51_mean = rowMeans(c_51, na.rm=TRUE)
c_51_sd = apply(c_51, 1, sd, na.rm=TRUE)
c_51_se = c_51_sd / rowSums(!is.na(c_51))
# 71. Herbaceous Vegetation
c_71 = extract(m ,info, '71', 'Landcover')
c_71_mean = rowMeans(c_71, na.rm=TRUE)
c_71_sd = apply(c_71, 1, sd, na.rm=TRUE)
c_71_se = c_71_sd / rowSums(!is.na(c_71))
# 81. Pasture
c_81 = extract(m ,info, '81', 'Landcover')
c_81_mean = rowMeans(c_81, na.rm=TRUE)
c_81_sd = apply(c_81, 1, sd, na.rm=TRUE)
c_81_se = c_81_sd / rowSums(!is.na(c_81))
# 90. Wetland
c_90 = extract(m ,info, '90', 'Landcover')
c_90_mean = rowMeans(c_90, na.rm=TRUE)
c_90_sd = apply(c_90, 1, sd, na.rm=TRUE)
c_90_se = c_90_sd / rowSums(!is.na(c_90))

## Raw profile data
#
# 11. Open Water
depth_11 <- mishradata$nodedepth[mishradata$idx[mishradata$Landcover == "11"]]
#c_data_11 <- mishradata$C.Density[mishradata$idx[mishradata$Landcover == "11"]]
c_data_11 <- mishradata$Ccontent[mishradata$idx[mishradata$Landcover == "11"]]
# 12. Perrenial Ice
depth_12 <- mishradata$nodedepth[mishradata$idx[mishradata$Landcover == "12"]]
#c_data_12 <- mishradata$C.Density[mishradata$idx[mishradata$Landcover == "12"]]
c_data_12 <- mishradata$Ccontent[mishradata$idx[mishradata$Landcover == "12"]]
# 21. Developed
depth_21 <- mishradata$nodedepth[mishradata$idx[mishradata$Landcover == "21"]]
#c_data_21 <- mishradata$C.Density[mishradata$idx[mishradata$Landcover == "21"]]
c_data_21 <- mishradata$Ccontent[mishradata$idx[mishradata$Landcover == "21"]]
# 31. Barren Land
depth_31 <- mishradata$nodedepth[mishradata$idx[mishradata$Landcover == "31"]]
#c_data_31 <- mishradata$C.Density[mishradata$idx[mishradata$Landcover == "31"]]
c_data_31 <- mishradata$Ccontent[mishradata$idx[mishradata$Landcover == "31"]]
# 41. Forrest
depth_41 <- mishradata$nodedepth[mishradata$idx[mishradata$Landcover == "41"]]
#c_data_41 <- mishradata$C.Density[mishradata$idx[mishradata$Landcover == "41"]]
c_data_41 <- mishradata$Ccontent[mishradata$idx[mishradata$Landcover == "41"]]
# 51. Scrub
depth_51 <- mishradata$nodedepth[mishradata$idx[mishradata$Landcover == "51"]]
#c_data_51 <- mishradata$C.Density[mishradata$idx[mishradata$Landcover == "51"]]
c_data_51 <- mishradata$Ccontent[mishradata$idx[mishradata$Landcover == "51"]]
# 71. Herbaceous vegetation
depth_71 <- mishradata$nodedepth[mishradata$idx[mishradata$Landcover == "71"]]
#c_data_71 <- mishradata$C.Density[mishradata$idx[mishradata$Landcover == "71"]]
c_data_71 <- mishradata$Ccontent[mishradata$idx[mishradata$Landcover == "71"]]
# 81. Pasture
depth_81 <- mishradata$nodedepth[mishradata$idx[mishradata$Landcover == "81"]]
#c_data_81 <- mishradata$C.Density[mishradata$idx[mishradata$Landcover == "81"]]
c_data_81 <- mishradata$Ccontent[mishradata$idx[mishradata$Landcover == "81"]]
# 90. Wetland
depth_90 <- mishradata$nodedepth[mishradata$idx[mishradata$Landcover == "90"]]
#c_data_90 <- mishradata$C.Density[mishradata$idx[mishradata$Landcover == "90"]]
c_data_90 <- mishradata$Ccontent[mishradata$idx[mishradata$Landcover == "90"]]

## Plot profiles
# Assign depths
depth <- seq(1,200,1)

# Open Water
# Open PDF file
if(plot_ensemble){
   pdf(file="Openwater_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Openwater.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_11)[2]), ncol=dim(c_11)[2], byrow=TRUE)
   matplot(c_11, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_11, depth_11, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
   lines(c_11_mean, depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_11_mean+c_11_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_11_mean-c_11_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,200.0,by=20),cex.axis=3)
axis(3,at=seq(0.0,200.0,by=20),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
legend(100, 150, legend = c("Fitted profile", "Fitted profile +- 1 STD"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
text(140, 120, "Open Water", cex = 3)
dev.off()

# Perrenial Ice
# Open PDF file
if(plot_ensemble){
   pdf(file="Perrenialice_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Perrenialice.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_12)[2]), ncol=dim(c_12)[2], byrow=TRUE)
   matplot(c_12, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_12, depth_12, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
   lines(c_12_mean, depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_12_mean+c_12_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_12_mean-c_12_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,40.0,by=5),cex.axis=3)
axis(3,at=seq(0.0,40.0,by=5),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
legend(92, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
legend(80, 150, legend = c("Fitted profile", "Fitted profile +- 1 STD"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
text(98, 120, "Perrenial Ice", cex = 3)
dev.off()

# Developed
# Open PDF file
if(plot_ensemble){
   pdf(file="Developed_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Developed.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_21)[2]), ncol=dim(c_21)[2], byrow=TRUE)
   matplot(c_21, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_21, depth_21, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
   lines(c_21_mean, depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_21_mean+c_21_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_21_mean-c_21_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,40.0,by=5),cex.axis=3)
axis(3,at=seq(0.0,40.0,by=5),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
legend(100, 150, legend = c("Fitted profile", "Fitted profile +- 1 STD"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
text(140, 120, "Developed", cex = 3)
dev.off()

# Barren land
# Open PDF file
if(plot_ensemble){
   pdf(file="Barrenland_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Barrenland.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_31)[2]), ncol=dim(c_31)[2], byrow=TRUE)
   matplot(c_31, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_31, depth_31, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
   lines(c_31_mean, depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_31_mean+c_31_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_31_mean-c_31_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,250.0,by=20),cex.axis=3)
axis(3,at=seq(0.0,250.0,by=20),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
legend(130, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
legend(120, 150, legend = c("Fitted profile", "Fitted profile +- 1 STD"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
text(155, 120, "Barren Land", cex = 3)
dev.off()

# Forest
# Open PDF file
if(plot_ensemble){
   pdf(file="Forest_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Forest.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_41)[2]), ncol=dim(c_41)[2], byrow=TRUE)
   matplot(c_41, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_41, depth_41, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
   lines(c_41_mean, depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_41_mean+c_41_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_41_mean-c_41_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,40.0,by=5),cex.axis=3)
axis(3,at=seq(0.0,40.0,by=5),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
legend(100, 150, legend = c("Fitted profile", "Fitted profile +- 1 STD"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
text(140, 120, "Forest", cex = 3)
dev.off()

# Scrub
# Open PDF file
if(plot_ensemble){
   pdf(file="Scrub_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Scrub.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_51)[2]), ncol=dim(c_51)[2], byrow=TRUE)
   matplot(c_51, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_51, depth_51, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
   lines(c_51_mean, depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_51_mean+c_51_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_51_mean-c_51_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,50.0,by=5),cex.axis=3)
axis(3,at=seq(0.0,50.0,by=5),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
legend(100, 150, legend = c("Fitted profile", "Fitted profile +- 1 STD"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
text(140, 120, "Scrub", cex = 3)
dev.off()

# Herbaceous
# Open PDF file
if(plot_ensemble){
   pdf(file="Herbaceous_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Herbaceous.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_71)[2]), ncol=dim(c_71)[2], byrow=TRUE)
   matplot(c_71, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_71, depth_71, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      xlim=c(0,70), ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
   lines(c_71_mean, depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_71_mean+c_71_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_71_mean-c_71_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,70.0,by=10),cex.axis=3)
axis(3,at=seq(0.0,70.0,by=10),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
legend(100, 150, legend = c("Fitted profile", "Fitted profile +- 1 STD"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
text(140, 120, "Herbaceous", cex = 3)
dev.off()

# Pasture
# Open PDF file
if(plot_ensemble){
   pdf(file="Pasture_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Pasture.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_81)[2]), ncol=dim(c_81)[2], byrow=TRUE)
   matplot(c_81, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_81, depth_81, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      xlim=c(0,100), ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
   lines(c_81_mean, depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_81_mean+c_81_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_81_mean-c_81_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,100.0,by=10),cex.axis=3)
axis(3,at=seq(0.0,100.0,by=10),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
legend(100, 150, legend = c("Fitted profile", "Fitted profile +- 1 STD"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
text(140, 120, "Pasture", cex = 3)
dev.off()

# Wetland
# Open PDF file
if(plot_ensemble){
   pdf(file="Wetland_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Wetland.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_90)[2]), ncol=dim(c_90)[2], byrow=TRUE)
   matplot(c_90, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_90, depth_90, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      xlim=c(0,60), ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
   lines(c_90_mean, depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_90_mean+c_90_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_90_mean-c_90_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,60.0,by=10),cex.axis=3)
axis(3,at=seq(0.0,60.0,by=10),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
legend(100, 150, legend = c("Fitted profile", "Fitted profile +- 1 STD"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
text(140, 120, "Wetland", cex = 3)
dev.off()


