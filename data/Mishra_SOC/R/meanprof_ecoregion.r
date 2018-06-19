# Mean profiles classified by Ecoregion
# must be called by MishraSOCspline.r 
# 9 different landcover types are involved
# 1 = Bering Taiga
# 2 = Aleutian Meadows
# 3 = Alaska Range Transition
# 4 = Coastal Rainforest
# 5 = Arctic Tundra
# 6 = Bering Tundra
# 7 = Pacific Mountain Transition
# 8 = Intermontane Boreal
# 9 = Coast Mountain Transition

# 1. Bering Taiga
c_1 = extract(m ,info, '1', 'Ecoregions')
c_1_mean = rowMeans(c_1, na.rm=TRUE)
c_1_sd = apply(c_1, 1, sd, na.rm=TRUE)
c_1_se = c_1_sd / rowSums(!is.na(c_1))
# 2. Aleutian Meadows
#c_2 = extract(m ,info, '2', 'Ecoregions')
#c_2_mean = rowMeans(c_2, na.rm=TRUE)
#c_2_sd = apply(c_2, 1, sd, na.rm=TRUE)
#c_2_se = c_2_sd / rowSums(!is.na(c_2))
# 3. Alaska Range Transition
c_3 = extract(m ,info, '3', 'Ecoregions')
c_3_mean = rowMeans(c_3, na.rm=TRUE)
c_3_sd = apply(c_3, 1, sd, na.rm=TRUE)
c_3_se = c_3_sd / rowSums(!is.na(c_3))
# 4. Coastal Rainforest
c_4 = extract(m ,info, '4', 'Ecoregions')
c_4_mean = rowMeans(c_4, na.rm=TRUE)
c_4_sd = apply(c_4, 1, sd, na.rm=TRUE)
c_4_se = c_4_sd / rowSums(!is.na(c_4))
# 5. Arctic Tundra
c_5 = extract(m ,info, '5', 'Ecoregions')
c_5_mean = rowMeans(c_5, na.rm=TRUE)
c_5_sd = apply(c_5, 1, sd, na.rm=TRUE)
c_5_se = c_5_sd / rowSums(!is.na(c_5))
# 6. Bering Tundra
c_6 = extract(m ,info, '6', 'Ecoregions')
c_6_mean = rowMeans(c_6, na.rm=TRUE)
c_6_sd = apply(c_6, 1, sd, na.rm=TRUE)
c_6_se = c_6_sd / rowSums(!is.na(c_6))
# 7. Pacific Mountain Transition
c_7 = extract(m ,info, '7', 'Ecoregions')
c_7_mean = rowMeans(c_7, na.rm=TRUE)
c_7_sd = apply(c_7, 1, sd, na.rm=TRUE)
c_7_se = c_7_sd / rowSums(!is.na(c_7))
# 8. Intermontane Boreal
c_8 = extract(m ,info, '8', 'Ecoregions')
c_8_mean = rowMeans(c_8, na.rm=TRUE)
c_8_sd = apply(c_8, 1, sd, na.rm=TRUE)
c_8_se = c_8_sd / rowSums(!is.na(c_8))
# 9. Coast Mountain Transition
#c_9 = extract(m ,info, '9', 'Ecoregions')
#c_9_mean = rowMeans(c_9, na.rm=TRUE)
#c_9_sd = apply(c_9, 1, sd, na.rm=TRUE)
#c_9_se = c_9_sd / rowSums(!is.na(c_9))

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

# 1. Bering Taiga
# Open PDF file
if(plot_ensemble){
   pdf(file="Beringtaiga_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Beringtaiga.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_1)[2]), ncol=dim(c_1)[2], byrow=TRUE)
   matplot(c_1, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_1, depth_1, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
   lines(c_1_mean, depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_1_mean+c_1_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_1_mean-c_1_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
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
text(140, 120, "Bering Taiga", cex = 3)
dev.off()

# 2. Aleutian Meadows
# Open PDF file
#if(plot_ensemble){
#   pdf(file="Aleutianmeadows_ensemble.pdf", width=17, height=10)
#} else{
#   pdf(file="Aleutianmeadows.pdf", width=17, height=10)
#}
## Set margin
#par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
#if(plot_ensemble){
#   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_2)[2]), ncol=dim(c_2)[2], byrow=TRUE)
#   matplot(c_2, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
#} else{
#   plot(c_data_2, depth_2, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
#      ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
#   lines(c_2_mean, depth, lwd = 4, xaxt='n', yaxt='n')
#   lines(c_2_mean+c_2_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
#   lines(c_2_mean-c_2_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
#}
## Move axis to top and left
#axis(3,at=seq(0.0,40.0,by=5),cex.axis=3)
#axis(3,at=seq(0.0,40.0,by=5),labels=FALSE,tcl=-0.25)
#axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
#axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
#grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
#legend(92, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
#legend(80, 150, legend = c("Fitted profile", "Fitted profile +- 1 STD"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
#     #  pch = "*", ncol = 4, cex = 0.8)
## Put Label for X and Y axis
#mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
#mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
#text(98, 120, "Aleutian Meadows", cex = 3)
#dev.off()

# 3. Alaska Range Transition
# Open PDF file
if(plot_ensemble){
   pdf(file="Alaskarangetransition_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Alaskarangetransition.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_3)[2]), ncol=dim(c_3)[2], byrow=TRUE)
   matplot(c_3, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_3, depth_3, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
   lines(c_3_mean, depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_3_mean+c_3_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_3_mean-c_3_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
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
text(140, 120, "Alaska Range Transition", cex = 3)
dev.off()

# 4. Coastal Rainforest
# Open PDF file
if(plot_ensemble){
   pdf(file="Coastalrainforest_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Coastalrainforest.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_4)[2]), ncol=dim(c_4)[2], byrow=TRUE)
   matplot(c_4, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_4, depth_4, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
   lines(c_4_mean, depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_4_mean+c_4_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_4_mean-c_4_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
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
text(155, 120, "Coastal Rainforest", cex = 3)
dev.off()

# 5. Arctic Tundra
# Open PDF file
if(plot_ensemble){
   pdf(file="Arctictundra_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Arctictundra.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_5)[2]), ncol=dim(c_5)[2], byrow=TRUE)
   matplot(c_5, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_5, depth_5, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
   lines(c_5_mean, depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_5_mean+c_5_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_5_mean-c_5_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
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
text(140, 120, "Arctic Tundra", cex = 3)
dev.off()

# 6. Bering Tundra
# Open PDF file
if(plot_ensemble){
   pdf(file="Beringtundra_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Beringtundra.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_6)[2]), ncol=dim(c_6)[2], byrow=TRUE)
   matplot(c_6, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_6, depth_6, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
   lines(c_6_mean, depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_6_mean+c_6_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_6_mean-c_6_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
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
text(140, 120, "Bering Tundra", cex = 3)
dev.off()

# 7. Pacific Mountain Transition
# Open PDF file
if(plot_ensemble){
   pdf(file="Pacificmountaintransition_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Pacificmountaintransition.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_7)[2]), ncol=dim(c_7)[2], byrow=TRUE)
   matplot(c_7, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_7, depth_7, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      xlim=c(0,70), ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
   lines(c_7_mean, depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_7_mean+c_7_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_7_mean-c_7_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
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
text(140, 120, "Pacific Mountain Transition", cex = 3)
dev.off()

# 8. Intermontane Boreal
# Open PDF file
if(plot_ensemble){
   pdf(file="Intermontaneboreal_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Intermontaneboreal.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_8)[2]), ncol=dim(c_8)[2], byrow=TRUE)
   matplot(c_8, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_8, depth_8, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      xlim=c(0,100), ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
   lines(c_8_mean, depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_8_mean+c_8_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_8_mean-c_8_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
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
text(140, 120, "Intermontane Boreal", cex = 3)
dev.off()

# 9. Coast Mountain Transition
# Open PDF file
#if(plot_ensemble){
#   pdf(file="Coastmountaintransition_ensemble.pdf", width=17, height=10)
#} else{
#   pdf(file="Coastmountaintransition.pdf", width=17, height=10)
#}
## Set margin
#par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
#if(plot_ensemble){
#   depth_ens = matrix(rep(seq(1,200,1),each=dim(c_9)[2]), ncol=dim(c_9)[2], byrow=TRUE)
#   matplot(c_9, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
#} else{
#   plot(c_data_9, depth_9, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
#      xlim=c(0,60), ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
#   lines(c_9_mean, depth, lwd = 4, xaxt='n', yaxt='n')
#   lines(c_9_mean+c_9_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
#   lines(c_9_mean-c_9_se, depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
#}
## Move axis to top and left
#axis(3,at=seq(0.0,60.0,by=10),cex.axis=3)
#axis(3,at=seq(0.0,60.0,by=10),labels=FALSE,tcl=-0.25)
#axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
#axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
#grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
#legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
#legend(100, 150, legend = c("Fitted profile", "Fitted profile +- 1 STD"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
#     #  pch = "*", ncol = 4, cex = 0.8)
## Put Label for X and Y axis
#mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
#mtext("Organic Carbon Density (kgC/m2)", side=3, outer=TRUE,line=-3.7,cex=3)
#text(140, 120, "Coast Mountain Transition", cex = 3)
#dev.off()

