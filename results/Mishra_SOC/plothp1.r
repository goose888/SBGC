# R plot using basic functions to generate plots
# Author: Shijie Shu

# Define reset function here
reset <- function() {
    par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
    plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
    }

# Load Library
library(ggplot2)

# User inputs
proj_1cycle_fname  = "proj_SOC_1cycle.csv"
proj_10cycle_fname = "proj_SOC_10cycle.csv"
#plot_path = ""

# Read in data frame
proj1_data  = read.csv(proj_1cycle_fname)
proj10_data = read.csv(proj_10cycle_fname)

# Separate data column to different variables
depth_model = proj1_data[,1]
tundra_contemp_1cycle = proj1_data[,3]
boreal_contemp_1cycle = proj1_data[,4]
tundra_forced_1cycle = proj1_data[,6]
boreal_forced_1cycle = proj1_data[,7]
tundra_free_1cycle = proj1_data[,9]
boreal_free_1cycle = proj1_data[,10]

tundra_contemp_10cycle = proj10_data[,3]
boreal_contemp_10cycle = proj10_data[,4]
tundra_forced_10cycle = proj10_data[,6]
boreal_forced_10cycle = proj10_data[,7]
tundra_free_10cycle = proj10_data[,9]
boreal_free_10cycle = proj10_data[,10]

# Plot figures: 4x3 curve charts
# Open PDF file
pdf(file="hp1_fig.pdf", width=20, height=12)
par(mfrow=c(1,2), omi=c(1.2,1.5,1.5,0.5), mar=c(1,1,1,1)) # all plots on one page
#par(mfrow=c(1,2), oma=c(4,1,1,1), mar=c(1,1,1,1)) # all plots on one page

## subplot one
# Set margin
plot(tundra_contemp_1cycle, depth_model, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
                       xlim=c(40,140), ylim=rev(c(0,160)), pch=17, cex= 0, xaxt='n', yaxt='n')

# Line plots
lines(tundra_contemp_1cycle, depth_model, lty = 3, lwd = 8, pch=17, cex=1, col="red", xaxt='n', yaxt='n')    # Ecoregion 5= Arctic Tundra
lines(boreal_contemp_1cycle, depth_model, lty = 3, lwd = 8, pch=17, cex=1, col="green", xaxt='n', yaxt='n')    # Ecoregion 8 =Intermontane Boreal
lines(tundra_forced_1cycle, depth_model, lty = 2, lwd = 4, pch=14, cex=1, col="red", xaxt='n', yaxt='n')    # Ecoregion 5= Arctic Tundra
lines(boreal_forced_1cycle, depth_model, lty = 2, lwd = 4, pch=14, cex=1, col="green", xaxt='n', yaxt='n')    # Ecoregion 8 =Intermontane Boreal
lines(tundra_free_1cycle, depth_model, lty = 1, lwd = 4, pch=11, cex=1, col="red", xaxt='n', yaxt='n')    # Ecoregion 5= Arctic Tundra
lines(boreal_free_1cycle, depth_model, lty = 1, lwd = 4, pch=11, cex=1, col="green", xaxt='n', yaxt='n')    # Ecoregion 8 =Intermontane Boreal

# Move axis to top and left
axis(3,at=seq(40.0,140.0,by=20),cex.axis=3)
axis(3,at=seq(40.0,140.0,by=20),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
text(45, 5, "(a)", cex=3)

# legend(80, 90, legend = c("Initial, Tundra", "Initial, Boreal"), cex=2.5, lty=3, col=c("red","green"), lwd=4, bty = "n")#, lty = 1:7,
# legend(80, 110, legend = c("Forced, Tundra", "Forced, Boreal"), cex=2.5, lty=2, col=c("red","green"), lwd=4, bty = "n")#, lty = 1:7,
# legend(80, 130, legend = c("Coupled, Tundra", "Coupled, Boreal"), cex=2.5, lty=1, col=c("red","green"), lwd=4, bty = "n")#, lty = 1:7,


#  legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
#   legend(60, 80, legend = c("Obs - Arctic Tundra", "Obs - Bering Tundra", "Obs - Intermontane Boreal"), cex=3, lty=1, col=c("red","blue","green"), lwd=3, bty = "n")#, lty = 1:7,
#   legend(60, 160, legend = c("Model - Arctic Tundra", "Model - Bering Tundra", "Model - Intermontane Boreal"), cex=3, lty=6, col=c("red","blue","green"), lwd=3, bty = "n")#, lty = 1:7,
#  pch = "*", ncol = 4, cex = 0.8)

## subplot two
par(mar=c(1,1,1,1))
plot(tundra_contemp_10cycle, depth_model, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
                       xlim=c(40,140), ylim=rev(c(0,160)), pch=17, cex= 0, xaxt='n', yaxt='n')
# Line plots
lines(tundra_contemp_10cycle, depth_model, lty = 3, lwd = 8, pch=17, cex=1, col="red", xaxt='n', yaxt='n')    # Ecoregion 5= Arctic Tundra
lines(boreal_contemp_10cycle, depth_model, lty = 3, lwd = 8, pch=17, cex=1, col="green", xaxt='n', yaxt='n')    # Ecoregion 8 =Intermontane Boreal
lines(tundra_forced_10cycle, depth_model, lty = 2, lwd = 4, pch=14, cex=1, col="red", xaxt='n', yaxt='n')    # Ecoregion 5= Arctic Tundra
lines(boreal_forced_10cycle, depth_model, lty = 2, lwd = 4, pch=14, cex=1, col="green", xaxt='n', yaxt='n')    # Ecoregion 8 =Intermontane Boreal
lines(tundra_free_10cycle, depth_model, lty = 1, lwd = 4, pch=11, cex=1, col="red", xaxt='n', yaxt='n')    # Ecoregion 5= Arctic Tundra
lines(boreal_free_10cycle, depth_model, lty = 1, lwd = 4, pch=11, cex=1, col="green", xaxt='n', yaxt='n')    # Ecoregion 8 =Intermontane Boreal

# Move axis to top and left
axis(3,at=seq(40.0,140.0,by=20),cex.axis=3)
axis(3,at=seq(40.0,140.0,by=20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)

#  legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
#   legend(60, 80, legend = c("Obs - Arctic Tundra", "Obs - Bering Tundra", "Obs - Intermontane Boreal"), cex=3, lty=1, col=c("red","blue","green"), lwd=3, bty = "n")#, lty = 1:7,
#   legend(60, 160, legend = c("Model - Arctic Tundra", "Model - Bering Tundra", "Model - Intermontane Boreal"), cex=3, lty=6, col=c("red","blue","green"), lwd=3, bty = "n")#, lty = 1:7,
#  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
text(45, 5, "(b)", cex=3)
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=4.7, cex=3)
mtext("Organic Carbon Density (kgC/m3)", side=3, outer=TRUE,line=3.5,cex=3)

#legend(40, 160, inset=c(-0.2,0), legend = c("Initial, Tundra", "Initial, Boreal"), cex=2.5, lty=3, col=c("red","green"), lwd=4, bty = "n")#, lty = 1:7,
#legend(-100, 180, inset=c(-0.2,0), legend = c("Forced, Tundra", "Forced, Boreal"), cex=2.5, lty=2, col=c("red","green"), lwd=4, bty = "n")#, lty = 1:7,
#legend(-100, 180, inset=c(-0.2,0), legend = c("Coupled, Tundra", "Coupled, Boreal"), cex=2.5, lty=1, col=c("red","green"), lwd=4, bty = "n")#, lty = 1:7,

dev.off()

