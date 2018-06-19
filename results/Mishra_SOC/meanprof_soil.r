# Mean profiles classified by Soil Order
# must be called by MishraSOCspline.r
# There are totally six soil orders and three soil suborders for gelisol here

# User control options
res <- 5

# Alfisol
c_alfisol = extract(m ,info, 'Alfisol', 'Soil.Order')
c_alfisol_mean = rowMeans(c_alfisol, na.rm=TRUE)
c_alfisol_sd = apply(c_alfisol, 1, sd, na.rm=TRUE)
# Filter the outliers
for (i in 1:dim(c_alfisol)[2]) {
  for (j in 1:200) {
     if(c_alfisol[j,i] > c_alfisol_mean[j] + 2*c_alfisol_sd[j]) {
        c_alfisol[j,i] <- NA
     }
  }
}
c_alfisol_sd = apply(c_alfisol, 1, sd, na.rm=TRUE)
c_alfisol_se = c_alfisol_sd / rowSums(!is.na(c_alfisol))

# Andisol
c_andisol = extract(m ,info, 'Andisol', 'Soil.Order')
c_andisol_mean = rowMeans(c_andisol, na.rm=TRUE)
c_andisol_sd = apply(c_andisol, 1, sd, na.rm=TRUE)
# Filter the outliers
for (i in 1:dim(c_andisol)[2]) {
  for (j in 1:200) {
     if(is.na(c_andisol_sd[j]) | is.na(c_andisol[j,i])){
        # Do nothing
     } else{
        if(c_andisol[j,i] > c_andisol_mean[j] + 2*c_andisol_sd[j]) {
           c_andisol[j,i] <- NA
        }
     }
  }
}
c_alfisol_sd = apply(c_alfisol, 1, sd, na.rm=TRUE)
c_andisol_se = c_andisol_sd / rowSums(!is.na(c_andisol))

# Entisol
c_entisol = extract(m ,info, 'Entisol', 'Soil.Order')
c_entisol_mean = rowMeans(c_entisol, na.rm=TRUE)
c_entisol_sd = apply(c_entisol, 1, sd, na.rm=TRUE)
# Filter the outliers
for (i in 1:dim(c_entisol)[2]) {
  for (j in 1:200) {
     if(is.na(c_entisol_sd[j]) | is.na(c_entisol[j,i])){
        # Do nothing
     } else{
        if(c_entisol[j,i] > c_entisol_mean[j] + 2*c_entisol_sd[j]) {
           c_entisol[j,i] <- NA
        }
     }
  }
}
c_entisol_sd = apply(c_entisol, 1, sd, na.rm=TRUE)
c_entisol_se = c_entisol_sd / rowSums(!is.na(c_entisol))

# Histosol
c_histosol = extract(m ,info, 'Histosol', 'Soil.Order')
c_histosol_mean = rowMeans(c_histosol, na.rm=TRUE)
c_histosol_sd = apply(c_histosol, 1, sd, na.rm=TRUE)
# Filter the outliers
for (i in 1:dim(c_histosol)[2]) {
  for (j in 1:200) {
     if(is.na(c_histosol_sd[j]) | is.na(c_histosol[j,i])){
        # Do nothing
     } else{
        if(c_histosol[j,i] > c_histosol_mean[j] + 2*c_histosol_sd[j]) {
           c_histosol[j,i] <- NA
        }
     }
  }
}
c_histosol_sd = apply(c_histosol, 1, sd, na.rm=TRUE)
c_histosol_se = c_histosol_sd / rowSums(!is.na(c_histosol))

# Inceptisol
c_inceptisol = extract(m ,info, 'Inceptisol', 'Soil.Order')
c_inceptisol_mean = rowMeans(c_inceptisol, na.rm=TRUE)
c_inceptisol_sd = apply(c_inceptisol, 1, sd, na.rm=TRUE)
# Filter the outliers
for (i in 1:dim(c_inceptisol)[2]) {
  for (j in 1:200) {
     if(is.na(c_inceptisol_sd[j]) | is.na(c_inceptisol[j,i])){
        # Do nothing
     } else{
        if(c_inceptisol[j,i] > c_inceptisol_mean[j] + 2*c_inceptisol_sd[j]) {
           c_inceptisol[j,i] <- NA
        }
     }
  }
}
c_inceptisol_sd = apply(c_inceptisol, 1, sd, na.rm=TRUE)
c_inceptisol_se = c_inceptisol_sd / rowSums(!is.na(c_inceptisol))

# Spodosol
c_spodosol = extract(m ,info, 'Spodosol', 'Soil.Order')
c_spodosol_mean = rowMeans(c_spodosol, na.rm=TRUE)
c_spodosol_sd = apply(c_spodosol, 1, sd, na.rm=TRUE)
# Filter the outliers
for (i in 1:dim(c_spodosol)[2]) {
  for (j in 1:200) {
     if(is.na(c_spodosol_sd[j]) | is.na(c_spodosol[j,i])){
        # Do nothing
     } else{
        if(c_spodosol[j,i] > c_spodosol_mean[j] + 2*c_spodosol_sd[j]) {
           c_spodosol[j,i] <- NA
        }
     }
  }
}
c_spodosol_sd = apply(c_spodosol, 1, sd, na.rm=TRUE)
c_spodosol_se = c_spodosol_sd / rowSums(!is.na(c_spodosol))

# Gelisol is seperated into three suborder
# Histel
c_gelisol_histel = extract(m ,info, 'Histel', 'Suborder')
c_gelisol_histel_mean = rowMeans(c_gelisol_histel, na.rm=TRUE)
c_gelisol_histel_sd = apply(c_gelisol_histel, 1, sd, na.rm=TRUE)
# Filter the outliers
for (i in 1:dim(c_gelisol_histel)[2]) {
  for (j in 1:200) {
     if(is.na(c_gelisol_histel_sd[j]) | is.na(c_gelisol_histel[j,i])){
        # Do nothing
     } else{
        if(c_gelisol_histel[j,i] > c_gelisol_histel_mean[j] + 2*c_gelisol_histel_sd[j]) {
           c_gelisol_histel[j,i] <- NA
        }
     }
  }
}
c_gelisol_histel_sd = apply(c_gelisol_histel, 1, sd, na.rm=TRUE)
c_gelisol_histel_se = c_gelisol_histel_sd / rowSums(!is.na(c_gelisol_histel))

# Orthel
c_gelisol_orthel = extract(m ,info, 'Orthel', 'Suborder')
c_gelisol_orthel_mean = rowMeans(c_gelisol_orthel, na.rm=TRUE)
c_gelisol_orthel_sd = apply(c_gelisol_orthel, 1, sd, na.rm=TRUE)
# Filter the outliers
for (i in 1:dim(c_gelisol_orthel)[2]) {
  for (j in 1:200) {
     if(is.na(c_gelisol_orthel_sd[j]) | is.na(c_gelisol_orthel[j,i])){
        # Do nothing
     } else{
        if(c_gelisol_orthel[j,i] > c_gelisol_orthel_mean[j] + 2*c_gelisol_orthel_sd[j]) {
           c_gelisol_orthel[j,i] <- NA
        }
     }
  }
}
c_gelisol_orthel_sd = apply(c_gelisol_orthel, 1, sd, na.rm=TRUE)
c_gelisol_orthel_se = c_gelisol_orthel_sd / rowSums(!is.na(c_gelisol_orthel))

# Turbel
c_gelisol_turbel = extract(m ,info, 'Turbel', 'Suborder')
c_gelisol_turbel_mean = rowMeans(c_gelisol_turbel, na.rm=TRUE)
c_gelisol_turbel_sd = apply(c_gelisol_turbel, 1, sd, na.rm=TRUE)
# Filter the outliers
for (i in 1:dim(c_gelisol_turbel)[2]) {
  for (j in 1:200) {
     if(is.na(c_gelisol_turbel_sd[j]) | is.na(c_gelisol_turbel[j,i])){
        # Do nothing
     } else{
        if(c_gelisol_turbel[j,i] > c_gelisol_turbel_mean[j] + 2*c_gelisol_turbel_sd[j]) {
           c_gelisol_turbel[j,i] <- NA
        }
     }
  }
}
c_gelisol_turbel_sd = apply(c_gelisol_turbel, 1, sd, na.rm=TRUE)
c_gelisol_turbel_se = c_gelisol_turbel_sd / rowSums(!is.na(c_gelisol_turbel))

## Raw profile data
# Alfisol
depth_alfisol <- mishradata$nodedepth[mishradata$idx[mishradata$Soil.Order == "Alfisol"]]
c_data_alfisol <- mishradata$C.Density[mishradata$idx[mishradata$Soil.Order == "Alfisol"]]
#c_data_alfisol <- mishradata$Ccontent[mishradata$idx[mishradata$Soil.Order == "Alfisol"]]
# Andisol
depth_andisol <- mishradata$nodedepth[mishradata$idx[mishradata$Soil.Order == "Andisol"]]
c_data_andisol <- mishradata$C.Density[mishradata$idx[mishradata$Soil.Order == "Andisol"]]
#c_data_andisol <- mishradata$Ccontent[mishradata$idx[mishradata$Soil.Order == "Andisol"]]
# Entisol
depth_entisol <- mishradata$nodedepth[mishradata$idx[mishradata$Soil.Order == "Entisol"]]
c_data_entisol <- mishradata$C.Density[mishradata$idx[mishradata$Soil.Order == "Entisol"]]
#c_data_entisol <- mishradata$Ccontent[mishradata$idx[mishradata$Soil.Order == "Entisol"]]
# Histosol
depth_histosol <- mishradata$nodedepth[mishradata$idx[mishradata$Soil.Order == "Histosol"]]
c_data_histosol <- mishradata$C.Density[mishradata$idx[mishradata$Soil.Order == "Histosol"]]
#c_data_histosol <- mishradata$Ccontent[mishradata$idx[mishradata$Soil.Order == "Histosol"]]
# Inceptisol
depth_inceptisol <- mishradata$nodedepth[mishradata$idx[mishradata$Soil.Order == "Inceptisol"]]
c_data_inceptisol <- mishradata$C.Density[mishradata$idx[mishradata$Soil.Order == "Inceptisol"]]
#c_data_inceptisol <- mishradata$Ccontent[mishradata$idx[mishradata$Soil.Order == "Inceptisol"]]
# Spodosol
depth_spodosol <- mishradata$nodedepth[mishradata$idx[mishradata$Soil.Order == "Spodosol"]]
c_data_spodosol <- mishradata$C.Density[mishradata$idx[mishradata$Soil.Order == "Spodosol"]]
#c_data_spodosol <- mishradata$Ccontent[mishradata$idx[mishradata$Soil.Order == "Spodosol"]]
# Gelisol is seperated into three suborder
# Histel
depth_gelisol_histel <- mishradata$nodedepth[mishradata$idx[mishradata$Suborder == "Histel"]]
c_data_gelisol_histel <- mishradata$C.Density[mishradata$idx[mishradata$Suborder == "Histel"]]
#c_data_gelisol_histel <- mishradata$Ccontent[mishradata$idx[mishradata$Suborder == "Histel"]]
# Orthel
depth_gelisol_orthel <- mishradata$nodedepth[mishradata$idx[mishradata$Suborder == "Orthel"]]
c_data_gelisol_orthel <- mishradata$C.Density[mishradata$idx[mishradata$Suborder == "Orthel"]]
#c_data_gelisol_orthel <- mishradata$Ccontent[mishradata$idx[mishradata$Suborder == "Orthel"]]
# Turbel
depth_gelisol_turbel <- mishradata$nodedepth[mishradata$idx[mishradata$Suborder == "Turbel"]]
c_data_gelisol_turbel <- mishradata$C.Density[mishradata$idx[mishradata$Suborder == "Turbel"]]
#c_data_gelisol_turbel <- mishradata$Ccontent[mishradata$idx[mishradata$Suborder == "Turbel"]]

## Plot profiles
# Assign depths
depth <- seq(1,200,res)

# Alfisol
# Open PDF file
if(plot_ensemble){
   pdf(file="Alfisol_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Alfisol_clean.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,res),each=dim(c_alfisol)[2]), ncol=dim(c_alfisol)[2], byrow=TRUE)
   matplot(c_alfisol, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_alfisol, depth_alfisol, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      xlim=c(0, 180), ylim=rev(c(0,200)), pch=17, cex= 0, xaxt='n', yaxt='n')
   lines(c_alfisol_mean[depth], depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_alfisol_mean[depth]+c_alfisol_sd[depth], depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_alfisol_mean[depth]-c_alfisol_sd[depth], depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,180.0,by=20),cex.axis=3)
axis(3,at=seq(0.0,180.0,by=20),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
#legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
legend(60, 130, legend = c("Fitted profile"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m3)", side=3, outer=TRUE,line=-3.7,cex=3)
text(70, 120, "Alfisol", cex = 4)
dev.off()

# Andisol
# Open PDF file
if(plot_ensemble){
   pdf(file="Andisol_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Andisol_clean.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,res),each=dim(c_andisol)[2]), ncol=dim(c_andisol)[2], byrow=TRUE)
   matplot(c_andisol, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_andisol, depth_andisol, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      xlim=c(0,180), ylim=rev(c(0,100)), pch=17, cex= 0, xaxt='n', yaxt='n')
   lines(c_andisol_mean[depth], depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_andisol_mean[depth]+c_andisol_sd[depth], depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_andisol_mean[depth]-c_andisol_sd[depth], depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,180.0,by=20),cex.axis=3)
axis(3,at=seq(0.0,180.0,by=20),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
#legend(92, 65, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
legend(80, 70, legend = c("Fitted profile"), cex=3, lty=1, col=c("black","gray20"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m3)", side=3, outer=TRUE,line=-3.7,cex=3)
text(105, 62, "Andisol", cex = 4)
dev.off()

# Entisol
# Open PDF file
if(plot_ensemble){
   pdf(file="Entisol_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Entisol_clean.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,res),each=dim(c_entisol)[2]), ncol=dim(c_entisol)[2], byrow=TRUE)
   matplot(c_entisol, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else {
   plot(c_data_entisol, depth_entisol, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      xlim=c(0,150), ylim=rev(c(0,200)), pch=17, cex= 0, xaxt='n', yaxt='n')
   lines(c_entisol_mean[depth], depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_entisol_mean[depth]+c_entisol_sd[depth], depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_entisol_mean[depth]-c_entisol_sd[depth], depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,200.0,by=20),cex.axis=3)
axis(3,at=seq(0.0,200.0,by=20),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
#legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
legend(70, 130, legend = c("Fitted profile"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m3)", side=3, outer=TRUE,line=-3.7,cex=3)
text(90, 120, "Entisol", cex = 4)
dev.off()

# Histosol
# Open PDF file
if(plot_ensemble){
   pdf(file="Histosol_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Histosol_clean.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,res),each=dim(c_histosol)[2]), ncol=dim(c_histosol)[2], byrow=TRUE)
   matplot(ens_Histosol, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_histosol, depth_histosol, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      xlim=c(0,160), ylim=rev(c(0,200)), pch=17, cex= 0, xaxt='n', yaxt='n')
   lines(c_histosol_mean[depth], depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_histosol_mean[depth]+c_histosol_sd[depth], depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_histosol_mean[depth]-c_histosol_sd[depth], depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,160.0,by=20),cex.axis=3)
axis(3,at=seq(0.0,160.0,by=20),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
#legend(130, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
legend(55, 140, legend = c("Fitted profile"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m3)", side=3, outer=TRUE,line=-3.7,cex=3)
text(100, 130, "Histosol", cex = 4)
dev.off()

# Inceptisol
# Open PDF file
if(plot_ensemble){
   pdf(file="Inceptisol_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Inceptisol_clean.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,res),each=dim(c_inceptisol)[2]), ncol=dim(c_inceptisol)[2], byrow=TRUE)
   matplot(c_inceptisol, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_inceptisol, depth_inceptisol, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      xlim=c(0,160), ylim=rev(c(0,200)), pch=17, cex= 0, xaxt='n', yaxt='n')
   lines(c_inceptisol_mean[depth], depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_inceptisol_mean[depth]+c_inceptisol_sd[depth], depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_inceptisol_mean[depth]-c_inceptisol_sd[depth], depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,160.0,by=20),cex.axis=3)
axis(3,at=seq(0.0,160.0,by=20),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
#legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
legend(60, 130, legend = c("Fitted profile"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m3)", side=3, outer=TRUE,line=-3.7,cex=3)
text(80, 120, "Inceptisol", cex = 4)
dev.off()

# Spodosol
# Open PDF file
if(plot_ensemble){
   pdf(file="Spodosol_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Spodosol_clean.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,res),each=dim(c_spodosol)[2]), ncol=dim(c_spodosol)[2], byrow=TRUE)
   matplot(c_spodosol, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_spodosol, depth_spodosol, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
      xlim=c(0,200), ylim=rev(c(0,200)), pch=17, cex= 0, xaxt='n', yaxt='n')
   lines(c_spodosol_mean[depth], depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_spodosol_mean[depth]+c_spodosol_sd[depth], depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_spodosol_mean[depth]-c_spodosol_sd[depth], depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,200.0,by=20),cex.axis=3)
axis(3,at=seq(0.0,200.0,by=20),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
#legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
legend(80, 130, legend = c("Fitted profile"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m3)", side=3, outer=TRUE,line=-3.7,cex=3)
text(110, 120, "Spodosol", cex = 4)
dev.off()

# Gelisol
# Histel
# Open PDF file
if(plot_ensemble){
   pdf(file="Histel_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Histel_clean.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,res),each=dim(c_gelisol_histel)[2]), ncol=dim(c_gelisol_histel)[2], byrow=TRUE)
   matplot(c_gelisol_histel, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_gelisol_histel, depth_gelisol_histel, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
   #   xlim=c(0,70), ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
     xlim=c(0,160), ylim=rev(c(0,200)), pch=17, cex= 0, xaxt='n', yaxt='n')
   lines(c_gelisol_histel_mean[depth], depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_gelisol_histel_mean[depth]+c_gelisol_histel_sd[depth], depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_gelisol_histel_mean[depth]-c_gelisol_histel_sd[depth], depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,200.0,by=10),cex.axis=3)
axis(3,at=seq(0.0,200.0,by=10),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
#legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
legend(60, 135, legend = c("Fitted profile"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m3)", side=3, outer=TRUE,line=-3.7,cex=3)
text(110, 130, "Gelisol - Histel", cex = 4)
dev.off()

# Orthel
# Open PDF file
if(plot_ensemble){
   pdf(file="Orthel_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Orthel_clean.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,res),each=dim(c_gelisol_orthel)[2]), ncol=dim(c_gelisol_orthel)[2], byrow=TRUE)
   matplot(c_gelisol_orthel, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_gelisol_orthel, depth_gelisol_orthel, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
    #  xlim=c(0,100), ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
      xlim=c(0,160), ylim=rev(c(0,200)), pch=17, cex= 0, xaxt='n', yaxt='n')
   lines(c_gelisol_orthel_mean[depth], depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_gelisol_orthel_mean[depth]+c_gelisol_orthel_sd[depth], depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_gelisol_orthel_mean[depth]-c_gelisol_orthel_sd[depth], depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,160.0,by=20),cex.axis=3)
axis(3,at=seq(0.0,160.0,by=20),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
#legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
legend(60, 130, legend = c("Fitted profile"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m3)", side=3, outer=TRUE,line=-3.7,cex=3)
text(100, 120, "Gelisol - Orthel", cex = 4)
dev.off()

# Turbel
# Open PDF file
if(plot_ensemble){
   pdf(file="Turbel_ensemble.pdf", width=17, height=10)
} else{
   pdf(file="Turbel_clean.pdf", width=17, height=10)
}
# Set margin
par(omi=c(0,0.3,0,0), mar=c(2,8,8,2)) # all plots on one page 
if(plot_ensemble){
   depth_ens = matrix(rep(seq(1,200,res),each=dim(c_gelisol_turbel)[2]), ncol=dim(c_gelisol_turbel)[2], byrow=TRUE)
   matplot(c_gelisol_turbel, depth_ens, type = 'l', lty = 1, lwd = 2, xaxt='n', yaxt='n', xlab="", ylab="", ylim=rev(c(0,200)))
} else{
   plot(c_data_gelisol_turbel, depth_gelisol_turbel, xlab="", ylab="", #xlab = "Cumulative root fraction", ylab = "Soil depth (cm)", 
    #  xlim=c(0,60), ylim=rev(c(0,200)), pch=17, cex= 1, xaxt='n', yaxt='n')
      xlim=c(0,160), ylim=rev(c(0,200)), pch=17, cex= 0, xaxt='n', yaxt='n')
   lines(c_gelisol_turbel_mean[depth], depth, lwd = 4, xaxt='n', yaxt='n')
   lines(c_gelisol_turbel_mean[depth]+c_gelisol_turbel_sd[depth], depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
   lines(c_gelisol_turbel_mean[depth]-c_gelisol_turbel_sd[depth], depth, lwd = 4, xaxt='n', yaxt='n', col="gray50")
}
# Move axis to top and left
axis(3,at=seq(0.0,160.0,by=20),cex.axis=3)
axis(3,at=seq(0.0,160.0,by=20),labels=FALSE,tcl=-0.25)
axis(2,at=seq(200,0,by=-20),cex.axis=3, las=1)
axis(2,at=seq(200,0,by=-20),labels=FALSE,tcl=-0.25)
grid(NULL, NULL, col = "lightgray", lty = "dotted", lwd = 3, equilogs = TRUE)
#legend(108, 135, legend = c("   Soil Samples"), cex=3, pch = 17, bty = "n")
legend(50, 130, legend = c("Fitted profile"), cex=3, lty=1, col=c("black","gray50"), lwd=3, bty = "n")#, lty = 1:7,
     #  pch = "*", ncol = 4, cex = 0.8)
# Put Label for X and Y axis
mtext("Soil Depth (cm)", side=2, outer=TRUE, line=-2, cex=3)
mtext("Organic Carbon Density (kgC/m3)", side=3, outer=TRUE,line=-3.7,cex=3)
text(100, 120, "Gelisol - Turbel", cex = 4)
dev.off()


