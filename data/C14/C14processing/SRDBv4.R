# srdb.R
#
# Script providing error checking and sample processing/analysis of the SRDB
# Ben Bond-Lamberty  December 2013 last updated November 2018 (by ORNL DAAC)

# User Note: run this script from the same directory as the files srdb-data-V4.csv, srdb-data_fields-V4.txt, and srdb-studies-V4.csv .

library('ggplot2')
library('ggExtra')
library('maps')
library('mapdata')

# -----------------------------------------------------------------------------
# Logging function
printlog <- function(msg, ..., ts = TRUE, cr = TRUE) {
  if(ts) cat(date(), " ")
  cat(msg, ...)
  if(cr) cat("\n")
}


# -----------------------------------------------------------------------------
# Main program

printlog("Welcome to srdb.R")

printlog("Reading data file...")
srdb <- read.csv("srdb-data-V4.csv")
printlog("Rows =", nrow(srdb), "columns =", ncol(srdb))

# replace fill-value "-9999" with NA
srdb[srdb == -9999] <- NA

printlog("Reading fields data file...")
srdb_info <- read.csv("srdb-data_fields-V4.txt", sep = ",", comment.char = "#")
printlog("Rows =", nrow(srdb_info), "columns =", ncol(srdb_info))

printlog("Reading studies data file...")
srdb_studies <- read.csv("srdb-studies-V4.csv")
printlog("Rows =", nrow(srdb_studies), "columns =", ncol(srdb_studies))

# -----------------------------------------------------------------------------
# Error checking

# Helper functions to check various field attributes
check_numeric <- function(d, dname = deparse(substitute(d))) { 
  printlog("Checking", dname, "is numeric")
  if(!is.numeric(d)) {
    message(paste(dname, "is not numeric"))
    d[ d == "" ] <- "0"	# don't want blanks to be listed below
    message(paste("- in records:", paste(which(is.na(as.numeric(d))), collapse = " ")))
  }
  invisible(is.numeric(d))
}
check_bounds <- function(d, lim, dname = deparse(substitute(d))) { 	# d should be within lim
  if(check_numeric(d, dname)) {
    printlog("Checking", dname, "in bounds (", lim[ 1 ], ",", lim[ 2 ], ")")
    oob <- d < lim[ 1 ] | d > lim[ 2 ]
    if(any(oob, na.rm = TRUE)) {
      message(paste(dname, "out of bounds (", lim[ 1 ], ",", lim[ 2 ], ")"))
      message(paste("- in records:", paste(which(oob), collapse = " ")))
    }
  }
}
check_lesseq <- function(d1, d2) { 	# d1 should be < =  d2
  greater <- d1 > d2
  if(any(greater, na.rm = TRUE)) {
    message(paste(deparse(substitute(d1)), "greater than", deparse(substitute(d2))))
    message(paste("- in records:", paste(which(greater), collapse = " ")))
  }
}
check_labels <- function(d, labs, dname = deparse(substitute(d))) {		# d should be ascending range
  printlog("Checking", dname, "in (", paste(labs, collapse = ", "), ")")
  inlabs <- d %in% labs
  if(any(!inlabs, na.rm = TRUE)) {
    message(paste(dname, "not in labels", paste(labs, collapse = " ")))
    message(paste("- in records:", paste(which(!inlabs), collapse = " ")))	
  }
}
check_fieldnames <- function(d, d_info) {
  fnames <- as.character(d_info[ d_info[, 2] !=  "", 2 ]) # Names in the srdb-data_fields.txt file
  ndb <- names(d)
  
  if(all(ndb == fnames)) {
    printlog("All names match!")
  } else {
    printlog("Following names do not match between field descriptions file and database:")
    mismatch <- ndb !=  fnames
    warning(c(rownumber = which(mismatch), data = ndb[ which(mismatch) ],
              descrip = fnames[ which(mismatch) ]))
  }
}
#_______________________________________________


printlog("-----------------------------------")
printlog("Error checking:")
check_fieldnames(srdb, srdb_info)

with(srdb, {
  check_bounds(Study_number, c(1, 19999))
  study_number_matches <- srdb$Study_number %in% srdb_studies$Study_number
  if(any(!study_number_matches)) {
    warning("Study numbers not in studies file: ", srdb$Study_number[!study_number_matches])    
  }
  check_bounds(Study_midyear, c(1960, 2018))
  check_bounds(YearsOfData, c(1, 99))
  check_bounds(Latitude, c(-90, 90))
  check_bounds(Longitude, c(-180, 180))
  check_bounds(Elevation, c(0, 7999))
  check_bounds(Age_ecosystem, c(0, 999))
  check_bounds(Age_disturbance, c(0, 999))
  unmanaged_ag <- srdb$Ecosystem_type == "Agriculture" & srdb$Ecosystem_state != "Managed"
  if(any(unmanaged_ag)) {
    warning("Non-managed agriculture in records: ", paste(srdb$Record_number[which(unmanaged_ag)], collapse = " "))    
  }
  check_labels(Soil_drainage, c("Dry", "Wet", "Medium", "Mixed", ""))
  check_bounds(Soil_BD, c(0.01, 99.9))
  check_bounds(Soil_CN, c(0.01, 99.9))
  check_bounds(Soil_sand, c(0.0, 999.9))
  check_bounds(Soil_silt, c(0.0, 999.9))
  check_bounds(Soil_clay, c(0.0, 999.9))
  check_bounds(MAT, c(-30, 40))
  check_bounds(MAP, c(0, 9999))
  check_bounds(PET, c(0, 9999))
  check_bounds(Study_temp, c(-30, 40))
  check_bounds(Study_precip, c(0, 9999))
  # TODO: meas_method one of a few values
  check_bounds(Meas_interval, c(0.01, 365))
  check_bounds(Annual_coverage, c(0.01, 1))
  # TODO: Partition_method one of a few values
  
  check_bounds(Rs_annual, c(-200, 5500))
  check_bounds(Rs_annual_err, c(0, 5500))
  check_bounds(Rs_interann_err, c(0, 5500))
  check_bounds(Rs_max, c(0, 100))
  check_bounds(Rs_maxday, c(1, 365))
  check_bounds(Rs_min, c(-1, 100))
  check_bounds(Rs_minday, c(1, 365))
  check_lesseq(Rs_min, Rs_max)
  # TODO: check following don't exceed Rs_annual
  check_bounds(Rlitter_annual, c(0, 5000))
  check_bounds(Ra_annual, c(0, 5000))
  check_bounds(Rh_annual, c(0, 5000))
  check_lesseq(Rlitter_annual, Rs_annual)
  check_lesseq(Ra_annual, Rs_annual)
  check_lesseq(Rh_annual, Rs_annual)
  check_bounds(RC_annual, c(0, 1))
  check_bounds(Rs_spring, c(0, 100))
  check_bounds(Rs_summer, c(0, 100))
  check_bounds(Rs_autumn, c(0, 100))
  check_bounds(Rs_winter, c(0, 100))
  check_bounds(Rs_growingseason, c(0, 100))
  check_bounds(Rs_wet, c(0, 100))
  check_bounds(Rs_dry, c(0, 100))
  check_bounds(RC_seasonal, c(0, 1))
  # TODO: RC_season one of a few values
  
  check_labels(Temp_effect, c("Positive", "Negative", "None", "Mixed", ""))
  check_numeric(Model_temp_min)
  check_numeric(Model_temp_max)
  check_lesseq(Model_temp_min, Model_temp_max)
  check_bounds(Model_N, c(2, 99999))
  check_bounds(Model_R2, c(0, 1))
  check_bounds(T_depth, c(-200, 100))
  check_numeric(Model_paramA)
  check_numeric(Model_paramB)
  check_numeric(Model_paramC)
  check_numeric(Model_paramD)
  check_numeric(Model_paramE)
  check_labels(WC_effect, c("Positive", "Negative", "None", "Mixed", ""))
  check_bounds(R10, c(0, 100))
  check_bounds(Q10_0_10, c(0.1, 200))
  check_bounds(Q10_5_15, c(0.1, 200))
  check_bounds(Q10_10_20, c(0.1, 200))
  check_bounds(Q10_0_20, c(0.1, 200))
  check_bounds(Q10_other1, c(0.1, 200))
  check_numeric(Q10_other1_temp_min)
  check_numeric(Q10_other1_temp_max)
  check_lesseq(Q10_other1_temp_min, Q10_other1_temp_max)
  check_bounds(Q10_other2, c(0.1, 200))
  check_numeric(Q10_other2_temp_min)
  check_numeric(Q10_other2_temp_max)
  check_lesseq(Q10_other2_temp_min, Q10_other2_temp_max)
  
  check_bounds(GPP, c(0, 9999))
  check_bounds(ER, c(0, 9999))
  check_bounds(NEP, c(-9999, 9999))
  check_bounds(NPP, c(0, 9999))
  check_bounds(ANPP, c(0, 9999))
  check_bounds(BNPP, c(0, 9999))
  check_bounds(NPP_FR, c(0, 9999))
  check_lesseq(ANPP, NPP)
  check_lesseq(BNPP, NPP)
  check_lesseq(NPP_FR, BNPP)
  
  check_bounds(TBCA, c(0, 9999))
  check_bounds(Litter_flux, c(0, 9999))
  check_bounds(Rootlitter_flux, c(0, 9999))
  check_bounds(TotDet_flux, c(0, 9999))
  check_lesseq(Litter_flux + Rootlitter_flux, TotDet_flux)
  
  check_bounds(Ndep, c(0, 999))
  check_bounds(LAI, c(0, 99))
  check_bounds(BA, c(0, 999))
  check_bounds(C_veg_total, c(0, 99999))
  check_bounds(C_AG, c(0, 999999))
  check_bounds(C_BG, c(0, 99999))
  check_bounds(C_CR, c(0, 99999))
  check_bounds(C_FR, c(0, 99999))
  check_lesseq(C_CR + C_FR, C_BG)
  check_lesseq(C_AG + C_BG, C_veg_total)	
  check_bounds(C_litter, c(0, 99999))
  check_bounds(C_soilmineral, c(0, 999999))
  check_bounds(C_soildepth, c(0, 200))
  
  # Outliers
})

srdb_spatial <- subset(srdb, !is.na(Latitude) & !is.na(Longitude))

#the following requires a file "fractional_land.0.5-deg.nc" which is provided on the SRDB github, but is not part of the DAAC's dataset package.

# # Check whether studies are all on land or not
# library(raster)
# sp <- SpatialPoints(cbind(srdb_spatial$Longitude, srdb_spatial$Latitude))
# landmask <- brick("fractional_land.0.5-deg.nc")
# srdb_spatial$landfrac <- extract(rotate(raster(landmask, 1)), sp)
# srdb_spatial$landfrac_cut <- cut(srdb_spatial$land, 4)
# 
# printlog("Checking points fall on land...")
# if(sum(srdb_spatial$landfrac <= 0.1, na.rm = TRUE)) {
#   notonland <- srdb_spatial$Record_number[srdb_spatial$landfrac < 0.1]
#   message(paste("- low-land records:", paste(notonland, collapse = " ")))
# }

printlog("All done with error checking. <RETURN>")
readline()

# -----------------------------------------------------------------------------
# Summary statistics

printlog("-----------------------------------")
printlog("Basic summary of SRDB and its data:")
printlog("Records  = ", nrow(srdb))
printlog("Fields  = ", ncol(srdb))
printlog("Rs_annual (gC/m2/yr) distribution:")
print(summary(as.numeric(srdb$Rs_annual)))
printlog("R10 (µmol/m2/s) distribution:")
print(summary(as.numeric(srdb$R10)))
printlog("Q10 (10-20 degC) distribution:")
print(summary(as.numeric(srdb$Q10_10_20)))
printlog("Biome distribution:")
print(summary(srdb$Biome))
printlog("Ecosystem distribution:")
print(summary(srdb$Ecosystem_type))
printlog("Ecosystem state distribution:")
print(summary(srdb$Ecosystem_state))
# etc.
printlog("-----------------------------------")


# -----------------------------------------------------------------------------
# Graphics
# -----------------------------------------------------------------------------

if(require(ggplot2) & require(ggExtra)) {
  
  theme_set(theme_bw())
  
  # -----------------------------------------------------------------------------
  # Points over time
  p_rs <- qplot(Study_midyear, Rs_annual, color = Biome,
                data = subset(srdb, !is.na(Study_midyear) & !is.na(Biome) & !is.na(Rs_annual) & Rs_annual >= 0 & Rs_annual <= 5000)) + 
    ylim(c(0, 5000)) + geom_vline(xintercept = 2008) + 
    xlab("Year") + ylab("Rs (gC/m2/yr)")
  png("rs_over_time.png", width = 720, height = 480)
  print(ggExtra::ggMarginal(p_rs, type = "histogram"))
  dev.off()
  p_rh <- qplot(Study_midyear, Rh_annual, color = Biome, 
                data = subset(srdb, !is.na(Study_midyear) & !is.na(Biome) & !is.na(Rh_annual) & Rh_annual >= 0 & Rh_annual <= 5000)) + 
    ylim(c(0, 5000)) + geom_vline(xintercept = 2008) + 
    xlab("Year") + ylab("Rh (gC/m2/yr)")
  png("rh_over_time.png", width = 720, height = 480)
  print(ggExtra::ggMarginal(p_rh, type = "histogram"))
  dev.off()
  
  

  # -----------------------------------------------------------------------------
  # Two simple maps
  
  if(require(maps) & require(mapdata)) {
    printlog("Making maps...")
    world <- map_data("world")
    srdb_spatial$long <- srdb_spatial$Longitude
    srdb_spatial$lat <- srdb_spatial$Latitude
    p1 <- ggplot(srdb_spatial, aes(x = long, y = lat)) + 
      geom_point(aes(color = Biome)) +
      geom_path(data = world, aes(group = group)) +
      scale_y_continuous(breaks = (-2:2) * 30) +
      scale_x_continuous(breaks = (-4:4) * 45) +
      coord_fixed(xlim = c(-180, 180), ylim = c(-90, 90)) 
    print(p1)
    ggsave("map1-world.pdf")
    
    p1na <- p1 + coord_fixed(xlim = c(-160, -50), ylim = c(15, 75))
    print(p1na)
    ggsave("map1-northamerica.pdf")
    
    p1europe <- p1 + coord_fixed(xlim = c(-10, 45), ylim = c(30, 70))
    print(p1europe)
    ggsave("map1-europe.pdf")
    
    p1china <- p1 + coord_fixed(xlim = c(75, 135), ylim = c(20, 55))
    print(p1china)
    ggsave("map1-china.pdf")
    
    # below does not work because landfrac is not found
    # p2 <- ggplot(subset(srdb_spatial, landfrac < 0.1), aes(x = long, y = lat)) + 
    #   geom_point(aes(color = Biome)) +
    #   geom_path(data = world, aes(group = group)) +
    #   scale_y_continuous(breaks = (-2:2) * 30) +
    #   scale_x_continuous(breaks = (-4:4) * 45) +
    #   ggtitle("Low land fraction records") +
    #   coord_fixed(xlim = c(-180, 180), ylim = c(-90, 90)) 
    # print(p2)
    # ggsave("map2-lowland-world.pdf")
    # 
    # p2na <- p2 + coord_fixed(xlim = c(-160, -50), ylim = c(15, 75))
    # print(p2na)
    # ggsave("map2-lowland-northamerica.pdf")
    # 
    # p2europe <- p2 + coord_fixed(xlim = c(-10, 45), ylim = c(30, 70))
    # print(p2europe)
    # ggsave("map2-lowland-europe.pdf")
    # 
    # p2china <- p2 + coord_fixed(xlim = c(75, 135), ylim = c(20, 55))
    # print(p2china)
    # ggsave("map2-lowland-china.pdf")
  } else printlog("No maps produced - maps not installed")
  
  # -----------------------------------------------------------------------------
  # Summary graphics
  
  # We limit the data (toss out negative values, etc.) to make figures more useful
  printlog("Limiting data set...")
  srdb1 <- subset(srdb, Rs_annual > 0 & Rs_annual < 4000)
  
  printlog("Generating and saving graphs...")
  print(qplot(Rs_annual, data = srdb1, fill = Biome))
  ggsave("Rs_annual-distribution-Biome.pdf")
  
  print(qplot(Rs_annual, data = srdb1, fill = Ecosystem_type))
  ggsave("Rs_annual-distribution-Ecosystem.pdf")
  
  print(qplot(Rs_interann_err, data = subset(srdb1, !is.na(Rs_interann_err)), fill = Ecosystem_type))
  ggsave("Rs_interann_err-distribution-Ecosystem.pdf")
  
  print(qplot(Rs_annual, Rh_annual, color = Biome,
              data = subset(srdb1, !is.na(Rs_annual) & !is.na(Rh_annual)))
        + geom_smooth(method = "lm"))
  ggsave("Rs_annual-vs-Rh_annual.pdf")
  
  srdb3 <- subset(srdb1, !is.na(TBCA) & Rs_annual > 0 & Ecosystem_state == "Natural")
  print(qplot(TBCA, Rs_annual, data = srdb3, color = Biome, main = "Natural ecosystems only") +
          geom_smooth(method = "lm", aes(group = 1))) # interesting!
  ggsave("TBCA-vs-Rs_annual.pdf")
} else printlog("No plots produced - ggplot2 not installed")

printlog("All done. Note that if any warnings occurred, these should be checked!")