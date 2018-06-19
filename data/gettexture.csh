#!/bin/tcsh

#----------------------------------------------------------------------------------
# Author: Shijie Shu
# Created on: Nov 2 2015
#----------------------------------------------------------------------------------
# User defined options and variables
set workdir=`pwd`
set textureinfo='mishratexture.csv'
set sitelist='gelisol_info.csv'

#----------------------------------------------------------------------------------
# DO NOT EDIT BELOW
#----------------------------------------------------------------------------------
# Extract site list
set sitelist=`awk -F , 'NR>1 {gsub(/\"/, "", $2); print $2}' $workdir/$sitelist`
# Extract columns and rows selected by site list
set i=1
while($i <= $#sitelist)
   awk -F , '$2 == "'"$sitelist[$i]"'" { print $2","$4","$6","$9","$10","$11}' $workdir/$textureinfo >> texture.csv
   @ i = $i + 1
end
