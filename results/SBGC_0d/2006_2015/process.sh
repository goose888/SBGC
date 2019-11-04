#!/bin/bash

flist=`ls *.nc`
ncea ${flist} mean_bgc2d_2006_2015.nc
