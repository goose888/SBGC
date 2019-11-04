#!/bin/bash

flist=`ls Global*.nc`
ncea ${flist} mean_bgp3d_1960_2009.nc
