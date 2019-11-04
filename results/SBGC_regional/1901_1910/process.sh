#!/bin/bash

flist=`ls Global*.nc`
ncea ${flist} mean_bgp3d_1901_1910.nc
