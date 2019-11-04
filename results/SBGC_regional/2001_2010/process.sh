#!/bin/bash

flist=`ls Global*.nc`
ncea ${flist} mean_bgp3d_2001_2010.nc
