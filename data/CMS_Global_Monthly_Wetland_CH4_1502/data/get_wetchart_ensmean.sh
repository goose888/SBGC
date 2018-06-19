#!/bin/bash

# ncks -d model,12,17 WetCHARTs_extended_ensemble.nc4 WetCHARTs_extended_high.nc4
# ncks -d time,0,119 WetCHARTs_extended_high.nc4 WetCHARTs_extended_high_0110.nc4
# # Get mean across different model.
# ncks -O --mk_rec_dmn model WetCHARTs_extended_high_0110.nc4 WetCHARTs_extended_high_0110.nc4
# # Attention, the following ncwa command requires a NCO version higher than 4.4.8, otherwise will report a bug
# ncwa -O -a model WetCHARTs_extended_high_0110.nc4 WetCHARTs_ext_mean_high_0110.nc4
# 
# Unit: mg m-2 day-1
# Get monthly total [gCH4 m-2 mon-1]
ncap2 -S mo2s.nco WetCHARTs_ext_mean_high_0110.nc4 WetCHARTs_ext_annual_high_0110.nc4
# Split into 10 years
ncks -d time,0,11 WetCHARTs_ext_annual_high_0110.nc4 WetCHARTs_ext_annual_high_2001.nc4
ncks -d time,12,23 WetCHARTs_ext_annual_high_0110.nc4 WetCHARTs_ext_annual_high_2002.nc4
ncks -d time,24,35 WetCHARTs_ext_annual_high_0110.nc4 WetCHARTs_ext_annual_high_2003.nc4
ncks -d time,36,47 WetCHARTs_ext_annual_high_0110.nc4 WetCHARTs_ext_annual_high_2004.nc4
ncks -d time,48,59 WetCHARTs_ext_annual_high_0110.nc4 WetCHARTs_ext_annual_high_2005.nc4
ncks -d time,60,71 WetCHARTs_ext_annual_high_0110.nc4 WetCHARTs_ext_annual_high_2006.nc4
ncks -d time,72,83 WetCHARTs_ext_annual_high_0110.nc4 WetCHARTs_ext_annual_high_2007.nc4
ncks -d time,84,95 WetCHARTs_ext_annual_high_0110.nc4 WetCHARTs_ext_annual_high_2008.nc4
ncks -d time,96,107 WetCHARTs_ext_annual_high_0110.nc4 WetCHARTs_ext_annual_high_2009.nc4
ncks -d time,108,119 WetCHARTs_ext_annual_high_0110.nc4 WetCHARTs_ext_annual_high_2010.nc4

# Accumulating data through each year
ncks -O --mk_rec_dmn time WetCHARTs_ext_annual_high_2001.nc4 WetCHARTs_ext_annual_high_2001.nc4
ncra -h -O -y ttl WetCHARTs_ext_annual_high_2001.nc4 WetCHARTs_ext_annualtot_high_2001.nc4
ncks -O --mk_rec_dmn time WetCHARTs_ext_annual_high_2002.nc4 WetCHARTs_ext_annual_high_2002.nc4
ncra -h -O -y ttl WetCHARTs_ext_annual_high_2002.nc4 WetCHARTs_ext_annualtot_high_2002.nc4
ncks -O --mk_rec_dmn time WetCHARTs_ext_annual_high_2003.nc4 WetCHARTs_ext_annual_high_2003.nc4
ncra -h -O -y ttl WetCHARTs_ext_annual_high_2003.nc4 WetCHARTs_ext_annualtot_high_2003.nc4
ncks -O --mk_rec_dmn time WetCHARTs_ext_annual_high_2004.nc4 WetCHARTs_ext_annual_high_2004.nc4
ncra -h -O -y ttl WetCHARTs_ext_annual_high_2004.nc4 WetCHARTs_ext_annualtot_high_2004.nc4
ncks -O --mk_rec_dmn time WetCHARTs_ext_annual_high_2005.nc4 WetCHARTs_ext_annual_high_2005.nc4
ncra -h -O -y ttl WetCHARTs_ext_annual_high_2005.nc4 WetCHARTs_ext_annualtot_high_2005.nc4
ncks -O --mk_rec_dmn time WetCHARTs_ext_annual_high_2006.nc4 WetCHARTs_ext_annual_high_2006.nc4
ncra -h -O -y ttl WetCHARTs_ext_annual_high_2006.nc4 WetCHARTs_ext_annualtot_high_2006.nc4
ncks -O --mk_rec_dmn time WetCHARTs_ext_annual_high_2007.nc4 WetCHARTs_ext_annual_high_2007.nc4
ncra -h -O -y ttl WetCHARTs_ext_annual_high_2007.nc4 WetCHARTs_ext_annualtot_high_2007.nc4
ncks -O --mk_rec_dmn time WetCHARTs_ext_annual_high_2008.nc4 WetCHARTs_ext_annual_high_2008.nc4
ncra -h -O -y ttl WetCHARTs_ext_annual_high_2008.nc4 WetCHARTs_ext_annualtot_high_2008.nc4
ncks -O --mk_rec_dmn time WetCHARTs_ext_annual_high_2009.nc4 WetCHARTs_ext_annual_high_2009.nc4
ncra -h -O -y ttl WetCHARTs_ext_annual_high_2009.nc4 WetCHARTs_ext_annualtot_high_2009.nc4
ncks -O --mk_rec_dmn time WetCHARTs_ext_annual_high_2010.nc4 WetCHARTs_ext_annual_high_2010.nc4
ncra -h -O -y ttl WetCHARTs_ext_annual_high_2010.nc4 WetCHARTs_ext_annualtot_high_2010.nc4

# Calculate decadal mean
ncea WetCHARTs_ext_annualtot_high_2001.nc4 WetCHARTs_ext_annualtot_high_2002.nc4 WetCHARTs_ext_annualtot_high_2003.nc4 WetCHARTs_ext_annualtot_high_2004.nc4 WetCHARTs_ext_annualtot_high_2005.nc4 WetCHARTs_ext_annualtot_high_2006.nc4 WetCHARTs_ext_annualtot_high_2007.nc4 WetCHARTs_ext_annualtot_high_2008.nc4 WetCHARTs_ext_annualtot_high_2009.nc4 WetCHARTs_ext_annualtot_high_2010.nc4 WetCHARTs_ch4_annualtot_high_0110.nc4
ncwa -O -a time WetCHARTs_ch4_annualtot_high_0110.nc4 WetCHARTs_ch4_annualtot_high_0110.nc4

# # Purge all temporary files
# rm *.tmp
# rm WetCHARTs_extended_high.nc4
# rm WetCHARTs_extended_high_0110.nc4
# rm WetCHARTs_ext_mean_high_0110.nc4
# rm WetCHARTs_ext_annual_high_2001.nc4
# rm WetCHARTs_ext_annual_high_2002.nc4
# rm WetCHARTs_ext_annual_high_2003.nc4
# rm WetCHARTs_ext_annual_high_2004.nc4
# rm WetCHARTs_ext_annual_high_2005.nc4
# rm WetCHARTs_ext_annual_high_2006.nc4
# rm WetCHARTs_ext_annual_high_2007.nc4
# rm WetCHARTs_ext_annual_high_2008.nc4
# rm WetCHARTs_ext_annual_high_2009.nc4
# rm WetCHARTs_ext_annual_high_2010.nc4
# rm WetCHARTs_ext_annual_high_0110.nc4
# rm WetCHARTs_ext_annualtot_high_2001.nc4
# rm WetCHARTs_ext_annualtot_high_2002.nc4
# rm WetCHARTs_ext_annualtot_high_2003.nc4
# rm WetCHARTs_ext_annualtot_high_2004.nc4
# rm WetCHARTs_ext_annualtot_high_2005.nc4
# rm WetCHARTs_ext_annualtot_high_2006.nc4
# rm WetCHARTs_ext_annualtot_high_2007.nc4
# rm WetCHARTs_ext_annualtot_high_2008.nc4
# rm WetCHARTs_ext_annualtot_high_2009.nc4
# rm WetCHARTs_ext_annualtot_high_2010.nc4
