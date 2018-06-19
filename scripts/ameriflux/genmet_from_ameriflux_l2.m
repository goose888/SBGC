clear;clc

%% Input needed in this code
inputfile = 'PHAC_forcing_h.nc';
sitename = 'PHAC';
tstep_sec = '1800';
site_timezone = '-6:00';
metfile = 'US-PHAC.nc';
height = 11;

%% Do not modify below
%% Open old NC file
ncid = netcdf.open(inputfile, 'NOWRITE');
% Get dimensions
tstpid = netcdf.inqDimID(ncid, 'tstep');
[~, steps] = netcdf.inqDim(ncid, tstpid);
% Get variables
latid = netcdf.inqVarID(ncid, 'nav_lat');
lat = netcdf.getVar(ncid, latid);
lonid = netcdf.inqVarID(ncid, 'nav_lon');
lon = netcdf.getVar(ncid, lonid);

yearid = netcdf.inqVarID(ncid, 'YEAR');
year = netcdf.getVar(ncid, yearid);
doyid = netcdf.inqVarID(ncid, 'DOY');
doy = netcdf.getVar(ncid, doyid);
hrminid = netcdf.inqVarID(ncid, 'HRMIN');
hrmin = netcdf.getVar(ncid, hrminid);
timeid = netcdf.inqVarID(ncid, 'TIME');
time = netcdf.getVar(ncid, timeid);
timestpid = netcdf.inqVarID(ncid, 'TIMEstp');
timestp = netcdf.getVar(ncid, timestpid);

rainfid = netcdf.inqVarID(ncid, 'Rainf');
rainf = netcdf.getVar(ncid, rainfid);
rainf = squeeze(rainf);
tairid = netcdf.inqVarID(ncid, 'Tair');
tair = netcdf.getVar(ncid, tairid);
tair = squeeze(tair);
rhid = netcdf.inqVarID(ncid, 'RH');
rh = netcdf.getVar(ncid, rhid);
rh = squeeze(rh);
vpdid = netcdf.inqVarID(ncid, 'VPD');
vpd = netcdf.getVar(ncid, vpdid);
vpd = squeeze(vpd);
qairid = netcdf.inqVarID(ncid, 'Qair');
qair = netcdf.getVar(ncid, qairid);
qair = squeeze(qair);
windid = netcdf.inqVarID(ncid, 'Wind');
wind = netcdf.getVar(ncid, windid);
wind = squeeze(wind);
swdownid = netcdf.inqVarID(ncid, 'SWdown');
swdown = netcdf.getVar(ncid, swdownid);
swdown = squeeze(swdown);
parid = netcdf.inqVarID(ncid, 'PAR');
par = netcdf.getVar(ncid, parid);
par = squeeze(par);
lwdownid = netcdf.inqVarID(ncid, 'LWdown');
lwdown = netcdf.getVar(ncid, lwdownid);
lwdown = squeeze(lwdown);
psurfid = netcdf.inqVarID(ncid, 'Psurf');
psurf = netcdf.getVar(ncid, psurfid);
psurf = squeeze(psurf);
aco2id = netcdf.inqVarID(ncid, 'aCO2');
aco2 = netcdf.getVar(ncid, aco2id);
aco2 = squeeze(aco2);
eco2id = netcdf.inqVarID(ncid, 'eCO2');
eco2 = netcdf.getVar(ncid, eco2id);
eco2 = squeeze(eco2);
ndepid = netcdf.inqVarID(ncid, 'Ndep');
ndep = netcdf.getVar(ncid, ndepid);
ndep = squeeze(ndep);

netcdf.close(ncid);


%% Create meteorology data
% Save Preprocessed data into a new netcdf file
ncid = netcdf.create(metfile, 'CLOBBER');
% Define dimensions
tid = netcdf.defDim(ncid, 't', steps);
xid = netcdf.defDim(ncid, 'x', 1);
yid = netcdf.defDim(ncid, 'y', 1);
zid = netcdf.defDim(ncid, 'z', 1);

% Some information needed from inputfile
startyear = num2str(year(1));
starttime = strcat(startyear, '-01-01 00:00:00 ', site_timezone);
createdate = strcat('File created from FACE site data on ', date);
sitelocation = strcat('Latitude:', num2str(lat(1)), ', Longitude:', num2str(lon(1)));

% Define variables
timeid = netcdf.defVar(ncid, 't', 'int', tid);
netcdf.putAtt(ncid, timeid, 't:long_name', 'Time axis');
netcdf.putAtt(ncid, timeid, 't:units', 'seconds since original timestamp');
netcdf.putAtt(ncid, timeid, 't:calendar', 'gregorian');
netcdf.putAtt(ncid, timeid, 't:time_origin', starttime);

tstpid = netcdf.defVar(ncid, 'timestp', 'int', tid);
netcdf.putAtt(ncid, tstpid, 'timestp:long_name', 'Time step axis');
netcdf.putAtt(ncid, tstpid, 'timestp:units', 'timesteps since original timestamp');
netcdf.putAtt(ncid, tstpid, 'timestp:tstep_sec', tstep_sec);
netcdf.putAtt(ncid, tstpid, 'timestp:time_origin', starttime);

lonid = netcdf.defVar(ncid, 'lon', 'double', xid);
netcdf.putAtt(ncid, lonid, 'lon:long_name', 'Longitude');
netcdf.putAtt(ncid, lonid, 'lon:units', 'degrees_east');

latid = netcdf.defVar(ncid, 'lat', 'double', yid);
netcdf.putAtt(ncid, latid, 'lat:long_name', 'Latitude');
netcdf.putAtt(ncid, latid, 'lat:units', 'degrees_north');

levelid = netcdf.defVar(ncid, 'level', 'double', zid);
netcdf.putAtt(ncid, levelid, 'level:long_name', 'Vertical levels');
netcdf.putAtt(ncid, levelid, 'level:units', 'm');

tairid = netcdf.defVar(ncid, 'Tair', 'double', [tid, zid, yid, xid]);
netcdf.putAtt(ncid, tairid, 'Tair:long_name', 'Near surface air temperature');
netcdf.putAtt(ncid, tairid, 'Tair:units', 'K');
netcdf.putAtt(ncid, tairid, 'Tair:_FillValue', -9999.);

% tairflagid = netcdf.defVar(ncid, 'Tair_flag', 'int', [xid, yid, zid, tid]);
% netcdf.putAtt(ncid, tairflagid, 'Tair_flag:long_name', 'Near surface air temperature');
% netcdf.putAtt(ncid, tairflagid, 'Tair_flag:flag_values', '0,1,2,3,4,6,7,8');
% netcdf.putAtt(ncid, tairflagid, 'Tair_flag:flag_meanings', strcat('0_Original, ', ...
%                          '1_Diurnal_mean_fill, 2_Daymet, 3_Daymet_and_dailyNCDC, ', ...
%                          '4_dailyNCDC, 6_hourlyNCDC, 7_nearby_tower, 8_multiple_var'));

qairid = netcdf.defVar(ncid, 'Qair', 'double', [tid, zid, yid, xid]);
netcdf.putAtt(ncid, qairid, 'Qair:long_name', 'Near surface specific humidity');
netcdf.putAtt(ncid, qairid, 'Qair:units', 'kg/kg');
netcdf.putAtt(ncid, qairid, 'Qair:_FillValue', -9999.);

% qairflagid = netcdf.defVar(ncid, 'Qair_flag', 'int', [xid, yid, zid, tid]);
% netcdf.putAtt(ncid, qairflagid, 'Qair_flag:long_name', 'Near surface specific humidity');
% netcdf.putAtt(ncid, qairflagid, 'Qair_flag:flag_values', '0,1,2,3,4,6,7,8');
% netcdf.putAtt(ncid, qairflagid, 'Qair_flag:flag_meanings', strcat('0_Original, ', ...
%                          '1_Diurnal_mean_fill, 2_Daymet, 3_Daymet_and_dailyNCDC, ', ...
%                          '4_dailyNCDC, 6_hourlyNCDC, 7_nearby_tower, 8_multiple_var'));

windid = netcdf.defVar(ncid, 'Wind', 'double', [tid, zid, yid, xid]);
netcdf.putAtt(ncid, windid, 'Wind:long_name', 'Near surface module of the wind');
netcdf.putAtt(ncid, windid, 'Wind:units', 'm/s');
netcdf.putAtt(ncid, windid, 'Wind:_FillValue', -9999.);

% windflagid = netcdf.defVar(ncid, 'Wind_flag', 'int', [xid, yid, zid, tid]);
% netcdf.putAtt(ncid, windflagid, 'Wind_flag:long_name', 'Near surface module of the wind');
% netcdf.putAtt(ncid, windflagid, 'Wind_flag:flag_values', '0,1,2,3,4,6,7,8');
% netcdf.putAtt(ncid, windflagid, 'Wind_flag:flag_meanings', strcat('0_Original, ', ...
%                          '1_Diurnal_mean_fill, 2_Daymet, 3_Daymet_and_dailyNCDC, ', ...
%                          '4_dailyNCDC, 6_hourlyNCDC, 7_nearby_tower, 8_multiple_var'));

rainfid = netcdf.defVar(ncid, 'Rainf', 'double', [tid, zid, yid, xid]);
netcdf.putAtt(ncid, rainfid, 'Rainf:long_name', 'Rainfall rate');
netcdf.putAtt(ncid, rainfid, 'Rainf:units', 'kg/m2/s');
netcdf.putAtt(ncid, rainfid, 'Rainf:_FillValue', -9999.);

% rainfflagid = netcdf.defVar(ncid, 'Rainf_flag', 'int', [xid, yid, zid, tid]);
% netcdf.putAtt(ncid, rainfflagid, 'Rainf_flag:long_name', 'Rainfall rate');
% netcdf.putAtt(ncid, rainfflagid, 'Rainf_flag:flag_values', '0,1,2,3,4,6,7,8');
% netcdf.putAtt(ncid, rainfflagid, 'Rainf_flag:flag_meanings', strcat('0_Original, ', ...
%                          '1_Diurnal_mean_fill, 2_Daymet, 3_Daymet_and_dailyNCDC, ', ...
%                          '4_dailyNCDC, 6_hourlyNCDC, 7_nearby_tower, 8_multiple_var'));

psurfid = netcdf.defVar(ncid, 'Psurf', 'double', [tid, zid, yid, xid]);
netcdf.putAtt(ncid, psurfid, 'Psurf:long_name', 'Surface pressure');
netcdf.putAtt(ncid, psurfid, 'Psurf:units', 'Pa');
netcdf.putAtt(ncid, psurfid, 'Psurf:_FillValue', -9999.);

% psurfflagid = netcdf.defVar(ncid, 'Psurf_flag', 'int', [xid, yid, zid, tid]);
% netcdf.putAtt(ncid, psurfflagid, 'Psurf_flag:long_name', 'Surface pressure');
% netcdf.putAtt(ncid, psurfflagid, 'Psurf_flag:flag_values', '0,1,2,3,4,6,7,8');
% netcdf.putAtt(ncid, psurfflagid, 'Psurf_flag:flag_meanings', strcat('0_Original, ', ...
%                          '1_Diurnal_mean_fill, 2_Daymet, 3_Daymet_and_dailyNCDC, ', ...
%                          '4_dailyNCDC, 6_hourlyNCDC, 7_nearby_tower, 8_multiple_var'));

swdownid = netcdf.defVar(ncid, 'SWdown', 'double', [tid, zid, yid, xid]);
netcdf.putAtt(ncid, swdownid, 'SWdown:long_name', 'Surface incident shortwave radiation');
netcdf.putAtt(ncid, swdownid, 'SWdown:units', 'W/m2');
netcdf.putAtt(ncid, swdownid, 'SWdown:_FillValue', -9999.);

% swdownflagid = netcdf.defVar(ncid, 'SWdown_flag', 'int', [xid, yid, zid, tid]);
% netcdf.putAtt(ncid, swdownflagid, 'SWdown_flag:long_name', 'Surface incident shortwave radiation');
% netcdf.putAtt(ncid, swdownflagid, 'SWdown_flag:flag_values', '0,1,2,3,4,6,7,8');
% netcdf.putAtt(ncid, swdownflagid, 'SWdown_flag:flag_meanings', strcat('0_Original, ', ...
%                          '1_Diurnal_mean_fill, 2_Daymet, 3_Daymet_and_dailyNCDC, ', ...
%                          '4_dailyNCDC, 6_hourlyNCDC, 7_nearby_tower, 8_multiple_var'));
  
lwdownid = netcdf.defVar(ncid, 'LWdown', 'double', [tid, zid, yid, xid]);
netcdf.putAtt(ncid, lwdownid, 'LWdown:long_name', 'Surface incident longwave radiation');
netcdf.putAtt(ncid, lwdownid, 'LWdown:units', 'W/m2');
netcdf.putAtt(ncid, lwdownid, 'LWdown:_FillValue', -9999.);

% lwdownflagid = netcdf.defVar(ncid, 'LWdown_flag', 'int', [xid, yid, zid, tid]);
% netcdf.putAtt(ncid, lwdownflagid, 'LWdown_flag:long_name', 'Surface incident longwave radiation');
% netcdf.putAtt(ncid, lwdownflagid, 'LWdown_flag:flag_values', '0,1,2,3,4,6,7,8');
% netcdf.putAtt(ncid, lwdownflagid, 'LWdown_flag:flag_meanings', strcat('0_Original, ', ...
%                          '1_Diurnal_mean_fill, 2_Daymet, 3_Daymet_and_dailyNCDC, ', ...
%                          '4_dailyNCDC, 6_hourlyNCDC, 7_nearby_tower, 8_multiple_var'));

vpdid = netcdf.defVar(ncid, 'VPD', 'double', [tid, zid, yid, xid]);
netcdf.putAtt(ncid, vpdid, 'VPD:long_name', 'Vapor Pressure Deficit');
netcdf.putAtt(ncid, vpdid, 'VPD:units', 'Pa');
netcdf.putAtt(ncid, vpdid, 'VPD:_FillValue', -9999.);

parid = netcdf.defVar(ncid, 'PAR', 'double', [tid, zid, yid, xid]);
netcdf.putAtt(ncid, parid, 'PAR:long_name', 'Incident or downward photosynthetically active radiation');
netcdf.putAtt(ncid, parid, 'PAR:units', 'umol/m2/s');
netcdf.putAtt(ncid, parid, 'PAR:_FillValue', -9999.);

co2airid = netcdf.defVar(ncid, 'CO2air', 'double', [tid, zid, yid, xid]);
netcdf.putAtt(ncid, co2airid, 'CO2air:long_name', 'Near surface CO2 concentration');
netcdf.putAtt(ncid, co2airid, 'CO2air:flask_site', sitename);
netcdf.putAtt(ncid, co2airid, 'CO2air:units', 'ppmv');
netcdf.putAtt(ncid, co2airid, 'CO2air:_FillValue', -9999.);

eco2airid = netcdf.defVar(ncid, 'ECO2air', 'double', [tid, zid, yid, xid]);
netcdf.putAtt(ncid, eco2airid, 'ECO2air:long_name', 'Near surface Elevated CO2 concentration');
netcdf.putAtt(ncid, eco2airid, 'ECO2air:flask_site', sitename);
netcdf.putAtt(ncid, eco2airid, 'ECO2air:units', 'ppmv');
netcdf.putAtt(ncid, eco2airid, 'ECO2air:_FillValue', -9999.);

ndepid = netcdf.defVar(ncid, 'Ndep', 'double', [tid, zid, yid, xid]);
netcdf.putAtt(ncid, ndepid, 'Ndep:long_name', 'Total N deposition over a time step of measurement (30 minutes)');
netcdf.putAtt(ncid, ndepid, 'Ndep:flask_site', sitename);
netcdf.putAtt(ncid, ndepid, 'Ndep:units', 'gm-2(30min)-1');
netcdf.putAtt(ncid, ndepid, 'Ndep:_FillValue', -9999.);

% Global attributes
gid = netcdf.getConstant('Global');
netcdf.putAtt(ncid, gid, 'institution', 'Oak Ridge National Laboratory');
netcdf.putAtt(ncid, gid, 'history', createdate);
netcdf.putAtt(ncid, gid, 'site_name', sitename);
netcdf.putAtt(ncid, gid, 'nearby_site', sitename);
netcdf.putAtt(ncid, gid, 'site_location', sitelocation);

% End definition
netcdf.endDef(ncid);
% Put values into file
netcdf.putVar(ncid, timeid, time);
netcdf.putVar(ncid, tstpid, timestp);
netcdf.putVar(ncid, lonid, lon);
netcdf.putVar(ncid, latid, lat);
netcdf.putVar(ncid, levelid, height);
netcdf.putVar(ncid, tairid, tair);
netcdf.putVar(ncid, qairid, qair);
netcdf.putVar(ncid, windid, wind);
netcdf.putVar(ncid, rainfid, rainf);
netcdf.putVar(ncid, psurfid, psurf);
netcdf.putVar(ncid, swdownid, swdown);
netcdf.putVar(ncid, lwdownid, lwdown);
netcdf.putVar(ncid, vpdid, vpd);
netcdf.putVar(ncid, parid, par);
netcdf.putVar(ncid, co2airid, aco2);
netcdf.putVar(ncid, eco2airid, eco2);
netcdf.putVar(ncid, ndepid, ndep);

% Close the file
netcdf.close(ncid);
