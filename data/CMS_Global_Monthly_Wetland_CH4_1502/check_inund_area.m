clear;clc

% WETX is the wetland area
load inund_extent.mat
% Load the wetland extent used in ISAM
fname = 'fw_frac_max.nc';
isam09 = permute(ncread(fname,'FW'),[2,1]);
lon = ncread(fname,'Lon');
lat = ncread(fname,'Lat');

fname = 'fw_frac_max.nc';
isam10 = permute(ncread(fname,'FW'),[2,1]);

fname = 'fw_frac_max_rcp45.nc';
isam45 = permute(ncread(fname,'FW'),[2,1]);

fname = 'fw_frac_max_rcp85.nc';
isam85 = permute(ncread(fname,'FW'),[2,1]);

tt=mean(WEXT,4);
%tt=WEXT(:,:,:,1);
%tt=WEXT(:,:,:,2);
%tt=WEXT(:,:,:,3);
%tt=WEXT(:,:,:,4);

wet09 = mean(tt(:,:,1:12),3);
wet10 = mean(tt(:,:,13:24),3);

temp = wet09(:,1:360);
wet09(:,1:360) = wet09(:,361:720);
wet09(:,361:720) = temp;

temp = wet10(:,1:360);
wet10(:,1:360) = wet10(:,361:720);
wet10(:,361:720) = temp;

% Store as netcdf file
ncid = netcdf.create('wet2009.nc', 'CLOBBER');
lonid = netcdf.defDim(ncid, 'Lon', 720);
latid = netcdf.defDim(ncid, 'Lat', 360);

longitudeid = netcdf.defVar(ncid, 'Lon', 'double', lonid);
latitudeid = netcdf.defVar(ncid, 'Lat', 'double', latid);
fwid = netcdf.defVar(ncid, 'FW', 'double', [lonid, latid]);

netcdf.endDef(ncid);

netcdf.putVar(ncid, longitudeid, lon);
netcdf.putVar(ncid, latitudeid, lat);
netcdf.putVar(ncid, fwid, wet09');

netcdf.close(ncid);

ncid = netcdf.create('wet2010.nc', 'CLOBBER');
lonid = netcdf.defDim(ncid, 'Lon', 720);
latid = netcdf.defDim(ncid, 'Lat', 360);

longitudeid = netcdf.defVar(ncid, 'Lon', 'double', lonid);
latitudeid = netcdf.defVar(ncid, 'Lat', 'double', latid);
fwid = netcdf.defVar(ncid, 'FW', 'double', [lonid, latid]);

netcdf.endDef(ncid);

netcdf.putVar(ncid, longitudeid, lon);
netcdf.putVar(ncid, latitudeid, lat);
netcdf.putVar(ncid, fwid, wet10');

netcdf.close(ncid);

% Check the area
% First open the nc file to obtain the land mask
ncid = netcdf.open('surfdata_05x05_13reg.nc', 'NOWRITE');
maskid = netcdf.inqVarID(ncid, 'REGION_MASK_CRU_NCEP');
glob_mask = netcdf.getVar(ncid, maskid);

netcdf.close(ncid);

glob_mask(glob_mask<12) = 0;
glob_mask(glob_mask>12) = 0;
glob_mask(glob_mask>0) = 1.0;
glob_mask_cal = double(glob_mask);

% Get the area
grid_area=zeros(720, 360);
EARTH_AREA=5.096e14;
lat=linspace(-89.75, 89.75, 360);
res = 0.5;
nlat = 360;
nlon = 720;

for i = 1:nlat
  for j = 1:nlon
     grid_area(j,i) = (EARTH_AREA/2)*abs(sin((lat(i) - res/2)*pi/180) - ...
              sin((lat(i) + res/2)*pi/180))/(360/res);
  end
end

conus_wet09 = glob_mask_cal .* wet09' .* grid_area;
conus_wet10 = glob_mask_cal .* wet10' .* grid_area;

conus_isam09 = glob_mask_cal .* isam09' .* grid_area;
conus_isam10 = glob_mask_cal .* isam10' .* grid_area;
conus_isam45 = glob_mask_cal .* isam45' .* grid_area;
conus_isam85 = glob_mask_cal .* isam85' .* grid_area;

