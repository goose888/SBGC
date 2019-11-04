function WCH4=WETCHARTS_CODE(EE)
%This code is used to generate the WetCHARTs wetland CH4 emission dataset 
%presented in Bloom et al., (2016).
%
%****USER INSTRUCTIONS****
% Before running the following code, make sure to:
%1. Update filepaths (marked as "[TO DO]" items)
%2. Obtain and rename ancillary dataset files  (marked as "[DOWNLOAD]" items)
%For additional information Contact: Anthony Bloom (abloom@jpl.nasa.gov)
%Copyright (c) California Institute of Technology, 2017
%
%Reference: Bloom, A. A., et al.: A global wetland methane emissions and uncertainty dataset for atmospheric chemical transport models, Geosci. Model Dev. Discuss., doi:10.5194/gmd-2016-224, in review, 2016. 
%
% Copyright 2017, by the California Institute of Technology. ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.
% This software may be subject to U.S. export control laws. By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. User has the responsibility to obtain export licenses, or other export authority as may be required before exporting such information to foreign countries or providing access to foreign persons.

if EE==1;    disp('Extended ensemble');end
if EE==0;    disp('Full ensemble');end


%Wetland CH4 emission derivation
WCH4=Wetland_CH4_emissions(EE);





end

%Wetland CH4 emission derivation
function WCH4=Wetland_CH4_emissions(EE)
%For each option define start:end years
styear=[2009,2001];sty=styear(EE+1);
endyear=[2010,2015];eny=endyear(EE+1);
%resolution
res=0.5;
%x-y mesh grid
[x,y,A]=loadworldxygrid(res);

%********************************************

% load heterotrophic respiration data
 RHET=load_heterotrophic_respiration_data(EE,sty,eny);

%***********************************
% derive & load wetland  inland water body extent parameterizations
 WEXT=load_extent_parameterizations(EE,sty,eny);

%********************************************
%derive & load CH4:C temperature dependence factor 
 TEMP=load_temp_ch4_c_dependence_factor(sty,eny);
%********************************************
TOTALCH4=[124.5,166,207.5];%Global totals for 2009-2010 emissions
%Configurations 
noii=size(WEXT,4);%Number of extent parameterizations
nott=size(TEMP,4);%Number of temperature parameterizations
nomm=size(RHET,4);%Number of heterotrophic respiration parameterizations
nosf=numel(TOTALCH4);%Number of scale factors

%Define WCH4 structure
WCH4.fluxname='Wetland CH4 emissions';
WCH4.units='mg CH4 m-2 day-1';
WCH4.month=mod([1:(eny-sty+1)*12]-1,12)+1;
WCH4.year=sty+floor(((1:(eny-sty+1)*12)-1)/12);
%Declare wetland emission array for efficiency
WCH4.data=zeros(360,720,numel(WCH4.month),nosf*nomm*nott*noii);
versions={'Full','Extended'};
WCH4.version=versions{EE+1};

save inund_extent.mat

%Step 3. Derive CH4 emissions
for s=1:nosf; %global scale factor loop
for mm=1:nomm; %heterotrophic respiration loop    
    for tt=1:nott;%3 temperature dependence loop
        for ii=1:noii; %extent parameterization loop
        %model index
        modix=(s-1)*nomm*nott*noii+(mm-1)*nott*noii+(tt-1)*noii+ii;
        %Wetland emissions
        WCH4.data(:,:,:,modix)=TEMP(:,:,:,tt).*RHET(:,:,:,mm).*WEXT(:,:,:,ii).*16/12*1e3;
        %Scale factor based on total 2009-2010 emissions
        SF=TOTALCH4(s)/sum(sum(mean(WCH4.data(:,:,WCH4.year>=2009 & WCH4.year<=2010,modix),3).*A*365.25/1e15));
        WCH4.data(:,:,:,modix)=WCH4.data(:,:,:,modix)*SF;
        WCH4.model_configuration(modix,:)=[modix,s,mm,tt,ii];
        WCH4.scalefactors(modix)=SF;
        end
    end    
end
end

%Storing additional information on ensemble CH4 emissions structure
if EE==0;
WCH4.configuration_info.configuration_matrix_columns={'Configuration number','Global Total','Heterotrophic Respiration','CH4:C Temperature dependence','Extent Parameterization'};
WCH4.configuration_info.globaltotals={'124.5 Tg/yr','166 Tg/yr','207.5 Tg/yr'};
WCH4.configuration_info.heterotrophic_respiration={'MsTMIP 1','MsTMIP 2','MsTMIP 3','MsTMIP 4','MsTMIP 5','MsTMIP 6','MsTMIP 7','MsTMIP 8','CARDAMOM'};
WCH4.configuration_info.ch4_temp_dependence={'CH4:C q10 = 1','CH4:C q10 = 2','CH4:C q10 = 3'};
WCH4.configuration_info.wetland_extent_parameterization={'SWAMPS & GLWD','SWAMPS & GLOBCOVER','PREC & GLWD','PREC & GLOBCOVER'};
elseif EE==1;
WCH4.configuration_info.configuration_matrix_columns={'Configuration number','Global Total','Heterotrophic Respiration','CH4:C Temperature dependence','Extent Parameterization'};
WCH4.configuration_info.globaltotals={'124.5 Tg/yr','166 Tg/yr','207.5 Tg/yr'};
WCH4.configuration_info.heterotrophic_respiration={'CARDAMOM'};
WCH4.configuration_info.ch4_temp_dependence={'CH4:C q10 = 1','CH4:C q10 = 2','CH4:C q10 = 3'};
WCH4.configuration_info.wetland_extent_parameterization={'PREC & GLWD','PREC & GLOBCOVER'};
end


end



%*********heterotrophic respiration datasets
function RHET=load_heterotrophic_respiration_data(EE,sty,eny)

    %Loading CARDAMOM 1x1 degree respiration dataset
    CARDAMOM=read_cardamom_respiration_dataset;
    %Nearest point interpolation to 0.5x0.5 degree resolution
    RhCAR05=zeros(360,720,numel(CARDAMOM.year));
    for r=1:2
        for c=1:2
           RhCAR05(r:2:end,c:2:end,:)=CARDAMOM.RHEmedian;
        end
    end

if EE==0
    %Loading 2001-2010 MsTMIP data
    MHR=read_mstmip_datasets; RHET=MHR.mstmipdata;
    %Adding CARDAMOM to the heterotrophic respiration datasets
    RHET(:,:,:,9)=RhCAR05(:,:,CARDAMOM.year<2011);
    %Years corresponding to each layer
    RHyears=MHR.year;
elseif EE==1;
    %CARDAMOM only
    RHET=RhCAR05;
    %Years corresponding to each layer
    RHyears=CARDAMOM.year;
end
    
%Final Heterotrophic respiration data ensemble
RHET=RHET(:,:,RHyears>=sty & RHyears<=eny,:);


end

%Wetland and inland water body extent parameterization
function INUNPREC=load_extent_parameterizations(EE,sty,eny)
%Wetland and inland water body extent parameterization
%Read precipitation dataset
EPREC=read_erai_precipitation;
%SWAMPS/MEaSUREs Surface inundation extent
MIA=read_swamps_surface_inundation_data;
%Load land-sea mask
LSM=loadlsmask;
%Load Globcover land cover data
GLOB=read_globcover_data;
%Load GLWD data
GLWD=read_glwd_data;

%Flooded, waterlogged soil and inland water body landcover extent.
WFRAC=sum(GLOB.frac(:,:,[16:18,21]),3);WFRAC(LSM==0)=0;
%wetlands and inland water body extent
MAXWETEXTENT=sum(GLWD.data(:,:,2:end),3);
%Scaling Iundated area extent maximum to maximum W&IWB extent
IAdata(:,:,:,1)=MIA.data.*repmat(MAXWETEXTENT./max(MIA.data,[],3),[1,1,size(MIA.data,3)]);
%Scaling mean wetland extent to Globcover W&IWB extent.
IAdata(:,:,:,2)=MIA.data.*repmat(WFRAC./mean(MIA.data,3),[1,1,size(MIA.data,3)]);
%Repeat for precipitation proxy
%Scaling maximum precipitation to maximum wetland extent
PREC(:,:,:,1)=EPREC.data.*repmat(MAXWETEXTENT./max(EPREC.data,[],3),[1,1,size(EPREC.data,3)]);
%IA.data(:,:,:,2)=MIA.data.*repmat(max(WFRAC./mean(MIA.data,3),1),[1,1,size(MIA.data,3)]);
%Scaling mean precipitation extent to Globcover mean wetland extent.
PREC(:,:,:,2)=EPREC.data.*repmat(WFRAC./mean(EPREC.data,3),[1,1,size(EPREC.data,3)]);
%For GLOBCOVER: setting maximum scaled extent to 1.
 PREC(PREC>1)=1;
 IAdata(IAdata>1)=1;
  IAdata(isnan(IAdata))=0;

if EE==0;
    for ii=1:2;INUNPREC(:,:,:,ii)=IAdata(:,:,MIA.year>=sty & MIA.year<=eny,ii);end
    INUNPREC(:,:,:,3:4)=PREC(:,:,EPREC.year>=sty & EPREC.year<=eny,:);
elseif EE==1
     INUNPREC(:,:,:,1:2)=PREC(:,:,EPREC.year>=sty & EPREC.year<=eny,:);
end
    
%Replacing empty values with zeros
INUNPREC(isnan(INUNPREC))=0;


end

%Temperature CH4:C dependence parameterization

function TEMP=load_temp_ch4_c_dependence_factor(sty,eny);
ESKT=read_erai_surface_skin_temp;
Q10ch4co2=[1,2,3];
eraitemp=ESKT.data(:,:,ESKT.year>=sty & ESKT.year<=eny);
for tt=1:3;    TEMP(:,:,:,tt)=(Q10ch4co2(tt)).^((eraitemp-10)/10);end

end






%read monthly ERA-interim precipitation
function ERAI=read_erai_precipitation
%[DOWNLOAD] How to obtain ERA-interim precipitation dataset:
%Step 1. go to http://apps.ecmwf.int/datasets/data/interim-full-mnth/levtype=sfc/
%Step 2. For each month and year, select the following: 
%parameter = "Total Precipitation"; step = "12"; time = "0:00:00" and "12:00:00"
%Step 3. "Retrieve Netcdf"
%Step 4. Change "grid" to "0.25x0.25"
%Step 5. Select "Retrieve now"
%Step 6. Rename each file in the following format for use here ERAI_025x025_prec_YYYY_MM.nc
%(YYYY = year; MM = month, e.g. "ERAI_025x025_prec_2001_01.nc")
%****Path name****
%[TO DO] Change pathname here:
pathname='./ERAI_PREC/';
quantity='prec';

fieldname='tp';
units='m/day';
conversion='*2';
minmax=0;
%find quantity 
resstr='025x025';

c=1;
for yr=2001:2015
    for m=1:12
        disp([yr,m]);
        filename=sprintf('%sERAI_%s_%s_%4i_%02i.nc',pathname,resstr,quantity,yr,m);
        ALL_DATA=grid_erai_at_half_degree_res(ncread(filename,fieldname));
        %Implement correction consistent with defined units
        eval(sprintf('ALL_DATA=ALL_DATA%s;',conversion));
        %Define output structure
        ERAI.fieldname=fieldname;
        ERAI.units=units;
        ERAI.data(:,:,c)=mean(ALL_DATA,3);
        ERAI.date(c)=datenum(sprintf('%02i/01/%i',m,yr))+365.25/24;
        ERAI.month(c)=m;
        ERAI.year(c)=yr;
        c=c+1;
    end
end
        
%Precipitation = positive-definite
%ERAI netcdf conversion leads to ~0 (including negligible -ve values)
%Setting all <0 values to 0
ERAI.data(ERAI.data<0)=0;



end

function ERAI=read_erai_surface_skin_temp
%[DOWNLOAD] How to obtain ERA-interim surface skin temperature dataset:
%Step 1. go to http://apps.ecmwf.int/datasets/data/interim-full-mnth/levtype=sfc/
%Step 2. For each month and year, select the following: 
%parameter = "Skin temperature"; step = "0"; time = "0:00:00" , "6:00:00" ,"12:00:00" and "18:00:00"
%Step 3. "Retrieve Netcdf"
%Step 4. Change "grid" to "0.25x0.25"
%Step 5. Select "Retrieve now"
%Step 6. Rename each file in the following format for use here ERAI_025x025_skt_YYYY_MM.nc
%(YYYY = year; MM = month, e.g. "ERAI_025x025_skt_2001_01.nc")
%[TO DO] Change pathname here:
pathname='./ERAI_SKT/';
quantity='skt';
fieldname='skt';
units='C';
conversion='-273.15';
%find quantity 
resstr='025x025';

%c = loop constant
c=1;
for yr=2001:2015
    for m=1:12
        disp([yr,m]);
        filename=sprintf('%sERAI_%s_%s_%4i_%02i.nc',pathname,resstr,quantity,yr,m);
        ALL_DATA=grid_erai_at_half_degree_res(ncread(filename,fieldname));
        %Implement correction consistent with defined units
        eval(sprintf('ALL_DATA=ALL_DATA%s;',conversion));
        %Define output structure
        ERAI.fieldname=fieldname;
        ERAI.units=units;
        ERAI.data(:,:,c)=mean(ALL_DATA,3);
        ERAI.date(c)=datenum(sprintf('%02i/01/%i',m,yr))+365.25/24;
        ERAI.month(c)=m;
        ERAI.year(c)=yr;
        c=c+1;
    end
end
        





end

function DOUT=grid_erai_at_half_degree_res(DIN)
%Re-grids era interim fields 
%Produces 0.5 x 0.5 dataset spanning LON=179.75W-179.75E, LAT=89.75S-89.75S
%current coordinates are LON=0W-359.5E, LAT=90N-90S
%fipping dimensions
DIN=permute(DIN,[2,1,3]);
%flipping upside down
DIN=DIN(end:-1:1,:,:);
%shifting to eurocentric
DIN=DIN(:,[end/2+1:end,1:end/2],:);
%repeating first column
DIN=[DIN,DIN(:,1,:)];
%Now coordinates are LON=180W-180E, LAT=90S-90N
%weighted averaging of dataset to 0.5x0.5 degrees
DIN2=DIN(1:2:end-2,:,:)+DIN(2:2:end-1,:,:)*2+DIN(3:2:end,:,:);
DIN3=DIN2(:,1:2:end-2,:)+DIN2(:,2:2:end-1,:)*2+DIN2(:,3:2:end,:);
DOUT=DIN3/16;
end


function GLOB=read_globcover_data
 %Original 300m GLOBCOVER dataset obtained from http://due.esrin.esa.int/page_globcover.php
 %The fractional cover for each type is calculated at 0.5 x 0.5 degree resolution.
 load WETCHARTS_AUXILIARY_DATA.mat GLOB
end

function GLWD=read_glwd_data
 %/GLWD-level3: Original dataset downloaded from http://gcmd.gsfc.nasa.gov
 %The fractional cover for each type calculated at 0.5 x 0.5 degree resolution.
 load WETCHARTS_AUXILIARY_DATA.mat GLWD
end


function MsTMIP=read_mstmip_datasets
%MsTMIP portal  http://nacp.ornl.gov/mstmipdata
%[DOWNLOAD] Download model monthly 2001-2010 "HeteroResp" datasets for models specified below
%[TO DO] Place files in directory specified by "pathname"
pathname='./MstMIP/';
MsTMIP.units='gC m-2 day-1';
MsTMIP.models={'BIOME-BGC','CLASS-CTEM-N','CLM4VIC','CLM4','DLEM','ISAM','TEM6','TRIPLEX-GHG'};
files=dir([pathname,'*','HeteroResp','.nc4']);

for m=1:numel(files)
    file=[pathname,files(m).name];
    %reading 'HeteroResp'
    %units: 'kg C m-2 s-1'
    %coverting to gC m-2 day-1
    DATA=permute(ncread(file,'HeteroResp'),[2,1,3]);
    DATA=DATA*1e3*3600*24;
    %Time: units = days since 01/01/1700
    %Monthly timesteps centered on ~14th of each month
    TIME=ncread(file,'time')+datenum('12/31/1699');
    %using 2001-2010 data only
    dts=TIME>=datenum('01/01/2001');
    MsTMIP.mstmipdata(:,:,:,m)=DATA(:,:,dts);
end
    MsTMIP.date=double(TIME(dts));
for n=1:numel(MsTMIP.date);
    MsTMIP.year(n)=year(MsTMIP.date(n));
    MsTMIP.month(n)=month(year(MsTMIP.date(n)));
end

MsTMIP.mstmipdata(isnan(MsTMIP.mstmipdata))=0;

end

function CARDAMOM=read_cardamom_respiration_dataset
%Load CARDAMOM (including 2011-2015 extension) 
load WETCHARTS_AUXILIARY_DATA.mat CARDAMOM
end

%Produces an x-y meshgrid and grid-cell area
function [x,y,A]=loadworldxygrid(res)
%Produces an x-y meshgrid 
hafresx=res/2;
hafresy=res/2;
[x,y]=meshgrid(-180+hafresx:res:180-hafresx,-90+hafresy:res:90-hafresy);
%Grid-cell area
A=111111.111^2*(res*res)*cos(y*pi/180);
 
end

function LSM=loadlsmask
 %Land sea mask
 load WETCHARTS_AUXILIARY_DATA.mat LSM
end


%SWAMPS
function WET=read_swamps_surface_inundation_data
%[DOWNLOAD] SWAMPS fw_28_swe dataset obtainable from wetlands.jpl.nasa.gov; 
%[TO DO] keep original file names and update directories
dpath='/data/jain1/c/sshu3/SBGC/runs/methane/inund_map/fw_14_swe/v2_ascat/ml/';

res=0.5;
[x,y,A]=loadworldxygrid(res);

WET.year=[];
WET.month=[];
WET.date=[];
WET.data=[];
WET.original=[];
WET.x=x;
WET.y=y;
[xe,ye]=ease_grid_coordinates;
WET.x_ease=xe;
WET.y_ease=ye;

k=1;
for y=2009:2012
    for m=1:12
    %read binary file gridth 1383 columns and 586 rows
    files=dir(sprintf('%s%4i%02i*.DAT',dpath,y,m));
    
    if isempty(files)==0;
    %declaring fields
    MONTHLY_DATA=zeros(1383,586)';
    MONTHLY_COUNT=zeros(1383,586)';

    for n=1:numel(files)
    fd=fopen([dpath,files(n).name]); 
    %**************extract date************
    %Dates only entered as such for compatibility and plotting purposes.
    %Data aggregated per calendar month
    %reading floating point data
    DATA=reshape(fread(fd,inf,'real*4'),1383,586)';    
    %aggregating data on a monthy timescale        
    MONTHLY_DATA=MONTHLY_DATA+DATA.*(DATA>=0);
    %Counting data on a monthly timescale
    MONTHLY_COUNT(DATA>=0)= MONTHLY_COUNT(DATA>=0)+1;
    %closing file
    fclose(fd);
    end
    
    %Mean monthly inundated area
    WETframe=MONTHLY_DATA./MONTHLY_COUNT;WETframe(MONTHLY_COUNT==0)=0;
    WET.IA_easegrid(:,:,k)=WETframe;
    WET.data(:,:,k)=regrid_ease_map(WETframe,res);        %regridding to 0.5x0.5
    WET.month(k)=m;
    WET.year(k)=y; 
    WET.date(k)=datenum(sprintf('%02i/01/%4i',m,y))+365.25/24;
    %loop constant
    k=k+1;
    end
    end
end

WET.data(isfinite(WET.data)==0)=0;



end


%********EASE GRID functions*****************
function [mapxy]=regrid_ease_map(mapease,res)
%Interpolates EASE grid map to cartesian 0.5 x 0.5 degree resolution
[x,y]=loadworldxygrid(res);
[xe,ye]=ease_grid_coordinates;
xe=[xe-360,xe,xe+360];
ye=[ye,ye,ye];
mapease=[mapease,mapease,mapease];
mapxy=interp2(xe,ye,mapease,x,y);
end

function [x,y]=ease_grid_coordinates;
%Produces grid of x and y EASE grid coordinates
xres=360/1383;
xv=[-180+xres/2:xres:180-xres/2];
for r=1:586/2
    A0=25*25*(r-1)*1383;
    A=25*25*r*1383;
    lat=asin(1- A/(2*pi*6371^2));
    lat0=asin(1-A0/(2*pi*6371^2));
    yv(r)=mean([lat,lat0])*180/pi;
end
yv(numel(yv)+1:numel(yv)*2)=-yv(end:-1:1);
[x,y]=meshgrid(xv,yv);
end





