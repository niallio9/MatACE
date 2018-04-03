function [ ] = make_ace_empty_netcdf_SDI( filename_out, gasname )
 %A function to create an empty ACE climatology netCDF file. The fields of
 %the file is in a standard format required for the SPARC Data Initiative
 %for a monthly mean zonal climatology time series. A serial monthly
 %climatology.
 %Data from ACE .mat climatology files can later be written to the .nc file
 %created here.
 
% *INPUT*
%           filename_out: STRING - the name of the .netcdf file that will
%           be created. This input can also contain a path to a directory
%           as well as the file name.
%
%           gasname: STRING - the name of the gas that the data is about.
%
% *OUTPUT*
%           A .nc file will be created
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 03/18

%% define some things
filename = filename_out;
if exist(filename,'file') == 2
    fprintf('\nThere is already a file called %s. It will be replaced.', filename)
    delete(filename);
end
gas = gasname;
% llev = 48;
llev = 28;
llat = 36;
ltime = 12;

%% create the file and add the variable information
fprintf('\ncreating the fields for %s...', gas)
nccreate(filename,gas,'Dimensions',{'lat', llat, 'plev', llev, 'time', ltime},'Datatype','double','Format','classic');
nccreate(filename,strcat(gas,'_STD'),'Dimensions',{'lat', llat, 'plev', llev, 'time', ltime},'Datatype','double','Format','classic');
nccreate(filename,strcat(gas,'_NR'),'Dimensions',{'lat', llat, 'plev', llev, 'time', ltime},'Datatype','int32','Format','classic');
nccreate(filename,'LST_MEAN','Dimensions',{'lat', llat, 'time', ltime},'Datatype','single','Format','classic');
nccreate(filename,'LST_MAX','Dimensions',{ 'lat', llat, 'time', ltime},'Datatype','single','Format','classic');
nccreate(filename,'LST_MIN','Dimensions',{ 'lat', llat, 'time', ltime},'Datatype','single','Format','classic');
nccreate(filename,'AVE_DOM','Dimensions',{ 'lat', llat, 'time', ltime},'Datatype','single','Format','classic');
nccreate(filename,'AVE_LAT','Dimensions',{ 'lat', llat, 'time', ltime},'Datatype','single','Format','classic');
nccreate(filename,'lat','Dimensions',{'lat', llat},'Datatype','double','Format','classic');
nccreate(filename,'plev','Dimensions',{'plev', llev},'Datatype','double','Format','classic');
nccreate(filename,'time','Dimensions',{'time', ltime},'Datatype','double','Format','classic');

% nccreate(filename,gas,'Dimensions',{'time', ltime, 'plev', llev,  'lat', llat},'Datatype','double','Format','classic');
% nccreate(filename,strcat(gas,'_STD'),'Dimensions',{'time', ltime, 'plev', llev,  'lat', llat},'Datatype','double','Format','classic');
% nccreate(filename,strcat(gas,'_NR'),'Dimensions',{'time', ltime, 'plev', llev,  'lat', llat},'Datatype','int32','Format','classic');
% nccreate(filename,'LST_MEAN','Dimensions',{'time', ltime, 'lat', llat},'Datatype','single','Format','classic');
% nccreate(filename,'LST_MAX','Dimensions',{'time', ltime, 'lat', llat},'Datatype','single','Format','classic');
% nccreate(filename,'LST_MIN','Dimensions',{'time', ltime, 'lat', llat},'Datatype','single','Format','classic');
% nccreate(filename,'AVE_DOM','Dimensions',{'time', ltime, 'lat', llat},'Datatype','single','Format','classic');
% nccreate(filename,'AVE_LAT','Dimensions',{'time', ltime, 'lat', llat},'Datatype','single','Format','classic');
% nccreate(filename,'lat','Dimensions',{'lat', llat},'Datatype','double','Format','classic');
% nccreate(filename,'plev','Dimensions',{'plev', llev},'Datatype','double','Format','classic');
% nccreate(filename,'time','Dimensions',{'time', ltime},'Datatype','double','Format','classic');

% write the attributes
fprintf('\nwriting the attributes...')
%gas
ncwriteatt(filename,gas,'long_name', sprintf('mixing ratio of %s in air', gas));
ncwriteatt(filename,gas,'cell_methods', 'lon: mean (zonal median), time: mean (of calendar month)');
ncwriteatt(filename,gas,'standard_name', sprintf('volume mixing ratio of %s in air', gas));
ncwriteatt(filename,gas,'units', 'ppv');
ncwriteatt(filename,gas,'FillValue', '-999');
%gas_STD
ncwriteatt(filename,strcat(gas,'_STD'),'long_name', sprintf('volume mixing ratio of %s in air standard deviation', gas));
ncwriteatt(filename,strcat(gas,'_STD'),'units', 'ppv');
ncwriteatt(filename,strcat(gas,'_STD'),'FillValue', '-999');
%gas_NR
ncwriteatt(filename,strcat(gas,'_NR'),'long_name', 'number of values');
ncwriteatt(filename,strcat(gas,'_NR'),'units', 'none');
ncwriteatt(filename,strcat(gas,'_NR'),'FillValue', '-999');
%LST_MEAN
ncwriteatt(filename,'LST_MEAN','long_name', 'mean of local solar time');
ncwriteatt(filename,'LST_MEAN','units', 'hours');
ncwriteatt(filename,'LST_MEAN','FillValue', '-999');
%LST_MAX
ncwriteatt(filename,'LST_MAX','long_name', 'maximum of local solar time');
ncwriteatt(filename,'LST_MAX','units', 'hours');
ncwriteatt(filename,'LST_MAX','FillValue', '-999');
%LST_MIN
ncwriteatt(filename,'LST_MIN','long_name', 'minimum of local solar time');
ncwriteatt(filename,'LST_MIN','units', 'hours');
ncwriteatt(filename,'LST_MIN','FillValue', '-999');
%AVE_DOM
ncwriteatt(filename,'AVE_DOM','long_name', 'average day of month');
ncwriteatt(filename,'AVE_DOM','units', 'days');
ncwriteatt(filename,'AVE_DOM','FillValue', '-999')
%AVE_LAT
ncwriteatt(filename,'AVE_LAT','long_name', 'average latitude');
ncwriteatt(filename,'AVE_LAT','units', 'degrees north');
ncwriteatt(filename,'AVE_LAT','FillValue', '-999')
%lat
ncwriteatt(filename,'lat','long_name', 'latitude');
ncwriteatt(filename,'lat','standard_name', 'latitude');
ncwriteatt(filename,'lat','units', 'degrees north');
% ncwriteatt(filename,'lat','_CoordinateAxisType', 'Lat');
ncwriteatt(filename,'lat','axis', 'Y')
%plev
ncwriteatt(filename,'plev','long_name', 'pressure');
ncwriteatt(filename,'plev','standard_name', 'air_pressure');
ncwriteatt(filename,'plev','positive', 'down');
ncwriteatt(filename,'plev','axis', 'Z');
ncwriteatt(filename,'plev','units', 'hPa');
% ncwriteatt(filename,'plev','_CoordinateAxisType', 'Pressure');
% ncwriteatt(filename,'plev','_CoordinateZisPositive', 'down');
%time 
ncwriteatt(filename,'time','long_name', 'time');
ncwriteatt(filename,'time','standard_name', 'time');
ncwriteatt(filename,'time','calendar', 'standard');
ncwriteatt(filename,'time','axis', 'T');
ncwriteatt(filename,'time','units', 'days since 1950-01-01 00:00:00');

%global attributes
ncwriteatt(filename,'/','Experiment', 'ACE-FTS');
ncwriteatt(filename,'/','Version', '3.5/3.6');
ncwriteatt(filename,'/','Organisation', 'CSA');
ncwriteatt(filename,'/','Type_of_Data', 'MEASUREMENTS');
ncwriteatt(filename,'/','Platform', 'SATELLITE');
ncwriteatt(filename,'/','Name_of_Platform', 'SCISAT-1');
ncwriteatt(filename,'/','PI_name', 'Prof. Kaley Walker');
ncwriteatt(filename,'/','File_Creation_Time', datestr(now));
ncwriteatt(filename,'/','File_Modification_Time', datestr(now));
ncwriteatt(filename,'/','Fields', 'T2Mz: Monthly-mean zonal mean 2-d atmosphere');
ncwriteatt(filename,'/','Climatology_version', 'i01');
ncwriteatt(filename,'/','Author', 'Niall Ryan');
ncwriteatt(filename,'/','Institute', 'University of Toronto');
ncwriteatt(filename,'/','Address', '60 St George Street, M5S 1A7, ON, Toronto, CA');
ncwriteatt(filename,'/','Email', 'kwalker@atmosp.physics.utoronto.ca,ajones@atmosp.physics.utoronto.ca');
ncwriteatt(filename,'/','Project_id', 'SPARC Data_Initiative');
ncwriteatt(filename,'/','Comment', '5 Points are required for an Monthly-Zonal-Average');
ncwriteatt(filename,'/','Scaling_of_data', 'None');
ncwriteatt(filename,'/','LST_level', '7 hPa'); %%%************************gotta check this********************************
ncwriteatt(filename,'/','AVE_DOM_Level', '7 hPa');
ncwriteatt(filename,'/','AVE_LAT_level', '7 hPa');
ncwriteatt(filename,'/','history', sprintf('%s:  Creating netCDF', datestr(now)));
ncwriteatt(filename,'/','Author', 'Niall Ryan');

fprintf('\nDone\n')

fileattrib(filename,'+w');
%
end

