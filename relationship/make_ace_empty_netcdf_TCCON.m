function [ ] = make_ace_empty_netcdf_TCCON( structin, filename_out )
 %A function to create an empty ACE climatology netCDF file. The fields of
 %the file is in a format for use by TCCON.
 %Data from ACE .mat TCCON files can later be written to the .nc file
 %created here. Use 'make_ace_data_for_tccon.m' to make the data that will
 %be written.
 
% *INPUT*
%           structin: STRING - the structure containing the gas and DMP
%           information. It can be created using
%           'make_ace_data_for_tccon.m'.
%
%           filename_out: STRING - the name of the .nc file that will be
%           created to hold the data.
%
% *OUTPUT*
%           A .nc file will be created
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 11/18

%% define some things
filename = filename_out;
if exist(filename,'file') == 2
    fprintf('\nThere is already a file called %s. It will be replaced.', filename)
    delete(filename);
end
data = structin;
gas1 = 'n2o';
gas2 ='ch4';
gas3 = 'hf';
gas4 = 'temperature';

% llev = 48;
llev = length(data.altitude_km);
ltime = length(data.occultation);
fillval = '-999';
clear data

%% create the file and add the variable information
disp('creating the fields for the netcdf file...')
nccreate(filename,'occultation','Dimensions',{'time', ltime},'Datatype','single','Format','classic');
nccreate(filename,'sr1ss0','Dimensions',{'time', ltime},'Datatype','single','Format','classic');
nccreate(filename,'time','Dimensions',{'time', ltime},'Datatype','double','Format','classic');
nccreate(filename,'zlev','Dimensions',{'zlev', llev},'Datatype','double','Format','classic');
nccreate(filename,'plev','Dimensions',{'zlev', llev, 'time', ltime},'Datatype','double','Format','classic');
nccreate(filename,'lat','Dimensions',{'zlev', llev, 'time', ltime},'Datatype','double','Format','classic');
nccreate(filename,'lon','Dimensions',{'zlev', llev, 'time', ltime},'Datatype','double','Format','classic');
nccreate(filename,'eql','Dimensions',{'zlev', llev, 'time', ltime},'Datatype','double','Format','classic');
nccreate(filename,'spv','Dimensions',{'zlev', llev, 'time', ltime},'Datatype','double','Format','classic');
nccreate(filename,gas1,'Dimensions',{'zlev', llev, 'time', ltime},'Datatype','double','Format','classic');
nccreate(filename,strcat(gas1,'_error'),'Dimensions',{'zlev', llev, 'time', ltime},'Datatype','double','Format','classic');
nccreate(filename,gas2,'Dimensions',{'zlev', llev, 'time', ltime},'Datatype','double','Format','classic');
nccreate(filename,strcat(gas2,'_error'),'Dimensions',{'zlev', llev, 'time', ltime},'Datatype','double','Format','classic');
nccreate(filename,gas3,'Dimensions',{'zlev', llev, 'time', ltime},'Datatype','double','Format','classic');
nccreate(filename,strcat(gas3,'_error'),'Dimensions',{'zlev', llev, 'time', ltime},'Datatype','double','Format','classic');
nccreate(filename,gas4,'Dimensions',{'zlev', llev, 'time', ltime},'Datatype','double','Format','classic');
% nccreate(filename,strcat(gas4,'_error'),'Dimensions',{'zlev', llev, 'time', ltime},'Datatype','double','Format','classic');


% write the attributes
fprintf('\nwriting the attributes...')
%n2o
ncwriteatt(filename,gas1,'long_name', sprintf('volume mixing ratio of %s in air', gas1));
ncwriteatt(filename,gas1,'units', 'ppv');
ncwriteatt(filename,gas1,'FillValue', fillval);
ncwriteatt(filename,strcat(gas1,'_error'),'long_name', sprintf('volume mixing ratio of %s in air error', gas1));
ncwriteatt(filename,strcat(gas1,'_error'),'units', 'ppv');
ncwriteatt(filename,strcat(gas1,'_error'),'FillValue', fillval);
%ch4
ncwriteatt(filename,gas2,'long_name', sprintf('volume mixing ratio of %s in air', gas2));
ncwriteatt(filename,gas2,'units', 'ppv');
ncwriteatt(filename,gas2,'FillValue', fillval);
ncwriteatt(filename,strcat(gas2,'_error'),'long_name', sprintf('volume mixing ratio of %s in air error', gas2));
ncwriteatt(filename,strcat(gas2,'_error'),'units', 'ppv');
ncwriteatt(filename,strcat(gas2,'_error'),'FillValue', fillval);
%hf
ncwriteatt(filename,gas3,'long_name', sprintf('volume mixing ratio of %s in air', gas3));
ncwriteatt(filename,gas3,'units', 'ppv');
ncwriteatt(filename,gas3,'FillValue', fillval);
ncwriteatt(filename,strcat(gas3,'_error'),'long_name', sprintf('volume mixing ratio of %s in air error', gas3));
ncwriteatt(filename,strcat(gas3,'_error'),'units', 'ppv');
ncwriteatt(filename,strcat(gas3,'_error'),'FillValue', fillval);
%temperature
ncwriteatt(filename,gas4,'long_name', 'air temperature');
ncwriteatt(filename,gas4,'units', 'K');
ncwriteatt(filename,gas4,'FillValue', fillval);
% ncwriteatt(filename,strcat(gas4,'_error'),'long_name', 'air temperature error');
% ncwriteatt(filename,strcat(gas4,'_error'),'units', 'K');
% ncwriteatt(filename,strcat(gas4,'_error'),'FillValue', fillval);
%lat
ncwriteatt(filename,'lat','long_name', 'latitude');
ncwriteatt(filename,'lat','units', 'degrees north');
% ncwriteatt(filename,'lat','_CoordinateAxisType', 'Lat');
% ncwriteatt(filename,'lat','axis', 'Y')
%lat
ncwriteatt(filename,'lon','long_name', 'longitude');
ncwriteatt(filename,'lon','units', 'degrees east');
%eql
ncwriteatt(filename,'eql','long_name', 'equivalent latitude');
ncwriteatt(filename,'eql','units', 'degrees north');
%spv
ncwriteatt(filename,'spv','long_name', 'scaled potential vorticity');
ncwriteatt(filename,'spv','units', 's-1');
%plev
ncwriteatt(filename,'plev','long_name', 'pressure');
ncwriteatt(filename,'plev','units', 'hPa');
% ncwriteatt(filename,'plev','standard_name', 'air_pressure');
% ncwriteatt(filename,'plev','positive', 'down');
% ncwriteatt(filename,'plev','axis', 'Z');
% ncwriteatt(filename,'plev','_CoordinateAxisType', 'Pressure');
% ncwriteatt(filename,'plev','_CoordinateZisPositive', 'down');
%zlev
ncwriteatt(filename,'zlev','long_name', 'altitude');
ncwriteatt(filename,'zlev','units', 'kilometres');
%time 
ncwriteatt(filename,'time','long_name', 'time');
% ncwriteatt(filename,'time','standard_name', 'time');
% ncwriteatt(filename,'time','calendar', 'standard');
% ncwriteatt(filename,'time','axis', 'T');
ncwriteatt(filename,'time','units', 'days since 1950-01-01 00:00:00');
%occultation
ncwriteatt(filename,'occultation','long_name', 'occultation number');
ncwriteatt(filename,'occultation','units', 'ID_number');
%sr1ss0
ncwriteatt(filename,'sr1ss0','long_name', 'sunrise/sunset indicator');
ncwriteatt(filename,'sr1ss0','units', '0/2:sunset, 1/3:sunrise');

%global attributes
ncwriteatt(filename,'/','Experiment', 'ACE-FTS');
ncwriteatt(filename,'/','Version', '3.5/3.6');
ncwriteatt(filename,'/','Organisation', 'CSA');
ncwriteatt(filename,'/','Type_of_Data', 'MEASUREMENTS AND DERIVED METEOROLOGICAL PRODUCTS (DMPs)');
ncwriteatt(filename,'/','Measurement Platform', 'SATELLITE');
ncwriteatt(filename,'/','Name_of_Platform', 'SCISAT-1');
ncwriteatt(filename,'/','Measurement_PI_name', 'Prof. Kaley Walker');
ncwriteatt(filename,'/','DMP_source', 'GEOS5MERRA2');
ncwriteatt(filename,'/','DMP_processing_centre', 'MLS_SFC');
ncwriteatt(filename,'/','DMP_author_name', 'Gloria Manney / Luis Millan');
ncwriteatt(filename,'/','File_Creation_Time', datestr(now));
ncwriteatt(filename,'/','File_Modification_Time', datestr(now));
ncwriteatt(filename,'/','Author', 'Niall Ryan');
ncwriteatt(filename,'/','Institute', 'University of Toronto');
ncwriteatt(filename,'/','Address', '60 St George Street, M5S 1A7, ON, Toronto, CA');
ncwriteatt(filename,'/','Email', 'kwalker@atmosp.physics.utoronto.ca, nryan@atmosp.physics.utoronto.ca');
ncwriteatt(filename,'/','Project_id', 'TCCON apriori gas relationships');
ncwriteatt(filename,'/','Important_Comment', 'Please sign up as an ACE data user and a DMP user if using the information in this file');
ncwriteatt(filename,'/','Author', 'Niall Ryan');

fprintf('\nDone\n')

fileattrib(filename,'+w');
%
end

