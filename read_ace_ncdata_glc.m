function [ glcstruct_out, gasname ] = read_ace_ncdata_glc( filename )
%A function to read the ACE GLC v3.5/6 .nc data and output a structure with
%relevant information.

% *INPUT*    
%           filename: STRING - the name of the .nc file that holds the ace
%           data
%
% *OUTPUT*
%           ace_out: MATLAB STRUCTURE with the following fields of ACE
%           data: filename, gas, occultation, beta_angle, srss,
%           altitude_km, date_mjd, tangent_lat, tangent_lon, lon, lat,
%           temperature, pressure_hPa.
%
%           gasname: STRING - this string is a constant: 'GLC'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 11/17

%% Define some things
filein = filename;
gasname = dir(filein); % gets the details of the .nc file in the directory
gasname = cell2mat({gasname.name}); % gets the name of the .nc file in a string format
gasname = gasname(16:end-3); % the filename is of the type: ACEFTS_L2_v3p6_O3.nc

%% set up the .mat file
out.source_file = filein;
out.occultation = ncread(filein,'orbit');
out.occultation = out.occultation';
out.sr1ss0 = ncread(filein,'sunset_sunrise');
out.sr1ss0 = out.sr1ss0';
out.beta_angle = ncread(filein,'beta_angle');
out.beta_angle = out.beta_angle';
%the date info is converted to a double so that the conversion to MJD works
out.date_mjd = date2mjd(double(ncread(filein,'year')),double(ncread(filein,'month')),double(ncread(filein,'day')),double(ncread(filein,'hour')));
out.date_mjd = out.date_mjd';
out.gas = gasname;
out.altitude_km = double(ncread(filein,'altitude'));
out.lon = double(ncread(filein,'GLC_longitude'));
out.lat = double(ncread(filein,'GLC_latitude'));
out.temperature = double(ncread(filein,'temperature'));
out.lat_tangent = double(ncread(filein,'latitude'));
out.lat_tangent = out.lat_tangent';
out.lon_tangent = double(ncread(filein,'longitude'));
out.lon_tangent = out.lon_tangent';
out.pressure_hPa = double(ncread(filein,'pressure')) * 1013.25; % convert from atm to hPa
out.quality_flag = ncread(filein,'quality_flag');

%% make the output file
glcstruct_out = out;

%
end

