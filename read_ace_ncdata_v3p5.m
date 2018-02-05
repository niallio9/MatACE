function [ ace_out, gasname ] = read_ace_ncdata_v3p5( filename )
%A function to read the ACE v3.5/6 .nc data and output a structure with
%the information.

% *INPUT*    
%           filename: STRING - the name of the .nc file that holds the ace
%           data
%
% *OUTPUT*
%           ace_out: MATLAB STRUCTURE with the following fields of ACE
%           data: filename, gas, occultation, beta_angle, sr1ss0,
%           altitude_km, date_mjd, lat_tangent, lon_tangent, vmr,
%           vmr_error, quality_flags, pressure_hPa. 
%
%           gasname: STRING - the name of the gas to which the data
%           corresponds 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 10/17

%% Define some things
filein = filename;
gasname = dir(filein); % gets the details of the .nc file in the directory
gasname = cell2mat({gasname.name}); % gets the name of the .nc file in a string format
gasname = gasname(16:end-3); % the filename is of the type: ACEFTS_L2_v3p5_O3.nc
gaserror = strcat(gasname,' error');

%% set up the .mat file
out.source_file = filein;
out.occultation = ncread(filein,'orbit');
out.occultation = out.occultation';
out.sr1ss0 = ncread(filein,'ss (0) or sr (1)');
out.sr1ss0 = out.sr1ss0';
out.beta_angle = ncread(filein,'beta angle');
out.beta_angle = out.beta_angle';
%the date info is converted to a double so that the conversion to MJD works
out.date_mjd = date2mjd(double(ncread(filein,'year')),double(ncread(filein,'month')),double(ncread(filein,'day')),double(ncread(filein,'hour')));
out.date_mjd = out.date_mjd';
out.gas = gasname;
out.altitude_km = double(ncread(filein,'altitude'));
out.vmr = double(ncread(filein,gasname));
out.vmr_error = double(ncread(filein,gaserror));
out.lat_tangent = double(ncread(filein,'latitude'));
out.lat_tangent = out.lat_tangent';
out.lon_tangent = double(ncread(filein,'longitude'));
out.lon_tangent = out.lon_tangent';
out.quality_flags = ncread(filein,'quality flag');
% there is no pressure information in the v3p5 files so it must be added
% later
%out.pressure_hPa = double(ncread(filein,'pressure')) * 1013.25; % convert from atm to hPa

%% make the output file
ace_out = out;

%
end

