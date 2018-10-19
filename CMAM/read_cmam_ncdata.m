function [cmam_out, gasname] = read_cmam_ncdata(filename_in)
%A function to read the CMAM30 .nc data and output a structure with
%the information.

% *INPUT*    
%           filename: STRING - the name of the .nc file that holds the cmam
%           data
%
% *OUTPUT*
%           cmam_out: MATLAB STRUCTURE with the following fields of cmam
%           data: source, gas, date_mjd, lat, lon, vmr, pressure_hPa. 
%
%           gasname: STRING - the name of the gas to which the data
%           corresponds 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 01/17
%   NJR - 04/18. changed the reding of 'gas' variable

cfile = filename_in;
[~, filename, extension] = fileparts(cfile);
cfile = strcat(filename,extension);
vmrgas = cfile(1:end-55); % this is for the cmam30 naming format of e.g., 'vmrclo_6hrChem_CMAM_CMAM30-SD_r1i1p1_2008010100-2008123118.nc'
disp('reading CMAM variables...')
cout.source = filename_in;
% cout.gas = vmrgas(4:end); %just take the gas name and not the 'vmr' preface
cout.gas = ncreadatt(cfile,vmrgas,'long_name');
cout.vmr = ncread(cfile,vmrgas);
cout.pressure_hPa = ncread(cfile, 'plev')/100;
cout.date_mjd = cmam2mjd(ncread(cfile,'time'));
cout.lat = ncread(cfile, 'lat');
cout.lon = ncread(cfile, 'lon');
disp('done')
%
cmam_out = cout;
gasname = cout.gas;

end









