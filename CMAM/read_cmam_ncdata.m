function [cmam_out, gasname] = read_cmam_ncdata(filename_in, latminmax, lonminmax)
%A function to read the CMAM30 .nc data and output a structure with
%the information.

% *INPUT*    
%           filename: STRING - the name of the .nc file that holds the cmam
%           data
%
%           latminmax: VECTOR - the minu=imum and maximum of the latitude
%           range you want to read. [latmin, latmax]
%
%           lonminmax: VECTOR - the minu=imum and maximum of the latitude
%           range you want to read. [lonmin, lonmax]
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
%latitude
latcmam = ncread(cfile,'lat')';
%longitude
loncmam = ncread(cfile,'lon')';
cout.source_file = filename_in;
cout.date_mjd = cmam2mjd(ncread(cfile,'time'))';
% cout.gas = vmrgas(4:end); %just take the gas name and not the 'vmr' preface
cout.gas = ncreadatt(cfile,vmrgas,'long_name');
if nargin == 1
    cout.vmr = ncread(cfile,vmrgas);
    cout.lat = latcmam;
    cout.lon = loncmam;
else
    %We will work from south to north, so make sure latend is larger than
    %latstart. Swap them around if not. Southern lats are negative.
    if latminmax(2) < latminmax(1)
        latminmax = fliplr(latminmax);
    end
    if lonminmax(2) < lonminmax(1)
        lonminmax = fliplr(lonminmax);
    end
    [~, ilatmin] = min(abs(latcmam - latminmax(1)))
    [~, ilatmax] = min(abs(latcmam - latminmax(2)))
    [~, ilonmin] = min(abs(loncmam - lonminmax(1)))
    [~, ilonmax] = min(abs(loncmam - lonminmax(2)))
    latout = latcmam(ilatmin:ilatmax);
    lonout = loncmam(ilonmin:ilonmax);
    cout.vmr = ncread(cfile,vmrgas,[ilonmin, ilatmin, 1, 1], [length(lonout), length(latout), Inf, Inf]);
    cout.lat = latout;
    cout.lon = lonout;
    
end
cout.pressure_hPa = ncread(cfile, 'plev')/100;
disp('done')
%
cmam_out = cout;
gasname = cout.gas;

end









