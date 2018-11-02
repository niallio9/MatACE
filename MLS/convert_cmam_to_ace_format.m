function [ tanstruct ] = convert_cmam_to_ace_format( cmamstruct )
%A function to read the ACE v3.5/6 .nc data and output a structure with
%the information.

% *INPUT*    
%           mlsstruct: STRUCTURE - the name of the .nc file that holds the
%           MLS data. Can be created with 'extract_mls_data.m'.
%
% *OUTPUT*
%           mls_out: STRUCTURE - with fields that correspond to the fields
%           of the standardised ACE structure ('made with
%           read_ace_ncdata.m'). Unnecessary 
%
%           gasname: STRING - the name of the gas to which the data
%           corresponds 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 04/18

%% Define some things
cmam = cmamstruct;
if isfield(cmam,'occultation')
   fprintf('The MLS structure already contains occultation information.\nNothing was done to the format.\n\n')
   tanstruct = cmam;
   return
end
disp('line 25')
lt = length(cmam.date_mjd);
lalt = length(cmam.pressure_hPa);
llon = length(cmam.lon);
llat = length(cmam.lat);
disp('line30')
vmrout = nan(lalt,llon*llat*lt);
latout = nan(1,llon*llat*lt); % this is a little more than one years worth of global data
lonout = nan(1,llon*llat*lt);
dateout = nan(1,llon*llat*lt);
disp('line 35')

%% set up the .mat file
out.source_file = cmam.source_file;
out.occultation = nan(1,lt);
out.sr1ss0 = nan(1,lt);
out.beta_angle = nan(1,lt);
out.date_mjd = nan; % a placeholder
out.gas = cmam.gas;
out.altitude_km = p2z_waccm(cmam.pressure_hPa * 100); % doing this for now, in case altitude_km is actualy needed. I think it is - NJR 04/18
disp('line45')
%% the vmr array
vmrcount = 0;
for k = 1:lt % go through times
    k
    %go through lons
    for i = 1:llon
        for j = 1:llat % go through lats
            vmrcount = vmrcount+1;
            vmrout(:,vmrcount) = cmam.vmr(i,j,:,k);
            dateout(vmrcount) = cmam.date_mjd(k);
            lonout(vmrcount) = cmam.lon(i);
            latout(vmrcount) = cmam.lat(j);    
        end
    end  
end
disp('loop done')
out.vmr = vmrout;
out.vmr_error = nan(size(out.vmr));
out.lat_tangent = latout;
out.lon_tangent = lonout;
out.date_mjd = dateout;
out.quality_flags = nan(lalt,lt);
out.pressure_hPa = repmat(cmam.pressure_hPa, 1, lt); % need the number or pressure profiles to match the number of occultations
out.lat = repmat(cmam.lat,lalt,1); % want these for now so that it mimics the inclusion of GLC information in the ACE data structures.
out.lon = repmat(cmam.lon,lalt,1);

%% make the output file
tanstruct = out; % the naming of the structure should be the same too, so that it can be used in the ACE functions.

%
end

