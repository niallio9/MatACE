function [ tanstruct ] = convert_mls_to_ace_format( mlsstruct )
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
mls = mlsstruct;
if isfield(mls,'occultation')
   fprintf('The MLS structure already contains occultation information.\nNothing was done to the format.\n\n')
   tanstruct = mls;
   return
end
lmls = length(mls.date_mjd);
lalt = length(mls.pressure_hPa);

%% set up the .mat file
out.source_file = mls.source_file;
out.occultation = nan(1,lmls);
out.sr1ss0 = nan(1,lmls);
out.beta_angle = nan(1,lmls);
out.date_mjd = mls.date_mjd;
out.gas = mls.gas;
out.altitude_km = p2z_waccm(mls.pressure_hPa * 100); % doing this for now, in case altitude_km is actualy needed. I think it is - NJR 04/18
out.vmr = mls.vmr;
out.vmr_error = mls.vmr_error;
out.lat_tangent = mls.lat;
out.lon_tangent = mls.lon;
out.quality_flags = nan(lalt,lmls);
out.pressure_hPa = repmat(mls.pressure_hPa, 1, lmls); % need the number or pressure profiles to match the number of occultations
out.lat = repmat(mls.lat,lalt,1); % want these for now so that it mimics the inclusion of GLC information in the ACE data structures.
out.lon = repmat(mls.lon,lalt,1);
if isfield(mls,'spv')
    out.spv = mls.spv;
end
if isfield(mls,'theta')
    out.theta = mls.theta;
end
if isfield(mls,'eql')
    out.eql = mls.eql;
end
%% make the output file
tanstruct = out; % the naming of the structure should be the same too, so that it can be used in the ACE functions.

%
end

