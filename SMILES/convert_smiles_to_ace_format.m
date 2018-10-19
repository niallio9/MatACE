function [ tanstruct ] = convert_smiles_to_ace_format( smilesstruct )
%A function to read the ACE v3.5/6 .nc data and output a structure with
%the information.

% *INPUT*    
%           smilesstruct: STRUCTURE - the name of the .nc file that holds the
%           MLS data. Can be created with 'extract_mls_data.m'.
%
% *OUTPUT*
%           smiles_out: STRUCTURE - with fields that correspond to the fields
%           of the standardised ACE structure ('made with
%           read_ace_ncdata.m'). Unnecessary 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 06/18

%% Define some things
smiles = smilesstruct;
if isfield(smiles,'occultation')
   fprintf('The MLS structure already contains occultation information.\nNothing was done to the format.\n\n')
   tanstruct = smiles;
   return
end
lsmiles = length(smiles.date_mjd);
lalt = length(smiles.altitude_km);

%% set up the .mat file
out.source_file = smiles.source_file;
out.occultation = nan(1,lsmiles);
out.sr1ss0 = nan(1,lsmiles);
out.beta_angle = nan(1,lsmiles);
out.date_mjd = smiles.date_mjd;
out.gas = smiles.gas;
out.altitude_km = smiles.altitude_km; % doing this for now, in case altitude_km is actualy needed. I think it is - NJR 04/18
out.vmr = smiles.vmr;
out.vmr_error = smiles.vmr_error;
out.lat_tangent = smiles.lat;
out.lon_tangent = smiles.lon;
out.quality_flags = nan(lalt,lsmiles);
out.pressure_hPa = smiles.pressure_hPa; % need the number or pressure profiles to match the number of occultations
out.lat = repmat(smiles.lat,lalt,1); % want these for now so that it mimics the inclusion of GLC information in the ACE data structures.
out.lon = repmat(smiles.lon,lalt,1);

%% add flags for when there are negative precisions
out.quality_flags(out.vmr_error < 0) = 2; % this will indicate a bad value

%% make the output file
tanstruct = out; % the naming of the structure should be the same too, so that it can be used in the ACE functions.

%
end

