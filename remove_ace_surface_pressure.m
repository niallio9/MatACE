function [ tanstruct_out ] = remove_ace_surface_pressure(tanstruct_in)
%A function to find ACE pressure measurements that are surface pressure
%values, and remove those measurement points for all data. This is because
%there seems to be multiple entries of surface pressure (1013.25 hPa). 

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same fields as the
%           input, but with the surface pressure values (1013.25 hPa)
%           replaced by NaNs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017

%% Define some things
datain = tanstruct_in;
out = tanstruct_in;

%% Change 1013.25hPa (1atm) pressure values to NaNs
[badI] = find(datain.pressure_hPa == 1013.25);
out.pressure_hPa(badI) = nan;
out.vmr(badI) = nan;
out.vmr_error(badI) = nan;
if length(out.altitude_km(1,:)) > 1 % for the data that has been interpolated to a pressure grid
    out.altitude_km(badI) = nan;
    fprintf('\nWarning: The input structure appears to have already been interpolated. Applying this filter here is a bit iffy...\n')
end
if isfield(out,'lon')
   out.lon(badI) = nan;
   out.lat(badI) = nan;
end

%%
tanstruct_out = out;
%
end

