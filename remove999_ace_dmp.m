function [ dmpstruct_out ] = remove999_ace_dmp( dmpstruct_in )
%A function to replace the -999 values in the ACE DMP data with NaNs.
%
% *INPUT*
%           dmpstruct_in: STRUCTURE - a .MAT structure containing ACE
%           dmp data. It is usually created using 'read_ace_dmp' or
%           'read_ace_dmp_for_mat'.          
%
% *OUTPUT*
%           dmpstruct_out: STRUCTURE - output has the same fields as the
%           input, but with the -999 data replaced with nans
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 10/2017

%% Define some things
dmpin = dmpstruct_in;
dmpout = dmpin;

%% replace the -999 values with NaNs
dmpout.T(dmpout.T == -999) = nan;
dmpout.pressure_hPa(dmpout.pressure_hPa == -999) = nan;
if isfield(dmpout,'altitude_km')
    dmpout.altitude_km(dmpout.altitude_km == -999) = nan;
end
dmpout.lon(dmpout.lon == -999) = nan;
dmpout.lat(dmpout.lat == -999) = nan;
dmpout.Theta(dmpout.Theta == -999) = nan;
dmpout.spv(dmpout.spv == -999) = nan;
dmpout.eql(dmpout.eql == -999) = nan;
%%
dmpstruct_out = dmpout;
end

