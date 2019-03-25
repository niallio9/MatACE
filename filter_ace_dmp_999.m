function [ dmpstruct_out ] = filter_ace_dmp_999( dmpstruct_in )
%A function to replace the -999 values in the ACE DMP data with NaNs. Will
%also replace the fill values of ... x 1e15.
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
% NJR - 10/2018 replaces new fill-values in the latest DMP data

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
dmpout.theta(dmpout.theta == -999) = nan;
dmpout.spv(dmpout.spv == -999) = nan;
dmpout.eql(dmpout.eql == -999) = nan;
%%
%% replace the 1e15 values with NaNs
dmpout.T(dmpout.T > 1e12) = nan;
dmpout.pressure_hPa(dmpout.pressure_hPa > 1e12) = nan;
if isfield(dmpout,'altitude_km')
    dmpout.altitude_km(dmpout.altitude_km > 1e12) = nan;
end
dmpout.lon(dmpout.lon > 1e12) = nan;
dmpout.lat(dmpout.lat > 1e12) = nan;
dmpout.theta(dmpout.theta> 1e12) = nan;
dmpout.spv(dmpout.spv > 1e12) = nan;
dmpout.eql(dmpout.eql > 1e12) = nan;
% tropopause data
dmpout.tropopauses_hPa_WMO(dmpout.tropopauses_hPa_WMO > 1e12) = nan;
dmpout.tropopauses_km_WMO(dmpout.tropopauses_km_WMO > 1e12) = nan;
dmpout.tropopauses_hPa_Dyn(dmpout.tropopauses_hPa_Dyn > 1e12) = nan;
dmpout.tropopauses_km_Dyn(dmpout.tropopauses_km_Dyn > 1e12) = nan;

%
dmpstruct_out = dmpout;
end

