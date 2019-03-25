function [ tanstruct_n2o_out ] = update_ace_n2o_flags( tanstruct_n2o_in, tanstruct_ch4_in, tanstruct_T_in )
%A function to catch some outliers in the N2O data that are likely caused
%by ACE looking through polar stratospheric clouds (from Patrick)

% *INPUT*
%           tanstruct_X_in: STRUCTURE - contains the gas specific ACE data,
%           for the gas 'X'.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
% *OUTPUT*
%           tanstruct_n2o_out: STRUCTURE - the same as the N2O input
%           structure, but with updated flags.

%% Define some things
n2o = tanstruct_n2o_in;
ch4 = tanstruct_ch4_in;
T = tanstruct_T_in;

if length(n2o.occultation) ~= length(ch4.occultation) || length(n2o.occultation) ~= length(T.occultation)
    error('the number of occultations in the input data structures down''t match. You should subset first')
end

%% find the bad data
vmr_ch4 = ch4.vmr; flag_ch4 = ch4.quality_flags;
vmr_ch4(flag_ch4>7) = nan; vmr_ch4(:, any(flag_ch4 >= 4 & flag_ch4 <= 6)) = nan;
vmr_ch4_20 = vmr_ch4(20,:);
vmr_T_20 = T.vmr(20,:);
vmr_n2o_20 = n2o.vmr(20,:);
lat = n2o.lat_tangent;

bad_occ = (vmr_ch4_20<10.6e-7 & vmr_n2o_20>1.6e-7) | (vmr_ch4_20<9.5e-7 & vmr_n2o_20>1.4e-7) & abs(lat) > 55 & vmr_T_20<196;

%% change the bad rows to have flag 7
tanstruct_n2o_out = n2o;
tanstruct_n2o_out.quality_flags(:, bad_occ) = 7;

%
end