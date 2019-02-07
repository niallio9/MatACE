function [tanstruct_inside, tanstruct_edge, tanstruct_outside] = ...
        subset_ace_by_vortex( tanstruct_in )
%A function to subset ace data corresponding to whether the location of a
%measurement point with respect to the polar vortex: inside, outside, or
%within the edge of the vortex

% *INPUT*
%           tanstruct_in: STRUCTURE - a .MAT structure containing ACE
%           data and metadata. It is usually created using
%           'read_ace_ncdata', 'read_ace_ncdata_for_mat'.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same
%           fields as the input, but with the data that coincides with the
%           dmp info.
%
%           dmpstruct_out: STRUCTURE - output has the same
%           fields as the input, but with the data that coincides with the
%           gas info.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 12/2018

%% Define some things
gas = tanstruct_in;
dmp_zlimit = 75; % the upper altitude limit of the dmp data

if ~isfield(gas,'spv')
    error('there is no sPV information int he structure\n you can use ''merge_ace_dmp.m''')
end
spv_in = 1.6e-4;
spv_out = 1.2e-4;


%% subset the ace data to match the altitude range of the dmp data
disp('subsetting the ACE data to the altitude range of the spv data (60km)...')
gas = subset_ace_by_alt(gas, gas.altitude_km(1), dmp_zlimit);
disp('done')

%% pick out the data that corresponds to the vortex relative location
ispv_in = find(abs(gas.spv) > spv_in); % get the indices of the data that lie in the vortex
ispv_out = find(abs(gas.spv) < spv_out); % get the indices of the latitudes that lie within the chosen range
ispv_edge = find(abs(gas.spv) >= spv_out & abs(gas.spv) <= spv_in); % get the indices of the latitudes that lie within the chosen range

tanstruct_inside = reduce_tanstruct_data_by_index(gas, ispv_in);
tanstruct_outside = reduce_tanstruct_data_by_index(gas, ispv_out);
tanstruct_edge = reduce_tanstruct_data_by_index(gas, ispv_edge);
%
end

