function [ vortex_pos ] = get_ace_vortex_position( tanstruct_in )
%A function to calculate the positions of ace measurement points with
%respect to a polar vortex.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
% *OUTPUT*
%           vortex_pos: ARRAY - an array of the same size as the
%           tanstruct_in.spv variable, containing numbers that indicate the
%           position of the measurement point with respect to a polar
%           vortex: 0 = inside, 1 = in the edge, 2 = outside.

%% Define some things
gas = tanstruct_in;
spv_in = 1.6e-4;
spv_out = 1.2e-4;

if ~isfield(gas,'spv')
    error('there is no sPV information in the structure. You can use merge_ace_dmp.m')
end

%% create the lst vector
vortex_pos = nan(size(gas.spv));
vortex_pos(abs(gas.spv) > spv_in) = 0;
vortex_pos(abs(gas.spv) < spv_out) = 2;
vortex_pos(abs(gas.spv) <= spv_in & abs(gas.spv) >= spv_out) = 1;

%
end