function [ tanstruct_out ] = include_ace_scaled_apriori(tanstruct_in)
%A function to change the flags of scaled a priori data from -888 to 0. The
%reason for this is so that scaled a priori data will then be treated as
%ordinary data points.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
% *OUTPUT*
%           tanstruct_out_flagged: STRUCTURE - output has the same
%           fields as the input, but with the -888 flags changed to 0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017

%% Alter the scaled a priori flags of 8
dataout = tanstruct_in;
dataout.quality_flags(dataout.quality_flags == 8) = 0;

%%
tanstruct_out = dataout;

end

