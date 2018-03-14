function [ tanstruct_before, tanstruct_after ] = subset_ace_by_lst( tanstruct_in, lst_break )
%A function to subset ace data into two parts, corresponding to before and
%after a chosen local solar time (LST). Empty arrays are produced if there
%is no data for an output. The longitude that has been added from the GLC
%files is used to calculate the LST.
 
% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           lst_break: FLOAT - the local solar time about which you want to
%           split the ace data. This is an optional input. The default LST
%           is 12.
%
% *OUTPUT*
%           tanstruct_before: STRUCTURE - output has the same
%           fields as the input, but with only data that corresponds to
%           local solar times before lst_break
%
%           tanstruct_after: STRUCTURE - output has the same
%           fields as the input, but with only data that corresponds to
%           local solar times after lst_break.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017

%% Define some things
gas = tanstruct_in;
if nargin == 1
    lstbreak = 12;
else
    lstbreak = lst_break;
end

%% pick out the data that corresponds to before and after the chosen LST
lst_ace = get_ace_lst(gas); % get an array of the LSTs of the occultations
ilst_before = find( lst_ace <= lstbreak ); % get the indices of the data points that are before the input LST
ilst_after = find( lst_ace > lstbreak ); % get the indices of the data points that are before the input LST

%Subset the data
gas_before = reduce_tanstruct_data_by_index(gas, ilst_before);
gas_after = reduce_tanstruct_data_by_index(gas, ilst_after);

tanstruct_before = gas_before;
tanstruct_after = gas_after;
%
end

