function [ tanstruct_out ] = subset_ace_by_lst_tangent( tanstruct_in, lst_start, lst_end )
%A function to subset ace data to that which has a local solar time (LST)
%that lies within a chosen range. Empty arrays are produced if there
%is no data for an output. The longitude that has been added from the GLC
%files is used to calculate the LST.
 
% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           lst_start: FLOAT - the start of the local solar time
%           range, to which you want to subset the ACE data.
%
%           lst_end: FLOAT - the end of the local solar time
%           range, to which you want to subset the ACE data.           
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same
%           fields as the input, but with only data that corresponds to
%           local solar times in the input range
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 04/2018

%% Define some things
gas = tanstruct_in;
lststart = lst_start;
lstend = lst_end;

%% pick out the data that corresponds to before and after the chosen LST
lst_ace = get_ace_lst_tangent(gas); % get an array of the LSTs of the occultations
if lststart < lstend
    ilst_good = find( lst_ace >= lststart & lst_ace <= lstend); % get the indices of the data points that are within the input LST range
else % lststart > lstend
    ilst_good = find( lst_ace >= lststart | lst_ace <= lstend); % get the indices of the data points that are within the input LST range
end
%Subset the data
gas_out = reduce_tanstruct_by_rowindex(gas, ilst_good);

tanstruct_out = gas_out;
%
end

