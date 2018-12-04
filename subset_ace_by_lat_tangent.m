function [ tanstruct_out ] = subset_ace_by_date( tanstruct_in, start_lat, end_lat )
%A function to subset ace data corresponding to a particular tangent
%latitude range. Empty arrays are produced if there are no data for that
%range.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           start_lat: FLOAT - The start point of the tangent latitude
%           range for which you want to extract ace data.
%
%           end_lat: FLOAT - The end point of the tnagent latitude range
%           for which you want to extract ace data.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same
%           fields as the input, but with only data in the chosen tangent
%           latitude range. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2018

%% Define some things
gas = tanstruct_in;
latace = gas.lat_tangent;
latstart = start_lat;
latend = end_lat;
%We will work from south to north, so make sure latend is larger than
%latstart. Swap them around if not. Southern lats are negative.
if latend < latstart
    latstart_old = latstart;
    latend_old = latend;
    latstart = latend_old;
    latend = latstart_old;
end

%% pick out the data that corresponds to the year
ilat_good = find(latace >= latstart & latace < latend); % get the indices of the dates which lie within the chosen range

%Subset the data
gasout = reduce_tanstruct_by_rowindex(gas,ilat_good);

tanstruct_out = gasout;

end

