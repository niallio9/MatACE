function [ tanstruct_out ] = subset_ace_by_lat( tanstruct_in, start_lat, end_lat )
%A function to subset ace data corresponding to a particular latitude
%range. Empty arrays are produced if there are no data for that range. The
%ace GLC data should be included in the input tanstruct. You can do this
%with 'merge_ace_glc.m'.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           start_lat: FLOAT - The start point of the latitude range for which
%           you want to extract ace data.
%
%           end_lat: FLOAT - The end point of the latitude range for which
%           you want to extract ace data.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same
%           fields as the input, but with only data in the chosen latitude
%           range.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017

%% Check if the structure has latitude and longitude information
if ~isfield(tanstruct_in,'lat')
    error('there is no ''lat'' field in the input. You can add the lat/lon information with ''add_ace_glc.m'' ')
end

%% Define some things
gas = tanstruct_in;
latace = gas.lat;
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

%% pick out the data that corresponds to the latitude bin
ilat = find(latace >= latstart & latace < latend); % get the indices of the latitudes that lie within the chosen range

%Subset the data
gasout = reduce_tanstruct_data_by_index(gas,ilat);

tanstruct_out = gasout;
%
end

