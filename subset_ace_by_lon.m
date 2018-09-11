function [ tanstruct_out ] = subset_ace_by_lon( tanstruct_in, start_lon, end_lon )
%A function to subset ace data corresponding to a particular latitude
%range. Empty arrays are produced if there are no data for that range. The
%ace GLC data should be included in the input tanstruct. You can do this
%with 'merge_ace_glc.m'.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           start_lon: FLOAT - The start point of the longitude range for which
%           you want to extract ace data.
%
%           end_lon: FLOAT - The end point of the longitude range for which
%           you want to extract ace data.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same
%           fields as the input, but with only data in the chosen longitude
%           range.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 07/2018

%% Check if the structure has latitude and longitude information
if ~isfield(tanstruct_in,'lon')
    error('there is no ''lon'' field in the input. You can add the lat/lon information with ''merge_ace_glc.m'' ')
end

%% Define some things
gas = tanstruct_in;
lonace = gas.lon;
lonstart = start_lon;
lonend = end_lon;
%We will work from west to east (-180 to 180), so make sure lonend is larger than
%lonstart. Swap them around if not.
if lonend < lonstart
    latstart_old = lonstart;
    latend_old = lonend;
    lonstart = latend_old;
    lonend = latstart_old;
end

%% pick out the data that corresponds to the latitude bin
ilon = find(lonace >= lonstart & lonace < lonend); % get the indices of the latitudes that lie within the chosen range

%Subset the data
gasout = reduce_tanstruct_data_by_index(gas,ilon);

tanstruct_out = gasout;
%
end

