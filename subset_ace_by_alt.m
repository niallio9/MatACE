function [ tanstruct_out, altace ] = subset_ace_by_alt( tanstruct_in, start_alt, end_alt )
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
% NJR - 09/2018

%% Check if the structure has latitude and longitude information
% if ~isfield(tanstruct_in,'lat')
%     error('there is no ''lat'' field in the input. You can add the lat/lon information with ''merge_ace_glc.m'' ')
% end

%% Define some things
gas = tanstruct_in;
altace = gas.altitude_km;
if length(altace(1,:)) == 1
    altace = repmat(altace, 1, length(gas.occultation));
end
altstart = start_alt;
altend = end_alt;
%We will work from low to high, so make sure altend is larger than
%altstart. Swap them around if not. Southern lats are negative.
if altend < altstart
    altstart_old = altstart;
    altend_old = altend;
    altstart = altend_old;
    altend = altstart_old;
end

%% pick out the data that corresponds to the altitude bin
ialt = find(altace >= altstart & altace < altend); % get the indices of the latitudes that lie within the chosen range

%Subset the data
gasout = reduce_tanstruct_data_by_index(gas,ialt);

tanstruct_out = gasout;
%
end

