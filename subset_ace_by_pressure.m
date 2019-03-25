function [ tanstruct_out, pres_ace ] = subset_ace_by_pressure( tanstruct_in, start_pressure, end_pressure )
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
pres_ace = gas.pressure_hPa;
if length(pres_ace(1,:)) == 1
    pres_ace = repmat(pres_ace, 1, length(gas.occultation));
end
pres_start = start_pressure;
pres_end = end_pressure;
%We will work from low to high, so make sure altend is smaller than
%altstart. Swap them around if not.
if pres_end > pres_start
    altstart_old = pres_start;
    altend_old = pres_end;
    pres_start = altend_old;
    pres_end = altstart_old;
end

%% pick out the data that corresponds to the altitude bin
ialt = find(pres_ace <= pres_start & pres_ace > pres_end); % get the indices of the altitudes that lie within the chosen range

%Subset the data
gasout = reduce_tanstruct_data_by_index(gas,ialt);

tanstruct_out = gasout;
%
end

