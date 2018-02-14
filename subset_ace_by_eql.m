function [ tanstruct_out ] = subset_ace_by_eql( tanstruct_in, start_eql, end_eql )
%A function to subset ace data corresponding to a particular equivalent
%latitude range. Empty arrays are produced if there are no data for that
%range. The ace GLC data should be included in the input tanstruct. You can
%do this with 'merge_ace_glc.m'.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           start_eql: FLOAT - The start point of the equivalent latitude
%           range for which you want to extract ace data.
%
%           end_eql: FLOAT - The end point of the equivalent latitude range
%           for which you want to extract ace data.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same
%           fields as the input, but with only data in the chosen
%           equivalent latitude range.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017

%% Check if the structure has latitude and longitude information
if ~isfield(tanstruct_in,'eql')
    error('there is no ''lat'' field in the input. You can add the eql information with ''merge_ace_dmp.m'' ')
end

%% Define some things
gas = tanstruct_in;
eqlace = gas.eql;
eqlstart = start_eql;
eqlend = end_eql;
%We will work from south to north, so make sure latend is larger than
%latstart. Swap them around if not. Southern lats are negative.
if eqlend < eqlstart
    eqlstart_old = eqlstart;
    eqlend_old = eqlend;
    eqlstart = eqlend_old;
    eqlend = eqlstart_old;
end

%% pick out the data that corresponds to the latitude bin
ieql = find(eqlace >= eqlstart & eqlace < eqlend); % get the indices of the latitudes that lie within the chosen range

%Subset the data
gasout = reduce_tanstruct_data_by_index(gas,ieql);

tanstruct_out = gasout;
%
end

