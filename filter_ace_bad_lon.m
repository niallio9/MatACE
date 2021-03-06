function [ tanstruct_out ] = filter_ace_bad_lon( tanstruct_in, lon_limit )
%A function to find ACE measurements that have erroneous GLC longitudes,
%and replace the latitudes for those measurements with the tangent
%longitude. 
%
% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'. 
%
%           lat_dif: FLOAT - the limit of acceptable difference in the glc
%           longitude and the tangent_longitude.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same
%           fields as the input, but with data filtered out according to
%           the limit above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017

%% Define some things
gasin = tanstruct_in;
gasout = gasin;
if nargin < 2
    lonlim = 40; % ad hoc value - NJR 10/18
else
    lonlim = lon_limit;
end
sizelon = size(gasout.lon(:,1)); % should be 150 levels for the ace measurement
lontangent = repmat(gasin.lon_tangent, [sizelon,1] ); % need this because matlab on deluge cant subtract a vector from a matrix

%% find the paces where the limit is broken and replace the data
diflon = gasin.lon - lontangent; % find the difference betweent he glc and tangent latitudes
[~,jbad] = find(abs(diflon) >= lonlim); % get the row indices of where the limit is broken
jbad = unique(jbad); % remove duplicates of the row indices
gasout.lon(:,jbad) = repmat(gasout.lon_tangent(jbad),sizelon); % replace the GLC latitudes with the tangent latitude at all altitude levels
%%
tanstruct_out = gasout;
end

