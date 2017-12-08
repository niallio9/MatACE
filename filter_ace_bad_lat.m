function [ tanstruct_out ] = filter_ace_bad_lat( tanstruct_in, lat_limit )
%A function to find ACE measurements that have erroneous GLC latitudes, and
%replace a;; the latitudes for those measurements with the tangent
%latitude.
%
% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'. 
%
%           lat_dif: FLOAT - the limit of acceptable difference in the glc
%           latitude and the tangent_latitude.
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
latlim = lat_limit;
sizelat = size(gasout.lat(:,1)); % should be 150 levels for the ace measurement

%% find the paces where the limit is broken and replace the data
diflat = gasin.lat - gasin.lat_tangent; % find the difference betweent he glc and tangent latitudes
[~,jbad] = find(diflat >= latlim); % get the row indices of where the limit is broken
jbad = unique(jbad); % remove duplicates of the row indices
gasout.lat(:,jbad) = repmat(gasout.lat_tangent(jbad),sizelat); % replace the GLC latitudes with the tangent latitude at all altitude levels
%%
tanstruct_out = gasout;
end

