function [ tanstruct_out ] = filter_ace_bad_lat( tanstruct_in, lat_limit )
%A function to find ACE measurements that have erroneous GLC latitudes, and
%replace the latitudes for those measurements with the tangent
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

%% Define some things
gasin = tanstruct_in;
gasout = gasin;
if nargin < 2
    latlim = 10; % ad hoc value - NJR 10/18
else
    latlim = lat_limit;
end
sizelat = size(gasout.lat(:,1)); % should be 150 levels for the ace measurement
lattangent = repmat(gasin.lat_tangent, [sizelat,1] ); % need this because matlab on deluge cant subtract a vector from a matrix

%% find the places where the limit is broken and replace the data
diflat = gasin.lat - lattangent; % find the difference between the glc and tangent latitudes
[~,jbad] = find(abs(diflat) >= latlim); % get the row indices of where the limit is broken
jbad = unique(jbad); % remove duplicates of the row indices
gasout.lat(:,jbad) = repmat(gasout.lat_tangent(jbad),sizelat); % replace the GLC latitudes with the tangent latitude at all altitude levels
%%
tanstruct_out = gasout;
end

