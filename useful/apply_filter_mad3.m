function [ array_out ] = apply_filter_mad3(array_in)
%A function to apply the a filtering whereby only data that lies within 3
%median absolute deviations of the median is kept. Each row of the data is
%filtered seperately.

% *INPUT*
%           array_in: ARRAY - array of data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
% *OUTPUT*
%           array_out: ARRAY - array of data with the filtered points
%           replaced by nans.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2018

%% Remove flagged values
datain = array_in;
% zlev = length(datain.altitude_km(:,1));
zlev = length(datain(:,1));
out = array_in; %this structure will be edited later on

medianvmr = nanmedian(datain,2); % get the median of the data at each level
madvmr = mad(datain,1,2); % get the median absolute deviation at each altitude
lim_upper = medianvmr + 3.*madvmr;
lim_lower = medianvmr - 3.*madvmr;

%% go through each altitude and remove values outside 3 mads
for n = 1:zlev
    ibad = find(datain(n,:) > lim_upper(n) | datain(n,:) < lim_lower(n));
    out(n,ibad) = nan; 
end

%% remove any columns that are only nans now


%%
array_out = out;

end

