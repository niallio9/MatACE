function [ climstruct_out ] = reduce_climstruct_data_by_obs_nr( climstruct_in, obs_min )
%A function to reduce the ace data according to the provided observation
% number minimum. Data that do not correspond to the indices are changed to
% NaNs.

% *INPUT*
%           climstruct_in: STRUCTURE - the ace climatology data. The
%           structure is created by running one of the
%           'make_ace_climatology...' matlab functions. 
%
%           obs_min: INTEGER - the minimum number of observations needed to
%           include a data point in the climatology.
%
% *OUTPUT*
%           climstruct_out: STRUCTURE - with the same fields as the
%           input but with data points removed that were the mean of fewer
%           measurements than 'obs_min'.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 03/2018

if floor(obs_min) ~= obs_min
    error('please provide and integer for the observation point minimum')
end

clim = climstruct_in;
obs = clim.obs_count;
bad_I = find(obs <= obs_min);

clim_out = clim;
clim_out.vmr_zonal(bad_I) = nan;
clim_out.vmr_zonal_var(bad_I) = nan;
clim_out.vmr_zonal_error(bad_I) = nan;
clim_out.vmr_zonal_standard_error(bad_I) = nan;
clim_out.lon_mean(bad_I) = nan;
clim_out.obs_count(bad_I) = nan;
for i = 1:length(bad_I)
    clim_out.obs_location{bad_I(i)} = [];
end

climstruct_out = clim_out;
%
end

