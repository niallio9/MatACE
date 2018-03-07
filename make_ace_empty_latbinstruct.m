function [ latbinstruct ] = make_ace_empty_latbinstruct( tanstruct_in, lat_bounds )
% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           lat_bounds: VECTOR - a vector containing the start and end
%           boundaries of the latitude bins. It is assumed that the end of
%           one bin corresponds to the start of another, meaning that the
%           binning is performed for all data within the first and last
%           entry of the vector. For the standard ACE climatology,
%           lat_bounds will be -90:5:90.
%
% *OUTPUT*
%           latbinstruct: STRUCTURE - output contains different
%           fields than the input. The fields relate to zonal mean values
%           of the input gas data. Differing output fields are
%           lat: the mid-point of the latitude bins
%           lat_bnds: the boundaries of the latitude bins
%           vmr_zonal: zonal means of the gas vmrs
%           vmr_zonal_var: the zonal standard deviation gas vmrs
%           vmr_zonal_error: the standard error on the zonal mean
%           obs_count: the number of observations per grid point
%           obs_location: the lon, lat, and alt of each data point in
%           each bin.
%           lst_median/mean/max/min/var: the median, mean, maximum,
%           minimum, and variance of the local solar times of the data
%           points in each bin.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017

%% Define some things
gas = tanstruct_in;

%% Check if the ace data structure is on a single altitude or pressure grid
if ~isfield(gas,'lat')
    error('there is no ''lat'' field in the input. You can add the lat/lon information with ''merge_ace_glc.m'' ')
end
if isfield(gas,'pressure_hPa') % exists if the ace data has been interpolated to a pressure grid, or is the new ace version that includes pressure
    if length(gas.pressure_hPa(1,:)) > 1 && length(gas.altitude_km(1,:)) > 1 % there should be either one altitude or one pressure profile for interpolated data
        error('There are multiple pressure and altitude grids. Binning requires that one of these be a common grid')
    end
else
    if length(gas.altitude_km(1,:)) > 1 % If pressure_hPa isn't part of the data structure, then there must be a common altitude grid
        error('There are multiple altitude grids in the ACE data. Binning requires at least one common grid (pressure or altitude)')
    end
end

%% Define some things
latbnds = lat_bounds;
llat = length(latbnds) - 1;
latout = nan(1,llat);
for i = 1:llat
    latout(i) = mean([latbnds(i),latbnds(i+1)]); %get the midpoints of the latitude bins
end
% lorbit = length(gas.occultation);
% lgrid = length(gas.altitude_km(:,1));
%maxl = lorbit*lgrid; % this is the maximum number of points that can go into a bin

%% Set up the output structure
gaslatbin.source_file = gas.source_file; % the same original source of the data
gaslatbin.beta_angle_mean = nan(1,llat); % one average per month
gaslatbin.beta_angle_var = nan(1,llat); % one average per month
% gaslatbin.date_mjd_mean = nan(1,llat); % make it a mean mjd **** this has been removed cos it's stupid
gaslatbin.start_date = []; % for the date of the earliest measurement
gaslatbin.end_date = []; % for the date of the latest measurement
gaslatbin.gas = gas.gas; % the same as the original file
gaslatbin.lon_tangent_mean = nan(1,llat);
gaslatbin.lat_tangent_mean = nan(1,llat);
%the data is either on a common altitude grid or pressure grid. It was
%checked above to make sure one of these is common to all data, i.e., one
%of either altitude or pressure is a vector and not an array
if length(gas.altitude_km(1,:)) == 1
    gaslatbin.altitude_km = gas.altitude_km;
    lalt = length(gas.altitude_km); % should be the same for altitude and pressure anyway
else
    gaslatbin.altitude_km_mean = nanmean(gas.altitude_km,2); %there should always be an altitude field in the ace data, whether interpolated or not
end
if length(gas.pressure_hPa(1,:)) == 1
    gaslatbin.pressure_hPa = gas.pressure_hPa; % there will only be a pressure field for data interpolated onto a pressure grid
    lalt = length(gas.pressure_hPa);% should be the same for altitude and pressure anyway
else
    gaslatbin.pressure_hPa_mean = nanmean(gas.pressure_hPa,2); %there should always be an altitude field in the ace data, whether interpolated or not   
end
gaslatbin.vmr_zonal = nan(lalt,llat);
gaslatbin.vmr_zonal_var = nan(lalt,llat);
gaslatbin.vmr_zonal_error = nan(lalt,llat);
gaslatbin.vmr_zonal_standard_error = nan(lalt,llat);
gaslatbin.lat = latout;
gaslatbin.lat_bounds = latbnds;
gaslatbin.lon_mean = nan(lalt,llat);
gaslatbin.obs_count = nan(lalt,llat);
gaslatbin.obs_location = cell(lalt,llat);% going to fill these cells with the time and location of each datapoint that goes in
% add in some local solar time stuff
gaslatbin.lst_median = nan(1,llat);
gaslatbin.lst_mean = nan(1,llat);
gaslatbin.lst_max = nan(1,llat);
gaslatbin.lst_min = nan(1,llat);
gaslatbin.lst_var = nan(1,llat);

%%
latbinstruct = gaslatbin;
end

