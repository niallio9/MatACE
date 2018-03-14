function [ eqlbinstruct ] = make_ace_empty_eqlbinstruct( tanstruct_in, eql_bounds )
% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           eql_bounds: VECTOR - a vector containing the start and end
%           boundaries of the equivalent latitude bins. It is assumed that
%           the end of one bin corresponds to the start of another, meaning
%           that the binning is performed for all data within the first and
%           last entry of the vector. For the standard ACE climatology,
%           eql_bounds will be -90:5:90.
%
% *OUTPUT*
%           eqlbinstruct: STRUCTURE - output contains different
%           fields than the input. The fields relate to zonal mean values
%           of the input gas data. Differing output fields are
%           eql: the mid-point of the equivalent latitude bins
%           eql_bnds: the boundaries of the latitude bins
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
% NJR - 02/2018

%% Define some things
gas = tanstruct_in;

%% Check if the ace data structure is on a single altitude or pressure grid
if ~isfield(gas,'eql')
    error('there is no ''eql'' field in the input. You can add the eql information with ''merge_ace_dmp.m'' ')
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
eqlbnds = eql_bounds;
leql = length(eqlbnds) - 1;
eqlout = nan(1,leql);
for i = 1:leql
    eqlout(i) = mean([eqlbnds(i),eqlbnds(i+1)]); %get the midpoints of the equivalent latitude bins
end
% lorbit = length(gas.occultation);
% lgrid = length(gas.altitude_km(:,1));
%maxl = lorbit*lgrid; % this is the maximum number of points that can go into a bin
%% Set up the output structure
gaseqlbin.source_file = gas.source_file; % the same original source of the data
gaseqlbin.beta_angle_mean = nan(1,leql); % one average per month  
gaseqlbin.beta_angle_var = nan(1,leql); % one average per month
gaseqlbin.doy_mean = nan(1,leql); % this has been added to get the mean day of the year. can be used to make mean day of month later
gaseqlbin.start_date = []; % for the date of the earliest measurement
gaseqlbin.end_date = []; % for the date of the latest measurement
gaseqlbin.gas = gas.gas; % the same as the original file
gaseqlbin.lon_tangent_mean = nan(1,leql);
gaseqlbin.lat_tangent_mean = nan(1,leql);
%the data is either on a common altitude grid or pressure grid. It was
%checked above to make sure one of these is common to all data, i.e., one
%of either altitude or pressure is a vector and not an array
if length(gas.altitude_km(1,:)) == 1
    gaseqlbin.altitude_km_mean = gas.altitude_km;
    lalt = length(gas.altitude_km); % should be the same for altitude and pressure anyway
else
    gaseqlbin.altitude_km_mean = nanmean(gas.altitude_km,2); % leaving like this for now because there are sometimes only one measurement in a given climatology
end
if length(gas.pressure_hPa(1,:)) == 1
    gaseqlbin.pressure_hPa = gas.pressure_hPa; % should be one pressure profile for pressure interpolated data
    lalt = length(gas.pressure_hPa);% should be the same for altitude and pressure anyway
else
    gaseqlbin.pressure_hPa_mean = nanmean(gas.pressure_hPa,2); %this is kind of unnecessary for pressure interpolateed profiles. will cause an error later if working with altitude interpolated profiles
end
gaseqlbin.vmr_zonal = nan(lalt,leql);
gaseqlbin.vmr_zonal_var = nan(lalt,leql);
gaseqlbin.vmr_zonal_error = nan(lalt,leql);
gaseqlbin.vmr_zonal_standard_error = nan(lalt,leql);
gaseqlbin.eql = eqlout;
gaseqlbin.eql_bounds = eqlbnds;
gaseqlbin.lon_mean = nan(lalt,leql);
gaseqlbin.obs_count = nan(lalt,leql);
gaseqlbin.obs_location = cell(lalt,leql);% going to fill these cells with the time and location of each datapoint that goes in
% add in some local solar time stuff
gaseqlbin.lst_median = nan(1,leql);
gaseqlbin.lst_mean = nan(1,leql);
gaseqlbin.lst_max = nan(1,leql);
gaseqlbin.lst_min = nan(1,leql);
gaseqlbin.lst_var = nan(1,leql);

%%
eqlbinstruct = gaseqlbin;
end

