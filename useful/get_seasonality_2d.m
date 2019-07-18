function [ seasonal_average, seasonal_timeseries, seasonal_timeseries_detrended, mid_dates_mjd  ] = get_seasonality_2d( data1, date_mjd )
%A funcion to compare a number of data sets with one other.

% *INPUT*
%           data1: 2D ARRAY - the data with shich the other data sets will
%           be compared
%
%           date_mjd: 1D ARRAY - the dates of the measurements vectors
%           (columns of data1) in Modified Julian Date
%
% *OUTPUT*
%           seasonal_average: 2D ARRAY - 4-column array with data averaged
%           by season. Columns are in order of [MAM, JJA, SON, DJF].
%
%           seasonal_timeseries: 2D ARRAY - time setries of seasonally
%           averaged data, beginning with MAM.
%
%           seasonal_timeseries_detrended: 2D ARRAY - time setries of
%           seasonally averaged data, beginning with MAM, with the total
%           average for each corresponding season subtracted. This is to
%           remove the seasonal cycle from the data.
%
%           start_dates_mjd: 1D ROW VECTOR - the start dates of the data
%           used in the averages for the seasonal time series.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% average by season to get seasonal trend
vdates = mjd2datevec(date_mjd);
ldata = length(data1(:, 1));
ndata = length(data1(1, :));
if ~isequal(ndata, length(date_mjd))
    error('columns of data and time arrays do no match in size')
end

seasonal_average = nan(ldata, 4); % empty array for data
i = find(vdates(:,2) == 3 | vdates(:,2) == 4 | vdates(:,2) == 5); % find data in MAM
seasonal_average(:, 1) = nanmean(data1(:, i), 2);
i = find(vdates(:,2) == 6 | vdates(:,2) == 7 | vdates(:,2) == 8); % JJA
seasonal_average(:, 2) = nanmean(data1(:, i), 2);
i = find(vdates(:,2) == 9 | vdates(:,2) == 10 | vdates(:,2) == 11); % SON
seasonal_average(:, 3) = nanmean(data1(:, i), 2);
i = find(vdates(:,2) == 12 | vdates(:,2) == 1 | vdates(:,2) == 2); % DJF
seasonal_average(:, 4) = nanmean(data1(:, i), 2);

%% make seasonal timeseries of Cly
years_unique = unique(vdates(:,1)); % the unique years in the time data
lyears_unique = length(years_unique);
lseasons = lyears_unique * 4 - 1;  % no djf in 1st year and last year
seasonal_timeseries = nan(ldata, lseasons);
seasonal_timeseries_detrended = nan(ldata, lseasons);
mid_dates_mjd = nan(1, lseasons);
season_mean_index = repmat(1:4, [1, lseasons]); % index representing the total seasonal average to be subtraced from each seasonal timeseries. starts with MAM
i = 0;
for j = 3: 3: lseasons * 3 % 3 months in every season. start in March
    i = i + 1;
%     i
%     j
    seasonal_timeseries(:,i) = nanmean(data1(:, j: j + 2), 2);
    
    seasonal_timeseries_detrended(:,i) = nanmean(data1(:, j: j + 2), 2) - seasonal_average(:, season_mean_index(i)); % remove seasonal average corresponding to the season
%     whos
    vdate_mid = vdates(j,:) + [0 1 0 0 0 0]; % add one month to the date to bring it to the mid-point of the season
    mid_dates_mjd(i) = datevec2mjd(vdate_mid); % the start times of the seasons
end
%
end

