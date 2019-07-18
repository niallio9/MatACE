function [ stat_out  ] = compare_data_sets_2d( data1, error1, varargin )
%A funcion to compare a number of data sets with one other.

% *INPUT*
%           data1: 2D ARRAY - the data with shich the other data sets will
%           be compared
%
%           error1: 2D ARRAY - the corresponding errors of the first data
%           array
%
%           varargin: STRING - The other inputs are the data and
%           error ARRAYs with which to compare data1. e.g., 'data2, error2,
%           ..., dataN, errorN'.
% *OUTPUT*
%           stat_out: STRUCTURE - the statistics of the data comparison.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isequal(size(data1), size(error1))
    error('data and error arrays are not the same size for the main dataset')
end

data_and_error = varargin;
% data_and_error = varargin(1: end - 1);
% do_plot = varargin{end};
% if ~isscalar(do_plot)
%     error('the last element of varargin must be a scalar value of 0 or 1')
% elseif do_plot ~= 0 && do_plot ~= 1
%     error('the last element of varargin must be a scalar value of 0 or 1')
% end
ndata = length(data_and_error) / 2; % to account for data and error entries

stat_out.data1_mean = nanmean(data1, 2);
stat_out.data1_std = nanstd(data1, [], 2);
stat_out.data1_mean_error = nanmean(error1, 2);
stat_out.data1_mean_error_percent = 100 * (stat_out.data1_mean_error ./ stat_out.data1_mean);
% lcol = length(data1(1, :));
lrow = length(data1(:, 1));
out_size = nan(lrow, ndata);
stat_out.dataN_mean = out_size;
stat_out.dataN_mean_error = out_size;
stat_out.dataN_mean_error_percent = out_size;
stat_out.dataN_std = out_size;
stat_out.difference_mean = out_size;
stat_out.difference_std = out_size;
stat_out.pdifference_mean = [];
stat_out.pdifference_std = [];
stat_out.correlation = out_size;
stat_out.regression = out_size;
stat_out.regression_error = out_size;
nfill = 0;

for n = 1 : 2 : ndata*2
    nfill = nfill + 1;
    data2 = data_and_error{n}; % a 2d array. data number n/2
    error2 = data_and_error{n + 1}; % a 2d array. error number n/2
    if ~isequal(size(data1), size(data2))
        error('data array %i is not the same size as the main dataset', nfill)
    end
    if ~isequal(size(error1), size(error2))
        error('error array %i is not the same size as the main dataset', nfill)
    end
    
    stat_out.dataN_mean(:, nfill) = nanmean(data2, 2);
    stat_out.dataN_std(:, nfill) = nanstd(data2, [], 2);
    stat_out.dataN_mean_error(:, nfill) = nanmean(error2, 2);
    stat_out.dataN_mean_error_percent(:, nfill) = 100 * (stat_out.dataN_mean_error(:, nfill) ./ stat_out.dataN_mean(:, nfill));
    datadif = data2 - data1;
    stat_out.difference_mean(:, nfill) = nanmean(datadif, 2);
    stat_out.difference_std(:, nfill) = nanstd(datadif, [], 2);
    
    for i = 1 : lrow
        Jnonan = find(~isnan(data1(i, :)) & ~isnan(data2(i, :))); % where there are not nans in both datasets
        r = corrcoef(data1(i, Jnonan), data2(i, Jnonan));
%         size(r)
        if isscalar(r)
            stat_out.correlation(i, nfill) = r(1, 1);
        else
            stat_out.correlation(i, nfill) = r(1, 2);
        end
        
        [~, stat_out.regression(i, nfill), ~, stat_out.regression_error(i, nfill)] = york_fit(data2(i, Jnonan), data1(i, Jnonan), error2(i, Jnonan), error1(i, Jnonan));
    end
end
stat_out.pdifference_mean = (100 * stat_out.difference_mean) ./ ((stat_out.data1_mean + stat_out.dataN_mean) ./ 2);
stat_out.pdifference_std = (100 * stat_out.difference_std) ./ ((stat_out.data1_mean + stat_out.dataN_mean) ./ 2);
%
end

