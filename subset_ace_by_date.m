function [ tanstruct_out ] = subset_ace_by_date( tanstruct_in, start_date, end_date )
%A function to subset ace data corresponding to a particular date ange.
%Empty arrays are produced if there are no data for the range. The start
%and end dates are included in the range.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           start_date: VECTOR - The start date of the time frame for which
%           you want to extract ace data. The format is [year,month,day].
%
%           end_date: VECTOR - The end date of the time frame for which
%           you want to extract ace data. The format is [year,month,day].
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same
%           fields as the input, but with only data in the chosen date
%           range. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017

%% Define some things
gas = tanstruct_in;
mjdace = gas.date_mjd;
mjdstart = date2mjd(start_date(1), start_date(2), start_date(3));
mjdend = date2mjd(end_date(1), end_date(2), end_date(3));

%% pick out the data that corresponds to the year
idate = find(mjdace >= mjdstart & mjdace <= mjdend); % get the indices of the dates which lie within the chosen range

%Subset the data
gasout = reduce_tanstruct_by_rowindex(gas,idate);

tanstruct_out = gasout;

end

