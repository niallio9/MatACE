function [ tanstruct_out ] = subset_ace_by_year( tanstruct_in, year_in )
%A function to subset ace data corresponding to a particular year. Empty
%arrays are produced if there are no data for that year.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           year_in: FLOAT - the year of the data that you want.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same fields as the
%           input, but with only data for the chosen year. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017

%% Define some things
gas = tanstruct_in;
yr = year_in;

%% pick out the data that corresponds to the year
vdate = datevec(mjd2datenum(gas.date_mjd)); % get a vector of the dates of the occultations
iyr = find( vdate(:,1) == yr ); % get the indices of the occulations from that year

%Subset the data
gasout = reduce_tanstruct_by_rowindex(gas,iyr);

tanstruct_out = gasout;

end

