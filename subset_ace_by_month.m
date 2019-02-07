function [ tanstruct_out, chosen_indices ] = subset_ace_by_month( tanstruct_in, month_in )
%A function to subset ace data corresponding to a particular month of the
%year. Empty arrays are produced if there is no data for that month.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           month_in: FLOAT - the month of the data that you want. Enter 1
%           for January, 2 for February, etc..
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same fields as the
%           input, but with only data for the chosen month. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017

%% Define some things
gas = tanstruct_in;
mn = month_in;

%% pick out the data that corresponds to the year
vdate = datevec(mjd2datenum(gas.date_mjd)); % get a vector of the dates of the occultations
imn = find( vdate(:,2) == mn ); % get the indices of the occulations from that month

%Subset the data
gasout = reduce_tanstruct_by_rowindex(gas,imn);

tanstruct_out = gasout;
chosen_indices = imn;

end

