function [ tanstruct_out ] = subset_ace_by_3month( tanstruct_in, months_in )
%A function to subset ace data corresponding to a particular 3-month period
%of the year. Empty arrays are produced if there is no data for that month.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           month_in: STRING - the possible inputs are 'djf', 'mam', 'jja',
%           and 'son'. They represent december-january-february, etc.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same fields as the
%           input, but with only data for the chosen month. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 01/2018

%% Define some things
gas = tanstruct_in;
switch months_in
    case 'DJF'
        mn = [12 1 2];
    case 'MAM'
        mn = 3:5;
    case 'JJA'
        mn = 6:8;
    case 'SON'
        mn = 9:11;
    otherwise
        error('invalid input for the 3-month period')
end

%% pick out the data that corresponds to the year
vdate = datevec(mjd2datenum(gas.date_mjd)); % get a vector of the dates of the occultations
imn = find( vdate(:,2) == mn(1) | vdate(:,2) == mn(2) | vdate(:,2) == mn(3)); % get the indices of the occulations from that year

%Subset the data
gasout = reduce_tanstruct_by_rowindex(gas,imn);

tanstruct_out = gasout;

end

