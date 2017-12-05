function [ dmpstruct_out ] = subset_ace_dmp_by_month( dmpstruct_in, month_in )
%A function to subset ACE DMP data corresponding to a particular month.
%Empty arrays are produced if there is no data for that month.

% *INPUT*
%           dmpstruct_in: STRUCTURE - containins the ACE dmp data. It is
%           usually created using 'read_ace_dmp' or 'read_ace_dmp_for_mat'.
%
%           month_in: FLOAT - the month of the data that you want. Enter 1
%           for January, 2 for February, etc..
%
% *OUTPUT*
%           dmpstruct_out: STRUCTURE - output has the same fields as the
%           input, but with only data that corresponds to the chosen
%           month.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017

%% Define some things
dmp = dmpstruct_in;
mn = month_in;

%% pick out the data that corresponds to the year
vdate = datevec(mjd2datenum(dmp.date_mjd)); % get a vector of the dates of the occultations
imn = find( vdate(:,2) == mn ); % get the indices of the occulations from that year

%Subset the data
dmpout = reduce_dmpstruct_by_rowindex(dmp,imn);

dmpstruct_out = dmpout;

end

