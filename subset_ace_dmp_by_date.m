function [ dmpstruct_out ] = subset_ace_dmp_by_date( dmpstruct_in, start_date, end_date )
%A function to subset ace data corresponding to a particular year. Empty
%arrays are produced if there is no data for that year.

% *INPUT*
%           dmpstruct_in: STRUCTURE - containins the ACE dmp data. It is
%           usually created using 'read_ace_dmp' or 'read_ace_dmp_for_mat'.
%
%           start_date: VECTOR - The start date of the time frame for which
%           you want to extract ace data. The format is [year,month,date].
%
%           end_date: VECTOR - The end date of the time frame for which
%           you want to extract ace data. The format is [year,month,date].
%
% *OUTPUT*
%           dmpstruct_out: STRUCTURE - output has the same fields as the
%           input, but with only data that corresponds to the chosen date
%           range.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017

%% Define some things
dmp = dmpstruct_in;
mjdace = dmp.date_mjd;
mjdstart = date2mjd(start_date(1), start_date(2), start_date(3));
mjdend = date2mjd(end_date(1), end_date(2), end_date(3));

%% pick out the data that corresponds to the year
idate = find(mjdace >= mjdstart & mjdace <= mjdend); % get the indices of the dates which lie within the chosen range

%Subset the data
dmpout = reduce_dmpstruct_by_rowindex(dmp,idate);

dmpstruct_out = dmpout;

end

