function [ dmpstruct_out ] = subset_ace_dmp_by_year( dmpstruct_in, year_in )
%A function to subset ace DMP data corresponding to a particular year.
%Empty arrays are produced if there are no data for that year.

% *INPUT*
%           dmpstruct_in: STRUCTURE - containins the ACE dmp data. It is
%           usually created using 'read_ace_dmp' or 'read_ace_dmp_for_mat'.
%
%           year_in: FLOAT - the year of the data that you want.
%
% *OUTPUT*
%           dmpstruct_out: STRUCTURE - output has the same fields as the
%           input, but with only data for the chosen year. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017

%% Define some things
dmp = dmpstruct_in;
yr = year_in;

%% pick out the data that corresponds to the year
vdate = datevec(mjd2datenum(dmp.date_mjd)); % get a vector of the dates of the occultations
iyr = find( vdate(:,1) == yr ); % get the indices of the occulations from that year

%Subset the data
dmpout = reduce_dmpstruct_by_rowindex(dmp,iyr);

dmpstruct_out = dmpout;

end

