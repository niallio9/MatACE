function [ mjd_out ] = datenum2mjd( datenum_data )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

dates = datevec(datenum_data);
mjd_out = date2mjd(dates(:,1),dates(:,2),dates(:,3),dates(:,4), dates(:,5), dates(:,6));

end

