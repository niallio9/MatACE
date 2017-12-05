function [ sdates_out ] = mjd2datenum( mjd_data )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

[dates(:,1), dates(:,2), dates(:,3), dates(:,4), dates(:,5), dates(:,6)] = mjd2date(mjd_data);
sdates_out = datenum(dates);

end

