function [ vdates_out ] = mjd2datevec( mjd_data )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

vdates_out = datevec(mjd2datenum(mjd_data));

end

