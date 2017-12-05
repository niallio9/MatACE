function [lst]=mjd2lst(mjd,longitude)
% slightly modified from Jaho's utc2lst
% NJR - 11/17 
[~,~,~,hour,minute,secs,ticks] = mjd2utc(mjd);

lst = mod(hour+minute/60+(secs+ticks)/60/60+longitude/15.0,24.0);
 


