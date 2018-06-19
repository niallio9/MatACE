function [lst]=mjd2lst(mjd,longitude)
% slightly modified from Jaho's utc2lst
% NJR - 11/17 

if isvector(mjd) && ismatrix(longitude) % for the code to run on stupid Deluge
    mjd = repmat(mjd, length(longitude(:,1)), 1);
elseif isvector(longitude) && ismatrix(mjd)
    longitude = repmat(longitude, 1, length(mjd(1,:)));
end

[~,~,~,hour,minute,secs,ticks] = mjd2utc(mjd);
% whos
lst = mod(hour+minute/60+(secs+ticks)/60/60+longitude/15.0,24.0);
 


