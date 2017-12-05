function [lst]=utc2lst(mjd,longitude)

[year,month,day,hour,minute,secs,ticks] = mjd2utc(mjd);

lst = mod(hour+minute/60+(secs+ticks)/60/60+longitude/15.0,24.0);
 


