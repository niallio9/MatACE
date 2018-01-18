function [ DOY ] = mjd2doy( MJD_in )
%%A short function to convert MJD to day of the year. It requires the
%%function 'datevec2doy.m'

DOY = datevec2doy(datevec(mjd2datenum(MJD_in)));

end

