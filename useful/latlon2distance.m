function [ distance_out ] = latlon2distance( lat_1, lon_1, lat_2, lon_2, altitude_in )
%A function to calculate the distance between two points on Earth. The
%answer is in units of km.

% *INPUT*
%           lat_1: VECTOR - the latitude(s) of the first location. Units
%           are degrees.
%
%           lon_1: VECTOR - the longitude(s) of the first location. Units
%           are degrees.
%
%           lat_2: VECTOR - the latitude(s) of the second location. Units
%           are degrees.
%
%           lon_2: VECTOR - the longitude(s) of the second location. Units
%           are degrees.
%
%           altitude_in: VECTOR - the altitude() of the locations, in km.
%
% *OUTPUT*
%           distance_out: VECTOR - the distance, in km, between the two
%           locations.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 04/18

%% define some things
Re = 6371; % radius of the earth
% lat1 = lat_1;
% lon1 = lon_1;
% lat2 = lat_2;
% lon2 = lon_2;
alt_km = altitude_in;
% change from degrees to radians
lat1 = (lat_1/180)*pi;
lon1 = (lon_1/180)*pi;
lat2 = (lat_2/180)*pi;
lon2 = (lon_2/180)*pi;
dlat = lat2 - lat1;
dlon = lon2 - lon1;
a = ((sin(dlat/2)).^2) + (cos(lat1).*cos(lat2)).*(sin(dlon/2)).^2;
c = 2*atan(sqrt(a)./sqrt(1-a));
% whos
d = (Re+alt_km).*c;

distance_out = d;

end

