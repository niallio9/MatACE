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
if isvector(lat1) && ismatrix(lat2) % for the code to run on stupid Deluge
    lat1 = repmat(lat1, 1, length(lat2(1,:)));
    lon1 = repmat(lon1, 1, length(lon2(1,:)));
elseif isvector(lat2) && ismatrix(lat1)
    lat2 = repmat(lat2, 1, length(lat1(1,:)));
    lon2 = repmat(lon2, 1, length(lon1(1,:)));
end
dlat = lat2 - lat1;
dlon = lon2 - lon1;
a = ((sin(dlat/2)).^2) + (cos(lat1).*cos(lat2)).*(sin(dlon/2)).^2;
c = 2*atan(sqrt(a)./sqrt(1-a));
if isvector(alt_km) && ismatrix(c) % for the code to run on stupid Deluge
    alt_km = repmat(alt_km, 1, length(c(1,:)));
end
% whos
d = (Re+alt_km).*c;

distance_out = d;
%
end

