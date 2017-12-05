function [ DOY ] = dayofyear( A )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

A(:,1) = m;
A(:,2) = d;
l = length(A)

for i = 1:l;
    if m(i) == 1;
        doy(i) = d(i);
    elseif m(i) == 2;
        doy(i) = d(i) + 31;
    elseif m(i) == 3;
        doy(i) = d(i) +59;
    elseif m(i) == 4;
        doy(i) = d(i) + 90;
    elseif m(i) == 5;
        doy(i) = d(i) + 120;
    elseif m(i) == 6;
        doy(i) = d(i) + 151;
    elseif m(i) == 7;
        doy(i) = d(i) + 181;
    elseif m(i) == 8;
        doy(i) = d(i) + 212;
    elseif m(i) == 9;
        doy(i) = d(i) + 243;
    elseif m(i) == 10;
        doy(i) = d(i) + 273;
    elseif m(i) == 11;
        doy(i) = d(i) + 304;
    elseif m(i) == 11;
        doy(i) = d(i) + 334;  
    end
end

DOY = doy;
end
    
