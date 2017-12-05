function [ mjd ] = datevec2mjd( datevec )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
l = length(datevec(:,1));
mjd = zeros(l);
for i = 1:l
mjd(i) = date2mjd(datevec(1),datevec(2),datevec(3), datevec(4),datevec(5),datevec(6));
end

end

