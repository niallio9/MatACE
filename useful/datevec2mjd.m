function [ mjd ] = datevec2mjd( datevec )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
lvec = length(datevec(:,1));
mjd = zeros(lvec,1);
for i = 1:lvec
mjd(i) = date2mjd(datevec(i,1),datevec(i,2),datevec(i,3), datevec(i,4),datevec(i,5),datevec(i,6));
end

end

