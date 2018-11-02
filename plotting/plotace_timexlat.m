function [ ] = plotace_timexlat( tanstruct_in )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

gas = tanstruct_in;

fs = 16;
figure
set(gcf,'Position', [246    95   719   570])
plot(mjd2datenum(gas.date_mjd), gas.lat_tangent, 'o')
title(sprintf('ACE-FTS %s measurement locations (30km tangent alt.)',gas.gas))
xlabel('date')
ylabel('latitude [degrees]')
% datetick('x','mmm')
datetick('x', 'mm/yy')
set(gca,'FontSize',fs)
end