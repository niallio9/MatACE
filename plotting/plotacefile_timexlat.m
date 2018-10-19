function [ ] = plotacefile_timexlat( gasname_in )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

gasname = gasname_in;
filein_pre = 'ACE_v3p6_';
filein_post = '.mat';
filein = strcat(filein_pre,gasname,filein_post);

temp = load(filein);
gas = temp.tanstruct;
clear temp
disp('applying flags...')
gas = apply_ace_flags(gas);
disp('done')

fs = 16;
figure
set(gcf,'Position', [265   128   719   577])
plot(mjd2datenum(gas.date_mjd), gas.lat_tangent, 'o')
title(sprintf('ACE-FTS %s measurement locations (30km tangent alt.)',gas.gas))
xlabel('date')
ylabel('latitude [degrees]')
% datetick('x','mmm')
% datetick('x', 'mm/yy')
dynamicDateTicks([],'x','mm/yy')
set(gca,'FontSize',fs)
end