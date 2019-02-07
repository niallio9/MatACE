function [ ] = plotacefile_timexlat_by_year( gasname_in, year_in )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% define some things
home_linux = '/home/niall/Dropbox/climatology/'; %#ok<NASGU>
home_mac = '/Users/niall/Dropbox/climatology/'; %#ok<NASGU>
home_windows = 'C:\Users\ryann\'; %#ok<NASGU>

gasname = gasname_in;
filein_pre = 'ACE_v3p6_';
filein_post = '.mat';
filein = strcat(home_windows,'/','ACE','/','matdata','/',filein_pre,gasname,filein_post);

temp = load(filein);
gas = temp.tanstruct;
clear temp
gas = subset_ace_by_year(gas,year_in);
disp('applying flags...')
gas = apply_ace_flags(gas);
disp('done')

n = length(gas.occultation);
sdate = mjd2datenum(gas.date_mjd);

fs = 16;
figure
set(gcf,'Position', [265   128   719   577])
plot(sdate, gas.lat_tangent, 'o')
title(sprintf('ACE-FTS %s measurement locations (30km tangent alt.)',gas.gas))
xlabel('date')
ylabel('latitude [degrees]')
% datetick('x','mmm')
dynamicDateTicks([],'x','mm/yy')
set(gca,'FontSize',fs)
text(mean(sdate),max(gas.lat_tangent),sprintf('N = %d', n),'Fontsize', 16)

end