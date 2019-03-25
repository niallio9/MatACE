function [ C, h ] = plot_on_ace_clim_latlev_contourf( climdata_in, title_in, c_grid )
%A function to plot the climatology data for ACE. The function uses
%'pcolor'.
%
% *INPUT*
%           climdata_in: ARRAY - the data to be plotted. It should be on a
%           grid that matches the lat x alt grid of the climatology data
%           (48x36), which is hard-coded below.
%
%
%           vmrzonvar_dif: STRING - the title of the figure  
%
% *OUTPUT*
%           makes a plot 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define some things
fs = 12;
% ytickspace = [10^-4 10^-2 1 10^2];
ytickspace = [10^-4 10^-3 10^-2 10^-1 1 10 30 60 10^2 200 500];
% yticknames = ['10^-4','10^-2',]
lat = flip(-87.5:5:87.5);
plev = [300 250 200 170 150 130 115 100 90 80 70 50 30 20 ...
    15 10 7 5 3 2 1.5 1 0.7 0.5 0.3 0.2 0.15 0.1]';

%% make the plot for the input variable
if nargin > 2
    cgrid = c_grid;
else
    colourmin = min(climdata_in(:));
    colourmax = max(climdata_in(:));
    numcolours = 10;
    cgrid = linspace(colourmin,colourmax,numcolours);
end
[C,h] = contourf(lat, plev, climdata_in, cgrid);
% h.LevelList = 10.^h.LevelList;
% h.ZData = 10.^h.ZData;
% clabel(C,h)
if nargin > 1
    title(title_in)
end

    
set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YMinorGrid','OFF', 'XGrid', 'ON');
yticks(ytickspace);
xlabel('lat bins [deg N]')
ylabel('pressure [hPa]')
%     cmap = lbmap(80,'RedBlue');
%     colormap(cmap);
%     caxis([-cmax,cmax]);
set(gca,'FontSize',fs)
%
end

