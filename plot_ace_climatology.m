function [ ] = plot_ace_climatology( climstruct_in )
%A function to plot the climatology data for ACE. The function uses
%'pcolor'.
%
% *INPUT*
%           climstruct_in: STRUCTURE - the data to be plotted. The structure
%           is created by running one of the 'make_ace_climatology...'
%           matlab functions.
%
% *OUTPUT*
%           makes a plot 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define some things
clim = climstruct_in;
gas = clim.gas;
climdata = clim.vmr_zonal;
if isfield(clim,'lat')
    lat = clim.lat;
%     lattitle = 'lat';
    latlabel = 'latitude [deg N]';
elseif isfield(clim, 'eql')
    lat = clim.eql;
%     lattitle = 'eql';
    latlabel = 'equivalent latitude [deg N]';
end
seasons = {'DJF', 'MAM', 'JJA', 'SON'};
if isfield(clim,'climatology_type')
    climtype = clim.climatology_type;
    climtime = clim.time;
    switch climtype
        case 'calendar_month'
            climtitle = sprintf('%s climatology, %02.0f', gas, climtime);
        case 'season'
            climtitle = sprintf('%s climatology, %s', gas, seasons{climtime});
        case 'serial_month'
            climtitle = sprintf('%s climatology, %i-%02.0f', gas, climtime(1), climtime(2));
    end
end
plev = clim.pressure_hPa;
fs = 16;
% yticknames = ['10^-4','10^-2',]
ytickspace = [10^-4 10^-2 1 10^2];

%% make the plot for the input variable
figure, pcolor(lat, plev, climdata);
if isfield(clim,'climatology_type')
    title(climtitle)
else
    title(gas,'interpreter','none')
end
set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on');
yticks(ytickspace);
xlabel(latlabel)
ylabel('pressure [hPa]')
c = colorbar;
% caxis([0,100]); % ONLY FOR SPECIFIC CASES***************************************************************
% c = colorbar('location','eastoutside','position',[0.923 0.309 0.016 0.35]);
if strcmp(gas,'T')
    title(c,'K [deg]','Position', [12 -23 0])
else
    title(c,'VMR','Position', [12 -23 0])
end
%     cmap = lbmap(80,'RedBlue');
%     colormap(cmap);
%     caxis([-cmax,cmax]);
set(gca,'FontSize',fs)
%
end

