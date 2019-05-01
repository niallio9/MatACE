function [ vmrzon, pace, lat ] = plot_ace_climatology_bygas_serialmonth_timeseries( gasname_in, years_in, lat_minmax, do_plot )
%A funcion to plot the serialmonth climatologies as a timeseries.

% *INPUT*
%           gasname_in: STRING - the name of the gas for which you want to
%           compare the climatologies.
%
%           years_in: STRING - the years over which you want to plot.
%
%           lat_minmax: STRING - the minimum and maximum of the latitude
%           range that you want to plot: [min, max].
%
% *OUTPUT*
%           vmr_out:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define which species you will be plotting
gasname = gasname_in;

%% define some things
home_linux = '/home/niall/Dropbox/climatology/'; %#ok<NASGU>
home_mac = '/Users/niall/Dropbox/climatology/'; %#ok<NASGU>
home_windows = 'C:\Users\ryann\jaho\'; %#ok<NASGU>
clim_dir = 'C:\Users\ryann\ACE\climdata\scaled_with_all_data\';
% clim_dir = '/Users/niall/ACE/climdata_pratmo_scaled/';

% clim_dir = '/Volumes/Seagate Backup Plus Drive/ACE/climdata/';
% clim_dir = '/Users/niall/Dropbox/'; %#ok<NASGU>

latmin = lat_minmax(1);
latmax = lat_minmax(2);
if mod(latmin,5) ~= 0 || mod(latmax,5) ~= 0
    error('you must choose latitude limits that are a multiple of 5 degrees, within [-90, 90]')
end
yearsin = years_in;
lyears = length(yearsin);
% months = [3,4,8,9];
months = 1:12;
% lmonths = length(months);
sdates = nan(1,12*lyears);

vmrzon = nan(48,36,12*lyears);
vmrzon_error = nan(48,36,12*lyears);
vmrzon_obscount = nan(48,36,12*lyears);

if nargin > 3
    yplot = do_plot;
else
    yplot = 1;
end

file_pre = 'ACEFTS_CLIM_v3_lat_'; % ACEFTS_CLIM_v3_O3_2008_12.mat
file_post = '.mat';
%         20km, 30km, 40km, 50km
%         50hPa, 10hPa, 2hPa, .5hPa
ipplot = [17, 21, 25, 29]; % the indexes of the pressure levels on which to plot
ipplot = [17, 21, 24, 27]; % the indexes of the pressure levels on which to plot
ipplot = flip(ipplot); % flip around the pressure levels for plotting

%% loop through the years and months and fill in the cells of arrays
ij = 0;
for j = 1:lyears
    %         j
    for i = 1:12 % go through all months
        ij = ij+1;
        %             ij
        sdates(ij) = datenum(yearsin(j),i,1); % the 15th of each month for each year;
        if  ismember(i,months)
            filein = strcat( clim_dir, gasname,'/','serial_month','/', file_pre, gasname, sprintf('_%i_%0.2i',yearsin(j), i), file_post)
            if exist(filein,'file') ~= 2
                fprintf('There is no file for %s_%i_%0.2i. Moving on...\n', gasname, yearsin(j), i)
                vmrzon(:,:,ij) = nan(48,36); % this is hardcoded here, for now.
                vmrzon_error(:,:,ij) = nan(48,36);
                vmrzon_obscount(:,:,ij) = nan(48,36);
            else
                clim = load(filein); clim = clim.climstruct; % the variable is called climstruct in the new data
%                 clim = reduce_climstruct_data_by_obs_nr(clim, 5);
                vmrzon(:,:,ij) = clim.vmr_zonal;
                vmrzon_error(:,:,ij) = sqrt(clim.vmr_zonal_var);
                vmrzon_obscount(:,:,ij) = clim.obs_count;
            end
        else
            vmrzon(:,:,ij) = nan(48,36); % this is hardcoded here, for now.
            vmrzon_error(:,:,ij) = nan(48,36);
            vmrzon_obscount(:,:,ij) = nan(48,36);
        end
    end
end
pace = clim.pressure_hPa;
zace = clim.altitude_km_mean;
% return
%% subset according to the chosen lat limits
% average over some latitude bins if needed
lat_bounds = clim.lat_bounds;
if latmax < latmin
    latmin_old = latmin;
    latmax_old = latmax;
    latmin = latmax_old;
    latmax = latmin_old;
end
ilatmin = find(lat_bounds == latmin);
ilatmax = find(lat_bounds == latmax) - 1;

vmrzon = squeeze(nanmean(vmrzon(:,ilatmin:ilatmax,:), 2 )) * 1e9;
vmrzon_error = squeeze(nanmean(vmrzon_error(:,ilatmin:ilatmax,:), 2 )) * 1e9; % this needs to be fixed for proper error propagation
vmrzon_obscount = squeeze(nansum(vmrzon_obscount(:,ilatmin:ilatmax,:), 2 ));

%% Make the plots if you want
lw = 1;
ms = 6;
% whos
%% line plot
nalt = length(ipplot);

if yplot == 1
    % Do for the zonal vmr
        figi = randi(100);
%     figi = 19;
    figure(figi), set(gcf,'Position', [97,49,852,630])
    axx = nan(nalt,1);
    disp('plotting by altitude')
    for i = 1:nalt
        figure(figi), axx(i) = subplot(nalt,1,i); hold(axx(i),'on')
        errorbar(axx(i),sdates, vmrzon(ipplot(i),:), vmrzon_error(ipplot(i),:),'bo-', 'Linewidth', lw, 'MarkerSize', ms )
        dynamicDateTicks([],'x','mm')
        %             ylim([0,4.0])
        ytext = max(vmrzon(ipplot(i),:))*1.1; % the position for the text in the plots, in ppbv
        text(sdates(5), ytext, sprintf('%0.1f hPa (~%0.1f km)', pace(ipplot(i)), zace(ipplot(i)) ));
        ylabel('VMR [ppb]')
%         set(gca,'XTickLabel',[]);
        if i == 1
            title(sprintf('%s zonal VMR, %i-%i%c latitude', gasname, latmin, latmax, char(176)))
        end
%         dynamicDateTicks([],'x','mm')
        if i == 4
            dynamicDateTicks([],'x','mm')
            %                 dynamicDateTicks
            xlabel('year')
            linkaxes(axx,'x')
        else
%                             datetick('x','mmyy')
            set(gca,'xticklabel',{[]})
        end
    end
    
    % Do for the zonal vmr
        figi = figi + 1;
%     figi = 19;
    figure(figi), set(gcf,'Position', [97,49,852,630])
    axx = nan(nalt,1);
    disp('plotting by altitude')
    for i = 1:nalt
        figure(figi), axx(i) = subplot(nalt,1,i); hold(axx(i),'on')
        stem(axx(i),sdates, vmrzon_obscount(ipplot(i),:),'bo-', 'Linewidth', lw, 'MarkerSize', ms )
        dynamicDateTicks([],'x','mm')
        %             ylim([0,4.0])
        ytext = max(vmrzon_obscount(ipplot(i),:))*1.1; % the position for the text in the plots, in ppbv
        text(sdates(5), ytext, sprintf('%0.1f hPa (~%0.1f km)', pace(ipplot(i)), zace(ipplot(i)) ));
        ylabel('# of obs''')
%         set(gca,'XTickLabel',[]);
        if i == 1
            title(sprintf('%s zonal observation count, %i-%i%c latitude', gasname, latmin, latmax, char(176)))
        end
%         dynamicDateTicks([],'x','mm')
        if i == 4
            dynamicDateTicks([],'x','mm')
            %                 dynamicDateTicks
            xlabel('year')
            linkaxes(axx,'x')
        else
%                             datetick('x','mmyy')
            set(gca,'xticklabel',{[]})
        end
    end
    
    %% contour plot
    ylim1 = 1;
    ylim2 = 10^3;
    fs = 12;
%     colourmin = min(vmrzon(:));
    colourmin = 0;
%     colourmax = 0.346;
    vmrzon_short = vmrzon(6:33,:);
    colourmax = max(vmrzon_short(:))
    numcolours = 10;
    cgrid = linspace(colourmin,colourmax,numcolours);
%     fs_title = 8;
    ytickspace = [10^-4 10^-2 1 10 10^2];
    %         cgrid = linspace(0,4,17);
    % %         cgrid = [linspace(0,1,5),linspace(1.2,4,5)] % the values corresponding to the contour colours
    % Do for the zonal vmr
    %     figi = randi(100);
    figii = figi + 2;
    figure(figii), set(gcf,'Position', [-1163 506 874 309])
    
    ylim2 = 300;
    ylim1 = 10^-1;
    
    figure(figii);
    disp('plotting contour series')
    box on;
    contourf(sdates, pace, vmrzon, cgrid); caxis([cgrid(1),cgrid(end)])
    set(gca,'color',[0 0 0]);
    title(sprintf('%s zonal VMR, %i-%i%c latitude', gasname, latmin, latmax, char(176)))
    colormap(parula(length(cgrid) - 1));
    c = colorbar('YTick', cgrid, 'Yticklabel', round2(cgrid,0.001), 'Location','eastoutside');
    title(c,'VMR [ppbv]')
%     c = colorbar('YTick', cgrid, 'location','eastoutside','position',[0.923 0.309 0.016 0.35]);
%     title(c,'VMR','Position', [8 -23 0])
    set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YGrid', 'ON','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);
    yticks(ytickspace);
    ylim([ylim1 ylim2]); % colormap(parula(length(cgrid) - 1));
    %                 colorbar('YTick', cgrid, 'Yticklabel', c, 'Location','southoutside')
    xlabel('year')
    ylabel('pressure [hPa]')
    dynamicDateTicks([],'x','mmyy')
end
%
end

