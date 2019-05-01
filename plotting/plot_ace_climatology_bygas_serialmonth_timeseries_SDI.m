function [ sdates, vmrzon, pace, lat ] = plot_ace_climatology_bygas_serialmonth_timeseries_SDI( gasname_in, years_in, lat_minmax, do_plot )
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
do_old_plot = 1;

%% define some things
home_linux = '/home/niall/Dropbox/climatology/'; %#ok<NASGU>
home_mac = '/Users/niall/Dropbox/climatology/'; %#ok<NASGU>
home_windows = 'C:\Users\ryann\jaho\'; %#ok<NASGU>
clim_dir = 'C:\Users\ryann\ACE\ACE_climdata_netcdf_SDI_20190417\';
clim_dir_old = 'C:\Users\ryann\ACE\SDI-ISSI\';
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
% months = 1:12;
% lmonths = length(months);
sdates = nan(1,12*lyears);

vmrzon = nan(36,28,12*lyears);
vmrzon_error = nan(36,28,12*lyears);
vmrzon_obscount = nan(36,28,12*lyears);

if nargin > 3
    yplot = do_plot;
else
    yplot = 1;
end

file_pre = 'SPARC_DI_T2Mz_'; % SPARC_DI_T2Mz_C2H2_2004_ACEFTS_v3.6_i01.nc
file_post = '_ACEFTS_v3.6_i01.nc';

file_pre_old = 'SPARC_DI_T2Mz_'; % SPARC_DI_T2Mz_SF6_2004_ACEFTS_v2.2_i01.nc
file_post_old = '_ACEFTS_v2.2_i02.nc';

%% get the name of the gas for the input folder name
if length(gasname) > 10 && strcmp(gasname(end-8:end-6), 'sap') % for '_sap_s10am' and '_sap_s10pm', etc.
    gasname_short = gasname(1:end-10);
    gasname_old = gasname_short;
    file_post_old = '_ACEFTS_v2.2_SCALED10AM_i02.nc';
elseif length(gasname) > 7 && strcmp(gasname(end-5:end), 'sap_am') % for '_sap_am' and '_sap_pm'
    gasname_short = gasname(1:end-7);
    gasname_old = strcat(gasname_short, '_am');
elseif length(gasname) > 7 && strcmp(gasname(end-5:end), 'sap_pm') % for '_sap_am' and '_sap_pm'
    gasname_short = gasname(1:end-7);
    gasname_old = strcat(gasname_short, '_pm');
elseif length(gasname) > 4 && strcmp(gasname(end-2:end), 'sap') % for '_sap'
    gasname_short = gasname(1:end-4);
    gasname_old = gasname_short;
elseif length(gasname) > 3 && strcmp(gasname(end-1:end), 'am') % for '_am'
    gasname_short = gasname(1:end-3);
    gasname_old = gasname;
elseif length(gasname) > 3 && strcmp(gasname(end-1:end), 'pm') % for '_pm'
    gasname_short = gasname(1:end-3);
    gasname_old = gasname;
else
    gasname_short = gasname;
    gasname_old = gasname;
end

%         20km, 30km, 40km, 50km
%         50hPa, 10hPa, 2hPa, .5hPa
ipplot = [17, 21, 25, 29]; % the indexes of the pressure levels on which to plot
ipplot = [12, 16, 19, 22]; % the indexes of the pressure levels on which to plot
ipplot = flip(ipplot); % flip around the pressure levels for plotting
nanout = nan(36, 28, 12);

%% loop through the years and months and fill in the cells of arrays
ij = 1:12;
for j = 1:lyears
    %         j
    %             ij
    %         sdates(ij) = datenum(yearsin(j),i,1); % the 15th of each month for each year;
    
    filein = strcat( clim_dir, gasname_short,'/', file_pre, gasname, sprintf('_%i',yearsin(j)), file_post)
    if exist(filein,'file') ~= 2
        fprintf('There is no file for %s_%i. Moving on...\n', gasname, yearsin(j))
        vmrzon(:,:,ij) = nanout; % this is hardcoded here, for now.
        vmrzon_error(:,:,ij) = nanout;
        vmrzon_obscount(:,:,ij) = nanout;
    else
        vmrzon(:,:,ij) = ncread(filein, gasname_short);
        vmrzon_error(:,:,ij) = ncread(filein, strcat(gasname_short, '_STD'));
        vmrzon_obscount(:,:,ij) = ncread(filein, strcat(gasname_short, '_NR'));
        sdates(ij) = ncread(filein, 'time');
    end
    if do_old_plot == 1
        filein_old = strcat( clim_dir_old, gasname_short,'/', file_pre_old, gasname_old, sprintf('_%i',yearsin(j)), file_post_old)
        if exist(filein_old,'file') ~= 2
            fprintf('There is no file for %s_%i. Moving on...\n', gasname_old, yearsin(j))
            vmrzon_old(:,:,ij) = nanout; % this is hardcoded here, for now.
            vmrzon_error_old(:,:,ij) = nanout;
            vmrzon_obscount_old(:,:,ij) = nanout;
        else
            vmrzon_old(:,:,ij) = ncread(filein_old, gasname_short);
            vmrzon_error_old(:,:,ij) = ncread(filein_old, strcat(gasname_short, '_STD'));
            vmrzon_obscount_old(:,:,ij) = ncread(filein_old, strcat(gasname_short, '_NR'));
        end
    end
    ij = ij + 12;
end
vmrzon(vmrzon == -999) = nan;
vmrzon_error(vmrzon_error == -999) = nan;
vmrzon_obscount(vmrzon_obscount == -999) = nan;
if do_old_plot == 1
    vmrzon_old(vmrzon_old == -999) = nan;
    vmrzon_error_old(vmrzon_error_old == -999) = nan;
    vmrzon_obscount_old(vmrzon_obscount_old == -999) = nan;
end

sdates = datenum(1950, 1, 1) + sdates;
pace = ncread(filein, 'plev');
% pace_old = ncread(filein_old, 'plev');
zace = p2z_waccm(pace * 100)/1000;
% pace = clim.pressure_hPa;
% zace = clim.altitude_km_mean;
% return
%% subset according to the chosen lat limits
% average over some latitude bins if needed
lat_bounds = fliplr(-90:5:90);
if latmax > latmin
    latmin_old = latmin;
    latmax_old = latmax;
    latmin = latmax_old;
    latmax = latmin_old;
end
ilatmin = find(lat_bounds == latmin);
ilatmax = find(lat_bounds == latmax) - 1;

vmrzon = squeeze(nanmean(vmrzon(ilatmin:ilatmax,:,:), 1 )) * 1e9;
vmrzon_error = squeeze(nanmean(vmrzon_error(ilatmin:ilatmax,:,:), 1 )) * 1e9; % this needs to be fixed for proper error propagation
vmrzon_obscount = squeeze(nansum(vmrzon_obscount(ilatmin:ilatmax,:,:), 1 ));
if do_old_plot == 1
    vmrzon_old = squeeze(nanmean(vmrzon_old(ilatmin:ilatmax,:,:), 1 )) * 1e9;
    vmrzon_error_old = squeeze(nanmean(vmrzon_error_old(ilatmin:ilatmax,:,:), 1 )) * 1e9; % this needs to be fixed for proper error propagation
    vmrzon_obscount_old = squeeze(nansum(vmrzon_obscount_old(ilatmin:ilatmax,:,:), 1 ));
end

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
%         whos sdates vmrzon vmrzon_error
        errorbar(axx(i),sdates, vmrzon(ipplot(i),:), vmrzon_error(ipplot(i),:),'bo', 'Linewidth', lw, 'MarkerSize', ms )
        if do_old_plot == 1
            hold on
            errorbar(axx(i),sdates, vmrzon_old(ipplot(i),:), vmrzon_error_old(ipplot(i),:),'ro', 'Linewidth', lw, 'MarkerSize', ms )
        end
        dynamicDateTicks([],'x','mm')
        xlim([sdates(1), sdates(end)])
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
        stem(axx(i),sdates, vmrzon_obscount(ipplot(i),:),'bo', 'Linewidth', lw, 'MarkerSize', ms )
        if do_old_plot == 1
            hold on
            stem(axx(i),sdates, vmrzon_obscount_old(ipplot(i),:),'ro', 'Linewidth', lw, 'MarkerSize', ms )
        end
        dynamicDateTicks([],'x','mm')
        xlim([sdates(1), sdates(end)])
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
    
    %     %% contour plot
    %     ylim1 = 1;
    %     ylim2 = 10^3;
    %     fs = 12;
    %     %     colourmin = min(vmrzon(:));
    %     colourmin = 0;
    %     %     colourmax = 0.346;
    %     vmrzon_short = vmrzon(6:33,:);
    %     colourmax = max(vmrzon_short(:))
    %     numcolours = 10;
    %     cgrid = linspace(colourmin,colourmax,numcolours);
    %     %     fs_title = 8;
    %     ytickspace = [10^-4 10^-2 1 10 10^2];
    %     %         cgrid = linspace(0,4,17);
    %     % %         cgrid = [linspace(0,1,5),linspace(1.2,4,5)] % the values corresponding to the contour colours
    %     % Do for the zonal vmr
    %     %     figi = randi(100);
    %     figii = figi + 2;
    %     figure(figii), set(gcf,'Position', [-1163 506 874 309])
    %
    %     ylim2 = 300;
    %     ylim1 = 10^-1;
    %
    %     figure(figii);
    %     disp('plotting contour series')
    %     box on;
    %     contourf(sdates, pace, vmrzon, cgrid); caxis([cgrid(1),cgrid(end)])
    %     set(gca,'color',[0 0 0]);
    %     title(sprintf('%s zonal VMR, %i-%i%c latitude', gasname, latmin, latmax, char(176)))
    %     colormap(parula(length(cgrid) - 1));
    %     c = colorbar('YTick', cgrid, 'Yticklabel', round2(cgrid,0.001), 'Location','eastoutside');
    %     title(c,'VMR [ppbv]')
    %     %     c = colorbar('YTick', cgrid, 'location','eastoutside','position',[0.923 0.309 0.016 0.35]);
    %     %     title(c,'VMR','Position', [8 -23 0])
    %     set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YGrid', 'ON','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);
    %     yticks(ytickspace);
    %     ylim([ylim1 ylim2]); % colormap(parula(length(cgrid) - 1));
    %     %                 colorbar('YTick', cgrid, 'Yticklabel', c, 'Location','southoutside')
    %     xlabel('year')
    %     ylabel('pressure [hPa]')
    %     dynamicDateTicks([],'x','mmyy')
end
%
end

